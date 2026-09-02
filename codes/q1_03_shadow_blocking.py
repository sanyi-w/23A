# ============================================================
# q1_03_shadow_blocking.py
#
# 2023A 问题一：
# 阴影遮挡效率 eta_sb
#
# 建模结构：
#
#   太阳方向 s
#       ↓
#   镜面法向 n_i
#       ↓
#   镜面真实三维矩形
#       ↓
#   ┌──────── Receiver Shadow
#   │
#   ├──────── Heliostat Shadowing
#   │
#   └──────── Heliostat Blocking
#               ↓
#        投影至目标镜局部 (u,v) 平面
#               ↓
#        多边形 intersection / union
#               ↓
#        eta_sb = 1 - A_loss / A_mirror
#
#
# 注意：
# 1. 当前采用太阳中心平行光模型。
# 2. Blocking 采用平面镜对应的平行反射方向 r_i。
# 3. 只考虑题目给出明确尺寸的圆柱集热器阴影；
#    不人为假定吸收塔支撑结构的横截面尺寸。
# 4. 当前仍未计算截断效率 eta_trunc。
# ============================================================

from pathlib import Path
from datetime import date
from time import perf_counter
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from shapely.geometry import Polygon, MultiPoint, box, GeometryCollection
from shapely.ops import unary_union


# ============================================================
# 0. 路径
# ============================================================

DATA_PATH = Path("data/附件.xlsx")

FIG_DIR = Path("figs")
RESULT_DIR = Path("results")

FIG_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. 题目常数
# ============================================================

# 镜场纬度
PHI = np.deg2rad(39.4)

# 镜面
MIRROR_WIDTH = 6.0
MIRROR_HEIGHT = 6.0
MIRROR_AREA = MIRROR_WIDTH * MIRROR_HEIGHT

# Q1 中镜中心安装高度
HELIOSTAT_Z = 4.0

# 镜面半对角线
# 用于候选镜保守包围球筛选
MIRROR_BOUND_RADIUS = 0.5 * np.sqrt(
    MIRROR_WIDTH ** 2 + MIRROR_HEIGHT ** 2
)

# 两块相同镜面的包围球半径和
PAIR_BOUND_RADIUS = 2.0 * MIRROR_BOUND_RADIUS

# 集热器
RECEIVER_CENTER = np.array([0.0, 0.0, 80.0])

RECEIVER_RADIUS = 3.5
RECEIVER_HEIGHT = 8.0

RECEIVER_Z_BOTTOM = 76.0
RECEIVER_Z_TOP = 84.0

# 圆柱集热器包围球半径
RECEIVER_BOUND_RADIUS = np.sqrt(
    RECEIVER_RADIUS ** 2
    + (RECEIVER_HEIGHT / 2.0) ** 2
)

# 圆柱边界离散数量
# 后续可做 72 / 144 / 288 收敛验证
RECEIVER_N_THETA = 144

# 目标镜在自身局部坐标中的固定区域
TARGET_RECT = box(
    -MIRROR_WIDTH / 2.0,
    -MIRROR_HEIGHT / 2.0,
    MIRROR_WIDTH / 2.0,
    MIRROR_HEIGHT / 2.0
)

# 数值容差
EPS = 1e-10
AREA_EPS = 1e-10

# 60 个题目工况
MONTHS = list(range(1, 13))
SOLAR_TIMES = [9.0, 10.5, 12.0, 13.5, 15.0]


# ============================================================
# 2. 数据读取
# ============================================================

def read_heliostat_coordinates(path):
    """
    读取附件中的 x、y 坐标。
    """

    xls = pd.ExcelFile(path)

    sheet_name = xls.sheet_names[0]

    df = pd.read_excel(
        path,
        sheet_name=sheet_name
    )

    df = (
        df
        .dropna(axis=0, how="all")
        .dropna(axis=1, how="all")
    )

    numeric_cols = df.select_dtypes(
        include=[np.number]
    ).columns.tolist()

    x_candidates = [
        c for c in numeric_cols
        if "x" in str(c).lower()
    ]

    y_candidates = [
        c for c in numeric_cols
        if "y" in str(c).lower()
    ]

    if x_candidates and y_candidates:
        x_col = x_candidates[0]
        y_col = y_candidates[0]

    elif len(numeric_cols) >= 2:
        x_col = numeric_cols[0]
        y_col = numeric_cols[1]

        print(
            "[提示] 未明确识别 x/y 列名，"
            "暂取前两个数值列。"
        )

    else:
        raise ValueError(
            "附件中无法识别 x、y 坐标列。"
        )

    x = df[x_col].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)

    x = x[mask]
    y = y[mask]

    H = np.column_stack([
        x,
        y,
        np.full(len(x), HELIOSTAT_Z)
    ])

    return H


# ============================================================
# 3. 太阳模型
# ============================================================

def day_from_spring_equinox(month, day=21):

    reference = date(2023, 3, 21)
    current = date(2023, month, day)

    return (current - reference).days


def solar_declination(D):

    sin_delta = (
        np.sin(2.0 * np.pi * D / 365.0)
        *
        np.sin(np.deg2rad(23.45))
    )

    sin_delta = np.clip(
        sin_delta,
        -1.0,
        1.0
    )

    return np.arcsin(sin_delta)


def solar_hour_angle(ST):

    return (
        np.pi / 12.0
        *
        (ST - 12.0)
    )


def solar_direction_enu(delta, omega):

    s = np.array([
        -np.cos(delta) * np.sin(omega),

        (
            np.cos(PHI) * np.sin(delta)
            -
            np.sin(PHI)
            * np.cos(delta)
            * np.cos(omega)
        ),

        (
            np.sin(PHI) * np.sin(delta)
            +
            np.cos(PHI)
            * np.cos(delta)
            * np.cos(omega)
        )
    ])

    return s / np.linalg.norm(s)


# ============================================================
# 4. 镜 -> 集热器方向
# ============================================================

def receiver_direction(H):

    vec = RECEIVER_CENTER - H

    d_hr = np.linalg.norm(
        vec,
        axis=1
    )

    r = vec / d_hr[:, None]

    return r, d_hr


# ============================================================
# 5. 定日镜法向
# ============================================================

def heliostat_normals(s, r):

    temp = r + s[None, :]

    temp_norm = np.linalg.norm(
        temp,
        axis=1,
        keepdims=True
    )

    return temp / temp_norm


# ============================================================
# 6. 镜面局部坐标
# ============================================================

def mirror_local_frames(normals):
    """
    构造：
        e_w：镜面水平边方向
        e_h：镜面高度方向

    题目条件：
        上下两边始终平行地面
    """

    N = len(normals)

    k = np.array([0.0, 0.0, 1.0])

    k_array = np.repeat(
        k[None, :],
        N,
        axis=0
    )

    # 水平且位于镜面内
    e_w = np.cross(
        k_array,
        normals
    )

    norms = np.linalg.norm(
        e_w,
        axis=1
    )

    # 极特殊情况：
    # 镜面法向恰好竖直
    # 此时任意水平轴都可作为宽度方向
    singular = norms < EPS

    e_w[~singular] /= (
        norms[~singular, None]
    )

    e_w[singular] = np.array(
        [1.0, 0.0, 0.0]
    )

    # 面内高度方向
    e_h = np.cross(
        normals,
        e_w
    )

    e_h /= np.linalg.norm(
        e_h,
        axis=1,
        keepdims=True
    )

    return e_w, e_h


# ============================================================
# 7. 镜面四角点
# ============================================================

def mirror_vertices(H, e_w, e_h):
    """
    顶点顺序：
        左下 -> 右下 -> 右上 -> 左上
    """

    half_w = MIRROR_WIDTH / 2.0
    half_h = MIRROR_HEIGHT / 2.0

    uv = np.array([
        [-half_w, -half_h],
        [ half_w, -half_h],
        [ half_w,  half_h],
        [-half_w,  half_h],
    ])

    vertices = (
        H[:, None, :]
        +
        uv[None, :, 0, None]
        * e_w[:, None, :]
        +
        uv[None, :, 1, None]
        * e_h[:, None, :]
    )

    return vertices


# ============================================================
# 8. 三维点 -> 目标镜局部二维坐标
# ============================================================

def to_local_uv(points, H_i, e_w_i, e_h_i):

    rel = points - H_i

    u = rel @ e_w_i
    v = rel @ e_h_i

    return np.column_stack([
        u,
        v
    ])


# ============================================================
# 9. 按 lambda >= 0 裁剪三维源多边形
# ============================================================

def clip_polygon_forward(points, lambdas):
    """
    保留沿指定投影方向能够真正到达目标平面的部分。

    lambda >= 0：
        说明投影方向正确。

    这里是在三维源多边形上进行一次线性半空间裁剪。
    """

    points = np.asarray(points)
    lambdas = np.asarray(lambdas)

    out = []

    n = len(points)

    for k in range(n):

        p1 = points[k]
        p2 = points[(k + 1) % n]

        l1 = lambdas[k]
        l2 = lambdas[(k + 1) % n]

        inside1 = l1 >= -EPS
        inside2 = l2 >= -EPS

        # inside -> inside
        if inside1 and inside2:

            out.append(p2)

        # inside -> outside
        elif inside1 and not inside2:

            t = l1 / (l1 - l2)

            p_int = (
                p1
                +
                t * (p2 - p1)
            )

            out.append(p_int)

        # outside -> inside
        elif (not inside1) and inside2:

            t = l1 / (l1 - l2)

            p_int = (
                p1
                +
                t * (p2 - p1)
            )

            out.append(p_int)
            out.append(p2)

        # outside -> outside
        else:
            pass

    if len(out) < 3:
        return None

    return np.asarray(out)


# ============================================================
# 10. 邻镜投影到目标镜
# ============================================================

def project_polygon_to_target(
    source_vertices,
    H_i,
    n_i,
    e_w_i,
    e_h_i,
    projection_direction
):
    """
    projection_direction：

    Shadowing:
        -s

    Blocking:
        -r_i

    返回：
        已经和目标镜矩形裁剪后的二维 Shapely geometry
    """

    d = projection_direction

    denom = np.dot(
        n_i,
        d
    )

    if abs(denom) < EPS:
        return None

    numer = (
        (source_vertices - H_i)
        @ n_i
    )

    lambdas = (
        -numer
        / denom
    )

    # 保留真实投影方向 lambda >= 0 的部分
    clipped_source = clip_polygon_forward(
        source_vertices,
        lambdas
    )

    if clipped_source is None:
        return None

    numer2 = (
        (clipped_source - H_i)
        @ n_i
    )

    lambda2 = (
        -numer2
        / denom
    )

    projected = (
        clipped_source
        +
        lambda2[:, None] * d
    )

    uv = to_local_uv(
        projected,
        H_i,
        e_w_i,
        e_h_i
    )

    # ----------------------------------------
    # Level 2：
    # 先做 AABB 快速筛选
    # ----------------------------------------

    u_min = np.min(uv[:, 0])
    u_max = np.max(uv[:, 0])

    v_min = np.min(uv[:, 1])
    v_max = np.max(uv[:, 1])

    half_w = MIRROR_WIDTH / 2.0
    half_h = MIRROR_HEIGHT / 2.0

    if (
        u_max < -half_w - EPS
        or
        u_min > half_w + EPS
        or
        v_max < -half_h - EPS
        or
        v_min > half_h + EPS
    ):
        return None

    # 源镜是凸多边形，
    # 平行投影仍然是凸多边形
    polygon = MultiPoint(
        uv.tolist()
    ).convex_hull

    if polygon.is_empty:
        return None

    if polygon.geom_type != "Polygon":
        return None

    overlap = polygon.intersection(
        TARGET_RECT
    )

    if (
        overlap.is_empty
        or
        overlap.area <= AREA_EPS
    ):
        return None

    return overlap


# ============================================================
# 11. 集热器圆柱边界点
# ============================================================

def build_receiver_boundary_points():

    theta = np.linspace(
        0.0,
        2.0 * np.pi,
        RECEIVER_N_THETA,
        endpoint=False
    )

    x = RECEIVER_RADIUS * np.cos(theta)
    y = RECEIVER_RADIUS * np.sin(theta)

    lower = np.column_stack([
        x,
        y,
        np.full(
            RECEIVER_N_THETA,
            RECEIVER_Z_BOTTOM
        )
    ])

    upper = np.column_stack([
        x,
        y,
        np.full(
            RECEIVER_N_THETA,
            RECEIVER_Z_TOP
        )
    ])

    return np.vstack([
        lower,
        upper
    ])


RECEIVER_BOUNDARY_POINTS = (
    build_receiver_boundary_points()
)


# ============================================================
# 12. 集热器阴影投影
# ============================================================

def receiver_shadow_polygon(
    H_i,
    n_i,
    e_w_i,
    e_h_i,
    s
):
    """
    将圆柱集热器沿 -s 投到目标镜 i。

    采用确定性的圆周边界离散，
    然后取凸包。

    不是 Monte Carlo。
    """

    d = -s

    denom = np.dot(
        n_i,
        d
    )

    if abs(denom) < EPS:
        return None

    points = RECEIVER_BOUNDARY_POINTS

    numer = (
        (points - H_i)
        @ n_i
    )

    lambdas = (
        -numer
        / denom
    )

    valid = lambdas >= -EPS

    if np.count_nonzero(valid) < 3:
        return None

    p = points[valid]
    lam = lambdas[valid]

    projected = (
        p
        +
        lam[:, None] * d
    )

    uv = to_local_uv(
        projected,
        H_i,
        e_w_i,
        e_h_i
    )

    # 快速 AABB
    if (
        np.max(uv[:, 0]) < -3.0
        or
        np.min(uv[:, 0]) > 3.0
        or
        np.max(uv[:, 1]) < -3.0
        or
        np.min(uv[:, 1]) > 3.0
    ):
        return None

    hull = MultiPoint(
        uv.tolist()
    ).convex_hull

    if hull.geom_type != "Polygon":
        return None

    overlap = hull.intersection(
        TARGET_RECT
    )

    if (
        overlap.is_empty
        or
        overlap.area <= AREA_EPS
    ):
        return None

    return overlap


# ============================================================
# 13. pairwise 中心距离平方
# ============================================================

def pairwise_distance_squared(H):
    """
    D2[i,j] = ||H_j - H_i||^2

    只存 N×N，
    避免保存 N×N×3 的三维差向量。
    """

    norm2 = np.sum(
        H * H,
        axis=1
    )

    D2 = (
        norm2[:, None]
        +
        norm2[None, :]
        -
        2.0 * H @ H.T
    )

    return np.maximum(
        D2,
        0.0
    )


# ============================================================
# 14. Shadowing 候选邻镜
# ============================================================

def shadow_candidate_lists(
    H,
    s,
    D2
):
    """
    Level 1：

    纵向：
        d_parallel

    横向：
        d_perp

    用镜面包围球构造保守光路圆柱。
    """

    # 每面镜中心沿太阳方向的坐标
    q = H @ s

    # d_parallel[i,j]
    # = (H_j - H_i) dot s
    d_parallel = (
        q[None, :]
        -
        q[:, None]
    )

    d_perp2 = (
        D2
        -
        d_parallel ** 2
    )

    d_perp2 = np.maximum(
        d_perp2,
        0.0
    )

    mask = (
        (d_perp2
         <= PAIR_BOUND_RADIUS ** 2 + EPS)
        &
        (d_parallel
         > -PAIR_BOUND_RADIUS)
    )

    np.fill_diagonal(
        mask,
        False
    )

    return [
        np.flatnonzero(mask[i])
        for i in range(len(H))
    ]


# ============================================================
# 15. Blocking 候选邻镜
# ============================================================

def blocking_candidate_lists(
    H,
    r,
    d_hr,
    D2
):
    """
    Blocking 方向 r_i 在 Q1 中与时间无关，
    所以候选列表可以预计算一次。
    """

    # M[i,j] = r_i dot H_j
    M = r @ H.T

    center_proj = np.sum(
        r * H,
        axis=1
    )

    d_parallel = (
        M
        -
        center_proj[:, None]
    )

    d_perp2 = (
        D2
        -
        d_parallel ** 2
    )

    d_perp2 = np.maximum(
        d_perp2,
        0.0
    )

    mask = (
        (
            d_perp2
            <= PAIR_BOUND_RADIUS ** 2 + EPS
        )
        &
        (
            d_parallel
            > -PAIR_BOUND_RADIUS
        )
        &
        (
            d_parallel
            <
            d_hr[:, None]
            + PAIR_BOUND_RADIUS
        )
    )

    np.fill_diagonal(
        mask,
        False
    )

    return [
        np.flatnonzero(mask[i])
        for i in range(len(H))
    ]


# ============================================================
# 16. 集热器阴影候选目标镜
# ============================================================

def receiver_shadow_candidate_mask(
    H,
    s
):
    """
    使用圆柱集热器的包围球进行保守筛选。

    集热器包围球半径：
        R_receiver_bound

    镜面包围球半径：
        R_mirror_bound

    构造沿 -s 的阴影管。
    """

    direction = -s

    delta = (
        H
        -
        RECEIVER_CENTER
    )

    d_parallel = (
        delta
        @ direction
    )

    perpendicular = (
        delta
        -
        d_parallel[:, None]
        * direction[None, :]
    )

    d_perp = np.linalg.norm(
        perpendicular,
        axis=1
    )

    tube_radius = (
        RECEIVER_BOUND_RADIUS
        +
        MIRROR_BOUND_RADIUS
    )

    mask = (
        (d_parallel > -tube_radius)
        &
        (d_perp <= tube_radius)
    )

    return mask


# ============================================================
# 17. 空 geometry
# ============================================================

def empty_geometry():

    return GeometryCollection()


# ============================================================
# 18. 单面目标镜 eta_sb
# ============================================================

def compute_eta_sb_for_target(
    i,
    H,
    normals,
    e_w,
    e_h,
    vertices,
    r,
    s,
    shadow_candidates,
    blocking_candidates,
    receiver_candidate,
    return_geometry=False
):
    """
    计算目标镜 i：

        receiver shadow
        + heliostat shadow
        + heliostat blocking

    最后统一求并集。
    """

    # --------------------------------------------------------
    # A. Receiver shadow
    # --------------------------------------------------------

    receiver_geom = None

    if receiver_candidate[i]:

        receiver_geom = receiver_shadow_polygon(
            H[i],
            normals[i],
            e_w[i],
            e_h[i],
            s
        )

    # --------------------------------------------------------
    # B. Heliostat shadowing
    # --------------------------------------------------------

    shadow_geoms = []

    for j in shadow_candidates[i]:

        geom = project_polygon_to_target(
            source_vertices=vertices[j],

            H_i=H[i],
            n_i=normals[i],

            e_w_i=e_w[i],
            e_h_i=e_h[i],

            projection_direction=-s
        )

        if geom is not None:
            shadow_geoms.append(geom)

    # --------------------------------------------------------
    # C. Blocking
    # --------------------------------------------------------

    blocking_geoms = []

    for j in blocking_candidates[i]:

        geom = project_polygon_to_target(
            source_vertices=vertices[j],

            H_i=H[i],
            n_i=normals[i],

            e_w_i=e_w[i],
            e_h_i=e_h[i],

            projection_direction=-r[i]
        )

        if geom is not None:
            blocking_geoms.append(geom)

    # --------------------------------------------------------
    # D. 分类求并
    # --------------------------------------------------------

    if shadow_geoms:
        shadow_union = unary_union(
            shadow_geoms
        )
    else:
        shadow_union = empty_geometry()

    if blocking_geoms:
        blocking_union = unary_union(
            blocking_geoms
        )
    else:
        blocking_union = empty_geometry()

    if receiver_geom is not None:
        receiver_union = receiver_geom
    else:
        receiver_union = empty_geometry()

    # --------------------------------------------------------
    # E. 全部无效区域统一求并
    # --------------------------------------------------------

    all_geoms = []

    if not receiver_union.is_empty:
        all_geoms.append(receiver_union)

    if not shadow_union.is_empty:
        all_geoms.append(shadow_union)

    if not blocking_union.is_empty:
        all_geoms.append(blocking_union)

    if all_geoms:

        loss_union = unary_union(
            all_geoms
        )

    else:

        loss_union = empty_geometry()

    # 理论上所有 geometry 已经裁到 TARGET_RECT
    loss_area = (
        0.0
        if loss_union.is_empty
        else loss_union.area
    )

    receiver_area = (
        0.0
        if receiver_union.is_empty
        else receiver_union.area
    )

    shadow_area = (
        0.0
        if shadow_union.is_empty
        else shadow_union.area
    )

    blocking_area = (
        0.0
        if blocking_union.is_empty
        else blocking_union.area
    )

    # 数值保护
    loss_area = np.clip(
        loss_area,
        0.0,
        MIRROR_AREA
    )

    eta_sb = (
        1.0
        -
        loss_area / MIRROR_AREA
    )

    result = {
        "eta_sb": eta_sb,

        "loss_area_m2": loss_area,

        "receiver_shadow_area_m2":
            receiver_area,

        "heliostat_shadow_area_m2":
            shadow_area,

        "blocking_area_m2":
            blocking_area,

        "shadow_candidate_count":
            len(shadow_candidates[i]),

        "blocking_candidate_count":
            len(blocking_candidates[i]),

        "actual_shadow_polygon_count":
            len(shadow_geoms),

        "actual_blocking_polygon_count":
            len(blocking_geoms),
    }

    if return_geometry:

        result["receiver_geom"] = (
            receiver_union
        )

        result["shadow_geom"] = (
            shadow_union
        )

        result["blocking_geom"] = (
            blocking_union
        )

        result["loss_geom"] = (
            loss_union
        )

    return result


# ============================================================
# 19. Shapely geometry 绘图
# ============================================================

def plot_shapely_geometry(
    ax,
    geom,
    facecolor,
    label,
    alpha=0.35
):

    if geom is None or geom.is_empty:
        return

    if geom.geom_type == "Polygon":

        polygons = [geom]

    elif geom.geom_type == "MultiPolygon":

        polygons = list(
            geom.geoms
        )

    else:
        return

    first = True

    for poly in polygons:

        xy = np.asarray(
            poly.exterior.coords
        )

        ax.fill(
            xy[:, 0],
            xy[:, 1],
            facecolor=facecolor,
            alpha=alpha,
            label=label if first else None
        )

        first = False


# ============================================================
# 20. 单工况计算
# ============================================================

def compute_condition(
    month,
    ST,
    H,
    r,
    d_hr,
    D2,
    blocking_candidates
):

    D = day_from_spring_equinox(
        month,
        21
    )

    delta = solar_declination(D)

    omega = solar_hour_angle(ST)

    s = solar_direction_enu(
        delta,
        omega
    )

    normals = heliostat_normals(
        s,
        r
    )

    e_w, e_h = mirror_local_frames(
        normals
    )

    vertices = mirror_vertices(
        H,
        e_w,
        e_h
    )

    shadow_candidates = shadow_candidate_lists(
        H,
        s,
        D2
    )

    receiver_candidate = (
        receiver_shadow_candidate_mask(
            H,
            s
        )
    )

    rows = []

    for i in range(len(H)):

        res = compute_eta_sb_for_target(
            i=i,

            H=H,
            normals=normals,
            e_w=e_w,
            e_h=e_h,
            vertices=vertices,

            r=r,
            s=s,

            shadow_candidates=
                shadow_candidates,

            blocking_candidates=
                blocking_candidates,

            receiver_candidate=
                receiver_candidate
        )

        rows.append({
            "month": month,
            "day": 21,
            "ST": ST,

            "mirror_id": i + 1,

            "x_m": H[i, 0],
            "y_m": H[i, 1],

            **res
        })

    df = pd.DataFrame(
        rows
    )

    state = {
        "D": D,
        "delta": delta,
        "omega": omega,
        "s": s,

        "normals": normals,
        "e_w": e_w,
        "e_h": e_h,
        "vertices": vertices,

        "shadow_candidates":
            shadow_candidates,

        "receiver_candidate":
            receiver_candidate
    }

    return df, state


# ============================================================
# 21. 工况结果图
# ============================================================

def plot_eta_sb_field(
    df,
    month,
    ST
):

    fig, ax = plt.subplots(
        figsize=(8, 8)
    )

    sc = ax.scatter(
        df["x_m"],
        df["y_m"],

        c=df["eta_sb"],

        s=18,

        vmin=0.0,
        vmax=1.0
    )

    cbar = fig.colorbar(
        sc,
        ax=ax
    )

    cbar.set_label(
        "Shadow-blocking efficiency"
    )

    ax.scatter(
        [0],
        [0],
        marker="*",
        s=180,
        label="Receiver tower"
    )

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    ax.set_xlabel(
        "x / m (East)"
    )

    ax.set_ylabel(
        "y / m (North)"
    )

    ax.set_title(
        f"Shadow-Blocking Efficiency\n"
        f"{month}/21  {ST:04.1f}"
    )

    ax.legend()

    fig.tight_layout()

    output = (
        FIG_DIR
        /
        (
            f"q1_result_eta_sb_"
            f"m{month:02d}_"
            f"t{ST:04.1f}.png"
        )
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

    return output


# ============================================================
# 22. 最差镜面的局部投影示意
# ============================================================

def plot_local_projection_example(
    df,
    state,
    H,
    r,
    blocking_candidates,
    month,
    ST
):

    # 选 eta_sb 最低的镜子
    row = df.loc[
        df["eta_sb"].idxmin()
    ]

    i = int(
        row["mirror_id"] - 1
    )

    result = compute_eta_sb_for_target(
        i=i,

        H=H,

        normals=state["normals"],
        e_w=state["e_w"],
        e_h=state["e_h"],
        vertices=state["vertices"],

        r=r,
        s=state["s"],

        shadow_candidates=
            state["shadow_candidates"],

        blocking_candidates=
            blocking_candidates,

        receiver_candidate=
            state["receiver_candidate"],

        return_geometry=True
    )

    fig, ax = plt.subplots(
        figsize=(7, 7)
    )

    # 目标镜
    x = [-3, 3, 3, -3, -3]
    y = [-3, -3, 3, 3, -3]

    ax.plot(
        x,
        y,
        linewidth=2,
        label="Target mirror"
    )

    plot_shapely_geometry(
        ax,
        result["receiver_geom"],
        facecolor="tab:orange",
        label="Receiver shadow"
    )

    plot_shapely_geometry(
        ax,
        result["shadow_geom"],
        facecolor="tab:blue",
        label="Heliostat shadow"
    )

    plot_shapely_geometry(
        ax,
        result["blocking_geom"],
        facecolor="tab:red",
        label="Blocking"
    )

    ax.set_xlim(
        -3.5,
        3.5
    )

    ax.set_ylim(
        -3.5,
        3.5
    )

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    ax.set_xlabel(
        "u / m"
    )

    ax.set_ylabel(
        "v / m"
    )

    ax.set_title(
        "Local Projection on Target Heliostat\n"
        f"Mirror {i + 1}, "
        f"{month}/21 {ST:04.1f}, "
        f"eta_sb={result['eta_sb']:.4f}"
    )

    ax.legend()

    ax.grid(
        alpha=0.25
    )

    fig.tight_layout()

    output = (
        FIG_DIR
        /
        "q1_verify_sb_local_projection.png"
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

    return output


# ============================================================
# 23. 摘要
# ============================================================

def summarize_condition(df):

    return {
        "eta_sb_mean":
            df["eta_sb"].mean(),

        "eta_sb_min":
            df["eta_sb"].min(),

        "eta_sb_max":
            df["eta_sb"].max(),

        "mean_receiver_shadow_area":
            df[
                "receiver_shadow_area_m2"
            ].mean(),

        "mean_heliostat_shadow_area":
            df[
                "heliostat_shadow_area_m2"
            ].mean(),

        "mean_blocking_area":
            df[
                "blocking_area_m2"
            ].mean(),

        "mean_loss_area":
            df[
                "loss_area_m2"
            ].mean(),

        "mean_shadow_candidates":
            df[
                "shadow_candidate_count"
            ].mean(),

        "mean_blocking_candidates":
            df[
                "blocking_candidate_count"
            ].mean(),

        "max_shadow_candidates":
            df[
                "shadow_candidate_count"
            ].max(),

        "max_blocking_candidates":
            df[
                "blocking_candidate_count"
            ].max(),
    }


# ============================================================
# 24. 主程序
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=["test", "full"],
        default="test"
    )

    args = parser.parse_args()

    start_time = perf_counter()

    # --------------------------------------------------------
    # A. 数据
    # --------------------------------------------------------

    H = read_heliostat_coordinates(
        DATA_PATH
    )

    N = len(H)

    print("=" * 70)
    print("Q1 Shadow / Blocking")
    print("=" * 70)

    print(
        f"Heliostat count N = {N}"
    )

    print(
        "Mirror bounding radius = "
        f"{MIRROR_BOUND_RADIUS:.6f} m"
    )

    print(
        "Pair screening radius = "
        f"{PAIR_BOUND_RADIUS:.6f} m"
    )

    # --------------------------------------------------------
    # B. 与时间无关的量
    # --------------------------------------------------------

    r, d_hr = receiver_direction(
        H
    )

    D2 = pairwise_distance_squared(
        H
    )

    print(
        "\nPrecomputing blocking candidates..."
    )

    blocking_candidates = (
        blocking_candidate_lists(
            H,
            r,
            d_hr,
            D2
        )
    )

    blocking_counts = np.array([
        len(x)
        for x in blocking_candidates
    ])

    print(
        "Mean blocking candidates = "
        f"{blocking_counts.mean():.3f}"
    )

    print(
        "Max blocking candidates = "
        f"{blocking_counts.max()}"
    )

    # --------------------------------------------------------
    # C. TEST 模式
    # --------------------------------------------------------

    if args.mode == "test":

        # 低太阳高度工况，
        # 更容易暴露 shadowing 问题
        month = 12
        ST = 9.0

        print(
            f"\nTEST condition: "
            f"{month}/21 {ST}"
        )

        df, state = compute_condition(
            month,
            ST,
            H,
            r,
            d_hr,
            D2,
            blocking_candidates
        )

        output_csv = (
            RESULT_DIR
            /
            "q1_sb_test.csv"
        )

        df.to_csv(
            output_csv,
            index=False,
            encoding="utf-8-sig"
        )

        summary = summarize_condition(
            df
        )

        print("\nTEST summary")

        for key, value in summary.items():

            print(
                f"{key:35s}"
                f"{value}"
            )

        # 数值约束检查
        assert (
            df["eta_sb"].between(
                -EPS,
                1.0 + EPS
            ).all()
        )

        assert (
            df["loss_area_m2"]
            <= MIRROR_AREA + 1e-8
        ).all()

        fig1 = plot_eta_sb_field(
            df,
            month,
            ST
        )

        fig2 = plot_local_projection_example(
            df=df,
            state=state,
            H=H,
            r=r,
            blocking_candidates=
                blocking_candidates,
            month=month,
            ST=ST
        )

        print("\nOutputs:")
        print(output_csv)
        print(fig1)
        print(fig2)

    # --------------------------------------------------------
    # D. FULL 60 工况
    # --------------------------------------------------------

    else:

        all_frames = []
        summary_rows = []

        for month in MONTHS:

            for ST in SOLAR_TIMES:

                t0 = perf_counter()

                print(
                    f"\nRunning "
                    f"{month}/21 {ST:04.1f}"
                )

                df, state = compute_condition(
                    month,
                    ST,
                    H,
                    r,
                    d_hr,
                    D2,
                    blocking_candidates
                )

                all_frames.append(
                    df
                )

                summary = summarize_condition(
                    df
                )

                summary_rows.append({
                    "month": month,
                    "day": 21,
                    "ST": ST,
                    **summary
                })

                elapsed = (
                    perf_counter()
                    -
                    t0
                )

                print(
                    "  eta_sb mean = "
                    f"{summary['eta_sb_mean']:.6f}"
                )

                print(
                    "  eta_sb min  = "
                    f"{summary['eta_sb_min']:.6f}"
                )

                print(
                    "  elapsed = "
                    f"{elapsed:.2f} s"
                )

        all_df = pd.concat(
            all_frames,
            ignore_index=True
        )

        summary_df = pd.DataFrame(
            summary_rows
        )

        all_path = (
            RESULT_DIR
            /
            "q1_sb_all.csv"
        )

        summary_path = (
            RESULT_DIR
            /
            "q1_sb_summary.csv"
        )

        all_df.to_csv(
            all_path,
            index=False,
            encoding="utf-8-sig"
        )

        summary_df.to_csv(
            summary_path,
            index=False,
            encoding="utf-8-sig"
        )

        # 选择冬季低太阳高度工况做展示
        selected = all_df[
            (all_df["month"] == 12)
            &
            (all_df["ST"] == 9.0)
        ]

        plot_eta_sb_field(
            selected,
            12,
            9.0
        )

        print("\nFull results:")
        print(all_path)
        print(summary_path)

    total_time = (
        perf_counter()
        -
        start_time
    )

    print(
        "\nTotal elapsed = "
        f"{total_time:.2f} s"
    )

    print(
        "\nNOTE:"
        "\neta_sb is now calculated."
        "\neta_trunc is NOT calculated yet."
        "\nTherefore this is still not the final Q1 optical efficiency."
    )


# ============================================================
# 25. 程序入口
# ============================================================

if __name__ == "__main__":
    main()