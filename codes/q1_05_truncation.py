# ============================================================
# q1_05_truncation.py
#
# 2023A Q1 截断效率
#
# 主模型：
#   有限太阳视盘
#       ↓
#   确定性角域求积
#       ↓
#   对每个太阳方向求实际反射方向
#       ↓
#   圆柱集热器沿反射方向逆投影至镜面
#       ↓
#   与有效镜面区域求交面积
#       ↓
#   eta_trunc
#
# 额外：
# 若 results/q1_basic_optical_all.csv 存在，
# 自动计算最终光学效率与镜场输出功率。
# ============================================================

from pathlib import Path
from time import perf_counter
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from numpy.polynomial.legendre import leggauss

from shapely.geometry import MultiPoint
from shapely.ops import unary_union

import q1_03_shadow_blocking as sb


# ============================================================
# 0. 输出目录
# ============================================================

RESULT_DIR = Path("results")
FIG_DIR = Path("figs")

RESULT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. 太阳视盘参数
# ============================================================

# 太阳视盘角半径
SUN_HALF_ANGLE = 4.65e-3       # rad

# 默认求积阶数
DEFAULT_N_MU = 3
DEFAULT_N_PSI = 12

# 集热器圆周离散
DEFAULT_N_THETA = 72

# 数值容差
EPS = 1e-12
AREA_EPS = 1e-10


# ============================================================
# 2. 在太阳中心方向周围建立局部正交基
# ============================================================

def build_sun_basis(s0):
    """
    构造：
        e1 ⟂ s0
        e2 ⟂ s0
        e1 ⟂ e2
    """

    ref = np.array([0.0, 0.0, 1.0])

    # 若 s0 太靠近竖直方向，
    # 换一个参考方向避免叉积退化
    if abs(np.dot(s0, ref)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])

    e1 = np.cross(ref, s0)

    e1 /= np.linalg.norm(e1)

    e2 = np.cross(s0, e1)

    e2 /= np.linalg.norm(e2)

    return e1, e2


# ============================================================
# 3. Uniform pillbox 太阳视盘确定性求积
# ============================================================

def solar_disk_quadrature(
    s0,
    n_mu=DEFAULT_N_MU,
    n_psi=DEFAULT_N_PSI,
):
    """
    太阳盘积分：

        dOmega = sin(beta) d beta d psi

    令：
        mu = cos(beta)

    则：
        dOmega = -dmu dpsi

    径向 mu：
        Gauss-Legendre

    周向 psi：
        等距周期求积

    返回：
        sun_dirs : (K,3)
        weights  : (K,)
    """

    e1, e2 = build_sun_basis(s0)

    # --------------------------------------------
    # Gauss-Legendre on [-1, 1]
    # --------------------------------------------

    xi, wi = leggauss(n_mu)

    mu_min = np.cos(SUN_HALF_ANGLE)

    # 映射到 [mu_min, 1]
    mu = (
        0.5 * (1.0 + mu_min)
        +
        0.5 * (1.0 - mu_min) * xi
    )

    w_mu = (
        0.5
        * (1.0 - mu_min)
        * wi
    )

    beta_sin = np.sqrt(
        np.maximum(
            0.0,
            1.0 - mu ** 2
        )
    )

    # --------------------------------------------
    # 周向
    # --------------------------------------------

    psi = (
        2.0 * np.pi
        *
        (
            np.arange(n_psi) + 0.5
        )
        / n_psi
    )

    w_psi = (
        2.0 * np.pi
        / n_psi
    )

    dirs = []
    weights = []

    for q in range(n_mu):

        for angle in psi:

            tangent_dir = (
                np.cos(angle) * e1
                +
                np.sin(angle) * e2
            )

            s = (
                mu[q] * s0
                +
                beta_sin[q]
                * tangent_dir
            )

            s /= np.linalg.norm(s)

            dirs.append(s)

            weights.append(
                w_mu[q] * w_psi
            )

    dirs = np.asarray(dirs)

    weights = np.asarray(weights)

    # 归一化，之后直接做加权平均
    weights /= np.sum(weights)

    return dirs, weights


# ============================================================
# 4. 镜面反射
# ============================================================

def reflect_sun_directions(
    sun_dirs,
    normal
):
    """
    入射传播方向：

        d_in = -s

    反射：

        d_out
        =
        d_in
        -
        2(d_in dot n)n
    """

    d_in = -sun_dirs

    dot_value = (
        d_in
        @ normal
    )

    d_out = (
        d_in
        -
        2.0
        * dot_value[:, None]
        * normal[None, :]
    )

    d_out /= np.linalg.norm(
        d_out,
        axis=1,
        keepdims=True
    )

    return d_out


# ============================================================
# 5. 有限圆柱边界
# ============================================================

def build_receiver_boundary(
    n_theta
):
    """
    取上下两个圆周。

    圆柱是凸集，
    平行投影后的区域可由这些边界点的凸包逼近。
    """

    theta = np.linspace(
        0.0,
        2.0 * np.pi,
        n_theta,
        endpoint=False
    )

    x = (
        sb.RECEIVER_RADIUS
        * np.cos(theta)
    )

    y = (
        sb.RECEIVER_RADIUS
        * np.sin(theta)
    )

    lower = np.column_stack([
        x,
        y,
        np.full(
            n_theta,
            sb.RECEIVER_Z_BOTTOM
        )
    ])

    upper = np.column_stack([
        x,
        y,
        np.full(
            n_theta,
            sb.RECEIVER_Z_TOP
        )
    ])

    return np.vstack([
        lower,
        upper
    ])


# ============================================================
# 6. 射线是否击中圆柱侧壁
# ============================================================

def ray_hits_receiver_side(
    points,
    direction
):
    """
    points:
        (M,3)

    direction:
        (3,)

    判断：
        P + t d
    是否击中：

        x^2 + y^2 = R^2
        z_bottom <= z <= z_top
    """

    points = np.asarray(points)

    dx, dy, dz = direction

    A = dx * dx + dy * dy

    if A < EPS:
        return np.zeros(
            len(points),
            dtype=bool
        )

    px = points[:, 0]
    py = points[:, 1]
    pz = points[:, 2]

    B = 2.0 * (
        px * dx
        +
        py * dy
    )

    C = (
        px ** 2
        +
        py ** 2
        -
        sb.RECEIVER_RADIUS ** 2
    )

    discriminant = (
        B ** 2
        -
        4.0 * A * C
    )

    hit = np.zeros(
        len(points),
        dtype=bool
    )

    valid = (
        discriminant >= 0.0
    )

    if not np.any(valid):
        return hit

    sqrt_disc = np.sqrt(
        np.maximum(
            discriminant[valid],
            0.0
        )
    )

    Bv = B[valid]

    t1 = (
        -Bv - sqrt_disc
    ) / (2.0 * A)

    t2 = (
        -Bv + sqrt_disc
    ) / (2.0 * A)

    pz_v = pz[valid]

    z1 = (
        pz_v
        +
        t1 * dz
    )

    z2 = (
        pz_v
        +
        t2 * dz
    )

    hit1 = (
        (t1 > EPS)
        &
        (z1 >= sb.RECEIVER_Z_BOTTOM - EPS)
        &
        (z1 <= sb.RECEIVER_Z_TOP + EPS)
    )

    hit2 = (
        (t2 > EPS)
        &
        (z2 >= sb.RECEIVER_Z_BOTTOM - EPS)
        &
        (z2 <= sb.RECEIVER_Z_TOP + EPS)
    )

    hit_valid = (
        hit1 | hit2
    )

    hit[
        np.flatnonzero(valid)
    ] = hit_valid

    return hit


# ============================================================
# 7. 实心圆柱投影与侧壁命中的等价条件
# ============================================================

def sidewall_equivalence_ratio(
    direction
):
    """
    若：

        horizontal / |vertical|
        >
        diameter / height
        =
        7/8
        =
        0.875

    则一条穿过实心有限圆柱的射线
    不可能只从一个底面进、另一个底面出，
    必然与圆柱侧壁相交。

    Q1 中理论上应始终满足。
    """

    horizontal = np.hypot(
        direction[0],
        direction[1]
    )

    vertical = abs(
        direction[2]
    )

    if vertical < EPS:
        return np.inf

    return (
        horizontal
        / vertical
    )


# ============================================================
# 8. 集热器逆投影到镜面
# ============================================================

def receiver_inverse_projection(
    receiver_boundary,
    H_i,
    n_i,
    e_w_i,
    e_h_i,
    direction,
):
    """
    光真实方向：
        +direction

    因此把集热器沿：
        -direction

    逆投回镜面。

    Q' =
        Q
        -
        [n·(Q-H)/(n·d)] d
    """

    denom = np.dot(
        n_i,
        direction
    )

    if abs(denom) < EPS:
        return None

    rel = (
        receiver_boundary
        -
        H_i
    )

    lam = (
        rel @ n_i
    ) / denom

    projected = (
        receiver_boundary
        -
        lam[:, None]
        * direction[None, :]
    )

    uv = sb.to_local_uv(
        projected,
        H_i,
        e_w_i,
        e_h_i
    )

    hull = MultiPoint(
        uv.tolist()
    ).convex_hull

    if (
        hull.is_empty
        or
        hull.geom_type != "Polygon"
    ):
        return None

    return hull


# ============================================================
# 9. 单工况的 SB 状态
# ============================================================

def build_condition_state(
    month,
    ST,
    H,
    r,
    D2
):
    """
    不调用 compute_condition()，
    避免先算一遍 eta_sb 后又重新计算 geometry。

    只构造状态。
    """

    D = sb.day_from_spring_equinox(
        month,
        21
    )

    delta = sb.solar_declination(
        D
    )

    omega = sb.solar_hour_angle(
        ST
    )

    s0 = sb.solar_direction_enu(
        delta,
        omega
    )

    normals = sb.heliostat_normals(
        s0,
        r
    )

    e_w, e_h = (
        sb.mirror_local_frames(
            normals
        )
    )

    vertices = sb.mirror_vertices(
        H,
        e_w,
        e_h
    )

    shadow_candidates = (
        sb.shadow_candidate_lists(
            H,
            s0,
            D2
        )
    )

    receiver_candidates = (
        sb.receiver_shadow_candidate_mask(
            H,
            s0
        )
    )

    return {
        "D": D,
        "delta": delta,
        "omega": omega,
        "s0": s0,

        "normals": normals,
        "e_w": e_w,
        "e_h": e_h,
        "vertices": vertices,

        "shadow_candidates":
            shadow_candidates,

        "receiver_candidates":
            receiver_candidates,
    }


# ============================================================
# 10. 单面镜截断效率
# ============================================================

def truncation_for_mirror(
    i,
    H,
    r,
    state,
    blocking_candidates,
    receiver_boundary,
    sun_dirs,
    sun_base_weights,
):
    """
    主结果：

        eta_trunc

    同时计算：

        eta_trunc_fullmirror

    后者忽略 SB 空间分布，
    用于诊断和文献比较。
    """

    # --------------------------------------------------------
    # A. 先重建该镜真正的 SB 无效区域
    # --------------------------------------------------------

    sb_result = (
        sb.compute_eta_sb_for_target(
            i=i,

            H=H,

            normals=
                state["normals"],

            e_w=
                state["e_w"],

            e_h=
                state["e_h"],

            vertices=
                state["vertices"],

            r=r,

            s=
                state["s0"],

            shadow_candidates=
                state[
                    "shadow_candidates"
                ],

            blocking_candidates=
                blocking_candidates,

            receiver_candidate=
                state[
                    "receiver_candidates"
                ],

            return_geometry=True
        )
    )

    loss_geom = (
        sb_result["loss_geom"]
    )

    if (
        loss_geom is None
        or
        loss_geom.is_empty
    ):
        effective_geom = (
            sb.TARGET_RECT
        )
    else:
        effective_geom = (
            sb.TARGET_RECT
            .difference(loss_geom)
        )

    if not effective_geom.is_valid:
        effective_geom = (
            effective_geom.buffer(0)
        )

    effective_area = (
        effective_geom.area
    )

    if effective_area <= AREA_EPS:

        return {
            "eta_sb":
                sb_result["eta_sb"],

            "eta_trunc":
                np.nan,

            "eta_trunc_fullmirror":
                np.nan,

            "quick_full_fraction":
                np.nan,

            "min_sidewall_ratio":
                np.nan,
        }

    # --------------------------------------------------------
    # B. 所有太阳方向经过该镜的真实反射方向
    # --------------------------------------------------------

    n_i = state["normals"][i]

    reflected_dirs = (
        reflect_sun_directions(
            sun_dirs,
            n_i
        )
    )

    # --------------------------------------------------------
    # C. 能量权重
    #
    # uniform sun radiance +
    # 实际镜面入射投影权重
    # --------------------------------------------------------

    incidence_cos = (
        sun_dirs
        @ n_i
    )

    incidence_cos = np.clip(
        incidence_cos,
        0.0,
        None
    )

    weights = (
        sun_base_weights
        *
        incidence_cos
    )

    weights_sum = np.sum(
        weights
    )

    if weights_sum <= EPS:
        raise RuntimeError(
            "Solar quadrature weight vanished."
        )

    weights /= weights_sum

    # --------------------------------------------------------
    # D. 中心反射验证
    # --------------------------------------------------------

    central_reflected = (
        reflect_sun_directions(
            state["s0"][None, :],
            n_i
        )[0]
    )

    center_error = np.linalg.norm(
        central_reflected
        -
        r[i]
    )

    # --------------------------------------------------------
    # E. 对太阳盘所有方向求面积
    # --------------------------------------------------------

    tau_effective = []

    tau_full = []

    quick_full_count = 0

    min_ratio = np.inf

    mirror_corners = (
        state["vertices"][i]
    )

    for k, direction in enumerate(
        reflected_dirs
    ):

        # ----------------------------------------
        # 几何安全检查：
        # 实心圆柱投影是否等价于侧壁命中
        # ----------------------------------------

        ratio = (
            sidewall_equivalence_ratio(
                direction
            )
        )

        min_ratio = min(
            min_ratio,
            ratio
        )

        if ratio <= (
            2.0
            * sb.RECEIVER_RADIUS
            / sb.RECEIVER_HEIGHT
            + 1e-8
        ):
            raise RuntimeError(
                "Sidewall equivalence failed: "
                f"mirror={i+1}, "
                f"ratio={ratio:.6f}"
            )

        # ----------------------------------------
        # 快速通道：
        # 如果目标镜4个角全部命中侧壁，
        # 由于接收区域是凸集，
        # 整块正方形都必然命中。
        # ----------------------------------------

        corner_hits = (
            ray_hits_receiver_side(
                mirror_corners,
                direction
            )
        )

        if np.all(corner_hits):

            tau_effective.append(1.0)
            tau_full.append(1.0)

            quick_full_count += 1

            continue

        # ----------------------------------------
        # 精确逆投影面积
        # ----------------------------------------

        acceptance = (
            receiver_inverse_projection(
                receiver_boundary=
                    receiver_boundary,

                H_i=H[i],
                n_i=n_i,

                e_w_i=
                    state["e_w"][i],

                e_h_i=
                    state["e_h"][i],

                direction=direction,
            )
        )

        if acceptance is None:

            tau_effective.append(0.0)
            tau_full.append(0.0)

            continue

        # 完整镜面
        full_hit = (
            sb.TARGET_RECT
            .intersection(
                acceptance
            )
            .area
        )

        tau_full_k = (
            full_hit
            / sb.MIRROR_AREA
        )

        # 去掉 SB 后的有效镜面
        effective_hit = (
            effective_geom
            .intersection(
                acceptance
            )
            .area
        )

        tau_effective_k = (
            effective_hit
            / effective_area
        )

        tau_full.append(
            np.clip(
                tau_full_k,
                0.0,
                1.0
            )
        )

        tau_effective.append(
            np.clip(
                tau_effective_k,
                0.0,
                1.0
            )
        )

    tau_effective = np.asarray(
        tau_effective
    )

    tau_full = np.asarray(
        tau_full
    )

    eta_trunc = np.sum(
        weights
        * tau_effective
    )

    eta_trunc_full = np.sum(
        weights
        * tau_full
    )

    return {
        "eta_sb":
            sb_result["eta_sb"],

        "eta_trunc":
            eta_trunc,

        "eta_trunc_fullmirror":
            eta_trunc_full,

        "center_reflection_error":
            center_error,

        "quick_full_fraction":
            (
                quick_full_count
                / len(reflected_dirs)
            ),

        "min_sidewall_ratio":
            min_ratio,
    }


# ============================================================
# 11. 单个太阳工况
# ============================================================

def compute_condition(
    month,
    ST,
    H,
    r,
    d_hr,
    D2,
    blocking_candidates,
    receiver_boundary,
    n_mu,
    n_psi,
):
    state = build_condition_state(
        month,
        ST,
        H,
        r,
        D2
    )

    sun_dirs, sun_weights = (
        solar_disk_quadrature(
            state["s0"],
            n_mu=n_mu,
            n_psi=n_psi,
        )
    )

    rows = []

    for i in range(len(H)):

        result = (
            truncation_for_mirror(
                i=i,

                H=H,
                r=r,

                state=state,

                blocking_candidates=
                    blocking_candidates,

                receiver_boundary=
                    receiver_boundary,

                sun_dirs=sun_dirs,

                sun_base_weights=
                    sun_weights,
            )
        )

        rows.append({
            "month": month,
            "day": 21,
            "ST": ST,

            "mirror_id": i + 1,

            "x_m": H[i, 0],
            "y_m": H[i, 1],

            **result
        })

    return pd.DataFrame(
        rows
    )


# ============================================================
# 12. 月平均图
# ============================================================

def plot_monthly_truncation(
    monthly_df
):
    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    ax.plot(
        monthly_df["month"],
        monthly_df["eta_trunc_monthly"],
        marker="o",
        label="Conditional truncation"
    )

    ax.plot(
        monthly_df["month"],
        monthly_df[
            "eta_trunc_fullmirror_monthly"
        ],
        marker="s",
        label="Full-mirror diagnostic"
    )

    ax.set_xlabel(
        "Month"
    )

    ax.set_ylabel(
        "Average truncation efficiency"
    )

    ax.set_title(
        "Monthly Average Truncation Efficiency"
    )

    ax.set_xticks(
        range(1, 13)
    )

    ax.grid(
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()

    output = (
        FIG_DIR
        /
        "q1_result_eta_trunc_monthly.png"
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    return output


# ============================================================
# 13. 典型工况空间分布
# ============================================================

def plot_spatial_truncation(
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
        c=df["eta_trunc"],
        s=18,
        vmin=0.0,
        vmax=1.0
    )

    cbar = fig.colorbar(
        sc,
        ax=ax
    )

    cbar.set_label(
        "Truncation efficiency"
    )

    ax.scatter(
        [0],
        [0],
        marker="*",
        s=160,
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
        "Spatial Distribution of "
        "Truncation Efficiency\n"
        f"{month}/21  {ST:04.1f}"
    )

    ax.legend()

    fig.tight_layout()

    output = (
        FIG_DIR
        /
        (
            f"q1_result_eta_trunc_"
            f"m{month:02d}_"
            f"t{ST:04.1f}.png"
        )
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    return output


# ============================================================
# 14. 如果基础光学结果已存在，自动合并最终 Q1
# ============================================================

def combine_final_optical(
    trunc_all
):
    basic_path = (
        RESULT_DIR
        /
        "q1_basic_optical_all.csv"
    )

    if not basic_path.exists():

        print(
            "\n[提示] 未找到 "
            "results/q1_basic_optical_all.csv"
        )

        print(
            "因此本次只输出 eta_trunc，"
            "暂不合并最终光学效率。"
        )

        return None

    basic = pd.read_csv(
        basic_path
    )

    keep_cols = [
        "month",
        "ST",
        "mirror_id",

        "DNI_kW_m2",

        "eta_cos",
        "eta_at",
        "eta_ref",
    ]

    basic = basic[
        keep_cols
    ]

    merged = basic.merge(
        trunc_all[
            [
                "month",
                "ST",
                "mirror_id",
                "eta_sb",
                "eta_trunc",
            ]
        ],

        on=[
            "month",
            "ST",
            "mirror_id"
        ],

        how="inner"
    )

    merged[
        "eta_optical"
    ] = (
        merged["eta_cos"]
        *
        merged["eta_at"]
        *
        merged["eta_ref"]
        *
        merged["eta_sb"]
        *
        merged["eta_trunc"]
    )

    merged[
        "power_kW"
    ] = (
        merged["DNI_kW_m2"]
        *
        sb.MIRROR_AREA
        *
        merged["eta_optical"]
    )

    all_out = (
        RESULT_DIR
        /
        "q1_final_optical_all.csv"
    )

    merged.to_csv(
        all_out,
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------
    # 60 工况
    # --------------------------------------------

    condition = (
        merged
        .groupby(
            ["month", "ST"],
            as_index=False
        )
        .agg(
            eta_optical_mean=(
                "eta_optical",
                "mean"
            ),

            eta_cos_mean=(
                "eta_cos",
                "mean"
            ),

            eta_sb_mean=(
                "eta_sb",
                "mean"
            ),

            eta_trunc_mean=(
                "eta_trunc",
                "mean"
            ),

            field_power_kW=(
                "power_kW",
                "sum"
            ),
        )
    )

    condition[
        "field_power_MW"
    ] = (
        condition[
            "field_power_kW"
        ]
        / 1000.0
    )

    total_area = (
        merged[
            "mirror_id"
        ].nunique()
        *
        sb.MIRROR_AREA
    )

    condition[
        "unit_area_power_kW_m2"
    ] = (
        condition[
            "field_power_kW"
        ]
        / total_area
    )

    condition_out = (
        RESULT_DIR
        /
        "q1_final_condition_summary.csv"
    )

    condition.to_csv(
        condition_out,
        index=False,
        encoding="utf-8-sig"
    )

    annual_optical = (
        condition[
            "eta_optical_mean"
        ].mean()
    )

    annual_power_MW = (
        condition[
            "field_power_MW"
        ].mean()
    )

    annual_unit_power = (
        condition[
            "unit_area_power_kW_m2"
        ].mean()
    )

    print()
    print("=" * 72)
    print("FINAL Q1 OPTICAL RESULT")
    print("=" * 72)

    print(
        f"Annual mean optical efficiency "
        f"= {annual_optical:.8f}"
    )

    print(
        f"Annual mean field power "
        f"= {annual_power_MW:.6f} MW"
    )

    print(
        f"Annual mean unit-area power "
        f"= {annual_unit_power:.6f} kW/m^2"
    )

    return {
        "annual_optical":
            annual_optical,

        "annual_power_MW":
            annual_power_MW,

        "annual_unit_power":
            annual_unit_power,
    }


# ============================================================
# 15. 主程序
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=[
            "test",
            "full"
        ],
        default="test"
    )

    parser.add_argument(
        "--n-mu",
        type=int,
        default=DEFAULT_N_MU
    )

    parser.add_argument(
        "--n-psi",
        type=int,
        default=DEFAULT_N_PSI
    )

    parser.add_argument(
        "--n-theta",
        type=int,
        default=DEFAULT_N_THETA
    )

    args = parser.parse_args()

    print("=" * 72)
    print("Q1 Truncation Efficiency")
    print("=" * 72)

    print(
        f"Sun half angle = "
        f"{SUN_HALF_ANGLE:.6e} rad"
    )

    print(
        f"Solar quadrature = "
        f"{args.n_mu} x "
        f"{args.n_psi} = "
        f"{args.n_mu * args.n_psi} "
        f"directions"
    )

    print(
        f"Receiver circle discretization "
        f"= {args.n_theta}"
    )

    # --------------------------------------------------------
    # A. 基础镜场
    # --------------------------------------------------------

    H = sb.read_heliostat_coordinates(
        sb.DATA_PATH
    )

    r, d_hr = sb.receiver_direction(
        H
    )

    D2 = sb.pairwise_distance_squared(
        H
    )

    blocking_candidates = (
        sb.blocking_candidate_lists(
            H,
            r,
            d_hr,
            D2
        )
    )

    receiver_boundary = (
        build_receiver_boundary(
            args.n_theta
        )
    )

    # --------------------------------------------------------
    # B. TEST
    # --------------------------------------------------------

    if args.mode == "test":

        month = 12
        ST = 9.0

        t0 = perf_counter()

        df = compute_condition(
            month=month,
            ST=ST,

            H=H,
            r=r,
            d_hr=d_hr,
            D2=D2,

            blocking_candidates=
                blocking_candidates,

            receiver_boundary=
                receiver_boundary,

            n_mu=args.n_mu,
            n_psi=args.n_psi,
        )

        elapsed = (
            perf_counter()
            -
            t0
        )

        out = (
            RESULT_DIR
            /
            "q1_trunc_test.csv"
        )

        df.to_csv(
            out,
            index=False,
            encoding="utf-8-sig"
        )

        print()
        print("=" * 72)
        print("TEST RESULT")
        print("=" * 72)

        print(
            f"eta_trunc mean = "
            f"{df['eta_trunc'].mean():.8f}"
        )

        print(
            f"eta_trunc min  = "
            f"{df['eta_trunc'].min():.8f}"
        )

        print(
            f"eta_trunc max  = "
            f"{df['eta_trunc'].max():.8f}"
        )

        print(
            f"full-mirror diagnostic mean = "
            f"{df['eta_trunc_fullmirror'].mean():.8f}"
        )

        print(
            f"max center-reflection error = "
            f"{df['center_reflection_error'].max():.3e}"
        )

        print(
            f"mean quick-full fraction = "
            f"{df['quick_full_fraction'].mean():.4f}"
        )

        print(
            f"minimum sidewall ratio = "
            f"{df['min_sidewall_ratio'].min():.6f}"
        )

        print(
            f"elapsed = "
            f"{elapsed:.2f} s"
        )

        fig = plot_spatial_truncation(
            df,
            month,
            ST
        )

        print()
        print(out)
        print(fig)

    # --------------------------------------------------------
    # C. FULL 60 CONDITIONS
    # --------------------------------------------------------

    else:

        all_frames = []
        condition_rows = []

        total_start = perf_counter()

        for month in sb.MONTHS:

            for ST in sb.SOLAR_TIMES:

                t0 = perf_counter()

                df = compute_condition(
                    month=month,
                    ST=ST,

                    H=H,
                    r=r,
                    d_hr=d_hr,
                    D2=D2,

                    blocking_candidates=
                        blocking_candidates,

                    receiver_boundary=
                        receiver_boundary,

                    n_mu=args.n_mu,
                    n_psi=args.n_psi,
                )

                all_frames.append(
                    df
                )

                elapsed = (
                    perf_counter()
                    -
                    t0
                )

                condition_rows.append({
                    "month": month,
                    "ST": ST,

                    "eta_sb_mean":
                        df["eta_sb"].mean(),

                    "eta_trunc_mean":
                        df[
                            "eta_trunc"
                        ].mean(),

                    "eta_trunc_fullmirror_mean":
                        df[
                            "eta_trunc_fullmirror"
                        ].mean(),

                    "center_reflection_error_max":
                        df[
                            "center_reflection_error"
                        ].max(),

                    "quick_full_fraction_mean":
                        df[
                            "quick_full_fraction"
                        ].mean(),

                    "min_sidewall_ratio":
                        df[
                            "min_sidewall_ratio"
                        ].min(),

                    "elapsed_s":
                        elapsed,
                })

                print(
                    f"{month:02d}/21 "
                    f"{ST:04.1f}   "
                    f"eta_trunc="
                    f"{df['eta_trunc'].mean():.6f}   "
                    f"time="
                    f"{elapsed:.2f}s"
                )

        all_df = pd.concat(
            all_frames,
            ignore_index=True
        )

        condition_df = pd.DataFrame(
            condition_rows
        )

        # --------------------------------------------
        # 月平均
        # --------------------------------------------

        monthly = (
            condition_df
            .groupby(
                "month",
                as_index=False
            )
            .agg(
                eta_sb_monthly=(
                    "eta_sb_mean",
                    "mean"
                ),

                eta_trunc_monthly=(
                    "eta_trunc_mean",
                    "mean"
                ),

                eta_trunc_fullmirror_monthly=(
                    "eta_trunc_fullmirror_mean",
                    "mean"
                ),
            )
        )

        annual_eta_trunc = (
            monthly[
                "eta_trunc_monthly"
            ].mean()
        )

        annual_eta_trunc_full = (
            monthly[
                "eta_trunc_fullmirror_monthly"
            ].mean()
        )

        annual_eta_sb = (
            monthly[
                "eta_sb_monthly"
            ].mean()
        )

        total_elapsed = (
            perf_counter()
            -
            total_start
        )

        # --------------------------------------------
        # 保存
        # --------------------------------------------

        all_path = (
            RESULT_DIR
            /
            "q1_trunc_all.csv"
        )

        condition_path = (
            RESULT_DIR
            /
            "q1_trunc_condition_summary.csv"
        )

        monthly_path = (
            RESULT_DIR
            /
            "q1_trunc_monthly_average.csv"
        )

        all_df.to_csv(
            all_path,
            index=False,
            encoding="utf-8-sig"
        )

        condition_df.to_csv(
            condition_path,
            index=False,
            encoding="utf-8-sig"
        )

        monthly.to_csv(
            monthly_path,
            index=False,
            encoding="utf-8-sig"
        )

        # --------------------------------------------
        # 输出
        # --------------------------------------------

        print()
        print("=" * 72)
        print("MONTHLY TRUNCATION EFFICIENCY")
        print("=" * 72)

        print(
            monthly.to_string(
                index=False,
                float_format=lambda x:
                    f"{x:.6f}"
            )
        )

        print()
        print("=" * 72)
        print("ANNUAL RESULT")
        print("=" * 72)

        print(
            f"Annual eta_sb "
            f"= {annual_eta_sb:.8f}"
        )

        print(
            f"Annual eta_trunc "
            f"= {annual_eta_trunc:.8f}"
        )

        print(
            f"Annual eta_trunc "
            f"(full mirror diagnostic) "
            f"= {annual_eta_trunc_full:.8f}"
        )

        print(
            f"Maximum center reflection error "
            f"= "
            f"{condition_df['center_reflection_error_max'].max():.3e}"
        )

        print(
            f"Minimum sidewall ratio "
            f"= "
            f"{condition_df['min_sidewall_ratio'].min():.6f}"
        )

        print(
            f"Total elapsed "
            f"= {total_elapsed:.2f} s"
        )

        # --------------------------------------------
        # 图
        # --------------------------------------------

        monthly_fig = (
            plot_monthly_truncation(
                monthly
            )
        )

        selected = all_df[
            (all_df["month"] == 6)
            &
            (all_df["ST"] == 12.0)
        ]

        spatial_fig = (
            plot_spatial_truncation(
                selected,
                6,
                12.0
            )
        )

        print()
        print("Outputs:")
        print(all_path)
        print(condition_path)
        print(monthly_path)
        print(monthly_fig)
        print(spatial_fig)

        # --------------------------------------------
        # 最终 Q1 光学效率
        # --------------------------------------------

        combine_final_optical(
            all_df
        )


if __name__ == "__main__":
    main()