# ============================================================
# q1_04_validate_and_average_sb.py
#
# 功能：
# 1. 验证候选镜筛选是否漏掉真实 Shadowing / Blocking
# 2. 验证四顶点平行投影面积是否符合解析面积公式
# 3. 计算 60 个题设工况下的月平均、年平均 eta_sb
# 4. 与公开资料结果自动比较
#
# 运行：
# python codes/q1_04_validate_and_average_sb.py
# ============================================================

from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from shapely.geometry import Polygon

import q1_03_shadow_blocking as sb


# ============================================================
# 0. 参数
# ============================================================

RESULT_DIR = Path("results")
RESULT_DIR.mkdir(parents=True, exist_ok=True)

# 每个验证工况暴力验证多少面目标镜
# 12 已经足够有代表性；
# 如果想更严格，可以改成 20 或 30
N_VALIDATE_TARGETS = 12

# 固定随机种子，保证结果可复现
RNG_SEED = 2026

# 两个验证工况：
# 冬季早晨 = Shadowing 较严重
# 夏季正午 = 高太阳高度的另一极端
VALIDATION_CONDITIONS = [
    (12, 9.0),
    (6, 12.0),
]

# 网上公开结果
WEB_REFERENCE_ANNUAL = {
    # 电子科技大学《实验科学与技术》
    # 基础镜场 1745 面、62820 m^2
    "UESTC_journal_MonteCarlo": 0.9738,

    # CN118428061A 中基础镜场计算
    "CN118428061A": 0.9240,
}

# 专利给出的基础镜场月平均 eta_sb
WEB_REFERENCE_MONTHLY_PATENT = {
    1: 0.9163,
    2: 0.9211,
    3: 0.9321,
    4: 0.9388,
    5: 0.9412,
    6: 0.9492,
    7: 0.9323,
    8: 0.9223,
    9: 0.9201,
    10: 0.9135,
    11: 0.9102,
    12: 0.8903,
}


# ============================================================
# 1. 构造“暴力全邻镜”候选列表
# ============================================================

def single_target_candidate_lists(N, i, neighbors):
    """
    compute_eta_sb_for_target() 只会读取第 i 项，
    所以只为目标镜 i 构造全邻镜列表。
    """

    empty = np.array([], dtype=int)

    lists = [
        empty
        for _ in range(N)
    ]

    lists[i] = neighbors

    return lists


# ============================================================
# 2. geometry 差异面积
# ============================================================

def symmetric_difference_area(g1, g2):

    if g1 is None and g2 is None:
        return 0.0

    if g1 is None:
        return 0.0 if g2.is_empty else g2.area

    if g2 is None:
        return 0.0 if g1.is_empty else g1.area

    return g1.symmetric_difference(g2).area


# ============================================================
# 3. 单个目标镜：筛选法 vs 暴力全搜索
# ============================================================

def validate_target_full_search(
    i,
    H,
    r,
    state,
    blocking_candidates,
):
    """
    当前快速模型：
        只检查 candidate_filter 得到的少量镜子

    暴力基准：
        对其余 N-1 面镜全部做投影判断

    如果二者 eta_sb、loss geometry 一致，
    说明候选筛选没有漏掉真正产生损失的邻镜。
    """

    N = len(H)

    # --------------------------------------------------------
    # A. 当前筛选模型
    # --------------------------------------------------------

    filtered = sb.compute_eta_sb_for_target(
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

        return_geometry=True,
    )

    # --------------------------------------------------------
    # B. 暴力全搜索
    # --------------------------------------------------------

    all_neighbors = np.delete(
        np.arange(N),
        i
    )

    full_shadow_lists = (
        single_target_candidate_lists(
            N,
            i,
            all_neighbors
        )
    )

    full_block_lists = (
        single_target_candidate_lists(
            N,
            i,
            all_neighbors
        )
    )

    # Receiver shadow 也强制计算，
    # 顺便验证 receiver candidate mask 是否漏检
    receiver_all = np.ones(
        N,
        dtype=bool
    )

    full = sb.compute_eta_sb_for_target(
        i=i,

        H=H,

        normals=state["normals"],
        e_w=state["e_w"],
        e_h=state["e_h"],
        vertices=state["vertices"],

        r=r,
        s=state["s"],

        shadow_candidates=
            full_shadow_lists,

        blocking_candidates=
            full_block_lists,

        receiver_candidate=
            receiver_all,

        return_geometry=True,
    )

    return {
        "mirror_id": i + 1,

        "eta_filter":
            filtered["eta_sb"],

        "eta_full":
            full["eta_sb"],

        "eta_abs_error":
            abs(
                filtered["eta_sb"]
                -
                full["eta_sb"]
            ),

        "loss_area_filter":
            filtered["loss_area_m2"],

        "loss_area_full":
            full["loss_area_m2"],

        "loss_area_abs_error":
            abs(
                filtered["loss_area_m2"]
                -
                full["loss_area_m2"]
            ),

        "loss_geometry_symdiff_area":
            symmetric_difference_area(
                filtered["loss_geom"],
                full["loss_geom"]
            ),

        "shadow_symdiff_area":
            symmetric_difference_area(
                filtered["shadow_geom"],
                full["shadow_geom"]
            ),

        "blocking_symdiff_area":
            symmetric_difference_area(
                filtered["blocking_geom"],
                full["blocking_geom"]
            ),

        "receiver_symdiff_area":
            symmetric_difference_area(
                filtered["receiver_geom"],
                full["receiver_geom"]
            ),

        "shadow_candidates_fast":
            len(
                state["shadow_candidates"][i]
            ),

        "blocking_candidates_fast":
            len(
                blocking_candidates[i]
            ),

        "full_neighbors":
            N - 1,
    }


# ============================================================
# 4. 平行投影面积解析验证
# ============================================================

def projection_area_validation(
    source_vertices,
    n_source,
    H_target,
    n_target,
    e_w_target,
    e_h_target,
    d,
):
    """
    平行投影面积解析式：

        A_proj
        =
        A_source *
        |n_source · d|
        / |n_target · d|

    与四顶点实际投影后的 polygon 面积比较。
    """

    denom = np.dot(
        n_target,
        d
    )

    if abs(denom) < 1e-10:
        return None

    numer = (
        (source_vertices - H_target)
        @ n_target
    )

    lam = (
        -numer / denom
    )

    projected = (
        source_vertices
        +
        lam[:, None] * d
    )

    uv = sb.to_local_uv(
        projected,
        H_target,
        e_w_target,
        e_h_target
    )

    poly = Polygon(
        uv
    )

    polygon_area = poly.area

    analytic_area = (
        sb.MIRROR_AREA
        *
        abs(
            np.dot(
                n_source,
                d
            )
        )
        /
        abs(
            np.dot(
                n_target,
                d
            )
        )
    )

    abs_error = abs(
        polygon_area
        -
        analytic_area
    )

    rel_error = (
        abs_error / analytic_area
        if analytic_area > 1e-14
        else np.nan
    )

    return {
        "polygon_area":
            polygon_area,

        "analytic_area":
            analytic_area,

        "abs_error":
            abs_error,

        "relative_error":
            rel_error,
    }


# ============================================================
# 5. 执行候选筛选验证
# ============================================================

def run_candidate_validation(
    H,
    r,
    d_hr,
    D2,
    blocking_candidates,
):
    rng = np.random.default_rng(
        RNG_SEED
    )

    validation_rows = []

    for month, ST in VALIDATION_CONDITIONS:

        print()
        print("=" * 72)
        print(
            f"Candidate validation: "
            f"{month}/21 {ST}"
        )
        print("=" * 72)

        # 先用快速模型算一次整个工况，
        # 找出最差镜面
        df, state = sb.compute_condition(
            month,
            ST,
            H,
            r,
            d_hr,
            D2,
            blocking_candidates
        )

        worst_i = int(
            df.loc[
                df["eta_sb"].idxmin(),
                "mirror_id"
            ]
            - 1
        )

        # 随机选其余目标镜
        pool = np.setdiff1d(
            np.arange(len(H)),
            [worst_i]
        )

        random_targets = rng.choice(
            pool,
            size=N_VALIDATE_TARGETS - 1,
            replace=False
        )

        targets = np.concatenate([
            [worst_i],
            random_targets
        ])

        print(
            "Selected target mirrors:",
            (targets + 1).tolist()
        )

        for count, i in enumerate(
            targets,
            start=1
        ):

            t0 = perf_counter()

            row = validate_target_full_search(
                i=i,
                H=H,
                r=r,
                state=state,
                blocking_candidates=
                    blocking_candidates
            )

            row["month"] = month
            row["ST"] = ST

            validation_rows.append(
                row
            )

            dt = perf_counter() - t0

            print(
                f"[{count:02d}/"
                f"{len(targets):02d}] "
                f"mirror={i+1:4d}  "
                f"eta_fast="
                f"{row['eta_filter']:.10f}  "
                f"eta_full="
                f"{row['eta_full']:.10f}  "
                f"error="
                f"{row['eta_abs_error']:.3e}  "
                f"time={dt:.2f}s"
            )

    validation_df = pd.DataFrame(
        validation_rows
    )

    out = (
        RESULT_DIR
        /
        "q1_sb_candidate_validation.csv"
    )

    validation_df.to_csv(
        out,
        index=False,
        encoding="utf-8-sig"
    )

    print()
    print("=" * 72)
    print("Candidate validation summary")
    print("=" * 72)

    print(
        "max eta error =",
        validation_df[
            "eta_abs_error"
        ].max()
    )

    print(
        "max loss-area error =",
        validation_df[
            "loss_area_abs_error"
        ].max()
    )

    print(
        "max geometry symdiff area =",
        validation_df[
            "loss_geometry_symdiff_area"
        ].max()
    )

    print(
        "mean fast shadow candidates =",
        validation_df[
            "shadow_candidates_fast"
        ].mean()
    )

    print(
        "mean fast blocking candidates =",
        validation_df[
            "blocking_candidates_fast"
        ].mean()
    )

    return validation_df


# ============================================================
# 6. 投影面积公式验证
# ============================================================

def run_projection_validation(
    H,
    r,
    d_hr,
    D2,
    blocking_candidates,
):
    """
    使用 12/21 9:00：
    分别检查 Shadowing 和 Blocking 投影面积。
    """

    month = 12
    ST = 9.0

    _, state = sb.compute_condition(
        month,
        ST,
        H,
        r,
        d_hr,
        D2,
        blocking_candidates
    )

    rows = []

    # ----------------------------------------
    # Shadowing 候选
    # ----------------------------------------

    for i in range(len(H)):

        for j in state[
            "shadow_candidates"
        ][i]:

            result = (
                projection_area_validation(
                    source_vertices=
                        state["vertices"][j],

                    n_source=
                        state["normals"][j],

                    H_target=H[i],

                    n_target=
                        state["normals"][i],

                    e_w_target=
                        state["e_w"][i],

                    e_h_target=
                        state["e_h"][i],

                    d=-state["s"],
                )
            )

            if result is not None:

                rows.append({
                    "type": "shadow",
                    "target_id": i + 1,
                    "source_id": j + 1,
                    **result
                })

    # ----------------------------------------
    # Blocking 候选
    # ----------------------------------------

    for i in range(len(H)):

        for j in blocking_candidates[i]:

            result = (
                projection_area_validation(
                    source_vertices=
                        state["vertices"][j],

                    n_source=
                        state["normals"][j],

                    H_target=H[i],

                    n_target=
                        state["normals"][i],

                    e_w_target=
                        state["e_w"][i],

                    e_h_target=
                        state["e_h"][i],

                    d=-r[i],
                )
            )

            if result is not None:

                rows.append({
                    "type": "blocking",
                    "target_id": i + 1,
                    "source_id": j + 1,
                    **result
                })

    df = pd.DataFrame(
        rows
    )

    out = (
        RESULT_DIR
        /
        "q1_sb_projection_area_validation.csv"
    )

    df.to_csv(
        out,
        index=False,
        encoding="utf-8-sig"
    )

    print()
    print("=" * 72)
    print("Projection-area validation")
    print("=" * 72)

    print(
        "cases =",
        len(df)
    )

    print(
        "max abs error =",
        df["abs_error"].max()
    )

    print(
        "max relative error =",
        df["relative_error"].max()
    )

    print(
        "mean relative error =",
        df["relative_error"].mean()
    )

    return df


# ============================================================
# 7. 计算完整 60 工况
# ============================================================

def run_full_60_conditions(
    H,
    r,
    d_hr,
    D2,
    blocking_candidates,
):

    condition_rows = []

    print()
    print("=" * 72)
    print("Running all 60 conditions")
    print("=" * 72)

    total_start = perf_counter()

    for month in sb.MONTHS:

        for ST in sb.SOLAR_TIMES:

            t0 = perf_counter()

            df, _ = sb.compute_condition(
                month,
                ST,
                H,
                r,
                d_hr,
                D2,
                blocking_candidates
            )

            eta_mean = (
                df["eta_sb"].mean()
            )

            elapsed = (
                perf_counter() - t0
            )

            condition_rows.append({
                "month": month,
                "ST": ST,
                "eta_sb_mean": eta_mean,
                "elapsed_s": elapsed,
            })

            print(
                f"{month:02d}/21 "
                f"{ST:04.1f}  "
                f"eta_sb="
                f"{eta_mean:.8f}  "
                f"time="
                f"{elapsed:.2f}s"
            )

    condition_df = pd.DataFrame(
        condition_rows
    )

    total_elapsed = (
        perf_counter() - total_start
    )

    # --------------------------------------------------------
    # 月平均：
    # 每个月 5 个时刻简单平均
    # --------------------------------------------------------

    monthly = (
        condition_df
        .groupby("month", as_index=False)
        ["eta_sb_mean"]
        .mean()
        .rename(
            columns={
                "eta_sb_mean":
                    "eta_sb_monthly"
            }
        )
    )

    # --------------------------------------------------------
    # 年平均：
    # 12个月月平均再平均
    # 等价于60个工况直接平均
    # --------------------------------------------------------

    annual_eta = (
        monthly[
            "eta_sb_monthly"
        ].mean()
    )

    # 加入网上专利结果
    monthly[
        "web_patent_eta_sb"
    ] = monthly["month"].map(
        WEB_REFERENCE_MONTHLY_PATENT
    )

    monthly[
        "difference_vs_patent"
    ] = (
        monthly["eta_sb_monthly"]
        -
        monthly["web_patent_eta_sb"]
    )

    monthly[
        "difference_percentage_point"
    ] = (
        100.0
        *
        monthly[
            "difference_vs_patent"
        ]
    )

    condition_out = (
        RESULT_DIR
        /
        "q1_sb_60_conditions.csv"
    )

    monthly_out = (
        RESULT_DIR
        /
        "q1_sb_monthly_average.csv"
    )

    condition_df.to_csv(
        condition_out,
        index=False,
        encoding="utf-8-sig"
    )

    monthly.to_csv(
        monthly_out,
        index=False,
        encoding="utf-8-sig"
    )

    print()
    print("=" * 72)
    print("Monthly average eta_sb")
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
    print("Annual average eta_sb")
    print("=" * 72)

    print(
        f"Our model = "
        f"{annual_eta:.8f}"
    )

    for name, reference in (
        WEB_REFERENCE_ANNUAL.items()
    ):

        diff = (
            annual_eta - reference
        )

        print(
            f"{name:30s}"
            f" reference={reference:.6f}"
            f"  diff={diff:+.6f}"
            f"  "
            f"({100*diff:+.3f} percentage points)"
        )

    print()
    print(
        f"Pure total compute time = "
        f"{total_elapsed:.2f}s"
    )

    return (
        condition_df,
        monthly,
        annual_eta
    )


# ============================================================
# 8. 主程序
# ============================================================

def main():

    print("=" * 72)
    print(
        "2023A Q1 - eta_sb "
        "Validation and Annual Average"
    )
    print("=" * 72)

    # --------------------------------------------------------
    # A. 基础数据
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

    # --------------------------------------------------------
    # B. 验证1：
    # 候选筛选 vs 暴力全搜索
    # --------------------------------------------------------

    run_candidate_validation(
        H,
        r,
        d_hr,
        D2,
        blocking_candidates
    )

    # --------------------------------------------------------
    # C. 验证2：
    # 四顶点面积 vs 解析面积
    # --------------------------------------------------------

    run_projection_validation(
        H,
        r,
        d_hr,
        D2,
        blocking_candidates
    )

    # --------------------------------------------------------
    # D. 完整 60 工况
    # --------------------------------------------------------

    _, monthly, annual_eta = (
        run_full_60_conditions(
            H,
            r,
            d_hr,
            D2,
            blocking_candidates
        )
    )

    # --------------------------------------------------------
    # E. 最终结果写文本
    # --------------------------------------------------------

    report_path = (
        RESULT_DIR
        /
        "q1_sb_validation_report.txt"
    )

    with open(
        report_path,
        "w",
        encoding="utf-8"
    ) as f:

        f.write(
            "2023A Q1 eta_sb validation\n"
        )

        f.write("=" * 60 + "\n\n")

        f.write(
            f"Our annual eta_sb = "
            f"{annual_eta:.10f}\n\n"
        )

        for name, reference in (
            WEB_REFERENCE_ANNUAL.items()
        ):

            diff = annual_eta - reference

            f.write(
                f"{name}: "
                f"{reference:.6f}, "
                f"diff={diff:+.6f}\n"
            )

        f.write("\n")

        f.write(
            monthly.to_string(
                index=False
            )
        )

    print()
    print(
        "Report saved to:",
        report_path
    )


if __name__ == "__main__":
    main()