# ============================================================
# q1_05a_check_sun_cone.py
#
# 目的：
# 检查 q1_05_truncation.py 中有限太阳视盘是否真正发挥作用
#
# 检查内容：
#
# 1. 太阳盘实际采样角 beta 是否非零
# 2. 反射后光锥角是否被保留
# 3. Point Sun 与 Finite Sun 的 eta_trunc 比较
# 4. quick-full 快速通道 vs 强制精确逆投影比较
# 5. 检查太阳锥角影响是否随镜塔距离增大
#
# 不修改原 q1_05_truncation.py
#
# 运行：
# python codes/q1_05a_check_sun_cone.py
#
# 或：
# python codes/q1_05a_check_sun_cone.py --month 12 --st 9
# ============================================================

from pathlib import Path
from time import perf_counter
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import q1_03_shadow_blocking as sb
import q1_05_truncation as tr


# ============================================================
# 0. 输出目录
# ============================================================

RESULT_DIR = Path("results")
FIG_DIR = Path("figs")

RESULT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. 基础工具
# ============================================================

def angle_between_vectors(a, b):
    """
    两个单位向量夹角，返回 rad。
    """

    value = np.dot(a, b)

    value = np.clip(
        value,
        -1.0,
        1.0
    )

    return np.arccos(value)


def angular_deviation(vectors, center):
    """
    多个单位向量相对中心方向的夹角。
    """

    values = vectors @ center

    values = np.clip(
        values,
        -1.0,
        1.0
    )

    return np.arccos(values)


# ============================================================
# 2. 重建目标镜的有效镜面 E_i
# ============================================================

def build_effective_geometry(
    i,
    H,
    r,
    state,
    blocking_candidates
):
    """
    E_i = R_i - L_i

    L_i：
        receiver shadow
        + heliostat shadow
        + blocking
    """

    sb_result = (
        sb.compute_eta_sb_for_target(
            i=i,

            H=H,

            normals=state["normals"],
            e_w=state["e_w"],
            e_h=state["e_h"],
            vertices=state["vertices"],

            r=r,
            s=state["s0"],

            shadow_candidates=
                state["shadow_candidates"],

            blocking_candidates=
                blocking_candidates,

            receiver_candidate=
                state["receiver_candidates"],

            return_geometry=True
        )
    )

    loss_geom = sb_result[
        "loss_geom"
    ]

    if (
        loss_geom is None
        or loss_geom.is_empty
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

    return (
        effective_geom,
        sb_result
    )


# ============================================================
# 3. 强制精确计算 eta_trunc
#
# 这里故意完全关闭 quick-full
# ============================================================

def exact_truncation_given_sun_dirs(
    i,
    H,
    r,
    state,
    receiver_boundary,
    effective_geom,
    sun_dirs,
    sun_base_weights
):
    """
    与主程序思想一致：

        每一个太阳方向
            ↓
        得到实际反射方向
            ↓
        集热器逆投影
            ↓
        与有效镜面求交

    关键区别：

        这里不允许 quick-full shortcut。

    所有太阳方向都真正执行：
        receiver_inverse_projection()
    """

    n_i = state["normals"][i]

    # --------------------------------------------------------
    # A. 实际反射方向
    # --------------------------------------------------------

    reflected_dirs = (
        tr.reflect_sun_directions(
            sun_dirs,
            n_i
        )
    )

    # --------------------------------------------------------
    # B. 太阳方向能量权重
    # --------------------------------------------------------

    incidence_cos = (
        sun_dirs @ n_i
    )

    incidence_cos = np.clip(
        incidence_cos,
        0.0,
        None
    )

    weights = (
        sun_base_weights
        * incidence_cos
    )

    weight_sum = np.sum(
        weights
    )

    if weight_sum <= tr.EPS:
        raise RuntimeError(
            "Solar weight vanished."
        )

    weights /= weight_sum

    # --------------------------------------------------------
    # C. 有效镜面面积
    # --------------------------------------------------------

    effective_area = (
        effective_geom.area
    )

    if effective_area <= tr.AREA_EPS:
        return np.nan, np.nan

    tau_effective = []
    tau_full = []

    # --------------------------------------------------------
    # D. 每个太阳方向精确逆投影
    # --------------------------------------------------------

    for direction in reflected_dirs:

        acceptance = (
            tr.receiver_inverse_projection(
                receiver_boundary=
                    receiver_boundary,

                H_i=H[i],
                n_i=n_i,

                e_w_i=
                    state["e_w"][i],

                e_h_i=
                    state["e_h"][i],

                direction=direction
            )
        )

        if acceptance is None:

            tau_effective.append(0.0)
            tau_full.append(0.0)

            continue

        # ------------------------------
        # 整块镜面诊断值
        # ------------------------------

        full_hit_area = (
            sb.TARGET_RECT
            .intersection(
                acceptance
            )
            .area
        )

        tau_full_k = (
            full_hit_area
            /
            sb.MIRROR_AREA
        )

        # ------------------------------
        # 真正有效镜面
        # ------------------------------

        effective_hit_area = (
            effective_geom
            .intersection(
                acceptance
            )
            .area
        )

        tau_effective_k = (
            effective_hit_area
            /
            effective_area
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

    eta_effective = np.sum(
        weights
        * tau_effective
    )

    eta_full = np.sum(
        weights
        * tau_full
    )

    return (
        eta_effective,
        eta_full
    )


# ============================================================
# 4. 单面镜完整诊断
# ============================================================

def diagnose_one_mirror(
    i,
    H,
    r,
    d_hr,
    state,
    blocking_candidates,
    receiver_boundary,
    finite_sun_dirs,
    finite_weights
):
    """
    比较：

    Point Sun
    Finite Sun exact
    Finite Sun quick
    """

    # --------------------------------------------------------
    # A. 有效镜面
    # --------------------------------------------------------

    effective_geom, sb_result = (
        build_effective_geometry(
            i=i,
            H=H,
            r=r,
            state=state,
            blocking_candidates=
                blocking_candidates
        )
    )

    # --------------------------------------------------------
    # B. Point Sun
    # --------------------------------------------------------

    point_dirs = (
        state["s0"][None, :]
    )

    point_weights = np.array(
        [1.0]
    )

    (
        eta_point,
        eta_point_full
    ) = exact_truncation_given_sun_dirs(
        i=i,

        H=H,
        r=r,

        state=state,

        receiver_boundary=
            receiver_boundary,

        effective_geom=
            effective_geom,

        sun_dirs=
            point_dirs,

        sun_base_weights=
            point_weights
    )

    # --------------------------------------------------------
    # C. Finite Sun
    #    强制精确逆投影
    # --------------------------------------------------------

    (
        eta_finite_exact,
        eta_finite_full
    ) = exact_truncation_given_sun_dirs(
        i=i,

        H=H,
        r=r,

        state=state,

        receiver_boundary=
            receiver_boundary,

        effective_geom=
            effective_geom,

        sun_dirs=
            finite_sun_dirs,

        sun_base_weights=
            finite_weights
    )

    # --------------------------------------------------------
    # D. Finite Sun
    #    原 q1_05 中的 quick-full 版本
    # --------------------------------------------------------

    quick_result = (
        tr.truncation_for_mirror(
            i=i,

            H=H,
            r=r,

            state=state,

            blocking_candidates=
                blocking_candidates,

            receiver_boundary=
                receiver_boundary,

            sun_dirs=
                finite_sun_dirs,

            sun_base_weights=
                finite_weights
        )
    )

    eta_finite_quick = (
        quick_result[
            "eta_trunc"
        ]
    )

    # --------------------------------------------------------
    # E. 反射锥角
    # --------------------------------------------------------

    reflected_dirs = (
        tr.reflect_sun_directions(
            finite_sun_dirs,
            state["normals"][i]
        )
    )

    reflected_angles = (
        angular_deviation(
            reflected_dirs,
            r[i]
        )
    )

    # --------------------------------------------------------
    # F. 太阳锥角空间展宽估计
    # --------------------------------------------------------

    max_reflected_angle = (
        np.max(
            reflected_angles
        )
    )

    # 单侧横向展宽
    spread_one_side = (
        d_hr[i]
        *
        np.tan(
            max_reflected_angle
        )
    )

    spread_full_width = (
        2.0
        *
        spread_one_side
    )

    return {
        "mirror_id":
            i + 1,

        "d_HR_m":
            d_hr[i],

        "eta_sb":
            sb_result["eta_sb"],

        "eta_trunc_point":
            eta_point,

        "eta_trunc_finite_exact":
            eta_finite_exact,

        "eta_trunc_finite_quick":
            eta_finite_quick,

        "finite_minus_point":
            (
                eta_finite_exact
                -
                eta_point
            ),

        "quick_minus_exact":
            (
                eta_finite_quick
                -
                eta_finite_exact
            ),

        "eta_point_fullmirror":
            eta_point_full,

        "eta_finite_fullmirror":
            eta_finite_full,

        "reflected_angle_min_mrad":
            (
                np.min(
                    reflected_angles
                )
                * 1000.0
            ),

        "reflected_angle_max_mrad":
            (
                np.max(
                    reflected_angles
                )
                * 1000.0
            ),

        "sun_spread_one_side_m":
            spread_one_side,

        "sun_spread_full_width_m":
            spread_full_width,

        "quick_full_fraction":
            quick_result[
                "quick_full_fraction"
            ],
    }


# ============================================================
# 5. 代表性镜面选择
# ============================================================

def representative_mirrors(
    d_hr
):
    """
    按镜塔距离选择：

        最近
        25%
        中位
        75%
        最远
    """

    order = np.argsort(
        d_hr
    )

    fractions = [
        0.00,
        0.25,
        0.50,
        0.75,
        1.00,
    ]

    indices = []

    for f in fractions:

        pos = int(
            round(
                f
                *
                (len(order) - 1)
            )
        )

        indices.append(
            order[pos]
        )

    return np.array(
        indices,
        dtype=int
    )


# ============================================================
# 6. 距离趋势采样
# ============================================================

def distance_sample_mirrors(
    d_hr,
    n_samples
):
    """
    不随机抽样。

    沿镜塔距离由近到远均匀选取镜面。
    """

    order = np.argsort(
        d_hr
    )

    positions = np.linspace(
        0,
        len(order) - 1,
        n_samples
    )

    positions = np.round(
        positions
    ).astype(int)

    indices = order[
        positions
    ]

    return np.unique(
        indices
    )


# ============================================================
# 7. 绘制锥角影响 vs 镜塔距离
# ============================================================

def plot_cone_effect(
    df
):
    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    ax.scatter(
        df["d_HR_m"],
        df[
            "finite_minus_point"
        ]
    )

    ax.axhline(
        0.0,
        linewidth=1
    )

    ax.set_xlabel(
        "Heliostat-receiver distance / m"
    )

    ax.set_ylabel(
        r"$\eta_{\mathrm{trunc}}^{finite}"
        r"-"
        r"\eta_{\mathrm{trunc}}^{point}$"
    )

    ax.set_title(
        "Finite-Sun Effect on "
        "Truncation Efficiency"
    )

    ax.grid(
        alpha=0.3
    )

    fig.tight_layout()

    output = (
        FIG_DIR
        /
        "q1_verify_sun_cone_effect_vs_distance.png"
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    return output


# ============================================================
# 8. 主程序
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--month",
        type=int,
        default=12
    )

    parser.add_argument(
        "--st",
        type=float,
        default=9.0
    )

    parser.add_argument(
        "--n-mu",
        type=int,
        default=3
    )

    parser.add_argument(
        "--n-psi",
        type=int,
        default=12
    )

    parser.add_argument(
        "--n-theta",
        type=int,
        default=72
    )

    parser.add_argument(
        "--samples",
        type=int,
        default=40
    )

    args = parser.parse_args()

    start_time = perf_counter()

    print("=" * 76)
    print(
        "Q1 Finite Sun / Solar Cone Diagnostic"
    )
    print("=" * 76)

    print(
        f"Condition = "
        f"{args.month}/21 "
        f"{args.st}"
    )

    print(
        f"Sun half angle = "
        f"{tr.SUN_HALF_ANGLE * 1000:.6f} mrad"
    )

    # ========================================================
    # A. 镜场基础量
    # ========================================================

    H = sb.read_heliostat_coordinates(
        sb.DATA_PATH
    )

    r, d_hr = (
        sb.receiver_direction(
            H
        )
    )

    D2 = (
        sb.pairwise_distance_squared(
            H
        )
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
        tr.build_receiver_boundary(
            args.n_theta
        )
    )

    state = (
        tr.build_condition_state(
            month=args.month,
            ST=args.st,
            H=H,
            r=r,
            D2=D2
        )
    )

    # ========================================================
    # B. 构造有限太阳视盘
    # ========================================================

    finite_sun_dirs, finite_weights = (
        tr.solar_disk_quadrature(
            state["s0"],
            n_mu=args.n_mu,
            n_psi=args.n_psi
        )
    )

    sun_angles = (
        angular_deviation(
            finite_sun_dirs,
            state["s0"]
        )
    )

    print()
    print("=" * 76)
    print("1. SOLAR DISK ANGULAR CHECK")
    print("=" * 76)

    print(
        f"Number of solar directions = "
        f"{len(finite_sun_dirs)}"
    )

    print(
        "Solar angular range = "
        f"{sun_angles.min() * 1000:.6f}"
        " ~ "
        f"{sun_angles.max() * 1000:.6f}"
        " mrad"
    )

    print(
        "Configured solar half angle = "
        f"{tr.SUN_HALF_ANGLE * 1000:.6f}"
        " mrad"
    )

    # --------------------------------------------------------
    # 基本断言
    # --------------------------------------------------------

    if (
        np.max(
            sun_angles
        )
        < 1e-6
    ):
        raise RuntimeError(
            "ERROR: Solar disk collapsed "
            "to a point."
        )

    print(
        "PASS: finite solar angular "
        "distribution exists."
    )

    # ========================================================
    # C. 检查反射是否保留角度
    # ========================================================

    rep_indices = (
        representative_mirrors(
            d_hr
        )
    )

    test_i = (
        rep_indices[-1]
    )

    reflected_dirs = (
        tr.reflect_sun_directions(
            finite_sun_dirs,
            state["normals"][test_i]
        )
    )

    reflected_angles = (
        angular_deviation(
            reflected_dirs,
            r[test_i]
        )
    )

    angular_error = (
        np.max(
            np.abs(
                reflected_angles
                -
                sun_angles
            )
        )
    )

    print()
    print("=" * 76)
    print("2. REFLECTED CONE ANGULAR CHECK")
    print("=" * 76)

    print(
        f"Test mirror = "
        f"{test_i + 1}"
    )

    print(
        f"d_HR = "
        f"{d_hr[test_i]:.6f} m"
    )

    print(
        "Incident cone range = "
        f"{sun_angles.min() * 1000:.6f}"
        " ~ "
        f"{sun_angles.max() * 1000:.6f}"
        " mrad"
    )

    print(
        "Reflected cone range = "
        f"{reflected_angles.min() * 1000:.6f}"
        " ~ "
        f"{reflected_angles.max() * 1000:.6f}"
        " mrad"
    )

    print(
        "Max angular-preservation error = "
        f"{angular_error:.3e} rad"
    )

    if angular_error < 1e-9:

        print(
            "PASS: mirror reflection "
            "preserves the solar cone angle."
        )

    else:

        print(
            "WARNING: angular preservation "
            "error is unexpectedly large."
        )

    # ========================================================
    # D. 五面代表镜：
    # Point / Finite / Quick 比较
    # ========================================================

    print()
    print("=" * 76)
    print(
        "3. POINT SUN vs FINITE SUN"
    )
    print("=" * 76)

    representative_rows = []

    for i in rep_indices:

        row = diagnose_one_mirror(
            i=i,

            H=H,
            r=r,
            d_hr=d_hr,

            state=state,

            blocking_candidates=
                blocking_candidates,

            receiver_boundary=
                receiver_boundary,

            finite_sun_dirs=
                finite_sun_dirs,

            finite_weights=
                finite_weights
        )

        representative_rows.append(
            row
        )

        print()
        print(
            f"Mirror {row['mirror_id']}"
            f"  d={row['d_HR_m']:.3f} m"
        )

        print(
            "  Point Sun eta_trunc       = "
            f"{row['eta_trunc_point']:.8f}"
        )

        print(
            "  Finite Sun exact eta      = "
            f"{row['eta_trunc_finite_exact']:.8f}"
        )

        print(
            "  Finite - Point            = "
            f"{row['finite_minus_point']:+.8f}"
        )

        print(
            "  Finite Sun quick eta      = "
            f"{row['eta_trunc_finite_quick']:.8f}"
        )

        print(
            "  Quick - Exact             = "
            f"{row['quick_minus_exact']:+.8e}"
        )

        print(
            "  quick-full fraction       = "
            f"{row['quick_full_fraction']:.4f}"
        )

        print(
            "  reflected cone max        = "
            f"{row['reflected_angle_max_mrad']:.4f}"
            " mrad"
        )

        print(
            "  estimated full sun spread = "
            f"{row['sun_spread_full_width_m']:.4f}"
            " m"
        )

    representative_df = (
        pd.DataFrame(
            representative_rows
        )
    )

    rep_output = (
        RESULT_DIR
        /
        "q1_sun_cone_representative.csv"
    )

    representative_df.to_csv(
        rep_output,
        index=False,
        encoding="utf-8-sig"
    )

    # ========================================================
    # E. 距离趋势
    # ========================================================

    print()
    print("=" * 76)
    print(
        "4. SUN-CONE EFFECT vs DISTANCE"
    )
    print("=" * 76)

    sample_indices = (
        distance_sample_mirrors(
            d_hr,
            args.samples
        )
    )

    trend_rows = []

    for count, i in enumerate(
        sample_indices,
        start=1
    ):

        (
            effective_geom,
            sb_result
        ) = build_effective_geometry(
            i=i,

            H=H,
            r=r,

            state=state,

            blocking_candidates=
                blocking_candidates
        )

        # ----------------------------------------
        # Point Sun
        # ----------------------------------------

        eta_point, _ = (
            exact_truncation_given_sun_dirs(
                i=i,

                H=H,
                r=r,

                state=state,

                receiver_boundary=
                    receiver_boundary,

                effective_geom=
                    effective_geom,

                sun_dirs=
                    state["s0"][None, :],

                sun_base_weights=
                    np.array([1.0])
            )
        )

        # ----------------------------------------
        # Finite Sun
        # ----------------------------------------

        eta_finite, _ = (
            exact_truncation_given_sun_dirs(
                i=i,

                H=H,
                r=r,

                state=state,

                receiver_boundary=
                    receiver_boundary,

                effective_geom=
                    effective_geom,

                sun_dirs=
                    finite_sun_dirs,

                sun_base_weights=
                    finite_weights
            )
        )

        cone_effect = (
            eta_finite
            -
            eta_point
        )

        estimated_spread = (
            2.0
            *
            d_hr[i]
            *
            np.tan(
                np.max(
                    sun_angles
                )
            )
        )

        trend_rows.append({
            "mirror_id":
                i + 1,

            "d_HR_m":
                d_hr[i],

            "eta_sb":
                sb_result["eta_sb"],

            "eta_point":
                eta_point,

            "eta_finite":
                eta_finite,

            "finite_minus_point":
                cone_effect,

            "estimated_full_spread_m":
                estimated_spread,
        })

        print(
            f"[{count:02d}/"
            f"{len(sample_indices):02d}] "
            f"id={i+1:4d}  "
            f"d={d_hr[i]:7.2f} m  "
            f"point={eta_point:.6f}  "
            f"finite={eta_finite:.6f}  "
            f"delta={cone_effect:+.6f}"
        )

    trend_df = pd.DataFrame(
        trend_rows
    )

    trend_output = (
        RESULT_DIR
        /
        "q1_sun_cone_distance_check.csv"
    )

    trend_df.to_csv(
        trend_output,
        index=False,
        encoding="utf-8-sig"
    )

    fig_output = (
        plot_cone_effect(
            trend_df
        )
    )

    # ========================================================
    # F. 汇总诊断
    # ========================================================

    max_abs_cone_effect = (
        trend_df[
            "finite_minus_point"
        ]
        .abs()
        .max()
    )

    mean_abs_cone_effect = (
        trend_df[
            "finite_minus_point"
        ]
        .abs()
        .mean()
    )

    max_quick_error = (
        representative_df[
            "quick_minus_exact"
        ]
        .abs()
        .max()
    )

    print()
    print("=" * 76)
    print("5. DIAGNOSTIC SUMMARY")
    print("=" * 76)

    print(
        "Maximum |Finite - Point| = "
        f"{max_abs_cone_effect:.8f}"
    )

    print(
        "Mean |Finite - Point| = "
        f"{mean_abs_cone_effect:.8f}"
    )

    print(
        "Maximum |Quick - Exact| = "
        f"{max_quick_error:.8e}"
    )

    if max_abs_cone_effect < 1e-4:

        print()
        print(
            "WARNING:"
        )

        print(
            "Finite Sun and Point Sun "
            "are almost identical."
        )

        print(
            "The solar cone exists in the "
            "direction calculation, but its "
            "effect may be disappearing in "
            "the receiver-acceptance geometry."
        )

    else:

        print()
        print(
            "PASS:"
        )

        print(
            "Finite solar angle produces "
            "a measurable truncation effect."
        )

    if max_quick_error > 1e-6:

        print()
        print(
            "WARNING:"
        )

        print(
            "The quick-full shortcut changes "
            "the result noticeably."
        )

        print(
            "Recommend disabling the shortcut "
            "in the final truncation solver."
        )

    else:

        print()
        print(
            "PASS:"
        )

        print(
            "quick-full shortcut agrees with "
            "the full inverse-projection result "
            "for the representative mirrors."
        )

    elapsed = (
        perf_counter()
        -
        start_time
    )

    print()
    print("=" * 76)
    print("OUTPUTS")
    print("=" * 76)

    print(rep_output)
    print(trend_output)
    print(fig_output)

    print(
        f"\nTotal elapsed = "
        f"{elapsed:.2f} s"
    )


# ============================================================
# 9. 程序入口
# ============================================================

if __name__ == "__main__":
    main()
    