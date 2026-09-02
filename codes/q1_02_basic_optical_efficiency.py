from pathlib import Path
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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

# 镜场纬度：北纬 39.4°
PHI = np.deg2rad(39.4)

# 场址海拔，DNI 公式要求单位 km
H_SITE = 3.0

# 太阳常数 kW/m^2
G0 = 1.366

# Q1 镜面安装高度 m
HELIOSTAT_HEIGHT = 4.0

# Q1 镜面尺寸 6 m × 6 m
MIRROR_WIDTH = 6.0
MIRROR_HEIGHT = 6.0
MIRROR_AREA = MIRROR_WIDTH * MIRROR_HEIGHT

# Q1 集热器中心：
# 吸收塔位于镜场中心，集热器中心高度为 80 m
RECEIVER = np.array([0.0, 0.0, 80.0])

# 镜面反射率
ETA_REF = 0.92

# 题目规定的五个当地时刻
SOLAR_TIMES = [9.0, 10.5, 12.0, 13.5, 15.0]

# 每月 21 日
MONTHS = list(range(1, 13))


# ============================================================
# 2. 工具函数：读取定日镜坐标
# ============================================================

def read_heliostat_coordinates(path):
    """
    读取附件中的定日镜 x, y 坐标。

    优先寻找列名中含 x / y 的数值列；
    若识别失败，则退化为前两列数值列。
    """

    xls = pd.ExcelFile(path)
    sheet_name = xls.sheet_names[0]

    df = pd.read_excel(path, sheet_name=sheet_name)

    # 去掉全空行、全空列
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")

    print("=" * 70)
    print("附件读取信息")
    print("=" * 70)
    print("工作表：", xls.sheet_names)
    print("当前读取：", sheet_name)
    print("数据形状：", df.shape)
    print("列名：", df.columns.tolist())

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    # 优先按列名寻找 x / y
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
        print(
            "\n[提示] 未能从列名明确识别 x、y，"
            "暂时使用前两个数值列。"
        )
        x_col = numeric_cols[0]
        y_col = numeric_cols[1]

    else:
        raise ValueError(
            "无法找到两列数值型坐标。"
            "请把 Excel 的列名和前几行截图发给我。"
        )

    print("识别 x 列：", x_col)
    print("识别 y 列：", y_col)

    x = df[x_col].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)

    # 去掉可能残留的 NaN
    mask = np.isfinite(x) & np.isfinite(y)

    x = x[mask]
    y = y[mask]

    return x, y


# ============================================================
# 3. 日期 D：以 3 月 21 日为第 0 天
# ============================================================

def day_from_spring_equinox(month, day=21):
    """
    使用一个非闰年 2023 作为日期计数载体。

    3 月 21 日：
        D = 0

    1、2 月会得到负 D，但赤纬公式是周期函数，
    因此与使用 365 周期内对应的正 D 等价。
    """

    reference = date(2023, 3, 21)
    current = date(2023, month, day)

    return (current - reference).days


# ============================================================
# 4. 太阳赤纬角 δ
# ============================================================

def solar_declination(D):
    """
    题目公式：

    sin(delta)
        = sin(2*pi*D/365) * sin(23.45°)

    返回：
        delta，单位 rad
    """

    sin_delta = (
        np.sin(2.0 * np.pi * D / 365.0)
        * np.sin(np.deg2rad(23.45))
    )

    # 防止浮点误差造成 arcsin 输入略超 [-1,1]
    sin_delta = np.clip(sin_delta, -1.0, 1.0)

    return np.arcsin(sin_delta)


# ============================================================
# 5. 太阳时角 ω
# ============================================================

def solar_hour_angle(ST):
    """
    题目公式：

    omega = pi/12 * (ST - 12)

    ST 单位：小时
    omega 单位：rad
    """

    return np.pi / 12.0 * (ST - 12.0)


# ============================================================
# 6. 太阳 ENU 单位方向
# ============================================================

def solar_direction_enu(delta, omega, phi=PHI):
    """
    ENU 坐标：
        x = East
        y = North
        z = Up

    s 的方向：
        镜场 -> 太阳中心
    """

    s_e = -np.cos(delta) * np.sin(omega)

    s_n = (
        np.cos(phi) * np.sin(delta)
        - np.sin(phi) * np.cos(delta) * np.cos(omega)
    )

    s_u = (
        np.sin(phi) * np.sin(delta)
        + np.cos(phi) * np.cos(delta) * np.cos(omega)
    )

    s = np.array([s_e, s_n, s_u], dtype=float)

    # 理论上已经是单位向量；
    # 这里再归一化只用于抵消微小浮点误差
    return s / np.linalg.norm(s)


# ============================================================
# 7. 题目太阳高度角公式
# ============================================================

def sin_solar_altitude_formula(delta, omega, phi=PHI):
    """
    用题目附录的标量公式独立计算 sin(alpha_s)，
    用来与太阳向量的 s_U 交叉验证。
    """

    return (
        np.cos(delta)
        * np.cos(phi)
        * np.cos(omega)
        +
        np.sin(delta)
        * np.sin(phi)
    )


# ============================================================
# 8. DNI
# ============================================================

def dni_from_sun_up(s_u, h_site=H_SITE):
    """
    题目给出的 DNI 经验公式。

    h_site 单位：km
    DNI 单位：kW/m^2
    """

    if s_u <= 0:
        # 本题规定工况太阳均应在地平线上方。
        # 留这个判断是为了让程序更健壮。
        return 0.0

    a = 0.4237 - 0.00821 * (6.0 - h_site) ** 2

    b = 0.5055 + 0.00595 * (6.5 - h_site) ** 2

    c = 0.2711 + 0.01858 * (2.5 - h_site) ** 2

    dni = G0 * (
        a
        +
        b * np.exp(-c / s_u)
    )

    return dni


# ============================================================
# 9. 镜子到集热器方向 r_i 和距离 d_HR
# ============================================================

def receiver_direction_and_distance(H):
    """
    H:
        shape = (..., 3)

    返回：
        r:
            镜面中心 -> 集热器中心的单位向量

        d:
            镜面中心 -> 集热器中心的距离
    """

    vec = RECEIVER - H

    d = np.linalg.norm(vec, axis=-1)

    r = vec / d[..., None]

    return r, d


# ============================================================
# 10. 大气透射率
# ============================================================

def atmospheric_efficiency(d_hr):
    """
    题目经验公式，d_hr 单位 m。
    """

    return (
        0.99321
        - 0.0001176 * d_hr
        + 1.97e-8 * d_hr ** 2
    )


# ============================================================
# 11. 镜面法向 n_i
# ============================================================

def heliostat_normal(s, r):
    """
    s:
        shape = (3,)

    r:
        shape = (N, 3)

    返回：
        n:
        shape = (N, 3)
    """

    temp = s[None, :] + r

    norm = np.linalg.norm(
        temp,
        axis=1,
        keepdims=True
    )

    return temp / norm


# ============================================================
# 12. 余弦效率
# ============================================================

def cosine_efficiency_method_1(s, n):
    """
    eta_cos = s dot n
    """

    return np.sum(
        n * s[None, :],
        axis=1
    )


def cosine_efficiency_method_2(s, r):
    """
    eta_cos
        = sqrt((1 + s dot r)/2)

    作为独立交叉验证公式。
    """

    sr = np.sum(
        r * s[None, :],
        axis=1
    )

    inside = (1.0 + sr) / 2.0

    # 消除极小浮点误差
    inside = np.clip(inside, 0.0, 1.0)

    return np.sqrt(inside)


# ============================================================
# 13. 反射定律回代
# ============================================================

def reflected_direction(s, n):
    """
    实际入射传播方向：
        d_in = -s

    反射公式：
        d_ref
        = d_in
          - 2(d_in dot n)n
    """

    d_in = -s

    dot_value = np.sum(
        n * d_in[None, :],
        axis=1
    )

    d_ref = (
        d_in[None, :]
        -
        2.0
        * dot_value[:, None]
        * n
    )

    return d_ref


# ============================================================
# 14. 主程序
# ============================================================

def main():

    # --------------------------------------------------------
    # A. 读取定日镜坐标
    # --------------------------------------------------------

    x, y = read_heliostat_coordinates(DATA_PATH)

    N = len(x)

    z = np.full(
        N,
        HELIOSTAT_HEIGHT
    )

    H = np.column_stack(
        [x, y, z]
    )

    print("\n定日镜总数 N =", N)


    # --------------------------------------------------------
    # B. Q1 中不随时间变化的空间量
    # --------------------------------------------------------

    r, d_hr = receiver_direction_and_distance(H)

    eta_at = atmospheric_efficiency(d_hr)

    # 检查题目经验公式适用范围
    if np.max(d_hr) > 1000.0:
        print(
            "\n[警告] 存在 d_HR > 1000 m 的镜子，"
            "超出题目给出的 eta_at 公式范围。"
        )

    print("\nd_HR 范围 / m:")
    print(
        f"{d_hr.min():.6f}"
        f" ~ "
        f"{d_hr.max():.6f}"
    )

    print("\neta_at 范围:")
    print(
        f"{eta_at.min():.8f}"
        f" ~ "
        f"{eta_at.max():.8f}"
    )


    # --------------------------------------------------------
    # C. 遍历题目 60 个工况
    # --------------------------------------------------------

    all_rows = []
    summary_rows = []

    validation = {
        "sun_unit_max_error": 0.0,
        "receiver_dir_unit_max_error": 0.0,
        "normal_unit_max_error": 0.0,
        "reflection_max_error": 0.0,
        "cosine_formula_max_error": 0.0,
        "solar_altitude_formula_max_error": 0.0,
    }

    # r 是静态量，可先验证
    r_norm = np.linalg.norm(r, axis=1)

    validation[
        "receiver_dir_unit_max_error"
    ] = np.max(
        np.abs(r_norm - 1.0)
    )


    for month in MONTHS:

        D = day_from_spring_equinox(
            month,
            21
        )

        delta = solar_declination(D)

        for ST in SOLAR_TIMES:

            omega = solar_hour_angle(ST)

            # --------------------
            # 太阳方向
            # --------------------

            s = solar_direction_enu(
                delta,
                omega
            )

            sun_norm_error = abs(
                np.linalg.norm(s) - 1.0
            )

            validation[
                "sun_unit_max_error"
            ] = max(
                validation[
                    "sun_unit_max_error"
                ],
                sun_norm_error
            )

            # --------------------
            # 与题目高度角公式验证
            # --------------------

            sin_alpha_formula = (
                sin_solar_altitude_formula(
                    delta,
                    omega
                )
            )

            altitude_error = abs(
                s[2] - sin_alpha_formula
            )

            validation[
                "solar_altitude_formula_max_error"
            ] = max(
                validation[
                    "solar_altitude_formula_max_error"
                ],
                altitude_error
            )

            # --------------------
            # DNI
            # --------------------

            dni = dni_from_sun_up(
                s[2]
            )

            # --------------------
            # 镜面法向
            # --------------------

            n = heliostat_normal(
                s,
                r
            )

            n_norm = np.linalg.norm(
                n,
                axis=1
            )

            normal_error = np.max(
                np.abs(n_norm - 1.0)
            )

            validation[
                "normal_unit_max_error"
            ] = max(
                validation[
                    "normal_unit_max_error"
                ],
                normal_error
            )

            # --------------------
            # 余弦效率：两种公式
            # --------------------

            eta_cos_1 = (
                cosine_efficiency_method_1(
                    s,
                    n
                )
            )

            eta_cos_2 = (
                cosine_efficiency_method_2(
                    s,
                    r
                )
            )

            cosine_error = np.max(
                np.abs(
                    eta_cos_1
                    -
                    eta_cos_2
                )
            )

            validation[
                "cosine_formula_max_error"
            ] = max(
                validation[
                    "cosine_formula_max_error"
                ],
                cosine_error
            )

            # --------------------
            # 反射定律验证
            # --------------------

            d_ref = reflected_direction(
                s,
                n
            )

            reflection_error = np.max(
                np.linalg.norm(
                    d_ref - r,
                    axis=1
                )
            )

            validation[
                "reflection_max_error"
            ] = max(
                validation[
                    "reflection_max_error"
                ],
                reflection_error
            )

            # --------------------
            # 基础诊断功率
            # 注意：不是最终功率
            # --------------------

            p_basic = (
                dni
                * MIRROR_AREA
                * eta_cos_1
                * eta_at
                * ETA_REF
            )

            # --------------------
            # 保存逐镜结果
            # --------------------

            for i in range(N):

                all_rows.append({
                    "month": month,
                    "day": 21,
                    "D": D,
                    "ST": ST,

                    "delta_deg": np.rad2deg(delta),
                    "omega_deg": np.rad2deg(omega),

                    "s_E": s[0],
                    "s_N": s[1],
                    "s_U": s[2],

                    "DNI_kW_m2": dni,

                    "mirror_id": i + 1,

                    "x_m": x[i],
                    "y_m": y[i],
                    "z_m": z[i],

                    "d_HR_m": d_hr[i],

                    "r_E": r[i, 0],
                    "r_N": r[i, 1],
                    "r_U": r[i, 2],

                    "n_E": n[i, 0],
                    "n_N": n[i, 1],
                    "n_U": n[i, 2],

                    "eta_cos": eta_cos_1[i],
                    "eta_at": eta_at[i],
                    "eta_ref": ETA_REF,

                    # 尚未计入：
                    # eta_sb
                    # eta_trunc
                    "P_basic_kW": p_basic[i],
                })

            # --------------------
            # 保存当前工况摘要
            # --------------------

            summary_rows.append({
                "month": month,
                "day": 21,
                "ST": ST,

                "D": D,
                "delta_deg": np.rad2deg(delta),
                "omega_deg": np.rad2deg(omega),

                "solar_E": s[0],
                "solar_N": s[1],
                "solar_U": s[2],

                "solar_altitude_deg":
                    np.rad2deg(
                        np.arcsin(s[2])
                    ),

                "DNI_kW_m2": dni,

                "eta_cos_mean":
                    np.mean(eta_cos_1),

                "eta_cos_min":
                    np.min(eta_cos_1),

                "eta_cos_max":
                    np.max(eta_cos_1),

                "eta_at_mean":
                    np.mean(eta_at),

                "P_basic_field_MW":
                    np.sum(p_basic) / 1000.0,
            })


    # --------------------------------------------------------
    # D. 转成 DataFrame
    # --------------------------------------------------------

    all_df = pd.DataFrame(
        all_rows
    )

    summary_df = pd.DataFrame(
        summary_rows
    )


    # --------------------------------------------------------
    # E. 上午—下午太阳几何对称性验证
    # --------------------------------------------------------

    symmetry_errors = []

    for month in MONTHS:

        month_df = summary_df[
            summary_df["month"] == month
        ]

        for t1, t2 in [
            (9.0, 15.0),
            (10.5, 13.5),
        ]:

            row1 = month_df[
                month_df["ST"] == t1
            ].iloc[0]

            row2 = month_df[
                month_df["ST"] == t2
            ].iloc[0]

            # 理论：
            # E 分量反号
            # N、U、DNI 相同

            err_e = abs(
                row1["solar_E"]
                +
                row2["solar_E"]
            )

            err_n = abs(
                row1["solar_N"]
                -
                row2["solar_N"]
            )

            err_u = abs(
                row1["solar_U"]
                -
                row2["solar_U"]
            )

            err_dni = abs(
                row1["DNI_kW_m2"]
                -
                row2["DNI_kW_m2"]
            )

            symmetry_errors.extend([
                err_e,
                err_n,
                err_u,
                err_dni,
            ])

    validation[
        "solar_symmetry_max_error"
    ] = max(symmetry_errors)


    # --------------------------------------------------------
    # F. 输出 CSV
    # --------------------------------------------------------

    all_path = (
        RESULT_DIR
        / "q1_basic_optical_all.csv"
    )

    summary_path = (
        RESULT_DIR
        / "q1_basic_optical_summary.csv"
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


    # --------------------------------------------------------
    # G. 输出验证日志
    # --------------------------------------------------------

    validation_path = (
        RESULT_DIR
        / "q1_basic_validation.txt"
    )

    with open(
        validation_path,
        "w",
        encoding="utf-8"
    ) as f:

        f.write(
            "2023A Q1 基础光学模块验证\n"
        )

        f.write("=" * 60 + "\n\n")

        f.write(
            f"heliostat_count = {N}\n\n"
        )

        for key, value in validation.items():

            f.write(
                f"{key} = "
                f"{value:.16e}\n"
            )

        f.write("\n")

        f.write(
            "说明：\n"
            "1. 这些是基础几何模块验证指标。\n"
            "2. 当前尚未计算 eta_sb 与 eta_trunc。\n"
            "3. P_basic 不是最终镜场输出功率。\n"
        )


    # ========================================================
    # H. 图 1：太阳几何验证
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    for month in MONTHS:

        month_df = summary_df[
            summary_df["month"] == month
        ]

        ax.plot(
            month_df["ST"],
            month_df["solar_altitude_deg"],
            marker="o",
            label=f"{month}月21日"
        )

    ax.set_xlabel(
        "Local time / h"
    )

    ax.set_ylabel(
        "Solar altitude / degree"
    )

    ax.set_title(
        "Solar Altitude at the 60 Prescribed Conditions"
    )

    ax.grid(
        alpha=0.3
    )

    ax.legend(
        ncol=3,
        fontsize=8
    )

    fig.tight_layout()

    solar_fig_path = (
        FIG_DIR
        / "q1_verify_solar_geometry.png"
    )

    fig.savefig(
        solar_fig_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


    # ========================================================
    # I. 图 2：典型时刻余弦效率空间分布
    # ========================================================

    # 选 6 月 21 日 12:00
    selected = all_df[
        (all_df["month"] == 6)
        &
        (all_df["ST"] == 12.0)
    ]

    fig, ax = plt.subplots(
        figsize=(8, 8)
    )

    sc = ax.scatter(
        selected["x_m"],
        selected["y_m"],
        c=selected["eta_cos"],
        s=16
    )

    cbar = fig.colorbar(
        sc,
        ax=ax
    )

    cbar.set_label(
        "Cosine efficiency"
    )

    ax.scatter(
        [0],
        [0],
        marker="*",
        s=150,
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
        "Spatial Distribution of Cosine Efficiency\n"
        "June 21, 12:00"
    )

    ax.legend()

    fig.tight_layout()

    cosine_fig_path = (
        FIG_DIR
        / "q1_mechanism_cosine_efficiency_field.png"
    )

    fig.savefig(
        cosine_fig_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


    # ========================================================
    # J. 终端摘要
    # ========================================================

    print("\n" + "=" * 70)
    print("基础光学模块计算完成")
    print("=" * 70)

    print("\n验证指标：")

    for key, value in validation.items():

        print(
            f"{key:40s}"
            f"{value:.6e}"
        )

    print("\n余弦效率总体范围：")

    print(
        f"min = "
        f"{all_df['eta_cos'].min():.6f}"
    )

    print(
        f"max = "
        f"{all_df['eta_cos'].max():.6f}"
    )

    print(
        f"mean = "
        f"{all_df['eta_cos'].mean():.6f}"
    )

    print("\n输出文件：")
    print(all_path)
    print(summary_path)
    print(validation_path)
    print(solar_fig_path)
    print(cosine_fig_path)

    print(
        "\n注意：当前 P_basic 尚未包含"
        "阴影遮挡效率 eta_sb 和"
        "截断效率 eta_trunc，"
        "不得作为问题一最终输出功率。"
    )


# ============================================================
# 15. 程序入口
# ============================================================

if __name__ == "__main__":
    main()