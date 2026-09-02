# 23A 模型推导记录

## Q1-01 太阳方向与定日镜法向量模型

### 1. 建模目标

问题一后续需要计算余弦效率、阴影遮挡效率及截断效率，因此首先必须确定：

1. 给定日期和时刻下，太阳在镜场坐标系中的单位方向；
2. 每一面定日镜指向集热器的目标反射方向；
3. 满足反射定律的定日镜镜面法向量。

最终希望建立计算链：

$$
(D,ST)
\rightarrow
(\delta,H)
\rightarrow
\mathbf s_h
\rightarrow
\mathbf s
\rightarrow
\mathbf r_i
\rightarrow
\mathbf n_i.
$$

其中：

- $D$：以春分日为第 0 天的日期序号；
- $ST$：题目规定的当地时间；
- $\delta$：太阳赤纬角；
- $H$：太阳时角；
- $\mathbf s_h$：太阳在时角赤道辅助坐标系中的单位方向；
- $\mathbf s$：太阳在镜场 ENU 坐标系中的单位方向；
- $\mathbf r_i$：第 $i$ 面定日镜指向集热器中心的单位向量；
- $\mathbf n_i$：第 $i$ 面定日镜的单位法向量。

---

## 2. 为什么不直接使用完整天球坐标系？

完整天文学中，可以首先利用太阳赤经 $\alpha$、赤纬 $\delta$ 描述太阳在固定赤道天球坐标系中的位置，再利用当地恒星时和经度计算当地时角：

$$
H=\mathrm{LST}-\alpha.
$$

但本题已经直接给出了太阳赤纬和太阳时角的简化计算方法，因此无需重新计算赤经、恒星时、岁差、章动等高精度天文量。

本题直接采用：

$$
\sin\delta
=
\sin\frac{2\pi D}{365}
\sin23.45^\circ,
$$

以及

$$
H
=
\omega
=
\frac{\pi}{12}(ST-12).
$$

因此模型从 $(H,\delta)$ 开始即可。

这种处理相当于保留与定日镜光学几何直接有关的太阳方向信息，而省略对本题结果贡献很小的完整天文坐标转换过程。

---

## 3. 建立时角赤道辅助坐标系

建立三维正交辅助坐标系

$$
O-X_hY_hZ_h,
$$

其原点取地心。

### 3.1 $Z_h$ 轴

规定

$$
Z_h
$$

沿地球自转轴指向北天极。

### 3.2 $X_h$ 轴

当地子午面是包含地心、观测点、当地天顶及北天极的平面。

取

$$
X_h
$$

为当地子午面与天球赤道面的交线方向。

因此 $X_h$ 与 $Z_h$ 均位于当地子午面内。

### 3.3 $Y_h$ 轴

按照右手坐标系定义

$$
\mathbf e_{Y_h}
=
\mathbf e_{Z_h}\times\mathbf e_{X_h}.
$$

由于 $Y_h$ 垂直于当地子午面，而当地东西方向也垂直于当地子午面，因此选择正方向后有

$$
\boxed{Y_h=E}
$$

即辅助坐标系的 $Y_h$ 轴已经与当地东方向一致。

因此：

$$
X_h-Z_h
$$

平面就是当地子午面，而

$$
Y_h
$$

即当地东方。

---

## 4. 太阳在辅助坐标系中的方向

太阳可视为位于单位天球上，因此只需描述其方向，不需要太阳到地球的实际距离。

太阳赤纬为 $\delta$，因此太阳方向在 $Z_h$ 轴上的分量为

$$
s_{Z_h}=\sin\delta.
$$

其在天球赤道面内的投影长度为

$$
\cos\delta.
$$

太阳时角 $H$ 表示太阳时圈相对于当地子午面的角距离，并规定向西为正。

因此赤道面内：

$$
s_{X_h}=\cos\delta\cos H.
$$

由于 $Y_h$ 规定向东为正，而时角向西为正，因此

$$
s_{Y_h}=-\cos\delta\sin H.
$$

于是太阳在辅助坐标系中的单位方向为

$$
\boxed{
\mathbf s_h=
\begin{pmatrix}
\cos\delta\cos H\\
-\cos\delta\sin H\\
\sin\delta
\end{pmatrix}
}
$$

并满足

$$
\|\mathbf s_h\|=1.
$$

---

## 5. 从辅助坐标系转换到当地 ENU 镜场坐标系

题目规定镜场坐标系为：

$$
x=E,\qquad y=N,\qquad z=U,
$$

其中：

- $E$：正东；
- $N$：正北；
- $U$：当地天顶方向。

所在地纬度记为

$$
\varphi.
$$

本题：

$$
\varphi=39.4^\circ.
$$

---

## 6. 当地天顶单位向量的确定

当地天顶 $U$ 位于当地子午面内，因此在时角辅助坐标系中：

$$
H_U=0.
$$

注意：

$$
H_U=0
$$

表示“天顶方向位于当地子午面”，并不表示任意时刻太阳的时角为零。

另一方面，位于纬度 $\varphi$ 的观测者，其当地天顶在赤道天球坐标中的赤纬恰好等于所在地纬度，即：

$$
\delta_U=\varphi.
$$

因此：

$$
(H_U,\delta_U)=(0,\varphi).
$$

代入辅助坐标系中的球面方向公式，得到当地天顶单位向量：

$$
\boxed{
\mathbf e_U=
\begin{pmatrix}
\cos\varphi\\
0\\
\sin\varphi
\end{pmatrix}
}
$$

其几何意义为：

$$
\angle(Z_h,U)=90^\circ-\varphi.
$$

即从北天极方向旋转一个余纬角即可得到当地天顶方向。

---

## 7. 当地东方和北方单位向量

辅助坐标系的 $Y_h$ 已经与东方重合，因此：

$$
\boxed{
\mathbf e_E=
\begin{pmatrix}
0\\
1\\
0
\end{pmatrix}
}
$$

当地北方向可以理解为沿地球表面纬度增加的方向。

对天顶向量关于纬度求导：

$$
\frac{\partial\mathbf e_U}{\partial\varphi}
=
\begin{pmatrix}
-\sin\varphi\\
0\\
\cos\varphi
\end{pmatrix}.
$$

该向量长度为 1，因此当地北方向单位向量为

$$
\boxed{
\mathbf e_N=
\begin{pmatrix}
-\sin\varphi\\
0\\
\cos\varphi
\end{pmatrix}
}
$$

于是：

$$
\mathbf e_E,\qquad
\mathbf e_N,\qquad
\mathbf e_U
$$

构成题目镜场坐标系的三个正交单位基向量。

---

## 8. 太阳在镜场 ENU 坐标系中的方向

将太阳方向 $\mathbf s_h$ 分别投影到三个当地基向量：

$$
s_E=\mathbf e_E\cdot\mathbf s_h,
$$

$$
s_N=\mathbf e_N\cdot\mathbf s_h,
$$

$$
s_U=\mathbf e_U\cdot\mathbf s_h.
$$

得到：

$$
\boxed{
s_E=-\cos\delta\sin H
}
$$

$$
\boxed{
s_N=
\cos\varphi\sin\delta
-
\sin\varphi\cos\delta\cos H
}
$$

$$
\boxed{
s_U=
\sin\varphi\sin\delta
+
\cos\varphi\cos\delta\cos H
}
$$

因此太阳在题目镜场坐标系中的单位方向为

$$
\boxed{
\mathbf s=
\begin{pmatrix}
-\cos\delta\sin H\\
\cos\varphi\sin\delta-\sin\varphi\cos\delta\cos H\\
\sin\varphi\sin\delta+\cos\varphi\cos\delta\cos H
\end{pmatrix}
}
$$

其中 $\mathbf s$ 的物理方向规定为：

$$
\boxed{\text{定日镜}\rightarrow\text{太阳中心}}
$$

而太阳光实际传播方向为：

$$
\boxed{
\mathbf d_{\mathrm{in}}=-\mathbf s
}
$$

---

## 9. 与题目太阳高度角公式进行验证

题目给出的太阳高度角满足：

$$
\sin\alpha_s
=
\cos\delta\cos\varphi\cos H
+
\sin\delta\sin\varphi.
$$

而模型得到：

$$
s_U
=
\sin\varphi\sin\delta
+
\cos\varphi\cos\delta\cos H.
$$

因此：

$$
\boxed{
s_U=\sin\alpha_s
}
$$

二者完全一致。

这说明所建立的太阳 ENU 向量模型与题目给出的太阳高度角公式一致。

太阳方位角若定义为从正北开始顺时针增加，则可以利用：

$$
\boxed{
\gamma_s
=
\operatorname{atan2}(s_E,s_N)
}
$$

统一判断上午和下午的方位象限，避免仅使用 $\arccos$ 时产生象限歧义。

---

## 10. 定日镜与集热器空间位置

对于第 $i$ 面定日镜，其中心坐标记为：

$$
\boxed{
\mathbf H_i=
\begin{pmatrix}
x_i\\
y_i\\
z_i
\end{pmatrix}
}
$$

问题一中所有定日镜安装高度均为：

$$
z_i=4\text{ m}.
$$

因此：

$$
\mathbf H_i=(x_i,y_i,4)^T.
$$

问题一中吸收塔位于镜场中心，集热器中心高度为 $80$ m，因此：

$$
\boxed{
\mathbf C=
\begin{pmatrix}
0\\
0\\
80
\end{pmatrix}
}
$$

---

## 11. 镜面到集热器的目标反射方向

从第 $i$ 面定日镜中心指向集热器中心：

$$
\mathbf C-\mathbf H_i.
$$

归一化得到目标反射方向：

$$
\boxed{
\mathbf r_i
=
\frac{\mathbf C-\mathbf H_i}
{\|\mathbf C-\mathbf H_i\|}
}
$$

其中：

$$
\boxed{
\mathbf r_i:
\text{定日镜中心}\rightarrow\text{集热器中心}
}
$$

---

## 12. 利用反射定律求定日镜法向

太阳光实际入射方向：

$$
\mathbf d_{\mathrm{in}}=-\mathbf s.
$$

要求中心反射光线指向集热器：

$$
\mathbf d_{\mathrm{out}}=\mathbf r_i.
$$

镜面反射公式为：

$$
\mathbf d_{\mathrm{out}}
=
\mathbf d_{\mathrm{in}}
-
2(\mathbf d_{\mathrm{in}}\cdot\mathbf n_i)\mathbf n_i.
$$

代入：

$$
\mathbf r_i
=
-\mathbf s
-
2(-\mathbf s\cdot\mathbf n_i)\mathbf n_i.
$$

整理：

$$
\mathbf r_i+\mathbf s
=
2(\mathbf s\cdot\mathbf n_i)\mathbf n_i.
$$

因此：

$$
\mathbf n_i
\parallel
\mathbf s+\mathbf r_i.
$$

归一化得到：

$$
\boxed{
\mathbf n_i
=
\frac{\mathbf s+\mathbf r_i}
{\|\mathbf s+\mathbf r_i\|}
}
$$

这也是“镜面法向为太阳方向与目标反射方向的角平分方向”的向量表达。

---

## 13. 法向量模型的物理意义

对于同一时刻，镜场尺度相对于日地距离极小，因此各定日镜可近似采用相同太阳中心方向：

$$
\mathbf s.
$$

但由于不同定日镜的位置 $\mathbf H_i$ 不同，其指向集热器的方向：

$$
\mathbf r_i
$$

不同，因此：

$$
\boxed{
\mathbf n_i
=
\mathbf n_i(\mathbf H_i,t)
}
$$

也随定日镜位置和时间变化。

换言之：

$$
\boxed{
\text{时间决定太阳方向}
+
\text{空间位置决定目标反射方向}
\rightarrow
\text{共同决定镜面姿态}
}
$$

---

## 14. 数值实现时必须进行的验证

### 验证 1：太阳方向单位化

应满足：

$$
\boxed{
\|\mathbf s\|=1
}
$$

数值计算中检查：

$$
\left|\|\mathbf s\|-1\right|<\varepsilon.
$$

### 验证 2：目标反射方向单位化

$$
\boxed{
\|\mathbf r_i\|=1
}
$$

### 验证 3：镜面法向单位化

$$
\boxed{
\|\mathbf n_i\|=1
}
$$

### 验证 4：反射定律回代

根据计算得到的 $\mathbf n_i$，重新计算：

$$
\mathbf d_{\mathrm{ref},i}
=
-\mathbf s
-
2[(-\mathbf s)\cdot\mathbf n_i]\mathbf n_i.
$$

理论上应满足：

$$
\boxed{
\mathbf d_{\mathrm{ref},i}
=
\mathbf r_i
}
$$

因此可以定义反射残差：

$$
\boxed{
\varepsilon_{\mathrm{ref},i}
=
\|\mathbf d_{\mathrm{ref},i}-\mathbf r_i\|
}
$$

要求：

$$
\varepsilon_{\mathrm{ref},i}\ll1.
$$

该残差可以作为后续程序正确性的第一项独立验证指标。

---

## 15. 当前模型中的简化

### 简化 1：太阳采用中心方向

镜面法向由太阳中心方向确定：

$$
\mathbf s.
$$

题目说明太阳具有有限锥角，但控制系统本身要求太阳中心光线经镜面中心反射后指向集热器中心，因此基础姿态模型使用太阳中心方向。

有限太阳角将在后续截断效率模型中进一步考虑。

### 简化 2：忽略镜场尺度导致的太阳视差

镜场半径仅为数百米，而日地距离为天文尺度，因此所有定日镜可采用同一太阳中心方向 $\mathbf s$。

---

## 16. 与后续模型的接口

本模块最终为每一个日期、时刻和定日镜输出：

$$
\boxed{
\mathbf s,\qquad
\mathbf r_i,\qquad
\mathbf n_i
}
$$

其中：

- $\mathbf s$ 用于余弦效率和太阳侧阴影计算；
- $\mathbf r_i$ 用于反射侧挡光计算；
- $\mathbf n_i$ 用于确定镜面真实空间姿态；
- $\mathbf H_i$ 与 $\mathbf n_i$ 将进一步用于构造定日镜矩形的四个顶点。

下一步首先建立：

$$
\boxed{\eta_{\cos,i}}
$$

余弦效率模型。

随后再进入：

$$
\eta_{sb,i},
\qquad
\eta_{at,i},
\qquad
\eta_{\mathrm{trunc},i}.
$$

---

## 17. 配套图件

计划配套以下推导图：

- `figs/q1_deriv_hour_declination_system.png`
  - 时角—赤纬辅助坐标系；
- `figs/q1_deriv_hourangle_to_enu.png`
  - 辅助坐标系向当地 ENU 镜场坐标系的转换；
- `figs/q1_deriv_heliostat_normal.png`
  - 太阳—定日镜—集热器反射几何与镜面法向。

图中使用的符号应与本节完全一致。

## Q1-02 太阳有限锥角的分层处理与模型简化依据

### 1. 问题来源

题目明确指出，太阳光并非严格平行光线，而是具有一定锥形角的一束锥形光线。

因此，从严格光学角度看，太阳不能完全表示为唯一的中心方向

$$
\mathbf s_0,
$$

而应理解为围绕太阳中心方向分布的一组入射方向：

$$
\boxed{
\mathbf s=\mathbf s(\beta,\psi)
}
$$

其中：

- $\beta$：相对于太阳中心方向的偏角；
- $\psi$：绕太阳中心轴的方位参数；
- $\mathbf s_0$：太阳视盘中心对应的单位方向。

因此，太阳有限角宽度理论上可能影响所有与光线方向有关的光学效率。

但是，各效率对太阳角宽度的敏感性明显不同。

若在所有效率计算中均进行完整太阳视盘积分或多光线追迹，将显著增加问题一的计算量，并进一步放大问题二、三优化过程中的计算成本。

因此，本模型采用

$$
\boxed{
\text{“按敏感程度分配模型精度”的分层处理策略}
}
$$

即：

$$
\boxed{
\text{弱敏感机制简化}
+
\text{中等敏感机制简化后验证}
+
\text{强敏感机制保留有限太阳角}
}
$$

而不是对所有效率统一采用同一精度等级。

---

## 2. 各项效率与太阳有限角宽度的关系

单面定日镜总光学效率为

$$
\eta_i
=
\eta_{\cos,i}
\eta_{sb,i}
\eta_{at,i}
\eta_{\mathrm{trunc},i}
\eta_{\mathrm{ref}}.
$$

分别考察太阳有限锥角对各项效率的影响。

---

## 3. 余弦效率：采用太阳中心方向近似

### 3.1 严格模型

对于单一太阳方向 $\mathbf s$，余弦效率为

$$
\eta_{\cos}
=
\mathbf n\cdot\mathbf s.
$$

若考虑有限太阳视盘，则不同入射光线具有不同方向

$$
\mathbf s(\beta,\psi),
$$

因此严格的平均余弦效率应写成对太阳视盘方向的辐射加权平均：

$$
\eta_{\cos}^{\mathrm{ext}}
=
\frac{
\displaystyle
\int_{\Omega_\odot}
L(\mathbf s)
\left(
\mathbf n\cdot\mathbf s
\right)
\,d\Omega
}{
\displaystyle
\int_{\Omega_\odot}
L(\mathbf s)\,d\Omega
}.
$$

其中 $L(\mathbf s)$ 表示太阳视盘不同方向的辐射强度。

---

### 3.2 为什么可以采用中心光线？

以太阳中心方向 $\mathbf s_0$ 为轴建立两个正交方向
$\mathbf e_1,\mathbf e_2$，则太阳视盘内任意方向可表示为

$$
\mathbf s(\beta,\psi)
=
\cos\beta\,\mathbf s_0
+
\sin\beta
\left(
\cos\psi\,\mathbf e_1
+
\sin\psi\,\mathbf e_2
\right).
$$

因此

$$
\mathbf n\cdot\mathbf s
=
\cos\beta
\left(
\mathbf n\cdot\mathbf s_0
\right)
+
\sin\beta
\left[
(\mathbf n\cdot\mathbf e_1)\cos\psi
+
(\mathbf n\cdot\mathbf e_2)\sin\psi
\right].
$$

由于太阳视盘关于中心方向近似轴对称，

$$
\int_0^{2\pi}\cos\psi\,d\psi=0,
$$

$$
\int_0^{2\pi}\sin\psi\,d\psi=0,
$$

故一阶方向扰动在太阳视盘平均中相互抵消。

又由于太阳角宽度较小，

$$
\cos\beta
=
1-\frac{\beta^2}{2}+O(\beta^4),
$$

因此：

$$
\boxed{
\eta_{\cos}^{\mathrm{ext}}
=
\eta_{\cos}^{\mathrm{center}}
+
O(\beta^2)
}
$$

即有限太阳角对余弦效率主要产生二阶小量修正。

因此问题一主模型采用太阳中心方向计算：

$$
\boxed{
\eta_{\cos,i}
\approx
\mathbf n_i\cdot\mathbf s_0
}
$$

该简化同时避免了对每个太阳视盘方向重复计算余弦投影。

### 3.3 模型决策

$$
\boxed{
\eta_{\cos}：
\text{采用中心太阳光线，不显式展开太阳锥角}
}
$$

理由：

1. 太阳视盘角宽度较小；
2. 太阳视盘近似轴对称，一阶方向扰动平均后抵消；
3. 余弦效率关于方向角连续、平滑；
4. 有限太阳角修正相对于其他光学损失较小；
5. 在问题二、三中，该简化可显著降低重复计算量。

---

## 4. 大气透射率：严格按照题设经验公式计算

题目直接给出：

$$
\eta_{at}
=
0.99321
-
0.0001176d_{HR}
+
1.97\times10^{-8}d_{HR}^{2},
$$

其中

$$
d_{HR}
=
\|\mathbf C-\mathbf H_i\|.
$$

因此本题的大气透射率被题目明确简化为：

$$
\boxed{
\eta_{at}=f(d_{HR})
}
$$

即仅由定日镜中心至集热器中心的距离确定。

严格来说，太阳视盘边缘光线与中心光线的实际传播距离存在极小差别，但这一效应：

1. 相对于数百米的传播距离极小；
2. 已超出题目给定经验公式的精度层级；
3. 若单独对此部分进行高精度光线积分，会造成模型各模块精度等级不一致。

因此不对题设公式继续精细化。

### 模型决策

$$
\boxed{
\eta_{at}：
\text{完全按照题设中心距离经验公式计算}
}
$$

不考虑太阳锥角引起的传播路径微小差异。

---

## 5. 镜面反射率：采用题设常数

题目允许镜面反射率取常数，例如：

$$
\boxed{
\eta_{\mathrm{ref}}=0.92
}
$$

现实情况下，反射率可能受到：

- 入射角；
- 光谱；
- 镜面粗糙度；
- 污染程度；

等因素影响。

但这些机制均已被题目主动简化。

因此本模型不再增加太阳入射角或太阳视盘方向对镜面反射率的修正。

### 模型决策

$$
\boxed{
\eta_{\mathrm{ref}}=0.92
}
$$

不考虑太阳锥角。

---

## 6. 阴影遮挡效率：主模型采用中心光线，保留有限太阳盘验证

阴影遮挡效率包含两种不同机制：

$$
\boxed{
\text{Shadowing}
}
$$

以及

$$
\boxed{
\text{Blocking}
}
$$

其中：

- Shadowing：邻近定日镜阻挡太阳到目标镜面的入射光；
- Blocking：邻近定日镜阻挡目标镜面反射至集热器的光线。

---

### 6.1 有限太阳角对 Shadowing 的影响

若将太阳视为点光源，则太阳方向唯一：

$$
\mathbf s_0.
$$

在该假设下，邻镜在目标镜面上的阴影具有确定的几何边界。

若考虑有限太阳盘，不同太阳盘方向对应：

$$
\mathbf s(\beta,\psi),
$$

不同方向产生的阴影边界略有差异。

因此真实阴影边缘会出现部分太阳盘可见、部分太阳盘被遮挡的区域，即：

$$
\boxed{
\text{半影区域（penumbra）}
}
$$

所以严格来说，太阳有限角宽度会影响阴影效率。

---

### 6.2 有限太阳角对 Blocking 的影响

类似地，不同入射方向

$$
\mathbf s(\beta,\psi)
$$

经过同一镜面后会产生略有差异的反射方向：

$$
\mathbf d_{\mathrm{out}}(\beta,\psi).
$$

因此邻镜可能只阻挡一部分太阳盘方向所对应的反射光。

所以有限太阳角理论上也影响挡光效率。

---

### 6.3 为什么主模型仍可采用中心太阳光线？

阴影遮挡效率计算本身已经需要处理：

$$
N
$$

面镜之间的空间几何干涉。

若进一步对太阳视盘离散为 $N_{\mathrm{sun}}$ 个方向，则计算规模近似增加为：

$$
N_{\mathrm{sun}}
\times
N
\times
N_{\mathrm{neighbor}}.
$$

这会显著增加问题一计算量，并在后续镜场优化中成为主要计算瓶颈。

已有定日镜阴影遮挡模型采用太阳中心光线进行投影和 ray tracing，并能够取得较高精度，因此中心光线模型本身是一种具有工程实践基础的快速近似。

因此主模型首先采用：

$$
\boxed{
\mathbf s(\beta,\psi)
\approx
\mathbf s_0
}
$$

建立几何阴影与挡光模型。

### 6.4 但不能直接认为该近似始终成立

阴影效率与几何边界有关，其对方向扰动的敏感性明显高于余弦效率。

特别是在：

- 太阳高度较低；
- 镜场排列较密；
- 邻镜投影接近镜面边界；

等工况下，有限太阳盘产生的半影可能变得更加明显。

因此本模型将中心太阳光线视为

$$
\boxed{
\text{阴影遮挡主模型的计算简化}
}
$$

而不是物理上的严格事实。

后续应选取典型高遮挡工况进行有限太阳盘敏感性验证。

定义：

$$
\eta_{sb}^{(0)}
=
\text{中心太阳光线模型结果},
$$

$$
\eta_{sb}^{(\mathrm{ext})}
=
\text{有限太阳盘模型结果}.
$$

比较：

$$
\boxed{
\Delta\eta_{sb}
=
\left|
\eta_{sb}^{(\mathrm{ext})}
-
\eta_{sb}^{(0)}
\right|
}
$$

若差异足够小，则保留中心光线简化；

若差异不可忽略，则将阴影遮挡模型升级为有限太阳盘模型。

### 模型决策

$$
\boxed{
\eta_{sb}：
\text{中心太阳光线主模型}
+
\text{有限太阳盘敏感性验证}
}
$$

---

## 7. 截断效率：不采用点太阳近似

截断效率定义为：

$$
\eta_{\mathrm{trunc}}
=
\frac{
\text{集热器实际接收能量}
}{
\text{扣除阴影遮挡后的镜面反射能量}
}.
$$

它所回答的问题本质上是：

$$
\boxed{
\text{定日镜形成的反射光斑有多少真正落在有限尺寸集热器上？}
}
$$

---

### 7.1 为什么太阳角宽度不能在这里忽略？

若只考虑太阳中心光线，定日镜中心光线在控制系统作用下被精确反射至集热器中心。

如果进一步把全部太阳光都视为与中心光线严格平行，则会严重低估反射光斑的空间扩展范围。

设太阳方向相对于中心方向存在角度偏差：

$$
\Delta\theta.
$$

光线传播距离约为：

$$
d_{HR}.
$$

则有限太阳角在接收器附近造成的横向光斑展宽尺度约为：

$$
\boxed{
\Delta l
\sim
d_{HR}\Delta\theta
}
$$

因此该影响会随着镜面至接收器距离增加而累计放大。

本题镜场空间尺度达到数百米，而集热器直径仅为有限的数米量级，因此太阳角宽度所造成的光斑扩展与集热器尺寸处于不可忽略的相对尺度。

此时有限太阳角直接决定：

$$
\boxed{
\text{光线命中}
\quad\text{或}\quad
\text{光线溢出}
}
$$

即直接影响截断效率。

因此，截断效率不能仅使用中心光线模型。

---

## 8. 截断效率的模型要求

在截断效率模块中，应将中心太阳方向

$$
\mathbf s_0
$$

扩展为太阳视盘方向集合：

$$
\boxed{
\mathbf s(\beta,\psi)
}
$$

再计算各入射方向经过有限尺寸镜面反射后与有限圆柱形集热器的空间关系。

该模块至少需要同时考虑：

1. 太阳有限角宽度；
2. 定日镜有限尺寸；
3. 定日镜空间姿态；
4. 镜面至集热器传播距离；
5. 集热器有限高度；
6. 集热器有限直径；
7. 圆柱外表受光几何。

因此截断效率是问题一中需要进行精细光学建模的核心模块之一。

### 模型决策

$$
\boxed{
\eta_{\mathrm{trunc}}：
\text{保留有限太阳角，不采用点太阳近似}
}
$$

后续再经过结构审查决定采用：

- 解析光斑模型；
- 确定性离散光线追迹；
- Monte Carlo 光线追迹；

或其他计算方式。

现阶段不提前确定数值算法。

---

## 9. 各效率的最终精度分配

| 效率 | 太阳有限角理论影响 | 主模型处理 | 后续验证 |
|---|---|---|---|
| $\eta_{\cos}$ | 弱 | 中心太阳方向 | 可选 |
| $\eta_{at}$ | 极弱，题目已经验化 | 完全按题设公式 | 无需 |
| $\eta_{\mathrm{ref}}$ | 题目已常数化 | $\eta_{\mathrm{ref}}=0.92$ | 无需 |
| $\eta_{sb}$ | 中等，边界处产生半影 | 中心太阳方向 | 必须选典型工况验证 |
| $\eta_{\mathrm{trunc}}$ | 强 | 保留有限太阳角 | 后续需做收敛验证 |

因此问题一不采用“所有模块统一精度”的建模方式，而采用：

$$
\boxed{
\text{余弦效率：简化}
}
$$

$$
\boxed{
\text{阴影遮挡：简化 + 验证}
}
$$

$$
\boxed{
\text{截断效率：精细建模}
}
$$

---

## 10. 本次简化的建模原则

本次处理体现一个可迁移的建模原则：

> 真实物理机制可能同时作用于多个子模块，但并不意味着必须在所有模块中以相同精度展开。应首先分析各子模块对该机制的敏感性，再将计算复杂度投入到对最终结果影响最大的部分。

可以概括为：

$$
\boxed{
\text{复杂物理效应}
\xrightarrow{\text{敏感性判断}}
\begin{cases}
\text{弱敏感} & \rightarrow \text{简化},\\
\text{中等敏感} & \rightarrow \text{简化并验证},\\
\text{强敏感} & \rightarrow \text{精细建模}.
\end{cases}
}
$$

这属于：

$$
\boxed{
\text{复杂机制局部化}
}
$$

而不是全局增加模型复杂度。

---

## 11. 当前模型状态

目前问题一各效率处理状态为：

$$
\eta_{\mathrm{ref}}
\quad\checkmark
$$

题设常数，已闭合；

$$
\eta_{at}
\quad\checkmark
$$

题设经验公式，已闭合；

$$
\eta_{\cos}
\quad\checkmark
$$

中心太阳光线模型，待数值实现；

$$
\eta_{sb}
\quad\triangle
$$

采用中心光线建立主模型，后续进行有限太阳盘敏感性验证；

$$
\eta_{\mathrm{trunc}}
\quad\triangle
$$

必须建立有限太阳角下的精细模型。

因此下一阶段首先完成：

$$
\boxed{
\eta_{\cos}
}
$$

随后重点进入：

$$
\boxed{
\eta_{sb}
}
$$

和

$$
\boxed{
\eta_{\mathrm{trunc}}
}
$$

两个尚未闭合的几何光学模块。

## Q1-03 余弦效率模型

### 1. 物理意义

定日镜接收的是法向直接辐射辐照度 $DNI$。

若太阳光线垂直入射镜面，则面积为 $A_i$ 的定日镜对太阳光的有效迎光面积等于其真实面积：

$$
A_{\mathrm{eff},i}=A_i.
$$

此时不存在余弦损失，故：

$$
\eta_{\cos,i}=1.
$$

若太阳光与镜面法向存在夹角 $\theta_i$，则从太阳入射方向观察，镜面的有效迎光面积为镜面在垂直太阳光方向平面上的正投影面积：

$$
\boxed{
A_{\mathrm{eff},i}
=
A_i\cos\theta_i
}
$$

因此该镜面能够截获的太阳辐射功率为：

$$
P_{\mathrm{inc},i}
=
DNI\cdot A_i\cos\theta_i.
$$

若镜面完全正对太阳，则理想接收功率为：

$$
P_{\mathrm{ideal},i}
=
DNI\cdot A_i.
$$

故余弦效率定义为：

$$
\eta_{\cos,i}
=
\frac{P_{\mathrm{inc},i}}
{P_{\mathrm{ideal},i}},
$$

即：

$$
\boxed{
\eta_{\cos,i}
=
\cos\theta_i
}
$$

因此，余弦损失的本质是：

$$
\boxed{
\text{镜面倾斜}
\rightarrow
\text{有效迎光投影面积减小}
\rightarrow
\text{可截获太阳辐射减少}
}
$$

---

### 2. 余弦效率的向量表达

前文定义太阳中心单位方向向量：

$$
\mathbf s
$$

其方向为：

$$
\boxed{
\text{定日镜中心}
\rightarrow
\text{太阳中心}
}
$$

第 $i$ 面定日镜的单位法向量记为：

$$
\mathbf n_i.
$$

设两者夹角为：

$$
\theta_i
=
\angle(\mathbf s,\mathbf n_i).
$$

由向量点积定义：

$$
\mathbf s\cdot\mathbf n_i
=
\|\mathbf s\|
\|\mathbf n_i\|
\cos\theta_i.
$$

由于：

$$
\|\mathbf s\|=1,
\qquad
\|\mathbf n_i\|=1,
$$

故：

$$
\boxed{
\eta_{\cos,i}
=
\mathbf s\cdot\mathbf n_i
}
$$

若在镜场 ENU 坐标系中写成分量形式：

$$
\mathbf s
=
\begin{pmatrix}
s_E\\
s_N\\
s_U
\end{pmatrix},
\qquad
\mathbf n_i
=
\begin{pmatrix}
n_{E,i}\\
n_{N,i}\\
n_{U,i}
\end{pmatrix},
$$

则：

$$
\boxed{
\eta_{\cos,i}
=
s_E n_{E,i}
+
s_N n_{N,i}
+
s_U n_{U,i}
}
$$

---

### 3. 入射传播方向的符号说明

太阳单位方向 $\mathbf s$ 的定义为：

$$
\text{定日镜}\rightarrow\text{太阳}.
$$

而太阳光实际传播方向恰好相反：

$$
\boxed{
\mathbf d_{\mathrm{in}}
=
-\mathbf s
}
$$

因此如果使用真实光线传播方向，则：

$$
\mathbf d_{\mathrm{in}}\cdot\mathbf n_i<0.
$$

此时余弦效率应表示为：

$$
\eta_{\cos,i}
=
-\mathbf d_{\mathrm{in}}\cdot\mathbf n_i.
$$

由于：

$$
-\mathbf d_{\mathrm{in}}
=
\mathbf s,
$$

仍可得到：

$$
\boxed{
\eta_{\cos,i}
=
\mathbf s\cdot\mathbf n_i.
}
$$

为避免后续公式频繁处理负号，本文统一采用从镜面指向太阳的单位向量 $\mathbf s$。

---

### 4. 利用反射几何进一步化简

前文已经建立第 $i$ 面定日镜指向集热器中心的目标反射方向：

$$
\boxed{
\mathbf r_i
=
\frac{
\mathbf C-\mathbf H_i
}{
\|\mathbf C-\mathbf H_i\|
}
}
$$

其中：

- $\mathbf H_i$：第 $i$ 面定日镜中心；
- $\mathbf C$：集热器中心。

根据反射定律，定日镜法向为太阳方向 $\mathbf s$ 与目标反射方向 $\mathbf r_i$ 的角平分方向：

$$
\boxed{
\mathbf n_i
=
\frac{
\mathbf s+\mathbf r_i
}{
\|\mathbf s+\mathbf r_i\|
}
}
$$

代入余弦效率：

$$
\eta_{\cos,i}
=
\mathbf s\cdot
\frac{
\mathbf s+\mathbf r_i
}{
\|\mathbf s+\mathbf r_i\|
}.
$$

展开分子：

$$
\eta_{\cos,i}
=
\frac{
\mathbf s\cdot\mathbf s
+
\mathbf s\cdot\mathbf r_i
}{
\|\mathbf s+\mathbf r_i\|
}.
$$

由于：

$$
\|\mathbf s\|=1,
$$

故：

$$
\mathbf s\cdot\mathbf s=1.
$$

于是：

$$
\eta_{\cos,i}
=
\frac{
1+\mathbf s\cdot\mathbf r_i
}{
\|\mathbf s+\mathbf r_i\|
}.
$$

另一方面：

$$
\|\mathbf s+\mathbf r_i\|^2
=
(\mathbf s+\mathbf r_i)
\cdot
(\mathbf s+\mathbf r_i).
$$

展开：

$$
\|\mathbf s+\mathbf r_i\|^2
=
\mathbf s\cdot\mathbf s
+
2\mathbf s\cdot\mathbf r_i
+
\mathbf r_i\cdot\mathbf r_i.
$$

由于：

$$
\|\mathbf s\|
=
\|\mathbf r_i\|
=
1,
$$

因此：

$$
\|\mathbf s+\mathbf r_i\|^2
=
2
+
2\mathbf s\cdot\mathbf r_i.
$$

即：

$$
\boxed{
\|\mathbf s+\mathbf r_i\|
=
\sqrt{
2(1+\mathbf s\cdot\mathbf r_i)
}
}
$$

代回余弦效率公式：

$$
\eta_{\cos,i}
=
\frac{
1+\mathbf s\cdot\mathbf r_i
}{
\sqrt{
2(1+\mathbf s\cdot\mathbf r_i)
}
}.
$$

最终得到：

$$
\boxed{
\eta_{\cos,i}
=
\sqrt{
\frac{
1+\mathbf s\cdot\mathbf r_i
}{
2
}
}
}
$$

---

### 5. 半角形式及几何意义

设太阳方向 $\mathbf s$ 与目标反射方向 $\mathbf r_i$ 的夹角为：

$$
\phi_i
=
\angle(\mathbf s,\mathbf r_i).
$$

由于二者均为单位向量：

$$
\mathbf s\cdot\mathbf r_i
=
\cos\phi_i.
$$

因此：

$$
\eta_{\cos,i}
=
\sqrt{
\frac{
1+\cos\phi_i
}{
2
}
}.
$$

利用半角公式：

$$
\cos\frac{\phi_i}{2}
=
\sqrt{
\frac{
1+\cos\phi_i
}{
2
}
},
$$

得到：

$$
\boxed{
\eta_{\cos,i}
=
\cos\frac{\phi_i}{2}
}
$$

另一方面，由于镜面法向 $\mathbf n_i$ 是 $\mathbf s$ 与 $\mathbf r_i$ 的角平分方向：

$$
\theta_i
=
\angle(\mathbf s,\mathbf n_i)
=
\frac{\phi_i}{2}.
$$

于是：

$$
\eta_{\cos,i}
=
\cos\theta_i
=
\cos\frac{\phi_i}{2},
$$

与面积投影模型完全一致。

因此余弦效率可以从两个等价角度理解：

$$
\boxed{
\eta_{\cos,i}
=
\mathbf s\cdot\mathbf n_i
}
$$

表示太阳方向与镜面法向之间的投影关系；

或：

$$
\boxed{
\eta_{\cos,i}
=
\sqrt{
\frac{
1+\mathbf s\cdot\mathbf r_i
}{2}
}
}
$$

表示太阳方向与镜面指向集热器方向之间的几何关系。

---

### 6. 余弦损失的几何本质

上述化简表明，余弦效率最终取决于：

$$
\boxed{
\text{太阳方向}
\quad\mathbf s
}
$$

与：

$$
\boxed{
\text{镜面到集热器方向}
\quad\mathbf r_i
}
$$

之间的夹角。

定日镜必须同时满足：

1. 接收太阳辐射；
2. 将太阳光反射至集热器。

因此镜面法向必须在太阳方向和目标反射方向之间进行几何折中。

当：

$$
\phi_i
$$

较小时，镜面能够较正地迎接太阳光：

$$
\eta_{\cos,i}\rightarrow1.
$$

当：

$$
\phi_i
$$

增大时，镜面倾斜程度增加，有效迎光投影面积降低：

$$
\eta_{\cos,i}\downarrow.
$$

因此可以概括为：

$$
\boxed{
\phi_i\uparrow
\Rightarrow
\theta_i=\frac{\phi_i}{2}\uparrow
\Rightarrow
A_{\mathrm{eff},i}\downarrow
\Rightarrow
\eta_{\cos,i}\downarrow
}
$$

---

### 7. 空间和时间依赖关系

在同一规定工况 $t$ 下，整个镜场可以近似采用相同太阳中心方向：

$$
\mathbf s(t).
$$

但不同定日镜中心位置 $\mathbf H_i$ 不同，因此：

$$
\mathbf r_i
=
\frac{
\mathbf C-\mathbf H_i
}{
\|\mathbf C-\mathbf H_i\|
}
$$

不同。

因此同一时刻不同定日镜的余弦效率不同：

$$
\boxed{
\eta_{\cos,i}
=
\eta_{\cos}(\mathbf H_i,t)
}
$$

同时，对于固定的定日镜：

$$
\mathbf r_i
$$

在问题一中不随时间改变，但太阳方向：

$$
\mathbf s(t)
$$

随日期和时刻改变，因此余弦效率也随时间变化。

这说明余弦效率同时具有：

$$
\boxed{
\text{空间非均匀性}
+
\text{时间变化性}
}
$$

---

### 8. 两种公式的代码交叉验证

数值实现中同时计算：

$$
\eta_{\cos,i}^{(1)}
=
\mathbf s\cdot\mathbf n_i
$$

和：

$$
\eta_{\cos,i}^{(2)}
=
\sqrt{
\frac{
1+\mathbf s\cdot\mathbf r_i
}{
2
}
}.
$$

理论上：

$$
\boxed{
\eta_{\cos,i}^{(1)}
=
\eta_{\cos,i}^{(2)}
}
$$

因此定义交叉验证残差：

$$
\boxed{
\varepsilon_{\cos,i}
=
\left|
\eta_{\cos,i}^{(1)}
-
\eta_{\cos,i}^{(2)}
\right|
}
$$

数值计算中应满足：

$$
\varepsilon_{\cos,i}\ll1.
$$

该验证不依赖最终镜场功率结果，而是直接验证：

- 太阳方向；
- 目标反射方向；
- 镜面法向；
- 余弦效率；

四个模块之间的内部一致性。

---

### 9. 当前结论

问题一中第 $i$ 面定日镜在工况 $t$ 下的余弦效率模型最终确定为：

$$
\boxed{
\eta_{\cos,i}(t)
=
\mathbf s(t)\cdot\mathbf n_i(t)
}
$$

并具有等价化简式：

$$
\boxed{
\eta_{\cos,i}(t)
=
\sqrt{
\frac{
1+\mathbf s(t)\cdot\mathbf r_i
}{
2
}
}
}
$$

其中第一式作为主模型表达，具有最直接的物理含义；

第二式作为几何化简式和程序交叉验证式。

余弦效率模型至此闭合。


## Q1-04 阴影遮挡效率模型

### 1. 模型目标

题目将阴影遮挡效率记为

$$
\eta_{sb,i},
$$

其本质是衡量第 $i$ 面定日镜中有多少有效镜面面积能够同时满足：

1. 能够接收到太阳直接辐射；
2. 反射后的光线不会被其他定日镜阻挡。

因此，阴影遮挡损失包含三类几何机制：

$$
\boxed{
\text{集热器阴影}
+
\text{邻镜 Shadowing}
+
\text{邻镜 Blocking}
}
$$

分别记为：

$$
T_i,\qquad
S_i,\qquad
B_i.
$$

其中：

- $T_i$：圆柱集热器在目标镜 $i$ 上形成的太阳阴影；
- $S_i$：其他定日镜阻挡太阳入射光形成的阴影区域；
- $B_i$：其他定日镜阻挡目标镜反射光形成的挡光区域。

最终将三类无效区域统一映射至目标镜自身二维平面，并通过集合并运算消除重复计数。

---

# 2. 点太阳近似

真实太阳具有有限视盘角宽度，不同太阳盘位置对应略有不同的入射方向。

在阴影遮挡主模型中，为降低大规模镜场计算复杂度，将有限太阳视盘退化为太阳中心方向

$$
\mathbf s.
$$

并认为在镜场尺度内，太阳中心对应的入射光线彼此平行。

因此实际传播方向为

$$
\boxed{
\mathbf d_{\mathrm{in}}
=
-\mathbf s
}
$$

需要强调：

> “点太阳”并不是将太阳看作有限距离的小型点光源，而是将有限太阳视盘的多方向平行光束压缩为太阳中心方向的一束平行光。

有限太阳视盘造成的半影效应暂不在阴影遮挡主模型中展开，后续可通过敏感性分析进行验证；太阳有限锥角将在截断效率模型中显式保留。

---

# 3. 有限尺寸定日镜的三维几何模型

前文已获得第 $i$ 面定日镜：

$$
\mathbf H_i
=
(x_i,y_i,z_i)^T
$$

及其单位法向：

$$
\mathbf n_i.
$$

但是“中心 + 法向”只能确定镜面所在平面，仍不能唯一确定矩形在平面内的转角。

题目规定定日镜上下两条边始终平行于地面。

令竖直单位向量为

$$
\mathbf k
=
\begin{pmatrix}
0\\
0\\
1
\end{pmatrix}.
$$

镜面水平边方向应同时：

1. 位于镜面内；
2. 垂直于竖直方向。

因此定义镜面宽度方向

$$
\boxed{
\mathbf e_{w,i}
=
\frac{
\mathbf k\times\mathbf n_i
}{
\|\mathbf k\times\mathbf n_i\|
}
}
$$

再由右手关系得到镜面高度方向

$$
\boxed{
\mathbf e_{h,i}
=
\mathbf n_i\times\mathbf e_{w,i}
}
$$

从而形成镜面局部正交基：

$$
\boxed{
\{
\mathbf e_{w,i},
\mathbf e_{h,i},
\mathbf n_i
\}
}
$$

对于宽度 $w_i$、高度 $h_i$ 的矩形定日镜，其任意一点可以表示为

$$
\boxed{
\mathbf P_i(u,v)
=
\mathbf H_i
+
u\mathbf e_{w,i}
+
v\mathbf e_{h,i}
}
$$

其中

$$
-\frac{w_i}{2}
\le u
\le
\frac{w_i}{2},
$$

$$
-\frac{h_i}{2}
\le v
\le
\frac{h_i}{2}.
$$

问题一中

$$
w_i=h_i=6\text{ m},
$$

因此目标镜在自身局部坐标系中恒表示为标准正方形

$$
\boxed{
R_i
=
[-3,3]\times[-3,3].
}
$$

镜面中心到角点的最大距离为半对角线

$$
\boxed{
\rho
=
\frac12\sqrt{6^2+6^2}
=
3\sqrt2
\approx4.243\text{ m}.
}
$$

该量后续用于候选邻镜的保守几何筛选。

---

# 4. 三维平行投影的一般公式

为统一 Shadowing 与 Blocking，引入一般平行投影模型。

设：

- $\mathbf Q$：待投影物体上的一点；
- $\mathbf d$：投影方向；
- $\mathbf H_i$：目标镜中心；
- $\mathbf n_i$：目标镜法向。

沿 $\mathbf d$ 作射线：

$$
\mathbf x(\lambda)
=
\mathbf Q+\lambda\mathbf d.
$$

目标镜平面满足

$$
\mathbf n_i\cdot
(\mathbf x-\mathbf H_i)
=
0.
$$

代入得到

$$
\mathbf n_i\cdot
(
\mathbf Q+\lambda\mathbf d-\mathbf H_i
)
=
0.
$$

因此

$$
\boxed{
\lambda
=
-
\frac{
\mathbf n_i\cdot(\mathbf Q-\mathbf H_i)
}{
\mathbf n_i\cdot\mathbf d
}
}
$$

对应投影点为

$$
\boxed{
\Pi_i(\mathbf Q;\mathbf d)
=
\mathbf Q
-
\frac{
\mathbf n_i\cdot(\mathbf Q-\mathbf H_i)
}{
\mathbf n_i\cdot\mathbf d
}
\mathbf d.
}
$$

仅当

$$
\boxed{
\lambda\ge0
}
$$

时，投影方向才与实际光路顺序一致。

---

# 5. 邻镜 Shadowing 模型

设 $j$ 为可能遮挡目标镜 $i$ 的邻镜。

太阳中心方向为

$$
\mathbf s,
$$

实际太阳光传播方向为

$$
-\mathbf s.
$$

因此 Shadowing 投影方向取

$$
\boxed{
\mathbf d_S
=
-\mathbf s.
}
$$

设镜 $j$ 的四个三维顶点为

$$
\mathbf Q_{j1},
\mathbf Q_{j2},
\mathbf Q_{j3},
\mathbf Q_{j4}.
$$

分别沿 $-\mathbf s$ 投影至目标镜平面：

$$
\boxed{
\mathbf Q'_{jk}
=
\Pi_i
(
\mathbf Q_{jk};
-\mathbf s
)
}
$$

得到镜 $j$ 在目标镜平面上的阴影投影四边形。

其物理含义为：

> 从挡光镜边界沿太阳中心光线传播方向作几何延拓，得到若无挡光镜存在时本应继续传播的光线位置，其与目标镜的交叠区域即为太阳阴影。

---

# 6. Blocking 模型

第 $i$ 面定日镜中心指向集热器中心的单位方向为

$$
\boxed{
\mathbf r_i
=
\frac{
\mathbf C-\mathbf H_i
}{
\|\mathbf C-\mathbf H_i\|
}.
}
$$

在当前

$$
\text{平面镜}
+
\text{太阳中心平行光}
$$

模型下，同一面定日镜上太阳中心光线的反射方向相同，均为

$$
\mathbf r_i.
$$

因此，若邻镜 $j$ 阻挡目标镜 $i$ 的反射光，可以沿反射方向的反方向

$$
\boxed{
\mathbf d_B
=
-\mathbf r_i
}
$$

将镜 $j$ 投影回目标镜平面。

于是镜 $j$ 顶点的挡光投影为

$$
\boxed{
\mathbf Q'_{jk}
=
\Pi_i
(
\mathbf Q_{jk};
-\mathbf r_i
)
}
$$

得到 Blocking 投影区域

$$
B_{ij}.
$$

需要注意：

> 定日镜跟踪控制保证的是“太阳中心光线经镜中心反射后指向集热器中心”，并不意味着整块平面镜所有光线严格汇聚于一点。对于理想平面镜，同一入射方向经过统一法向反射后仍保持平行。

因此本模型采用平行反射方向计算 Blocking。

---

# 7. 三维投影向目标镜二维局部坐标的降维

所有 Shadowing / Blocking 投影点均已经位于目标镜 $i$ 所在平面。

因此无需继续进行三维多边形运算。

对于任意投影点 $\mathbf P$，定义：

$$
\boxed{
u
=
(\mathbf P-\mathbf H_i)
\cdot
\mathbf e_{w,i}
}
$$

$$
\boxed{
v
=
(\mathbf P-\mathbf H_i)
\cdot
\mathbf e_{h,i}.
}
$$

于是完成映射：

$$
\boxed{
(x,y,z)
\longrightarrow
(u,v).
}
$$

此时目标镜始终为

$$
R_i=[-3,3]^2,
$$

而邻镜投影成为目标镜局部二维平面中的凸四边形。

因此原始三维镜间遮挡问题被转化为：

$$
\boxed{
\text{二维凸多边形相交与集合并问题}.
}
$$

---

# 8. 圆柱集热器阴影

问题一中集热器中心为

$$
\mathbf C=(0,0,80)^T,
$$

其半径为

$$
R_c=3.5\text{ m},
$$

高度为

$$
h_c=8\text{ m}.
$$

因此圆柱集热器区域为

$$
\boxed{
\mathcal C
=
\{
(x,y,z):
x^2+y^2\le3.5^2,\;
76\le z\le84
\}.
}
$$

由于题目未给出吸收塔支撑结构的横截面尺寸，本模型不人为增加塔身尺寸，仅对题目明确给定尺寸的圆柱集热器计算太阳阴影。

圆柱集热器属于凸集，平行投影保持凸性。

因此将其上下圆周分别进行确定性离散：

$$
\theta_k
=
\frac{2\pi k}{N_\theta},
$$

$$
\mathbf Q_k^{-}
=
\begin{pmatrix}
R_c\cos\theta_k\\
R_c\sin\theta_k\\
76
\end{pmatrix},
$$

$$
\mathbf Q_k^{+}
=
\begin{pmatrix}
R_c\cos\theta_k\\
R_c\sin\theta_k\\
84
\end{pmatrix}.
$$

沿太阳传播方向

$$
-\mathbf s
$$

投影至目标镜平面，并转化为 $(u,v)$ 坐标。

对所有投影边界点取二维凸包，即得到集热器阴影区域

$$
\boxed{
T_i.
}
$$

当前程序采用

$$
N_\theta=144.
$$

该方法属于确定性边界离散，而非 Monte Carlo 随机抽样。

---

# 9. 候选邻镜筛选

若对 $N$ 面镜逐对进行精确多边形运算，则每个工况需处理约

$$
N(N-1)
$$

组镜对。

问题一中

$$
N=1745,
$$

单工况理论镜对数量超过 $3\times10^6$。

因此在精确投影前进行保守候选筛选。

---

## 9.1 光路方向分解

对于目标镜 $i$ 和邻镜 $j$：

$$
\Delta\mathbf H_{ji}
=
\mathbf H_j-\mathbf H_i.
$$

给定单位主方向 $\mathbf d$，定义纵向分量：

$$
\boxed{
d_{\parallel}
=
\Delta\mathbf H_{ji}
\cdot
\mathbf d
}
$$

以及横向分量：

$$
\boxed{
\Delta\mathbf H_{\perp}
=
\Delta\mathbf H_{ji}
-
d_{\parallel}\mathbf d
}
$$

其横向距离为

$$
\boxed{
d_{\perp}
=
\|
\Delta\mathbf H_{\perp}
\|.
}
$$

由于每面镜完全包含于半径

$$
\rho=3\sqrt2
$$

的包围球内，若两镜在垂直主光路方向上的中心距离满足

$$
\boxed{
d_{\perp}>2\rho,
}
$$

则不存在一条平行于该光路方向的射线能够同时通过两面镜，因此该邻镜可安全排除。

---

## 9.2 Shadowing 候选

Shadowing 以朝太阳方向

$$
\mathbf d=\mathbf s
$$

进行镜面前后排序。

第一层保守候选条件为：

$$
\boxed{
d_{\perp}^{S}
\le2\rho,
}
$$

$$
\boxed{
d_{\parallel}^{S}
>-2\rho.
}
$$

通过该筛选的邻镜才进入实际太阳投影计算。

---

## 9.3 Blocking 候选

Blocking 以目标镜反射方向

$$
\mathbf d=\mathbf r_i
$$

作为主方向。

同时挡光镜必须位于目标镜与接收器之间的反射光路附近，因此采用：

$$
\boxed{
d_{\perp}^{B}
\le2\rho,
}
$$

$$
\boxed{
-2\rho
<
d_{\parallel}^{B}
<
d_{HR,i}+2\rho.
}
$$

其中：

$$
d_{HR,i}
=
\|\mathbf C-\mathbf H_i\|.
$$

上述条件仅用于排除绝不可能造成干涉的镜子，并不直接判断是否真正发生 Shadowing 或 Blocking。

---

# 10. 第二层 AABB 快速筛选

对第一层候选镜进行四顶点投影后，首先计算二维投影多边形：

$$
u_{\min},
\quad
u_{\max},
\quad
v_{\min},
\quad
v_{\max}.
$$

若满足：

$$
u_{\max}<-3,
$$

或

$$
u_{\min}>3,
$$

或

$$
v_{\max}<-3,
$$

或

$$
v_{\min}>3,
$$

则投影多边形的轴对齐包围盒与目标镜完全不相交，可直接排除。

只有通过该筛选的镜子才进行精确二维多边形相交计算。

因此整体算法采用：

$$
\boxed{
\text{光路管筛选}
\rightarrow
\text{AABB筛选}
\rightarrow
\text{精确多边形运算}.
}
$$

---

# 11. 多边形集合运算与阴影遮挡效率

定义：

$$
S_i
=
\bigcup_{j\in\mathcal S_i}
S_{ij},
$$

$$
B_i
=
\bigcup_{j\in\mathcal B_i}
B_{ij}.
$$

其中 $\mathcal S_i$、$\mathcal B_i$ 分别表示目标镜 $i$ 的 Shadowing 与 Blocking 候选集合。

三类无效区域统一表示为

$$
\boxed{
L_i
=
R_i
\cap
\left[
T_i
\cup
S_i
\cup
B_i
\right].
}
$$

这里必须采用集合并，而不能简单将各损失面积相加，因为：

- 不同邻镜阴影可能相互重叠；
- 不同 Blocking 区域可能重叠；
- Shadowing 与 Blocking 也可能发生空间重叠；
- 集热器阴影也可能与邻镜遮挡区域重合。

最终阴影遮挡效率定义为

$$
\boxed{
\eta_{sb,i}
=
1-
\frac{
\operatorname{Area}(L_i)
}{
A_i
}
}
$$

问题一中

$$
A_i=36\text{ m}^2.
$$

---

# 12. 数值实现

程序实现结构为：

$$
\boxed{
\text{镜面三维几何}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{候选邻镜筛选}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Shadowing / Blocking / Receiver Shadow 投影}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{局部二维 }(u,v)\text{ 坐标}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Polygon Intersection / Union}
}
$$

$$
\downarrow
$$

$$
\boxed{
A_{\mathrm{loss},i}
\rightarrow
\eta_{sb,i}.
}
$$

三维投影、局部坐标变换和候选镜筛选均由本文模型显式建立。

二维多边形交、并及面积作为标准计算几何运算，由成熟几何库执行。

---

# 13. 模型交叉验证

为了避免“程序能够运行”被误认为“几何模型正确”，本模型进行了两层独立验证。

## 13.1 候选筛选与全邻镜搜索验证

典型高遮挡工况下，对抽样目标镜分别采用：

1. 候选镜快速筛选模型；
2. 对其余全部 $1744$ 面镜进行暴力精确投影。

得到：

$$
\boxed{
\max
|
\eta_{sb}^{\mathrm{fast}}
-
\eta_{sb}^{\mathrm{full}}
|
=
0
}
$$

同时：

$$
\boxed{
\max
|
A_{\mathrm{loss}}^{\mathrm{fast}}
-
A_{\mathrm{loss}}^{\mathrm{full}}
|
=
0
}
$$

且两个损失区域的二维几何对称差面积为

$$
\boxed{0}.
$$

说明抽样验证中候选筛选未遗漏真实遮挡镜。

---

## 13.2 平行投影面积解析验证

平面区域沿单位方向 $\mathbf d$ 投影至另一平面时，完整投影面积应满足

$$
\boxed{
A_{\mathrm{proj}}
=
A_j
\frac{
|\mathbf n_j\cdot\mathbf d|
}{
|\mathbf n_i\cdot\mathbf d|
}.
}
$$

对程序产生的四顶点投影多边形面积与该解析公式进行交叉验证。

共验证

$$
5211
$$

组投影关系。

最大绝对误差：

$$
\boxed{
8.24\times10^{-13}\text{ m}^2
}
$$

最大相对误差：

$$
\boxed{
2.32\times10^{-14}
}
$$

平均相对误差：

$$
\boxed{
3.04\times10^{-15}.
}
$$

误差处于双精度浮点计算误差量级。

---

# 14. 阴影遮挡效率结果

对题目规定的

$$
12\times5=60
$$

个太阳工况计算全部 $1745$ 面定日镜阴影遮挡效率。

最终年平均阴影遮挡效率为：

$$
\boxed{
\bar\eta_{sb}
=
0.927943.
}
$$

其中12月21日9:00低太阳高度典型工况下：

$$
\bar\eta_{sb}=0.858600.
$$

该工况下平均 Shadowing 候选镜数约为

$$
1.68,
$$

平均 Blocking 候选镜数约为

$$
1.30,
$$

最大候选数均仅为

$$
3.
$$

说明基于光路方向的候选筛选显著降低了精确几何运算数量。

完整60工况纯计算时间约为

$$
17.54\text{ s}.
$$

阴影遮挡模型至此闭合。

---

# 15. 当前 Q1 模型状态

目前已经完成：

$$
\boxed{
\mathbf s
}
\quad\checkmark
$$

$$
\boxed{
\mathbf r_i
}
\quad\checkmark
$$

$$
\boxed{
\mathbf n_i
}
\quad\checkmark
$$

$$
\boxed{
DNI
}
\quad\checkmark
$$

$$
\boxed{
\eta_{\cos}
}
\quad\checkmark
$$

$$
\boxed{
\eta_{at}
}
\quad\checkmark
$$

$$
\boxed{
\eta_{\mathrm{ref}}
}
\quad\checkmark
$$

$$
\boxed{
\eta_{sb}
}
\quad\checkmark
$$

尚未闭合：

$$
\boxed{
\eta_{\mathrm{trunc}}
}
$$

下一阶段进入有限太阳视盘及有限尺寸集热器下的截断效率模型。