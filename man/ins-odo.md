# 轮式里程计与运动约束

本文说明本工程 INS / GNSS 组合里三类陆地车辆约束：

1. **轮速里程计**（`src/ins-gnss/ins-odo.cc`）
2. **非完整性约束 NHC**（`src/ins-gnss/ins-nhc.cc`）
3. **零速 / 零角速度**（`ins-zvu.cc`、`ins-zaru.cc`、`ins-static-detect.cc`）

视觉里程计（VO）是另一条链路，本文不覆盖。

船体侧壁吸附清洁机器人（相对船体惯导 / 轮速 / NHC）见 [`ins-hull-rel.md`](ins-hull-rel.md)，**不能**直接使用本文的 ECEF `odo()` / `nhc()` / `zvu()`。

参考文献（源文件头注释）：Groves *Principles of GNSS, Inertial, and Multisensor Integrated Navigation*；Shin *Accuracy Improvement of Low Cost INS/GPS for Land Applications*。

---

## 1. 在组合导航里起什么作用

IMU 机械编排给出 ECEF 位置 $r^e$、速度 $v^e$、姿态 $C_b^e$。GNSS 更新在遮挡、周跳、固定失败时变弱，速度和航向会漂。

陆地车体近似满足：

$$
v_b = C_e^b v^e
\approx
\begin{bmatrix} v_{\mathrm{fwd}} \\ 0 \\ 0 \end{bmatrix}
$$

也就是 **只能沿车头滚、不能侧向滑、不能离地**。工程上拆成两路虚测量，再加停车约束：

| 约束 | 代码 | 测什么 | 主要防什么 |
|--|--|--|--|
| 轮速 | `odo()` | $v_{b,x}$（或 $v^r_x$）等于轮速 | GNSS 中断后前向速度漂 |
| NHC | `nhc()` | $v_{b,y}\approx 0$，$v_{b,z}\approx 0$ | 侧向滑移式漂移、高度乱跳 |
| ZVU | `zvu()` | 静止时 $v^e\approx 0$ | 停车位置/速度随机游走 |
| ZARU | `zaru()` | 静止时 $\omega\approx 0$ | 停车陀螺零偏漂 |

**没有**单独的车轮打滑状态（不估滑移率、不比较 IMU 加速度与轮速导数）。打滑相关处理是：**大转弯时不用 odo/NHC**，以及 **侧向已经很大时不再强行 NHC**。详见第 8 节。

`carvig-tc.conf` 里默认 `ins-nhc=0`、`ins-zvu=0`、`ins-zaru=0`，代码有、当前配置关着。

---

## 2. 坐标系与安装参数

| 符号（代码） | 含义 |
|--|--|
| $b$ | IMU 本体坐标系（车辆：$x$ 前，$y$ 右，$z$ 下，取决于 `ins-imucoors`） |
| $e$ | ECEF |
| $r$ | 后轴/轮速坐标系（默认与 $b$ 重合） |
| $C_b^e$（`ins->Cbe`） | $b\to e$ |
| $C_b^r$（`ins->Cbr`） | $b\to r$。`initodo()` 初始化为单位阵 |
| $l_{br}$（`ins->rbl`，`odopt.lever`） | 轮速原点相对 IMU 的杆臂，**在 $b$ 系** |
| $s$（`ins->os`，`odopt.s`） | 比例因子（轮胎半径误差、充气、磨损） |

预测轮速时（`ve2vr`）把 INS 的 ECEF 速度变到轮速点、再变到 $r$ 系：

$$
v_b^{\mathrm{odo}} = C_e^b\, v^e + \bigl[\omega_{ib}^b - C_e^b\,\omega_{ie}^e\bigr]\times l_{br}
$$

$$
\hat{v}^r = (1+s)\,C_b^r\, v_b^{\mathrm{odo}}
$$

其中 $\omega_{ib}^b$ 用当前 IMU 陀螺（`imu->gyro`），地球自转 $\omega_{ie}^e = (0,0,\Omega)^{\mathrm{T}}$。杠杆臂项把 IMU 处速度换到轮心；转弯时这项不可忽略。

NHC **不用**杠杆臂，直接用 IMU 点速度 $v_b=C_e^b v^e$。

---

## 3. 原始脉冲如何变成速度

GI310 一类 IMU 报文里带 16 位脉冲计数（`rcv/imu.cc`：`decode_odo_data`）。相邻历元差 $\Delta c$ 先按溢出折回，再：

$$
\Delta r = \frac{\Delta c}{N_{\mathrm{res}}}\,\pi D
$$

- $N_{\mathrm{res}}$：`ins-odo_res`（默认 2048）
- $D$：`ins-odo_d`（默认 0.73 m，轮径）

结果写入 `odod_t`：`dt`、`dr`。IMU 初始化时 `initlc()` 会调用 `odores()` / `odod()` 把分辨率和轮径灌进解码器。

滤波并不对每个 IMU 拍都更新。`odomeas()` 把 $\Delta t,\Delta r$ 累加，直到跨度超过 `ins-odo_odt`（默认 1 s，未设则 0.5 s）才输出：

$$
v_x^r = \frac{\sum\Delta r}{\sum\Delta t}
$$

默认 `ins-odo_all=0`，**只用这一个前向速度**；`all=1` 时用 `vr[0..2]` 三维（侧向/天向一般不可信，需轮速能给出这些分量）。

`ins-odo_dir` 在配置里存在，**当前 `ins-odo.cc` 未使用**。

---

## 4. 里程计量测与 EKF

误差状态是闭环 EKF：预测后 $x=0$，量测更新得到 $\delta x$，再反馈到 $C_b^e,v^e,r^e$ 以及 $s,C_b^r,l_{br}$。

### 4.1 新息

$$
v = v_{\mathrm{odo}}^r - \hat{v}^r
$$

单轴时 $v = v_x^{\mathrm{odo}} - \hat{v}_x^r$。  
$|v|>30\,\mathrm{m/s}$ 只打日志，**不拒绝更新**（`MAXINOV`）。

另有硬门槛：$\lVert\omega\rVert > 30^\circ/\mathrm{s}$ 时整次更新放弃（大转弯轮速易滑移）。

### 4.2 观测噪声

`odoHVR()` 里对角 $R = (0.1\,\mathrm{m/s})^2$（宏 `VARVEL`）。  
配置项 `ins-odo_ostd` **没有接到 $R$**，改 conf 目前不会改变量测标准差。

### 4.3 雅可比（量测对误差状态）

$H$ 只填打开了的状态列（见第 5 节）：

| 状态 | 函数 | 作用 |
|--|--|--|
| INS 速度 $\delta v^e$ | `jacobian_ov_dv` | $\partial\hat{v}^r/\partial v^e \propto C_b^r C_e^b$ |
| 姿态 $\delta\phi$ | `jacobian_ov_da` | 速度投影 + 杠杆臂随姿态变化 |
| 陀螺零偏 $\delta b_g$ | `jacobian_ov_bg` | $\omega\times l$ |
| 比例因子 $\delta s$ | `jacobian_ov_ds` | $\approx -\hat{v}^r$ |
| 安装角 $\delta\alpha_{br}$ | `jacobian_ov_dma` | $C_b^r$ 小角 |
| 杆臂 $\delta l$ | `jacobian_ov_lever` | 角速度叉乘 |

`UPD_IN_EULER` 打开时，姿态列再乘扰动欧拉角雅可比 `jacobian_prot_pang`。

标准卡尔曼：`filter(x, ins->P, H, v, R, nx, nv)`。成功则 `ins->stat = INSS_ODO`，`odoclp()` 反馈：

- $s \leftarrow s - \delta s$
- $C_b^r$ 用 `corratt` 修正
- $l_{br} \leftarrow l_{br} + \delta l$
- 位置/速度/姿态/IMU 零偏走公共 `clp()`（`ins-nhc.cc`）

---

## 5. 哪些状态会进滤波器

里程计相关维数由开关决定（`ins-gnss-state.cc`），排在 GNSS 天线杆臂之后：

| 开关 | 维数 | 初值协方差 |
|--|--|--|
| `ins-estodos` | 1（$s$） | `ins-uncos`，缺省 $0.01^2$ |
| `ins-estodol` | 3（杆臂） | `ins-uncol` |
| `ins-estodoa` | 3（安装角） | `ins-uncoa` |

全关时，量测仍更新 **速度、姿态、陀螺零偏**（后两者须 `ins-bgopt` 等已打开）。这是最常见用法：把轮速当速度观测，不估轮胎半径。

时间更新里 `stochasticPhi` 按 `ins-osproopt` / `olproopt` / `oaproopt` 处理这几维：

- `1` 随机游走、`2` 随机常数：$\Phi_{ii}=1$
- `3` 高斯–马尔可夫：$\Phi_{ii}=\mathrm{e}^{-|\Delta t|/\tau}$

`getQ()` **没有**给 $s$、杆臂、安装角加过程噪声，只有陀螺/加计/零偏/钟差。这几维主要靠 $\Phi$ 和量测，而不是 $Q$ 驱动。

NHC / ZVU / ZARU **不增加状态**，只对已有速度、姿态（ZARU 还对陀螺零偏）做虚测量。

---

## 6. 调用时机（和 GNSS 的关系）

运动约束挂在 **INS 机械编排 + GNSS 更新之后**，不是 GNSS 观测方程里的一项。

**事后紧组合** `tc-postpos`（`tcpostpos.cc`）：

1. `tcigpos(..., INSUPD_MEAS 或 TIME)`：传播 +（有观测时）GNSS 紧组合
2. `ins-nhc=1` 时按节拍 `nhc()`
3. `insopt.odo != 0` 时 `odo(...)`
4. 静止检测：IMU 检测 **或** `detstatic_ODO`
5. `zvu` / `zaru`

**实时** `carvigsvr.cc`：`lcigpos`/`tcigpos` → `motion()`（NHC + ZVU/ZARU）→ `odo()`。

纯 GNSS 程序 `rnx2rtkp` 会链入这些 `.cc`，但 GNSS-only 模式不会走这条更新。

门控变量是 `insopt.odo`。当前 `options.cc` **没有** `ins-odo` 这一项，conf 里的 `ins-odo_*` 只配分辨率、轮径、杆臂等。要用轮速融合，需在代码里置 `insopt.odo=1`，或自行把开关加进 sysopts。NHC / ZVU / ZARU 的开关已经在 conf 里。

---

## 7. 非完整性约束（NHC）

源码 `ins-nhc.cc`。陆地车轮与地面点接触、无侧向/法向自由度时：

$$
v_{b,y} = 0,\qquad v_{b,z} = 0
$$

这是两条 **虚速度观测**，不是传感器。新息就是当前 INS 投影出来的侧向/天向速度：

$$
v_b = C_e^b v^e,\qquad
\tilde{z}
=
\begin{bmatrix} v_{b,y} \\ v_{b,z} \end{bmatrix}
-
\begin{bmatrix} 0 \\ 0 \end{bmatrix}
$$

对误差状态的雅可比（`bldnhc`）：

$$
\delta v_b
\approx
C_e^b\,\delta v^e
+
C_e^b [v^e\times]\,\delta\phi
$$

只取 $y,z$ 两行写入 $H$ 的速度列和姿态列。观测噪声对角：

$$
R = \mathrm{diag}(0.05^2,\,0.05^2)\quad (\mathrm{m/s})^2
$$

滤波成功则 `ins->stat = INSS_NHC`，同样走 `clp()` 反馈速度/姿态。

### 7.1 什么时候不用 NHC

对 $y$、$z$ 分别判断，不满足就 **跳过这一轴**（另一轴仍可用）：

$$
|v_{b,i}| > 0.5\,\mathrm{m/s}
\quad\text{或}\quad
\lVert\omega\rVert > 30^\circ/\mathrm{s}
$$

侧向已经大于 $0.5\,\mathrm{m/s}$ 时再把测量当成 0，会把真实侧滑硬拉成直线，所以放弃。大转弯同样放弃（见第 8 节）。

### 7.2 更新节拍

`motion()` / `tcpostpos`：

```text
if (nhc && (++nc > nhz)) { nc = 0; nhc(ins, opt, imu); }
```

`ins-nhz` 注释写成“更新频率”，代码里是 **IMU 计数阈值**。`carvig-tc.conf` 中 `ins-hz=100`、`ins-nhz=100` 时，大约每 101 个 IMU 拍做一次 NHC（约 1 Hz），不是 100 Hz。

---

## 8. 打滑：有防护，没有滑移估计

工程里 **没有** 滑移率状态、没有“IMU 加速度 vs 轮速微分”一致性检验、没有左右轮差速打滑检测。和打滑/漂移有关的只有下面几道门。

### 8.1 大转弯不用轮速、不用 NHC

轮速模型按纯滚。横摆太大时轮胎更容易滑，量测会把 INS 速度拉错。两边用同一门槛：

$$
\lVert\omega_{ib}^b\rVert > 30^\circ/\mathrm{s}
\quad\Rightarrow\quad
\text{本次 } \mathtt{odo}/\mathtt{nhc}\text{ 整次或该轴放弃}
$$

- `odofilt()`：整次 `return 0`
- `bldnhc()`：该轴 `continue`

这是 **回避** 滑移，不是估计滑移。直线急加速抱死、低附侧滑但横摆不大时，轮速仍会进滤波。

### 8.2 NHC 侧向过大则不约束

$|v_{b,y}|$ 或 $|v_{b,z}|$ 已超过 $0.5\,\mathrm{m/s}$ 时不再施加 NHC。含义是：滤波器已经看到明显侧向速度，宁可让它漂，也不用错误的“必须为 0”把轨迹掰直。

### 8.3 比例因子 $s$ 代不了打滑

$s$ 吸收的是缓慢的轮胎半径/充气误差。突发打滑是有偏速度，不是常数比例。打开 `ins-estodos` 也不会当滑移观测器。

### 8.4 停车才是真正的“防漂移”

GNSS 没了、车停着，INS 速度会随机游走、位置线性漂。这时靠 **判静 + ZVU**（第 9 节），比轮速更有效：轮速为 0 只能说明轮子没转，车可能仍在滑（冰面），所以 ZVU 还要求 IMU 窗口也像静止。

---

## 9. 零速 / 零角速度（防停车漂移）

### 9.1 先判静，再更新

`detstc()` 按 `ins-detst` 选检测器（`ins-static-detect.cc`），窗口长度 `ins-zv_ws`（默认 5 拍）：

| `ins-detst` | 名称 | 大意 |
|--|--|--|
| 0 / 1 | GLRT | 加计相对重力方向 + 陀螺能量，Skog 广义似然 |
| 2 | MV | 加计滑动方差 |
| 3 | MAG | $\|f\|$ 与 $\|g\|$ 之差 |
| 4 | ARE | 陀螺能量 |
| 5 | ALL | 以上同时成立 |

GLRT 统计量（$\sigma_a=$ `ins-zv_siga`，$\sigma_g=$ `ins-zv_sigg`）：

$$
T=\frac{1}{N}\sum_{k=1}^{N}
\left(
\frac{\lVert\omega_k\rVert^2}{\sigma_g^2}
+
\frac{\lVert f_k - \frac{\|g\|}{\|\bar f\|}\bar f\rVert^2}{\sigma_a^2}
\right)
$$

$T<\gamma$（`ins-zv_gamma1` 等）判静止。

里程计另开一路 `detstatic_ODO`：窗口 `ins-zv_odt`（默认 0.1 s）内

$$
\frac{\sum\Delta r}{\sum\Delta t}=0
\quad\Rightarrow\quad \text{轮速认为静止}
$$

与 IMU 检测 **或** 运算：`zf = detstc(...) | detstatic_ODO(...)`。  
轮子没转但车在滑时，ODO 会误判静；IMU 检测可挡一部分。两者都开更稳。

`detstatic_ODO` **只改标志，不更新状态**。

### 9.2 ZVU：速度拉回 0

`ins-zvu.cc`，需 `ins-zvu=1` 且 `zf=1`。再加两道门槛，且每累计 `MINZC=15` 次静止才真正滤一次：

$$
\tilde{z} = v^e - 0,\qquad
H = \bigl[0_{3\times 3}\ \ I_3\ \ \cdots\bigr]
$$

$$
R = (0.05\,\mathrm{m/s})^2 I_3
$$

仅当

$$
\lVert v^e\rVert < 0.1\,\mathrm{m/s}
\quad\text{且}\quad
\lVert\omega\rVert < 10^\circ/\mathrm{s}
$$

才更新。成功则 `INSS_ZVU`。这是停车位置不乱爬的主手段。

### 9.3 ZARU：陀螺当 0，估零偏

`ins-zaru.cc`，需 `ins-zaru=1`、`ins-bgopt` 在估零偏、且 `zf=1`。每 `MINZAC=100` 次静止滤一次：

$$
\tilde{z} = 0 - \omega_{\mathrm{meas}},\qquad
H \text{ 对准 } \delta b_g
$$

$$
R = (1^\circ/\mathrm{s})^2 I_3
$$

仅当 $\lVert\omega\rVert<5^\circ/\mathrm{s}$ 且 $\lVert v^e\rVert<0.1\,\mathrm{m/s}$。成功则 `INSS_ZARU`。停车给陀螺零偏可观性，减轻之后再走时的航向漂。

---

## 10. 配置项对照

| conf | 代码 | 说明 |
|--|--|--|
| （缺开关） | `insopt.odo` | 是否做轮速更新 |
| `ins-odo_res` / `ins-odo_d` | `odopt.res` / `d` | 脉冲→弧长 |
| `ins-odo_odt` | `odopt.odt` | 速度平滑窗口（s） |
| `ins-odo_s` | `odopt.s` → `ins->os` | 比例因子初值 |
| `ins-odo_lever1..3` | `odopt.lever` | $b$ 系杆臂（m） |
| `ins-odo_all` | `odopt.all` | 0：仅前向；1：三维 |
| `ins-odo_ostd` | `odopt.ostd` | **未接入 $R$** |
| `ins-odo_dir` | `odopt.dir` | **未使用** |
| `ins-estodos/a/l` | `estodos/a/l` | 是否估计 $s$ / 安装角 / 杆臂 |
| `ins-uncos/oa/ol` | `unc.os/oa/ol` | 上述状态初值标准差 |
| `ins-osproopt` 等 | `osproopt` 等 | $1$ RW / $2$ 常数 / $3$ GM |
| `ins-nhc` | `nhc` | 非完整性约束开关（当前 conf 默认 0） |
| `ins-nhz` | `nhz` | NHC 间隔：每 `nhz+1` 个 IMU 拍一次 |
| `ins-zvu` / `ins-zaru` | `zvu` / `zaru` | 零速 / 零角速度（当前 conf 默认 0） |
| `ins-detst` | `detst` | 静止检测器类型（0/1=GLRT） |
| `ins-zv_ws` | `zvopt.ws` | IMU 判静窗口 |
| `ins-zv_odt` | `zvopt.odt` | 里程计判静窗口 |
| `ins-zv_siga/sigg/gamma*` | `zvopt` | GLRT 等阈值 |

建议陆地车同时开：`nhc` 压侧向，`odo` 压前向，停车再开 `zvu`（有轮速则 ODO 判静作或）。不要指望 NHC/odo 处理大侧滑。

---

## 11. 适用与限制

**适合**

- 轮式车辆、轮速与 IMU 已时间对齐（odo 塞在 `imud_t` 同拍）
- GNSS 中断或半空时抑制速度发散（前向靠 odo，侧向靠 NHC）
- 停车：判静 + ZVU 防止位置爬走；ZARU 稳住陀螺零偏

**限制**

- **无滑移状态**：直线打滑、冰面侧滑、轮速虚高，只要横摆 $<30^\circ/\mathrm{s}$，轮速仍会进滤波
- NHC 在已有明显侧向速度时主动放弃，不能把侧滑“纠正”回来
- 轮速 $R$ 固定 $0.1\,\mathrm{m/s}$；NHC $R$ 固定 $0.05\,\mathrm{m/s}$
- `ins-odo_ostd`、`ins-odo_dir` 未接入；`insopt.odo` 使能开关未进 conf
- `carvig-tc.conf` 默认 NHC/ZVU/ZARU 全关，开组合导航后需自己打开
- 与 GNSS 是 **序贯更新**（先 GNSS 再 NHC/odo/ZVU），不是同一历元联合 $H$

**和 GNSS 半空问题的关系**：这些约束的是 **速度 / 停车**，不能修正墙沿 NLOS 伪距相位偏差。GNSS 仍差时，主要让速度合理、停车不漂；水平位置仍可能被错误固定解拉偏。
