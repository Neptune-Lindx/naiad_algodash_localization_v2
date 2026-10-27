# 半边天 RTK：正确固定解

本文整理园区集装箱墙 / 半边天场景下，如何让 `rnx2rtkp` 出**正确固定解**（尤其高度）。  
只保留已在本仓库验证过的结论，和文献里对得上的做法。细节讨论见 `dev_log/20260819_2253_…`、`20260820_0005_…`、`20260907_1254_…`。

INS / 轮速 / NHC 见 [`ins-odo.md`](ins-odo.md)。C/N0 标定工具在 `src/tools/C_N0_param_calibrate/`。

---

## 1. 问题是什么

空旷 RTK 正常；半边天经常出现：

- 报了固定解（`Q=1`），平面还凑合，**高度偏几米并钉住**；
- 或 ratio 看起来不低，整数其实是错的。

这不是「高度角模型还差 0.5 dB」。主因两层：

1. **墙沿衍射 / NLOS**：相位仍能锁，观测带**几米级偏差**（非色散几何路程 $\delta\rho$），把浮点模糊度带偏。  
2. **好星挤在同一半球**：VDOP 天生差；PAR 门槛过松、`fix-and-hold` 会把错高度焊死。

手剔墙沿 11 颗（G05 G15 G24 E29 J03 C01 C06 C07 C08 C10 C22）后，高程才对得上。`rnx2rtkp` 还曾在剔星后对了 2 秒，第 3 秒 PAR 用 **ratio=1.4** 换成 $U\approx +4.6\,\mathrm{m}$（`20260819_2253`）。Qt 用经典 LAMBDA、`val thres=2.0`，不会收 1.4。

下面这张表来自**基站开阔、两端天线可对比**的那份园区数据（`20260819_2253`），说明墙沿不是低高度角：

| | 墙沿组 | 健康星 |
|--|--|--|
| 流动站 SNR | ~29（G24/C06 仍有 36～37） | ~45 |
| $\Delta\mathrm{SNR}$（rover−base） | **−7～−19 dB** | +2～+8 |
| 缺相位 | ~50% | ~18% |

G24/C06 仰角可过 15°，相对空旷 $C/N0^*$ 的 S1 常只有 2～4，过不了 `cn0-t0=8`。

**现实约束（本文按这个设计，不按上表当默认）**

1. **基站和流动站可以同时贴同一面墙。** 公共遮挡下 $\Delta\mathrm{SNR}\approx 0$，墙沿星看不出来。  
2. **两端天线可以不同**（增益、方向图、罩）。绝对载噪比之差混进天线差，不是遮挡差。  
3. 因此 **rover−base 的 $\Delta\mathrm{SNR}$ 只是「基站开阔且天线可比」时的对照实验，不能当主方案。**

同一面墙、短基线时，公共绕路 $\delta\rho$ 双差会消掉一截。要打的是**两端不一致**的那部分，以及**方位上天际线里不该进 AR 的星**。乱剔两边共视的墙沿星，几何更差。

---

## 2. 观测里什么会被组合消掉

$$
P = \rho + c\Delta t + I + T + \delta\rho + M_P + \varepsilon_P
$$

$$
\Phi = \rho + c\Delta t - I + T + \delta\rho + \lambda N + M_\Phi + \varepsilon_\Phi
$$

| 名称 | 符号 | 色散？ | GF / MW / CMC 能否看见 |
|--|--|--|--|
| 几何距离 | $\rho$ | 否 | 消掉 |
| 对流层 | $T$ | 否（非色散，**不是**几何） | 消掉 |
| 绕路路程 | $\delta\rho$（衍射/NLOS） | 近似否 | **消掉** |
| 电离层 | $I\propto 1/f^2$ | 是 | 留得住 |

墙沿高度偏差主要是 $\delta\rho$ 或错整数，不是电离层。因此 **TurboEdit / MP / CMC 不能当半边天主方案**。

DOP（尤其 VDOP）半边天会升高，但它是**整场一个数**，不能写成 $\sigma_i^2\propto\mathrm{VDOP}$。剔 NLOS 后 PDOP 往往略升，ADOP 反而更好，固定更稳（Hsu 等 3DMA RTK）。

---

## 3. 本仓库现状（2026-09）

| 模块 | 配置 / 代码 | 半边天评价 |
|--|--|--|
| 手剔墙沿星 | `pos1-exclsats` 11 颗 | 高程对得上；不能泛化 |
| 高度角权 | `varerr` $a^2+b^2/\sin^2 E$ | 挡不住墙沿 |
| S1 γ | $T_0=8$，$S1=\max(S1_r,S1_b)$ | 两端都贴墙时公共星被双倍放大；墙沿常 S1&lt;8 |
| S3 | 已算，**未进权** | 只覆盖 C/N0 跳，稳衍射无效 |
| SNR 门 | `snrmask=off` | 36 dB 墙沿过得去 |
| 周跳 | LLI / GF / MW 两历元差 / SD | 周跳用，不是 NLOS 选星 |
| 钟跳 | `repair_clkjmp` 现 `#if 1` | 无明确 ms 跳时应关（`20260729_2205`） |
| AR | PAR + `fix-and-hold`；`FindMinRatio≈1.3` | 假固定主通道 |
| `dynamics` | off | 保持 off |
| IGG3 | `stats-rkf=on` | 先验 $R$ 偏小时打不中自洽 NLOS |
| `rejgdop=30` | 只卡 SPP | RTK 固定几乎不管 VDOP |
| C/N0 模板 | 空旷拟合 Ave；BDS F2/F3 抄 F1；eq(10) θ≥30° | S1 尺度漂，才把 $T_0$ 抬到 8 |

`rnx2rtkp.cc` 本身不做预处理，清洗在 `rtkpos.cc`。

---

## 4. 方案：最优先 → 备选

原则：**先挡带偏的相位，再收紧固定，最后才精细化 $R$。**  
降权不够，墙沿要 **exclusion**。不要为压 PDOP 而留脏星。

默认工况：**两端都可能贴墙，天线可以不同。** 选星不能假设基站是干净参考，也不能比绝对载噪比。

### P0 — 方位天际线（不依赖天线、不依赖谁开阔）

墙是一条方位带，不是 15° 截止。两端贴同一面墙时，临界星挤在**同一方位扇区**。

- 天空图、粗楼块或一次静态天空扫描，定扇区；扇区内的星不进 AR（可仍进浮点）。  
- 有 3D 模型则用地理截止 / 3DMA（文献 1、2、4）。  
- 连续 $N$ 历元在扇区内再确认，避免一颗星刚擦边就被永久踢掉。

这比单历元 S1 稳，也比 $\Delta\mathrm{SNR}$ 稳：方向图不同的天线，遮挡方位仍然一致。

### P1 — 各站相对自己的空旷模板（天线分开标定）

比较对象是 **该天线的 $C/N0^*(\theta)$**，不是对方的 SNR。流动站、基站各拟合一套 $a_k,b_k$（或分 bin 模板），不要 Base/Rover 平均，也不要抄对方天线。

记 $\delta$ 或 S1 为「相对本机模板的亏缺」。短基线共用星：

| 两端相对各自模板 | 含义 | 处理 |
|--|--|--|
| 都不大 | 相对各自天线都干净 | 高度角权，可进 AR |
| **都大且接近** | 公共墙沿 | **保留双差**，不要 $\max(S1_r,S1_b)$ 再乘两次 $\gamma$ |
| **只有一端大** | 一边探出墙角 / 贴墙距离不同 | 这条 DD 不进 AR |

无相位、只有伪距：不进模糊度（与天线无关）。连续 $N$ 历元「单端大」再踢，连续 $M$ 历元对称或都干净再放回。

现状 $S1=\max(S1_r,S1_b)$ 且共用一套 Ave 系数：天线不同时两端 S1 不可比，两端都贴墙时公共星被过度降权。

### P2 — AR 纪律（选星之后仍可能假固定）

1. 半边天关掉 `fix-and-hold`，改 `continuous` / `instantaneous`。  
2. PAR 最低 ratio ≥ 2.0，或半边天走经典 LAMBDA（Qt 已验证）。`pos2-arthres` 管不到 `FindMinRatio`。  
3. 选星后 ns / VDOP（或 ADOP）不达标 → **宁浮点、不固定**。  
4. `pos1-dynamics` 保持 off。

### P3 — 模板做对（给 P1 用，本身不是剔星）

1. 空旷、**分天线**拟合；工作仰角收到 `elmask`（15°），避免 15°–30° 外推。  
2. 频点对上实测类型；BDS 无 B2a 不要抄 F1。  
3. $T_0$ 按该天线空旷残差的 $q_{95}$ 定，不要沿用 8，也不要直接抄 Alloy 的 2（模板尺度不对时 2 会误伤）。  
4. 分 bin 中位数比三次式稳（SIGMA-δ）。  
5. S3 只抓 C/N0 跳；两端都稳、只是都比各自模板低时，S3 帮不上，交给 P0/P1。

S1 假阳性（打低仰角开敞星）+ 假阴性（留墙沿）会比纯高度角更差。

### P4 — 周跳 / 码多径窗口（次要）

把 GF/MW 从「两历元差」升级为 TurboEdit 窗口（均值 + $k\sigma$）。  
MP1 / CMC 滑动 RMS **只降伪距权**。  
无双频码相则不算 MW，直接该频不进 AR。  
**不要上 Hatch** 救固定高度。这些组合消掉 $\delta\rho$，看不见公共绕路。

### P5 — 高度约束（不等价于 GNSS 自己固定对）

VDOP 差是几何事实。紧组合 INS / NHC / 轮速管短时 $U$。GNSS 仍按 P0+P1+P2。错的 GNSS 高度会把 INS 带偏。

### 条件可用，不能当默认：$\Delta\mathrm{SNR}$

仅当**同时**满足：基站确认开阔、两端天线型号与安装接近（或事先用空旷共视星把天线差标成常数）。  
此时 $\mathrm{SNR}_\mathrm{rov}-\mathrm{SNR}_\mathrm{base}\lesssim -6\sim-8\,\mathrm{dB}$ 可复现手剔 11 颗（及 E33、C39、R14、R17）。  
天线一换、或基站也贴墙，这条失效，退回 P0+P1。

### 不建议当主方案

| 做法 | 原因 |
|--|--|
| 默认用 rover−base $\Delta\mathrm{SNR}$ | 同墙 $\approx 0$；异天线混进增益差 |
| $S1=\max(S1_r,S1_b)$ 且共用一套系数 | 公共墙沿被双倍放大；天线不同时不可比 |
| 再降 $T_0$、加满 eq(12) S3 | 打不中「仰角够、C/N0 还不低」的墙沿 |
| 只开 `snrmask` | 绝对门限跨天线失效；36 dB 墙沿也过得去 |
| 按 $\sigma^2\propto\mathrm{DOP}$ 逐星加权 | 无此主流模型；DOP 是整场门控 |
| 最小 PDOP 选星 | 会留下 NLOS（文献 5） |
| GNSS 层 RANSAC | 野值可自洽（`20260825_1219`） |
| 无 ms 跳还开 `repair_clkjmp` | 掉固定率 |
| IGG3 单独补先验 | $R$ 按高度角估小时，自洽 NLOS 不野 |

---

## 5. 建议落地顺序

按「两端都可能贴墙、天线可以不同」：

1. **方位扇区**：墙的方位不进 AR。这是不看载噪比的第一刀。  
2. **分天线空旷模板**：各站 S1/$\delta$ 只和自己的 $C/N0^*$ 比。两端都大且接近则留双差；只有一端大则不进 AR。缺相位不进模糊度。去掉共用 Ave 和 $\max(S1_r,S1_b)$。  
3. 关 hold；PAR ratio ≥ 2 或改 LAMBDA。看是否还出现「对 2 秒再跳到 +几米」。  
4. 几何不够（VDOP/ADOP）就出浮点。  
5. 需要走车且 $U$ 不能漂：INS/轮速，GNSS 侧仍 1+2+3。

$\Delta\mathrm{SNR}$ 只在基站开阔且天线可比时，当作和第 1 节那张表的对照，不写进默认配置。

验证时先分清：浮点 $U$ 是否已偏（选星问题）还是一固定才跳（AR 问题）。顺手核对基站坐标和天线 ARP。

---

## 6. 参考文献（最优 → 最次）

编号与第 4 节优先级对齐：越靠前越该先读、越贴「半边天正确固定」。每条都给可打开的链接。

本仓库已验证（比外文更贴你们数据）：

- [`dev_log/20260819_2253_为什么除去这些卫星，我rtkpost_q.md`](../dev_log/20260819_2253_为什么除去这些卫星，我rtkpost_q.md) — 手剔、$\Delta\mathrm{SNR}$、PAR ratio=1.4 假固定、同一环境  
- [`dev_log/20260820_0005_G05_G15_G24_E29_J03_.md`](../dev_log/20260820_0005_G05_G15_G24_E29_J03_.md) — 衍射 / NLOS 检索词  
- [`dev_log/20260825_1219_这个用ransac剔除效果好吗，有相关论.md`](../dev_log/20260825_1219_这个用ransac剔除效果好吗，有相关论.md) — RANSAC 不宜当主方案  
- [`dev_log/20260729_2205_repair_clkjmp时钟跳变修复移植.md`](../dev_log/20260729_2205_repair_clkjmp时钟跳变修复移植.md) — 钟跳默认应关  
- [`dev_log/20260904_0012_将伽利略的和gps的合并，不算glona.md`](../dev_log/20260904_0012_将伽利略的和gps的合并，不算glona.md) — C/N0 分组与抄系数  

### 外文（P0 / P1 / P2：天际线、本机模板、固定几何）

不依赖「基站比流动站干净」，也不比两端绝对载噪比。

1. Zhang, Z., et al. (2022). A composite stochastic model considering the terrain topography for real-time GNSS monitoring in canyon environments. *Journal of Geodesy*. 方位地理截止，两端同墙时第一优先。  
   <https://doi.org/10.1007/s00190-022-01660-7>

2. Hsu, L.-T., et al. 3D mapping database-aided GNSS RTK and its assessments in urban canyons. 天际线剔 NLOS 后 PDOP 略升、ADOP 降，固定变好。  
   [PDF](https://www.polyu.edu.hk/aae/ipn-lab/us/publications/Fullpaper/3D%20Mapping%20Database-Aided%20GNSS%20RTK%20and%20Its%20Assessments%20in%20Urban%20Canyons.pdf)

3. Brunner, F. K., Hartinger, H., & Troyer, L. (1999). GPS signal diffraction modelling: the stochastic SIGMA-δ model. *Journal of Geodesy*, 73, 259–267. 相对**本天线**空旷模板的 $\Delta$C/N0，不是相对基站。  
   <https://doi.org/10.1007/s001900050242>

4. Geng, J., Xin, S., & Hsu, L.-T. (2023). 3D-mapping-aided PPP-RTK for deep urban canyons. 3DMA 能打中高 C/N0 的 NLOS。  
   <https://doi.org/10.57757/iugg23-0273>

5. Hotta, K., et al. (2025). A novel GNSS satellite and pseudorange selection method for robust positioning in urban canyons. 写明纯 PDOP 选星会留下 NLOS。  
   <https://doi.org/10.1109/CANDARW68385.2025.00060>

6. Teunissen, P. J. G. (1997). A canonical theory for short GPS baselines. Part IV: precision versus reliability. *Journal of Geodesy*, 71, 513–525. ADOP：该不该固定看模糊度几何，不是只看 PDOP。  
   <https://doi.org/10.1007/s001900050119>  
   同主题会议稿：<https://www.ion.org/publications/abstract.cfm?articleID=2811>

7. Wen, W., et al. (2022). 3D LiDAR aided GNSS NLOS mitigation for reliable GNSS-RTK positioning in urban canyons. 遮挡降低 AR 可观测；补街向几何。  
   <https://arxiv.org/abs/2212.05477>

### 外文（P3：C/N0 / 仰角权，辅助）

8. Tay, S., & Marais, J. Weighting models for GPS pseudorange observations for land transportation in urban canyons. C/N0 × 高度角 × LOS/NLOS。  
   <https://hal.science/hal-00942180/PDF/doc00016838.pdf>

9. Hartinger, H., & Brunner, F. K. (1999). Variances of GPS phase observations: the SIGMA-ε model. *GPS Solutions*, 2(4), 35–43. 按 C/N0 定相位方差。  
   <https://doi.org/10.1007/PL00012765>

10. Sairo, H., Akopian, D., & Takala, J. (2003). Weighted dilution of precision as quality measure in satellite positioning. *IEE Proc. Radar Sonar Navig.*, 150(6), 430–436. WDOP：权来自观测质量，DOP 只是选子集准则。  
    <https://doi.org/10.1049/ip-rsn:20031008>

11. Won, D. H., et al. (2012). Weighted DOP with consideration on elevation-dependent range errors of GNSS satellites. *IEEE Trans. Instrum. Meas.*  
    <https://doi.org/10.1109/TIM.2012.2205512>

12. Park, C. S., Kim, I., Lee, J. G., & Jee, G.-I. (1996). A satellite selection criterion incorporating the effect of elevation angle in GPS positioning. *Control Engineering Practice*, 4(12), 1741–1746.  
    <https://doi.org/10.1016/S0967-0661(96)00192-X>

### 外文（P4：周跳 / 码多径；消几何，对墙沿 $\delta\rho$ 弱）

13. Blewitt, G. (1990). An automatic editing algorithm for GPS data. *Geophysical Research Letters*, 17(3), 199–202. TurboEdit（MW + GF 窗口）。  
    <https://doi.org/10.1029/GL017i003p00199>  
    开放 PDF：<https://nbmg.unr.edu/Staff/pdfs/Blewitt_GL017i003p00199.pdf>

14. Liu, Z. (2011). A new automated cycle slip detection and repair method for a single dual-frequency GPS receiver. *Journal of Geodesy*, 85, 171–183.  
    <https://doi.org/10.1007/s00190-010-0426-y>  
    作者页 PDF：<http://www.lsgi.polyu.edu.hk/staff/george.liu/Reprint_Liu2010(A%20new%20automated%20cycle%20slip%20detection%20and%20repair%20method%20for%20a%20single%20dual-frequency%20GPS%20receiver).pdf>

15. Estey, L. H., & Meertens, C. M. (1999). TEQC: The multi-purpose toolkit for GPS/GLONASS data. *GPS Solutions*, 3(1), 42–49. MP1/MP2。  
    <https://doi.org/10.1007/PL00012778>

16. Melbourne, W. G. (1985). The case for ranging in GPS-based geodetic systems. *Proc. 1st Int. Symp. on Precise Positioning with GPS*, Rockville.  
    Wübbena, G. (1985). Software developments for geodetic positioning with GPS using TI 4100 code and carrier measurements. 同上。MW 组合定义（会议录，无统一 DOI）。  
    检索：<https://scholar.google.com/scholar?q=Melbourne+1985+The+case+for+ranging+in+GPS-based+geodetic+systems>  
    <https://scholar.google.com/scholar?q=W%C3%BCbbena+1985+Software+developments+for+geodetic+positioning+TI+4100>

### 外文（最次：本场景不要当主方案）

17. Hatch, R. (1982). The synergism of GPS code and carrier measurements. *Proc. 3rd Int. Symp. on Satellite Doppler Positioning*, Las Cruces, 1213–1231. 码平滑；半边天不建议靠它救固定高度。  
    <https://scholar.google.com/scholar?q=Hatch+1982+The+synergism+of+GPS+code+and+carrier+measurements>

18. 武大等. A realistic stochastic model considering the multipath effects based on C/N0 for GNSS positioning in harsh environments. *IEEE Trans. Instrum. Meas.*（UM-ELE 五段 $\gamma$，S1/S3）。未公开 $a_k,b_k$；Alloy $T_0=2,\,T_1=6$。**解决不了沿墙几何塌和 PAR 假固定**，本仓库已对过。  
    <https://doi.org/10.1109/TIM.2026.3655902>

