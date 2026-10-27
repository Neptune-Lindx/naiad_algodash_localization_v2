# naiad_localization_v2 框架说明

> 算法核拷自 `/home/neptune-lindx/src/C/gnss_imu_couple` 的 `src/ins-gnss`。
> 船体相对融合的公式以同目录 [`ins-hull-rel.md`](ins-hull-rel.md) 为准。
> 陆地里程计对照见 [`ins-odo.md`](ins-odo.md)。

## 1. 这个包要做成什么

第一版只做大船坐标系下的 IMU 与轮速融合，对应 `ins-hull-rel.md` 的模式 A。节点不订阅 `sensor_rtk/data` 和 `sensor_xcar/data`。紧组合和壁虎模式 B 留到后续版本。

对准之后的目标仍然是：

| 场景 | 算法 | 观测 |
|--|--|--|
| 水面，RTK 固定 | RTK 与 IMU 紧组合 | 载波/伪距 + 机器人 IMU。本版未接 |
| 水下，没有 RTK | 大船坐标系下 IMU 与轮速融合 | 机器人 IMU + 轮速。本版走模式 A |
| 刚出水面或水面失锁 | 保持惯导，快速搜星，固定后回到紧组合 | 本版未接 |

## 2. 数据怎么走

订阅与 `naiad_localization` 的话题输入对齐，QoS 也对齐：除 `robot/profile` 用 `QoS(1).transient_local()` 外，队列深度都是 10。节点启动时把解析后的话题名打到日志，格式是 `subscribe <topic> (<param>)`。

| 话题 | 消息 | 参数名 | 缓存给谁 |
|--|--|--|--|
| `motion/data` | `MotionData` | `motion_topic` | 轮速、转向反馈 |
| `motion/odometer` | `MotionOdometer` | `wheel_topic` | 轮速增量，准静态和船体融合 |
| `motor/info` | `DeviceInfo` | `motor_info_topic` | 电机设备信息 |
| `fog/data` | `FogData` | `fog_topic` | 机器人 IMU |
| `fog/state` | `FogState` | `fog_state_topic` | FOG 对准与故障状态 |
| `sensor_depth/data` | `DepthData` | `depth_topic` | 压力，水面 / 水下 |
| `ais/data` | `AisData` | `ais_topic` | 大船 AIS 航向 |
| `robot/profile` | `std_msgs/String` | `profile_topic` | 机器人剖面 |
| `task/parameter` | `TaskParameter` | `task_parameter_topic` | 作业参数，含作业舷 |
| `pumpbox/state` | `PumpboxState` | `pumpbox_topic` | 泵盒状态 |
| `pushbox/state` | `PushboxState` | `pushbox_topic` | 顶出盒状态 |

`sensor_rtk/data` 和 `sensor_xcar/data` 本版不订阅。

```text
motion/data  motion/odometer  motor/info
fog/data  fog/state  sensor_depth/data
ais/data  robot/profile  task/parameter
pumpbox/state  pushbox/state
        │
        ▼
NodeLocalizationV2
        │
        ├─ DepthGate        压力滞回，得到 submerged
        ├─ StaticGate       水面准静态保持，锁存 align_ready
        └─ HullFusion       模式 A：IMU + 轮速，大船坐标系
```

节点名 `naiad_localization_v2`。参数在 `config/localization_v2.yaml`，启动文件 `launch/localization_v2.launch.py`。

输出：

| 话题 | 消息 | 何时发 |
|--|--|--|
| `localization_v2/status` | `std_msgs/String` | 每拍。字段含 backend、filter_wired、submerged、align_ready |
| `localization_v2/odometry` | `naiad_interfaces/RobotOdometry` | 仅 `solution_valid=1`。坐标系 `ship`，子坐标系 `imu` |

不往 `robot/odometry` 发。现有 `naiad_localization` 仍占用那条话题；本包的滤波接上之前，把空解发到业务话题会把下游带偏。

没有有效深度样本时按水面处理。压力无效时 `DepthGate` 保持上一状态。本版水面和水下都走船体模式 A。

准静态对准看两件事同时成立并保持 `align_hold_sec`（默认 30 s）：三轮距离增量绝对值之和不超过 `wheel_static_m`，陀螺模长不超过 `gyro_static`。陀螺阈值与 `fog/data.gyroscope` 同一单位。这只锁存“可以开始对准”，对准方程仍在 `ins-align.cc` / `ins-init-rt.cc`，还没有从桥接层调用。

## 3. 四个相位

`ModeManager` 的转移：

1. **align**。起点。水下则转入 underwater_hull（对准未完成也会转，船体融合侧以后标成降级）。水面且对准保持已满足、RTK 已固定，转入 surface_tight。对准保持已满足但还没有固定解，转入 resurface_search。
2. **surface_tight**。下潜转入 underwater_hull。RTK 丢失但仍在水面，转入 resurface_search，不改走船体融合。
3. **underwater_hull**。压力回到水面阈值以下，转入 resurface_search，并调用 `arm_fast_search()`。
4. **resurface_search**。再次下潜回到 underwater_hull。对准保持仍有效且 RTK 重新固定，回到 surface_tight。

`arm_fast_search()` 目前只置位。接到 `ins-gnss-tc.cc` 时，这个标志表示保留星历和模糊度，走重新搜星，而不是把滤波器清零。

## 4. 两套算法分别放哪

紧组合继续用拷贝来的库，不改陆地接口的含义：

| 步骤 | 现有文件 | 以后由谁调用 |
|--|--|--|
| 机械编排、圆锥/划桨 | `src/ins-gnss/ins.cc` 的 `updateins` | 紧组合与船体融合都用机器人 IMU 这一路 |
| RTK 位置 | `src/ins-gnss/rtkpos.cc` | 紧组合 |
| 紧组合量测更新 | `src/ins-gnss/ins-gnss-tc.cc` | `GnssTcBridge::step_tight` |
| 序贯顺序：机械编排 → GNSS → 运动约束 | `src/ins-gnss/tcpostpos.cc` | 紧组合主循环参照这个顺序 |
| 初始对准 | `src/ins-gnss/ins-align.cc`、`ins-init-rt.cc` | `GnssTcBridge::step_align` |
| 紧组合配置 | `config/carvig-tc.conf`（`pos1-posmode = ins-tightly-coupled`） | 读入选项，不再走文件后处理主程序 |

船体融合是新模块，公式在 [`ins-hull-rel.md`](ins-hull-rel.md)。实现时加独立分支，不要改 `cpt_tc` 那条陆地路径。

| 陆地函数 | 文件 | 船上为什么不能直接开 |
|--|--|--|
| `odo()` / `ve2vr` | `ins-odo.cc` | 把 ECEF 速度往轮速上拉。轮子停转时对地速度仍是船速 |
| `nhc()` | `ins-nhc.cc` | 约束的是对地机体系侧向/法向速度 |
| `zvu()` | `ins-zvu.cc` | 新息是 `ins->ve`，吸附停车时对地速度不为 0 |
| `zaru()` | `ins-zaru.cc` | 约束的是对地角速度，不是相对船体角速度 |

本版 `HullFusion` 固定为模式 A（只有艏向，`ω^s = 0`，NHC 噪声要放宽）。模式 B 依赖壁虎，等 `sensor_xcar/data` 接回来再用。定义见 `ins-hull-rel.md` 第 3 节。

## 5. 和 gnss_imu_couple 的差别

`gnss_imu_couple` 是 CMake 可执行文件 `tc-postpos` / `rnx2rtkp`，从 RINEX、IMU 文件读入，把解写到 `bin/*.pos`。

本包把同一套源码放在 `src/ins-gnss/`，输入输出改成 ROS2 话题。`src/ins-gnss/app/` 里的 `main` 保留作后处理对照，不链进节点。

默认 **不编译** `ins-gnss`（`NAIAD_BUILD_INS_GNSS=OFF`）。打开后会编静态库 `ins_gnss`，需要本机已有的 pthread，以及这套源码原本依赖的 lapack/blas。节点在桥接层真正调用 `ins-gnss-tc` / 船体更新之前不链接它。

没有拷入本包的内容：

| 原路径 | 原因 |
|--|--|
| `src/relative-pose`、`src/libviso`、`src/msckf-python` | 视觉里程计，不是水面紧组合或水下轮速融合 |
| `include/relative-pose.h` | 同上 |
| `src/tools`、`plot` | 离线画图 |
| `src/ins-gnss/ins-vo.zip` | 打包的视觉代码 |
| `rcv/io_png.*`、`rcv/jpeg.cc`、`rcv/m39-mix.cc` | 图像解码，原 CMake 也没有编进 `tc-postpos` |
| `man/theory.pdf` | RTKLIB 理论手册，约 860 KB，仍在 `gnss_imu_couple/man/theory.pdf` |
| `build/`、`bin/`、`.git` | 构建产物和原仓库元数据 |

`src/ins-gnss` 里仍留着相机、PnP、光流等文件。它们和紧组合编在同一个库里，删掉可能缺符号，所以一并拷了。水面/水下主路径不调用它们。

## 6. 现在已经接通的，和还没有接通的

已经接通：

- ROS2 订阅下面 11 个话题。不订阅 `sensor_rtk/data`、`sensor_xcar/data`
- 水面/水下压力滞回、准静态保持
- 船体模式 A 入口（IMU + 轮速，不算状态）
- `localization_v2/status`

还没有接通，所以 `filter_wired=0`，并且不发里程：

- `ins-align.cc` 的准静态对准
- `ins-gnss-tc.cc` 的紧组合更新
- `ins-hull-rel.md` 里的相对运动学、轮速观测、NHC、相对零速
- 把 `NavStatus.p_s` / `v_s` 填成船体坐标系的位置和速度

下一步是在 `GnssTcBridge` 和 `HullFusion` 里调用上述源码，置 `filter_wired` 和 `solution_valid`。置位之后节点才会发布 `localization_v2/odometry`。

## 7. 目录

```text
naiad_algodash_localization_v2/
  package.xml                  ROS2 包清单，包名 naiad_localization_v2
  CMakeLists.txt               编 ROS2 节点；可选编 ins_gnss
  README.md                    入口说明
  rules/ws_common.cmake        与 naiad_localization 相同的工作空间路径规则
  tools/package_version.py     从 package.xml 读版本号
  launch/localization_v2.launch.py
  config/localization_v2.yaml  节点参数
  config/carvig-tc.conf        紧组合选项，拷自 app/carvig-tc.conf
  data/geomag/                 IGRF/WMM 系数，geomag.cc 使用
  man/                         本说明、船体设计、陆地里程计对照、RTK 笔记
  include/localization_v2/     本包新接口
  include/ins-gnss/            carvig 公共头，编译 ins-gnss 时加入 include 路径
  src/ros/                     节点与 main
  src/mode/                    深度滞回、准静态保持、相位机
  src/bridge/                  紧组合调用边界
  src/hull/                    船体融合调用边界
  src/ins-gnss/                拷贝的算法实现
```

## 8. 每个文件存放什么

下面按职责列出。`src/ins-gnss` 与 `include/ins-gnss` 的内容来自 `gnss_imu_couple`，职责根据文件头注释和原 `CMakeLists.txt` 的 `tc-postpos` 源文件列表整理。

### 8.1 本包新增

| 文件 | 存放什么 |
|--|--|
| `package.xml` | ROS2 包名、版本、对 rclcpp / naiad_interfaces 的依赖 |
| `CMakeLists.txt` | 节点 `naiad_localization_v2` 与核心库的编译安装规则 |
| `README.md` | 如何启动、文档入口 |
| `.gitignore` | 忽略 build、install、log 和后处理解文件 |
| `rules/ws_common.cmake` | 解析 ROS 发行版、工作空间安装前缀、应用版本宏 |
| `tools/package_version.py` | 读取 `package.xml` 的 `<version>` |
| `launch/localization_v2.launch.py` | 带 yaml 参数启动节点 |
| `config/localization_v2.yaml` | 调度周期、潜水压力阈值、对准保持、话题名、RTK 固定类型 |
| `config/carvig-tc.conf` | 紧组合运行选项：频点、模糊度、输出、INS 噪声 |
| `include/localization_v2/work_mode.hpp` | 四个工作相位枚举和名字 |
| `include/localization_v2/nav_status.hpp` | 算法入口返回值。`filter_wired` 为假时禁止当导航解发布 |
| `include/localization_v2/depth_gate.hpp` | 压力滞回接口 |
| `include/localization_v2/static_gate.hpp` | 准静态保持接口 |
| `include/localization_v2/mode_manager.hpp` | 相位机输入和接口 |
| `include/localization_v2/gnss_tc_bridge.hpp` | 对准、滑行、紧组合、快速搜星标志 |
| `include/localization_v2/hull_fusion.hpp` | 船体模式 A/B/auto 和 `step` |
| `include/localization_v2/node_localization_v2.hpp` | ROS2 节点成员：订阅、缓存、参数 |
| `src/mode/depth_gate.cpp` | 高于 `dive_on_bar` 进水下，低于 `dive_off_bar` 出水面 |
| `src/mode/static_gate.cpp` | 连续静止达到保持时间后锁存对准就绪 |
| `src/mode/mode_manager.cpp` | 第 3 节的相位转移 |
| `src/bridge/gnss_tc_bridge.cpp` | 紧组合三个入口的空实现，返回 `filter_wired=0` |
| `src/hull/hull_fusion.cpp` | 按壁虎是否新鲜选择 `hull_a` 或 `hull_b`，不算状态 |
| `src/ros/node_localization_v2.cpp` | 订阅第 2 节的 11 个话题，调度船体模式 A，发 status；有解才发里程 |
| `src/ros/main.cpp` | `rclcpp::init` 后 spin 节点 |
| `src/ins-gnss/CMakeLists.txt` | 可选静态库源文件列表，不含 `app/` |
| `man/framework.md` | 本说明 |
| `man/ins-hull-rel.md` | 船体相对系运动学、模式 A/B、轮速与 NHC 设计 |
| `man/ins-odo.md` | 陆地 `odo()` 现状，说明它为什么不能直接用于船体 |
| `man/rtk-halfsky-fix.md` | 半边天场景下的 RTK 固定笔记，出水搜星时参考 |
| `man/rtk-qt-par-minsat.md` | 部分模糊度与最少卫星数笔记 |
| `data/geomag/IGRF12.COF` | IGRF 地磁系数 |
| `data/geomag/WMM2015.COF` | 世界磁场模型系数 |

### 8.2 `include/ins-gnss` 公共头

| 文件 | 存放什么 |
|--|--|
| `carvig.h` | INS/GNSS 组合库总头：状态、选项、对外函数声明 |
| `egm9615.h` | EGM96 大地水准面系数表 |
| `geomag.h` | 地磁模型接口 |
| `queue.h` | 环形队列 |
| `r8lib.h` | 双精度数组工具声明 |
| `sse2neon.h` | x86 SSE 到 ARM NEON 的映射 |
| `sse2neonbinding.h` | 上述映射的绑定声明 |
| `triangle.h` | 三角剖分接口。视觉几何用，主路径不调用 |
| `structure/uthash.h` | uthash 哈希表 |
| `structure/utlist.h` | 链表宏 |
| `structure/utarray.h` | 动态数组宏 |
| `structure/utstack.h` | 栈宏 |
| `structure/utstring.h` | 动态字符串宏 |
| `structure/utringbuffer.h` | 环形缓冲宏 |

### 8.3 惯导、紧组合和运动约束

| 文件 | 存放什么 |
|--|--|
| `ins.cc` | 机械编排、圆锥误差和划桨补偿、`updateins` |
| `ins-gnss.cc` | INS/GNSS 松组合公共过程 |
| `ins-gnss-tc.cc` | 紧组合量测更新。水面 RTK 相位的算法主体 |
| `ins-gnss-state.cc` | 组合状态向量的索引和接口 |
| `ins-gnss-iekf.cc` | 松组合迭代卡尔曼 |
| `ins-gnss-rb.cc` | 以机器人为中心的松组合表述 |
| `ins-gnss-rts.cc` | RTS 平滑，后处理用 |
| `ins-gnss-fbsm.cc` | 前向/反向固定区间平滑 |
| `ins-align.cc` | 静基座初始对准 |
| `ins-init-rt.cc` | 实时初始化 |
| `ins-static-detect.cc` | 用 IMU 判断静止 |
| `ins-back-mech.cc` | 反向机械编排 |
| `ins-odo.cc` | 陆地里程计辅助。船体模式不要调用 |
| `ins-nhc.cc` | 陆地非完整约束。船体模式不要调用 |
| `ins-zvu.cc` | 陆地零速修正。船体模式不要调用 |
| `ins-zaru.cc` | 陆地零角速度修正。船体模式不要调用 |
| `ins-doppler.cc` | 多普勒辅助 |
| `ins-dualant-calib.cc` | 双天线安装标定 |
| `ins-geomag.cc` | 地磁量测辅助 |
| `ins-sim.cc` | 惯导仿真 |
| `quaternion.cc` | 四元数 |
| `so3.cc` | SO(3) 李群 |
| `se3.cc` | SE(3) 李群 |
| `exp-matrix.cc` | 矩阵指数 |
| `postpos.cc` | 后处理定位流程 |
| `tcpostpos.cc` | 紧组合后处理流程：先 GNSS，再运动约束 |
| `solution.cc` | 把解格式化成文件输出 |
| `options.cc` | 解析 `carvig-*.conf` |
| `null.cc` | 关闭可选功能时的空实现 |

### 8.4 GNSS 定位、模糊度和产品

| 文件 | 存放什么 |
|--|--|
| `rtkcmn.cc` | RTKLIB 公共：时间、坐标、矩阵、卫星号 |
| `rtkpos.cc` | RTK 解算 |
| `pntpos.cc` | 单点定位 |
| `ppp.cc` | 精密单点定位 |
| `ppp_ar.cc` | PPP 模糊度固定 |
| `ppp_corr.cc` | PPP 改正数 |
| `ar.cc` | 模糊度固定策略 |
| `lambda.cc` | LAMBDA 整数最小二乘 |
| `parlambda.cc` | 部分模糊度搜索 |
| `parlambda.h` | 部分模糊度搜索声明 |
| `PAR_LAMBDA.cc` | 另一套部分模糊度实现 |
| `PAR_LAMBDA.h` | 上一文件的声明 |
| `refsat.cc` | 按观测质量选基准星 |
| `cn0_varerr.cc` | 用载噪比相对高度角模板估计伪距方差 |
| `ephemeris.cc` | 广播星历和钟差 |
| `preceph.cc` | 精密星历和钟差 |
| `rinex.cc` | RINEX 观测和导航电文 |
| `ionex.cc` | 电离层格网 |
| `tides.cc` | 潮汐改正 |
| `geoid.cc` | 大地水准面高 |
| `datum.cc` | 基准转换 |
| `sbas.cc` | SBAS 改正 |
| `qzslex.cc` | QZSS LEX |
| `tle.cc` | NORAD 两行根数 |
| `gnss-simobs.cc` | 仿真 GNSS 观测 |

### 8.5 差分电文、数据流和接收机

这些文件保留原库的电文解码。本版节点不订阅 RTK 话题，紧组合接上后再把原始观测送进来。

| 文件 | 存放什么 |
|--|--|
| `rtcm.cc` | RTCM 公共 |
| `rtcm2.cc` | RTCM 2 解码 |
| `rtcm3.cc` | RTCM 3 解码 |
| `rtcm3e.cc` | RTCM 3 编码 |
| `stream.cc` | 串口、TCP、文件流 |
| `streamsvr.cc` | 流服务器 |
| `rcvraw.cc` | 接收机原始报文分发 |
| `convrnx.cc` | 原始报文转 RINEX |
| `convkml.cc` | 解转 KML |
| `convgpx.cc` | 解转 GPX |
| `gis.cc` | Shapefile 读写 |
| `vt.cc` | 文本虚拟控制台 |
| `carvigsvr.cc` | 库的服务端循环，原工程在线入口之一 |
| `rcv/ublox.cc` | u-blox 二进制 |
| `rcv/novatel.cc` | NovAtel OEM |
| `rcv/novatel-oem6.cc` | NovAtel OEM6 |
| `rcv/septentrio.cc` | Septentrio SBF |
| `rcv/comnav.cc` | ComNav OEM |
| `rcv/tersus.cc` | Tersus |
| `rcv/switnav.cc` | Swift Navigation 二进制 |
| `rcv/javad.cc` | Javad GREIS |
| `rcv/skytraq.cc` | SkyTraq |
| `rcv/nvs.cc` | NVS BINR |
| `rcv/crescent.cc` | Hemisphere Crescent |
| `rcv/ss2.cc` | NovAtel Superstar II |
| `rcv/rt17.cc` | Trimble RT-17 |
| `rcv/gsof.cc` | Trimble GSOF |
| `rcv/cmr.cc` | Trimble CMR |
| `rcv/binex.cc` | BINEX |
| `rcv/gw10.cc` | Furuno GW-10 |
| `rcv/rcvlex.cc` | QZSS LEX 接收机 |
| `rcv/rinex-rt.cc` | 实时 RINEX 流 |
| `rcv/imu.cc` | IMU 原始报文 |
| `rcv/stim300.cc` | STIM300 IMU |
| `rcv/imucsv.cc` | IMU CSV |
| `rcv/igv-sim.cc` | INS/GNSS/视觉仿真数据 |
| `rcv/ground-truth.cc` | 读真值，用于离线对比 |
| `rcv/pgm.cc` | PGM 图像。视觉用，主路径不调用 |
| `rcv/serialisation.h` | 序列化声明 |
| `rcv/serialisation.cc` | 序列化实现 |

### 8.6 数值、地磁和视觉残留

这些文件随库拷贝。紧组合和水下轮速融合的主路径不依赖视觉文件；可选静态库若链接失败，先查它们的外部依赖，而不是船体公式。

| 文件 | 存放什么 |
|--|--|
| `r8lib.cc` | 双精度工具函数 |
| `qrdcmp.cc` | QR 分解 |
| `svd.cc` | 奇异值分解 |
| `geomag.cc` | IGRF 计算，读 `data/geomag` |
| `sse2neon-binding.cc` | SSE/NEON 绑定实现 |
| `triangle.cc` | 三角剖分 |
| `ins-camera.cc` | 相机模型 |
| `ins-pose.cc` | 位姿量测辅助 |
| `ins-pnp.cc` | PnP |
| `ins-flow-track.cc` | KLT 光流 |
| `ins-track.cc` | 特征跟踪 |
| `ins-track-visualizer.cc` | 特征跟踪显示 |
| `ins-libfast.cc` | FAST 角点 |
| `kalibr-rosbag-creator.cc` | 导出 Kalibr 标定包 |

### 8.7 `src/ins-gnss/app` 文件型程序

ROS 节点替代它们的输入输出。需要离线复现时，仍可在原 `gnss_imu_couple` 工程里编 `tc-postpos`。

| 文件 | 存放什么 |
|--|--|
| `app/tc-postpos.cc` | 紧组合后处理 `main` |
| `app/rnx2rtkp.cc` | RINEX 定位 `main` |
| `app/lc-rts.cc` | 松组合 RTS 平滑 `main` |
| `app/lc-fbsm.cc` | 松组合前后向平滑 `main` |
| `app/igvonav.cc` | 带视觉的控制台 `main`。本包主路径不使用 |
| `app/carvig-tc.conf` | 紧组合配置原件，`config/carvig-tc.conf` 是它的副本 |
| `app/carvig.conf` | 通用配置 |
| `app/carvig-rtk.conf` | 纯 RTK 配置 |
| `app/carvig-lc.conf` | 松组合配置 |
| `app/carvig-sim.conf` | 仿真配置 |
