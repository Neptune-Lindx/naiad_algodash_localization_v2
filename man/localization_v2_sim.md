# localization_v2 仿真启动流程

World：`neptune-sim/sim/worlds/neptune_cleaning_floating_dual_demo.sdf`  
机器人：`N1_robot`（右舷）  
三个终端，按顺序启动，不要关前面的进程。

本机路径：

- 仿真仓库：`~/src/C/neptune-sim`
- 定位仓库：`~/src/C/naiad_workspace/src/naiad_algodash_localization_v2`
- 定位二进制：`~/src/C/naiad_workspace/staging/install.x86_64/lib/naiad_localization_v2/naiad_localization_v2`

节点参数：`config/localization_v2.yaml`（脚本直接读源文件，不读安装目录里的那份）  
仿真参数：`config/neptune_sim.yaml`  
脚本：`scripts/start_neptune_sim_localization_v2.sh`

和 `naiad_algodash_localization_sim` 的差别：v2 没有 `algo/transform`。设位服务是 `algo/set_position`，姿态取最近一帧 `fog/data` 的 roll、pitch、yaw，四元数不参与设位。不要同时再跑旧的 `naiad_localization`，两个节点都会占用 `algo/set_position`。

---

## 终端 1 — Gazebo + ros_gz_bridge

```bash
cd ~/src/C/neptune-sim
./sim/scripts/run_sim_floating_dual_ship_gps_ros2_bridge.sh
```

等 GUI 起来、world 初始化完成。不要关。

---

## 终端 2 — naiad 消息桥

```bash
cd ~/src/C/neptune-sim
./sim/scripts/run_naiad_sim_msg_bridge.sh
```

看到 `Pipeline summary: enabled=6 active=6` 即可。不要关。

桥提供：`/fog/data`、`/sensor_depth/data`、`/motion/odometer`、`/ais/data`。

---

## 终端 3 — localization_v2

```bash
cd ~/src/C/naiad_workspace/src/naiad_algodash_localization_v2
./scripts/start_neptune_sim_localization_v2.sh start
```

脚本会：

1. 检查 `/fog/data` 已经在。
2. 启动 `naiad_localization_v2`，参数用 `config/localization_v2.yaml`。这时还没有设位，不发里程。
3. 等到 `/sim/N1_robot/rel_odo` 有一帧。这一帧表示 Gazebo 里的黑线起点已经记下。然后再设位，红线的第一帧才能和黑线起点用同一次偏移对齐。
4. 补仿真里没有的话题：
   - `robot/profile`：`naiad_launch/profiles/robot_n1.yaml`
   - `fog/state`：`skip_fog_align: true` 时状态是 `ALIGNMENT_DONE`，跳过真实 FOG 对准
   - `task/parameter`：右舷、吃水见 `neptune_sim.yaml`
   - `motion/data`：空消息。设位不看它
5. 调用 `algo/set_position`。位置是 `initial_x` / `initial_y` / `initial_z`，大船右-前-上，单位米，默认 `0 0 0`。`orientation` 不使用。姿态取对准后的 `fog/data.orientation`，换成前-右-下到北-东-地，再用 `ais/data.heading` 得到相对大船的初始姿态。没对准或没有船首向时设位失败，脚本会重试。
6. 等到 `/localization_v2/odometry` 出一帧。

这次 Gazebo 里如果红线已经按旧的第一帧锁过偏移，只重启本脚本不会改那次偏移。要红线和黑线重新对齐，先关掉 Gazebo 再开。

v2 的设位不等 `align_hold_sec`。`align_hold_sec` 只影响状态字里的 `align_ready`，不挡住里程输出。设位要求 `fog/state` 已对准，并且 `ais/data.heading` 在 0 到 360 度。没收到 `fog/data` 时返回「没有fog数据」，没对准返回「fog未对准」，没有船首向返回「没有大船航向」。

改初值后先停再起：

```bash
./scripts/start_neptune_sim_localization_v2.sh stop
# 编辑 config/neptune_sim.yaml 的 initial_x / initial_y / initial_z
./scripts/start_neptune_sim_localization_v2.sh start
```

融合频率、杆臂、IGG3 门槛在 `config/localization_v2.yaml`。这个脚本读的是源文件，改完重新 `start` 即可，不用为了这几个参数再编译。

---

## 验证

```bash
cd ~/src/C/naiad_workspace/src/naiad_algodash_localization_v2
./scripts/start_neptune_sim_localization_v2.sh status
```

另开终端：

```bash
cd ~/src/C/naiad_workspace
source /opt/ros/humble/setup.bash
source staging/install.x86_64/setup.bash
ros2 topic hz /localization_v2/odometry
ros2 topic echo /localization_v2/odometry --once
ros2 topic echo /localization_v2/status --field data
```

`/localization_v2/odometry` 是 `RobotOdometry`。位置和速度是大船右-前-上。`pose.euler` 仍是滤波内部的机体相对大船欧拉角，航向绕朝下的 z 轴、顺时针为正，和 status 里取反后的航向符号相反。

`/localization_v2/status` 是一行字符串。设位成功后才有姿态：

- `rpy_s=(横滚,俯仰,航向)`：横滚、俯仰仍是机体相对大船的前-右-下，度。航向在打印时取反，从上往下看逆时针为正
- `rpy_n=(横滚,俯仰,航向)`：横滚、俯仰仍是机体相对北-东-地，度。航向同样只在打印时取反，逆时针为正
- `t=`：最近一帧 `fog/data` 的时间戳
- `ps` / `vs`：大船右-前-上的位置和速度
- `llh` / `vn`：经纬高和北东地速度

只看姿态：

```bash
ros2 topic echo /localization_v2/status --field data | grep -oE 'rpy_[sn]=\([^)]*\)'
```

日志：`logs/neptune_sim_<时间戳>/`

| 文件 | 内容 |
|------|------|
| `naiad_localization_v2.log` | 定位节点 |
| `sim_init.log` | 补话题和 `algo/set_position` |
| `init_done` | 设位用的 `x y z`。有这个文件表示设位成功 |

`sim_init.log` 里应有 `set_position OK`。

---

## 停止（逆序）

```bash
# 终端 3
cd ~/src/C/naiad_workspace/src/naiad_algodash_localization_v2
./scripts/start_neptune_sim_localization_v2.sh stop

# 终端 2、终端 1：Ctrl+C
```

`stop` 只停 localization_v2 和补话题的进程，不停 Gazebo 和消息桥。
