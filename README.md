# naiad_localization_v2

第一版节点只做大船坐标系下的 IMU 与轮速融合（船体模式 A）。不订阅 `sensor_rtk/data` 和 `sensor_xcar/data`。

包名 `naiad_localization_v2`，可执行文件同名。框架说明在 `man/framework.md`。

## 环境

- ROS 2 Humble
- 工作空间 `/home/neptune-lindx/src/C/naiad_workspace` 已能编过 `naiad_interfaces`，安装目录是 `staging/install.x86_64`
- 在工作空间根目录执行 `./build.py`。若提示缺少 `config.json`，先做一次：

```sh
cd /home/neptune-lindx/src/C/naiad_workspace
ln -s config/config_naiad_v2.json config.json
```

## 编译

```sh
cd /home/neptune-lindx/src/C/naiad_workspace
source /opt/ros/humble/setup.bash
./build.py -b naiad_localization_v2
```

产物在 `staging/install.x86_64`。默认不编译 `src/ins-gnss`。要编那份静态库时，在上述命令之外给 CMake 加 `-DNAIAD_BUILD_INS_GNSS=ON`；当前节点还不链接它。

## 运行

```sh
cd /home/neptune-lindx/src/C/naiad_workspace
source /opt/ros/humble/setup.bash
source env_x86_64.sh
ros2 launch naiad_localization_v2 localization_v2.launch.py
```

`env_x86_64.sh` 要在工作空间根目录执行，它会把 `staging/install.x86_64` 加进环境。

## 设位

节点起来并且已经收到 `fog/data`、`fog/state`（对准完成）和 `ais/data` 之后，另开一个终端设初始位置。`position` 是后轮中心在对外大船坐标下的位置，右 $x$、前 $y$、上 $z$，单位米。固定锚点经纬高也是后轮中心。`orientation` 当前不用。初始姿态用最近一帧 FOG 四元数：前-上-右到北-天-东换成前-右-下到北-东-地，再用 AIS 船首向（从北向东为正）转到相对大船的姿态。没有 RTK 时，全局锚点固定为经度 114.0°、纬度 30.0°、高度 0，速度为 0。

```sh
cd /home/neptune-lindx/src/C/naiad_workspace
source /opt/ros/humble/setup.bash
source env_x86_64.sh
ros2 service call /algo/set_position naiad_interfaces/srv/RobotSetPosition "{frame_id: '', position: {x: 0.0, y: 0.0, z: 0.0}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 1.0}}"
```

成功后控制台先打一行 `algo/set_position ship=...`，之后每一帧 IMU 打融合位置 `fusion ship p=... llh=...`。还没收到 `fog/data` 时服务返回「没有fog数据」。

不用 launch、直接跑节点：

```sh
ros2 run naiad_localization_v2 naiad_localization_v2 --ros-args \
  --params-file src/naiad_algodash_localization_v2/config/localization_v2.yaml
```

启动后日志里每个订阅打一行 `subscribe <话题> (<参数名>)`。状态在 `localization_v2/status`。滤波还没接入时 `filter_wired=0`，节点不发布 `localization_v2/odometry`。

## 当前订阅

| 话题 | 消息 |
|--|--|
| `motion/data` | `MotionData` |
| `motion/odometer` | `MotionOdometer` |
| `motor/info` | `DeviceInfo` |
| `fog/data` | `FogData` |
| `fog/state` | `FogState` |
| `sensor_depth/data` | `DepthData` |
| `ais/data` | `AisData` |
| `robot/profile` | `std_msgs/String` |
| `task/parameter` | `TaskParameter` |
| `pumpbox/state` | `PumpboxState` |
| `pushbox/state` | `PushboxState` |

`robot/profile` 的 QoS 是 `transient_local`、深度 1。其余队列深度 10。话题名可以在 `config/localization_v2.yaml` 里改。
