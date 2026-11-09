# test_odo

`tools/test_odo.cpp` 编出来的节点。左右后轮大轮里程取平均，用 FOG 姿态和 AIS 船首向转到大船系后积分，不进卡尔曼。

消息用现有的 `RobotOdometry`。位置和速度是大船右 $x$、前 $y$、上 $z$，与 `localization_v2/odometry` 同一套坐标。起点是进程起来后的第一帧有效轮速。

## 编译

在工作空间根目录，和主节点一起编。`test_odo` 没有单独的包。

```sh
cd /home/neptune-lindx/src/C/naiad_workspace
source /opt/ros/humble/setup.bash
./build.py -b naiad_localization_v2
```

产物：

```text
staging/install.x86_64/lib/naiad_localization_v2/test_odo
```

## 运行

`env_x86_64.sh` 要在工作空间根目录执行。

```sh
cd /home/neptune-lindx/src/C/naiad_workspace
source /opt/ros/humble/setup.bash
source env_x86_64.sh
ros2 run naiad_localization_v2 test_odo
```

启动后控制台打一行实际话题名：`test_odo fog=... wheel=... ais=... publish=...`。

三路都到了才开始积分并发布：`fog/data` 的姿态四元数（前-上-右到北-天-东），`ais/data` 的船首向（度，0 到 360），以及 `motion/odometer` 从第二帧起的左右大轮速度。缺任何一路都不发位置。不用 roll、pitch、yaw 重建姿态。

看输出：

```sh
ros2 topic echo /localization_v2/wheel_odometry
```

## 参数

默认话题如下。要改名时用 `--ros-args -p`。

| 参数 | 默认 | 消息 |
|--|--|--|
| `fog_topic` | `fog/data` | `FogData` |
| `wheel_topic` | `motion/odometer` | `MotionOdometer` |
| `ais_topic` | `ais/data` | `AisData` |
| `odometry_topic` | `localization_v2/wheel_odometry` | `RobotOdometry` |

```sh
ros2 run naiad_localization_v2 test_odo --ros-args \
  -p odometry_topic:=localization_v2/wheel_odometry
```
