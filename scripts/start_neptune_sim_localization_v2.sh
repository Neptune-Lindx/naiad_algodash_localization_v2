#!/usr/bin/env bash
# 在 neptune-sim 已出传感器话题的前提下，启动 naiad_localization_v2：
#   跳过 FOG 对准（发布 ALIGNMENT_DONE）→ algo/set_position 写入位置初值 → /localization_v2/odometry
#
# 用法：
#   ./scripts/start_neptune_sim_localization_v2.sh start
#   ./scripts/start_neptune_sim_localization_v2.sh status
#   ./scripts/start_neptune_sim_localization_v2.sh stop
#
# 前置（不要关）：
#   cd /home/neptune-lindx/src/C/neptune-sim
#   ./sim/scripts/run_sim_floating_dual_ship_gps_ros2_bridge.sh
#   ./sim/scripts/run_naiad_sim_msg_bridge.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
NODE_PARAMS="${PKG_ROOT}/config/localization_v2.yaml"
SIM_PARAMS="${PKG_ROOT}/config/neptune_sim.yaml"
WS_ROOT="${NAIAD_WS:-/home/neptune-lindx/src/C/naiad_workspace}"
INSTALL_SETUP="${WS_ROOT}/staging/install.x86_64/setup.bash"
LOC_BIN="${LOC_BIN:-${WS_ROOT}/staging/install.x86_64/lib/naiad_localization_v2/naiad_localization_v2}"
STATE_FILE="${PKG_ROOT}/logs/.neptune_sim_loc_v2.env"
ODOM_TOPIC="/localization_v2/odometry"

source_ros() {
  if [[ ! -f "${INSTALL_SETUP}" ]]; then
    echo "Missing ${INSTALL_SETUP}" >&2
    exit 1
  fi
  set +u
  # shellcheck disable=SC1090
  source /opt/ros/humble/setup.bash
  # shellcheck disable=SC1090
  source "${INSTALL_SETUP}"
  set -u
}

load_state() {
  if [[ ! -f "${STATE_FILE}" ]]; then
    echo "No session. Run: $0 start" >&2
    exit 1
  fi
  # shellcheck disable=SC1090
  source "${STATE_FILE}"
}

wait_for_topic() {
  local topic="$1" timeout_s="${2:-30}" elapsed=0
  while (( elapsed < timeout_s )); do
    if ros2 topic list 2>/dev/null | grep -qx "${topic}"; then
      return 0
    fi
    sleep 1
    elapsed=$((elapsed + 1))
  done
  echo "ERROR: ${topic} not available after ${timeout_s}s" >&2
  return 1
}

# 黑线起点由 Gazebo 插件锁在第一帧船体系真值上。红线偏移用「黑线第一个点 − 红线第一帧」。
# 红线第一帧必须在黑线已经有点之后才发出，否则那一帧只被记下来、不画，起点会错开。
# rel_odo 由仿真桥以 best_effort 发出，先按这个可靠性等。
wait_for_black_origin() {
  local topic="/sim/N1_robot/rel_odo" timeout_s="${1:-30}"
  echo "Waiting for one ${topic} message before set_position..."
  if timeout "${timeout_s}" ros2 topic echo "${topic}" --qos-reliability best_effort --once >/dev/null 2>&1; then
    return 0
  fi
  if timeout 8 ros2 topic echo "${topic}" --once >/dev/null 2>&1; then
    return 0
  fi
  echo "ERROR: ${topic} had no message. Black trajectory origin is not locked yet." >&2
  echo "Start Gazebo with the ship-relative trajectory plugin, then run this script again." >&2
  return 1
}

stop_all() {
  if [[ -f "${STATE_FILE}" ]]; then
    # shellcheck disable=SC1090
    source "${STATE_FILE}" || true
  fi
  for p in "${HELPER_PID:-}" "${LOC_PID:-}"; do
    if [[ -n "${p}" ]] && kill -0 "${p}" 2>/dev/null; then
      kill "${p}" 2>/dev/null || true
      wait "${p}" 2>/dev/null || true
    fi
  done
  local orphan
  while read -r orphan; do
    [[ -n "${orphan}" ]] && kill "${orphan}" 2>/dev/null || true
  done < <(pgrep -f "naiad_localization_v2 --ros-args --params-file ${NODE_PARAMS}" 2>/dev/null || true)
  rm -f "${STATE_FILE}"
  echo "Stopped neptune-sim localization_v2."
}

cmd_start() {
  if [[ -f "${STATE_FILE}" ]]; then
    echo "Session exists. Run: $0 stop" >&2
    exit 1
  fi
  if [[ ! -x "${LOC_BIN}" ]]; then
    echo "ERROR: ${LOC_BIN} not found. Build naiad_localization_v2 or set LOC_BIN." >&2
    exit 1
  fi
  if [[ ! -f "${NODE_PARAMS}" || ! -f "${SIM_PARAMS}" ]]; then
    echo "ERROR: missing ${NODE_PARAMS} or ${SIM_PARAMS}" >&2
    exit 1
  fi

  LOG_DIR="${PKG_ROOT}/logs/neptune_sim_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "${LOG_DIR}"
  source_ros

  echo "Checking neptune-sim bridge topics..."
  wait_for_topic "/fog/data" 30 || {
    echo "Start Gazebo and naiad_sim_msg_bridge first. /fog/data is required." >&2
    exit 1
  }
  for topic in /motion/odometer /sensor_depth/data /ais/data; do
    if ! ros2 topic list 2>/dev/null | grep -qx "${topic}"; then
      echo "WARN: ${topic} not found." >&2
    fi
  done

  echo "Starting naiad_localization_v2 with ${NODE_PARAMS}"
  "${LOC_BIN}" --ros-args --params-file "${NODE_PARAMS}" \
    >"${LOG_DIR}/naiad_localization_v2.log" 2>&1 &
  LOC_PID=$!
  sleep 1
  if ! kill -0 "${LOC_PID}" 2>/dev/null; then
    echo "ERROR: naiad_localization_v2 exited. See ${LOG_DIR}/naiad_localization_v2.log" >&2
    tail -30 "${LOG_DIR}/naiad_localization_v2.log" >&2 || true
    exit 1
  fi

  if ! wait_for_black_origin 30; then
    kill "${LOC_PID}" 2>/dev/null || true
    wait "${LOC_PID}" 2>/dev/null || true
    exit 1
  fi

  python3 - "${SIM_PARAMS}" "${LOG_DIR}/init_done" >"${LOG_DIR}/sim_init.log" 2>&1 <<'PY' &
import sys

import rclpy
import yaml
from geometry_msgs.msg import Point, Quaternion
from rclpy.node import Node
from rclpy.qos import DurabilityPolicy, HistoryPolicy, QoSProfile
from std_msgs.msg import String

from naiad_interfaces.msg import FogState, MotionData, ShipInfo, TaskParameter
from naiad_interfaces.srv import RobotSetPosition

params_path, done_path = sys.argv[1], sys.argv[2]
with open(params_path, encoding="utf-8") as fh:
    doc = yaml.safe_load(fh)
cfg = doc["neptune_sim"]["ros__parameters"]

skip_align = bool(cfg.get("skip_fog_align", True))
mark = (
    float(cfg.get("initial_x", 0.0)),
    float(cfg.get("initial_y", 0.0)),
    float(cfg.get("initial_z", 0.0)),
)
board = int(cfg.get("ship_work_board", 0))
draft = float(cfg.get("ship_draft", 0.0))
profile_path = cfg["robot_profile"]
profile_text = open(profile_path, encoding="utf-8").read()


class NeptuneSimInit(Node):
    def __init__(self):
        super().__init__("neptune_sim_loc_v2_init")
        self._done = False
        self._busy = False
        self._attempts = 0
        qos = QoSProfile(
            depth=1,
            durability=DurabilityPolicy.TRANSIENT_LOCAL,
            history=HistoryPolicy.KEEP_LAST,
        )
        self._profile = self.create_publisher(String, "robot/profile", qos)
        self._fog = self.create_publisher(FogState, "fog/state", 10)
        self._task = self.create_publisher(TaskParameter, "task/parameter", 10)
        self._motion = self.create_publisher(MotionData, "motion/data", 10)
        self._set_position = self.create_client(RobotSetPosition, "algo/set_position")
        self.create_timer(0.1, self._tick_fog)
        self.create_timer(1.0, self._tick_task)
        self.create_timer(0.05, self._tick_motion)
        self.create_timer(0.5, self._tick_profile)
        self.create_timer(1.0, self._try_init)
        self.get_logger().info(
            f"skip_fog_align={skip_align} set_position={mark} ship_work_board={board}"
        )

    def _tick_profile(self):
        self._profile.publish(String(data=profile_text))

    def _tick_fog(self):
        msg = FogState()
        msg.header.stamp = self.get_clock().now().to_msg()
        msg.header.frame_id = "fog"
        msg.actived = True
        if skip_align:
            msg.state = FogState.STATE_ALIGNMENT_DONE
            msg.alignment_count = 1
            msg.aligning_elapsed_time = 0
            msg.aligning_remaining_time = 0
        else:
            msg.state = FogState.STATE_ALIGNING
        msg.model_id = FogState.MODEL_FOGTJ_FOG70D_NB
        msg.model = "TJ-CMP-70D"
        msg.version = "sim"
        msg.serial_number = "neptune-sim"
        msg.work_mode = 3
        msg.self_checking = 2047
        self._fog.publish(msg)

    def _tick_task(self):
        msg = TaskParameter()
        msg.header.stamp = self.get_clock().now().to_msg()
        msg.header.frame_id = "map"
        msg.ship_info.type = ShipInfo.TYPE_UNKNOWN
        msg.ship_draft = draft
        msg.ship_work_board = board
        msg.ship_outer_board = (
            ShipInfo.PORT_SIDE if board == ShipInfo.STARBOARD_SIDE else ShipInfo.STARBOARD_SIDE
        )
        msg.operation_scene = TaskParameter.SCENE_ANCHOR
        msg.bilge_protect = False
        msg.track_source = TaskParameter.SMALL_FUSION
        msg.crawl_ability = TaskParameter.CRAWL_NORMAL
        self._task.publish(msg)

    def _tick_motion(self):
        msg = MotionData()
        msg.header.stamp = self.get_clock().now().to_msg()
        msg.header.frame_id = "base_link"
        msg.period_us = 20000
        self._motion.publish(msg)

    def _try_init(self):
        if self._done or self._busy:
            return
        self._attempts += 1
        if self._attempts > 45:
            self.get_logger().error("gave up waiting for algo/set_position")
            self._done = True
            return
        if not self._set_position.service_is_ready():
            self.get_logger().info(f"waiting for algo/set_position ({self._attempts})")
            return
        req = RobotSetPosition.Request()
        req.frame_id = "ship"
        req.position = Point(x=mark[0], y=mark[1], z=mark[2])
        req.orientation = Quaternion(x=0.0, y=0.0, z=0.0, w=1.0)
        self._busy = True
        future = self._set_position.call_async(req)
        future.add_done_callback(self._on_position)

    def _on_position(self, future):
        self._busy = False
        try:
            resp = future.result()
        except Exception as exc:
            self.get_logger().warning(f"algo/set_position failed: {exc}")
            return
        if resp.result.status != 0:
            self.get_logger().warning(
                f"algo/set_position rejected: status={resp.result.status} message={resp.result.message!r}"
            )
            return
        self._done = True
        with open(done_path, "w", encoding="utf-8") as fh:
            fh.write(f"{mark[0]} {mark[1]} {mark[2]}\n")
        self.get_logger().info(
            f"set_position OK {mark}; /localization_v2/odometry is live"
        )


def main():
    rclpy.init()
    node = NeptuneSimInit()
    try:
        rclpy.spin(node)
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
PY
  HELPER_PID=$!

  cat >"${STATE_FILE}" <<EOF
LOG_DIR=${LOG_DIR}
LOC_PID=${LOC_PID}
HELPER_PID=${HELPER_PID}
NODE_PARAMS=${NODE_PARAMS}
SIM_PARAMS=${SIM_PARAMS}
EOF

  echo "Waiting for set_position, then ${ODOM_TOPIC} ..."
  local elapsed=0
  local ready=0
  while (( elapsed < 90 )); do
    if [[ -f "${LOG_DIR}/init_done" ]] && ros2 topic list 2>/dev/null | grep -qx "${ODOM_TOPIC}"; then
      if timeout 8 ros2 topic echo "${ODOM_TOPIC}" --once >/dev/null 2>&1; then
        ready=1
        break
      fi
    fi
    if ! kill -0 "${LOC_PID}" 2>/dev/null; then
      echo "ERROR: naiad_localization_v2 exited. See ${LOG_DIR}/naiad_localization_v2.log" >&2
      tail -30 "${LOG_DIR}/naiad_localization_v2.log" >&2 || true
      exit 1
    fi
    sleep 2
    elapsed=$((elapsed + 2))
  done

  echo
  echo "localization_v2 sim started."
  echo "  node params     : ${NODE_PARAMS}"
  echo "  sim params      : ${SIM_PARAMS}"
  echo "  localization PID: ${LOC_PID}"
  echo "  init PID        : ${HELPER_PID}"
  echo "  output topic    : ${ODOM_TOPIC}"
  echo "  status topic    : /localization_v2/status"
  echo "  logs            : ${LOG_DIR}"
  if (( ready == 1 )); then
    echo "  init            : set_position done ($(cat "${LOG_DIR}/init_done"))"
  else
    echo "  init            : not confirmed within 90s. See ${LOG_DIR}/sim_init.log"
    echo "                    set_position needs /fog/data. Do not run naiad_localization at the same time."
  fi
  echo
  echo "Check:  ros2 topic echo ${ODOM_TOPIC} --once"
  echo "Stop:   $0 stop"
}

cmd_status() {
  load_state
  source_ros
  echo "LOG_DIR=${LOG_DIR}"
  if [[ -n "${LOC_PID:-}" ]] && kill -0 "${LOC_PID}" 2>/dev/null; then
    echo "  localization=${LOC_PID} running"
  else
    echo "  localization=${LOC_PID:-} not running"
  fi
  if [[ -n "${HELPER_PID:-}" ]] && kill -0 "${HELPER_PID}" 2>/dev/null; then
    echo "  init=${HELPER_PID} running"
  else
    echo "  init=${HELPER_PID:-} not running"
  fi
  if [[ -f "${LOG_DIR}/init_done" ]]; then
    echo "  set_position: $(cat "${LOG_DIR}/init_done")"
  else
    echo "  set_position: not done"
  fi
  if ros2 topic list 2>/dev/null | grep -qx "${ODOM_TOPIC}"; then
    timeout 5 ros2 topic hz "${ODOM_TOPIC}" --window 3 2>/dev/null || echo "  ${ODOM_TOPIC}: no messages in 5s"
  else
    echo "  ${ODOM_TOPIC}: not found"
  fi
}

usage() {
  sed -n '3,13p' "$0"
}

CMD="${1:-start}"
shift || true
case "${CMD}" in
  start) cmd_start "$@" ;;
  stop) stop_all ;;
  status) cmd_status ;;
  *) usage; exit 1 ;;
esac
