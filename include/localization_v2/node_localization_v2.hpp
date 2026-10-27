#ifndef LOCALIZATION_V2_NODE_LOCALIZATION_V2_HPP_
#define LOCALIZATION_V2_NODE_LOCALIZATION_V2_HPP_

#include <string>

#include "localization_v2/depth_gate.hpp"
#include "localization_v2/gnss_tc_bridge.hpp"
#include "localization_v2/hull_fusion.hpp"
#include "localization_v2/mode_manager.hpp"
#include "localization_v2/static_gate.hpp"
#include "localization_v2/wheel_speed.hpp"
#include "naiad_interfaces/msg/api_state.hpp"
#include "naiad_interfaces/srv/robot_set_position.hpp"
#include "naiad_interfaces/msg/ais_data.hpp"
#include "naiad_interfaces/msg/depth_data.hpp"
#include "naiad_interfaces/msg/device_info.hpp"
#include "naiad_interfaces/msg/fog_data.hpp"
#include "naiad_interfaces/msg/fog_state.hpp"
#include "naiad_interfaces/msg/motion_data.hpp"
#include "naiad_interfaces/msg/motion_odometer.hpp"
#include "naiad_interfaces/msg/pumpbox_state.hpp"
#include "naiad_interfaces/msg/pushbox_state.hpp"
#include "naiad_interfaces/msg/robot_odometry.hpp"
#include "naiad_interfaces/msg/task_parameter.hpp"
#include "rclcpp/rclcpp.hpp"
#include "std_msgs/msg/string.hpp"

namespace naiad
{
namespace localization_v2
{

class NodeLocalizationV2 : public rclcpp::Node
{
public:
  explicit NodeLocalizationV2(const std::string & name);

private:
  template <typename MsgT>
  typename rclcpp::Subscription<MsgT>::SharedPtr subscribe_input(
    const std::string & param_name, const rclcpp::QoS & qos, MsgT & slot, bool & have);

  void on_tick();
  void on_set_position(
    const std::shared_ptr<naiad_interfaces::srv::RobotSetPosition::Request> request,
    std::shared_ptr<naiad_interfaces::srv::RobotSetPosition::Response> response);

  bool wheel_speed_is_zero() const;
  void load_lever(const char * name, double out[3]);
  bool rtk_available() const { return false; }
  double stamp_sec(const builtin_interfaces::msg::Time & stamp) const;

  DepthGate depth_gate_;
  StaticGate static_gate_;
  ModeManager mode_manager_;
  GnssTcBridge gnss_tc_;
  HullFusion hull_;
  WheelOdomSpeed wheel_odom_speed_;
  WheelLeverArms wheel_levers_;
  WheelSpeedAtFog wheel_speed_at_fog_;

  rclcpp::Subscription<naiad_interfaces::msg::MotionData>::SharedPtr motion_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::MotionOdometer>::SharedPtr wheel_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::DeviceInfo>::SharedPtr motor_info_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::FogData>::SharedPtr fog_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::FogState>::SharedPtr fog_state_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::DepthData>::SharedPtr depth_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::AisData>::SharedPtr ais_sub_;
  rclcpp::Subscription<std_msgs::msg::String>::SharedPtr profile_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::TaskParameter>::SharedPtr task_parameter_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::PumpboxState>::SharedPtr pumpbox_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::PushboxState>::SharedPtr pushbox_sub_;
  rclcpp::Service<naiad_interfaces::srv::RobotSetPosition>::SharedPtr set_position_service_;
  rclcpp::Publisher<std_msgs::msg::String>::SharedPtr status_pub_;
  rclcpp::Publisher<naiad_interfaces::msg::RobotOdometry>::SharedPtr odom_pub_;
  rclcpp::TimerBase::SharedPtr timer_;

  naiad_interfaces::msg::MotionData latest_motion_;
  naiad_interfaces::msg::MotionOdometer latest_wheel_;
  naiad_interfaces::msg::DeviceInfo latest_motor_info_;
  naiad_interfaces::msg::FogData latest_fog_;
  naiad_interfaces::msg::FogState latest_fog_state_;
  naiad_interfaces::msg::DepthData latest_depth_;
  naiad_interfaces::msg::AisData latest_ais_;
  std_msgs::msg::String latest_profile_;
  naiad_interfaces::msg::TaskParameter latest_task_parameter_;
  naiad_interfaces::msg::PumpboxState latest_pumpbox_;
  naiad_interfaces::msg::PushboxState latest_pushbox_;
  bool have_motion_{false};
  bool have_wheel_{false};
  bool have_motor_info_{false};
  bool have_fog_{false};
  bool have_fog_state_{false};
  bool have_depth_{false};
  bool have_ais_{false};
  bool have_profile_{false};
  bool have_task_parameter_{false};
  bool have_pumpbox_{false};
  bool have_pushbox_{false};

  double dive_on_bar_{1.20};
  double dive_off_bar_{1.12};
  double align_hold_sec_{30.0};
  double wheel_static_m_{1.0e-3};
  double gyro_static_{0.05};
  double wheel_speed_zero_{1.0e-3};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
