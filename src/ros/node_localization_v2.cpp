#include "localization_v2/node_localization_v2.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <sstream>
#include <vector>

namespace naiad
{
namespace localization_v2
{
namespace
{

/// 算法船体系是前 x、右 y、下 z。对外大船坐标是右 x、前 y、上 z。
void ship_frd_to_rfu(const double frd[3], double rfu[3])
{
  rfu[0] = frd[1];
  rfu[1] = frd[0];
  rfu[2] = -frd[2];
}

void ship_rfu_to_frd(const double rfu[3], double frd[3])
{
  frd[0] = rfu[1];
  frd[1] = rfu[0];
  frd[2] = -rfu[2];
}

/// 内部航向绕朝下的 z 轴，从上往下看顺时针为正。status 只在打印时改成逆时针为正。
double status_yaw_ccw_deg(double yaw_deg)
{
  double yaw = -yaw_deg;
  if (yaw > 180.0) {
    yaw -= 360.0;
  } else if (yaw <= -180.0) {
    yaw += 360.0;
  }
  return yaw;
}

}  // namespace

template <typename MsgT>
typename rclcpp::Subscription<MsgT>::SharedPtr NodeLocalizationV2::subscribe_input(
  const std::string & param_name, const rclcpp::QoS & qos, MsgT & slot, bool & have)
{
  const std::string topic = get_parameter(param_name).as_string();
  auto subscription = create_subscription<MsgT>(
    topic, qos, [this, &slot, &have](const typename MsgT::SharedPtr msg) {
      slot = *msg;
      have = true;
    });
  RCLCPP_INFO(
    get_logger(), "subscribe %s (%s)", subscription->get_topic_name(), param_name.c_str());
  return subscription;
}

NodeLocalizationV2::NodeLocalizationV2(const std::string & name)
: rclcpp::Node(name)
{
  declare_parameter<int>("schedule_period_ms", 20);
  declare_parameter<double>("wheel_fusion_hz", 50.0);
  declare_parameter<double>("nhc_hz", 100.0);
  declare_parameter<double>("wheel_igg3_k0", 5.0);
  declare_parameter<double>("wheel_igg3_k1", 15.0);
  declare_parameter<double>("dive_on_bar", dive_on_bar_);
  declare_parameter<double>("dive_off_bar", dive_off_bar_);
  declare_parameter<double>("align_hold_sec", align_hold_sec_);
  declare_parameter<double>("wheel_static_m", wheel_static_m_);
  declare_parameter<double>("gyro_static", gyro_static_);
  declare_parameter<double>("wheel_speed_zero", wheel_speed_zero_);
  declare_parameter<std::vector<double>>("lever_fog_to_rear_center", {0.0, 0.0, 0.0});
  declare_parameter<std::vector<double>>("lever_rear_to_motor_right", {0.0, 0.0, 0.4235});
  declare_parameter<std::vector<double>>("lever_rear_to_motor_left", {0.0, 0.0, -0.4235});
  declare_parameter<std::vector<double>>("lever_rear_to_motor_front", {0.68, 0.0, 0.0});
  declare_parameter<std::vector<double>>("lever_rear_to_small_right", {0.0, -0.075, 0.24});
  declare_parameter<std::vector<double>>("lever_rear_to_small_left", {0.0, -0.075, -0.24});
  declare_parameter<std::string>("motion_topic", "motion/data");
  declare_parameter<std::string>("wheel_topic", "motion/odometer");
  declare_parameter<std::string>("motor_info_topic", "motor/info");
  declare_parameter<std::string>("fog_topic", "fog/data");
  declare_parameter<std::string>("fog_state_topic", "fog/state");
  declare_parameter<std::string>("depth_topic", "sensor_depth/data");
  declare_parameter<std::string>("ais_topic", "ais/data");
  declare_parameter<std::string>("profile_topic", "robot/profile");
  declare_parameter<std::string>("task_parameter_topic", "task/parameter");
  declare_parameter<std::string>("pumpbox_topic", "pumpbox/state");
  declare_parameter<std::string>("pushbox_topic", "pushbox/state");
  declare_parameter<std::string>("status_topic", "localization_v2/status");
  declare_parameter<std::string>("odometry_topic", "localization_v2/odometry");

  const double wheel_hz = get_parameter("wheel_fusion_hz").as_double();
  const double nhc_hz = get_parameter("nhc_hz").as_double();
  if (!(wheel_hz > 0.0) || !(nhc_hz > 0.0)) {
    RCLCPP_ERROR(
      get_logger(), "wheel_fusion_hz and nhc_hz must be positive, got %.3f and %.3f", wheel_hz,
      nhc_hz);
  } else {
    hull_.set_update_hz(wheel_hz, nhc_hz);
    RCLCPP_INFO(get_logger(), "wheel fusion %.3f Hz, nhc %.3f Hz", wheel_hz, nhc_hz);
  }
  const double igg3_k0 = get_parameter("wheel_igg3_k0").as_double();
  const double igg3_k1 = get_parameter("wheel_igg3_k1").as_double();
  if (!(igg3_k0 > 0.0) || !(igg3_k1 > igg3_k0)) {
    RCLCPP_ERROR(
      get_logger(), "wheel IGG3 needs 0 < k0 < k1, got %.3f and %.3f", igg3_k0, igg3_k1);
  } else {
    hull_.set_wheel_igg3(igg3_k0, igg3_k1);
    RCLCPP_INFO(get_logger(), "wheel IGG3 k0=%.3f k1=%.3f", igg3_k0, igg3_k1);
  }
  dive_on_bar_ = get_parameter("dive_on_bar").as_double();
  dive_off_bar_ = get_parameter("dive_off_bar").as_double();
  align_hold_sec_ = get_parameter("align_hold_sec").as_double();
  wheel_static_m_ = get_parameter("wheel_static_m").as_double();
  gyro_static_ = get_parameter("gyro_static").as_double();
  wheel_speed_zero_ = get_parameter("wheel_speed_zero").as_double();
  load_lever("lever_fog_to_rear_center", wheel_levers_.fog_to_rear);
  hull_.set_fog_to_rear_fur(wheel_levers_.fog_to_rear);
  load_lever("lever_rear_to_motor_right", wheel_levers_.motor_right);
  load_lever("lever_rear_to_motor_left", wheel_levers_.motor_left);
  load_lever("lever_rear_to_motor_front", wheel_levers_.motor_front);
  load_lever("lever_rear_to_small_right", wheel_levers_.small_right);
  load_lever("lever_rear_to_small_left", wheel_levers_.small_left);
  hull_.set_variant(HullVariant::HeadingOnly);
  if (!HullFusion::level_specific_force_holds()) {
    RCLCPP_ERROR(get_logger(), "mode A level specific-force check failed");
  }
  if (!HullFusion::wheel_observation_holds()) {
    RCLCPP_ERROR(get_logger(), "four-wheel observation H check failed");
  }
  if (!HullFusion::wheel_igg3_holds()) {
    RCLCPP_ERROR(get_logger(), "wheel IGG3 weight check failed");
  }
  if (!HullFusion::nhc_observation_holds()) {
    RCLCPP_ERROR(get_logger(), "mode A NHC observation H check failed");
  }

  const auto period = std::chrono::milliseconds(get_parameter("schedule_period_ms").as_int());
  const rclcpp::QoS sensor_qos(10);
  motion_sub_ = subscribe_input<naiad_interfaces::msg::MotionData>(
    "motion_topic", sensor_qos, latest_motion_, have_motion_);
  {
    const std::string topic = get_parameter("wheel_topic").as_string();
    wheel_sub_ = create_subscription<naiad_interfaces::msg::MotionOdometer>(
      topic, sensor_qos, [this](const naiad_interfaces::msg::MotionOdometer::SharedPtr msg) {
        latest_wheel_ = *msg;
        have_wheel_ = true;
        const auto & motor = msg->motor_odometer;
        const auto & small = msg->wheel_odometer;
        const WheelSpeed speed = wheel_odom_speed_.update(
          stamp_sec(msg->header.stamp), motor[0], motor[1], motor[2], small[0], small[1]);
        // if (speed.valid) {
        //   RCLCPP_INFO(
        //     get_logger(),
        //     "motor wheel right=%.4f left=%.4f front=%.4f m/s small wheel right=%.4f left=%.4f m/s "
        //     "dt=%.4f",
        //     speed.motor_right_mps, speed.motor_left_mps, speed.motor_front_mps,
        //     speed.small_right_mps, speed.small_left_mps, speed.dt_sec);
        // }
        if (!speed.valid || !have_fog_) {
          return;
        }
        const auto & gyro = latest_fog_.gyroscope;
        wheel_speed_at_fog_ = project_wheel_speed_to_fog(
          speed, wheel_levers_, gyro.x, gyro.y, gyro.z);
        hull_.update_wheels(stamp_sec(msg->header.stamp), wheel_speed_at_fog_);
      });
    RCLCPP_INFO(get_logger(), "subscribe %s (wheel_topic)", wheel_sub_->get_topic_name());
  }
  motor_info_sub_ = subscribe_input<naiad_interfaces::msg::DeviceInfo>(
    "motor_info_topic", sensor_qos, latest_motor_info_, have_motor_info_);
  {
    const std::string topic = get_parameter("fog_topic").as_string();
    fog_sub_ = create_subscription<naiad_interfaces::msg::FogData>(
      topic, sensor_qos, [this](const naiad_interfaces::msg::FogData::SharedPtr msg) {
        latest_fog_ = *msg;
        have_fog_ = true;
        if (!hull_.integrate_fog(
              stamp_sec(msg->header.stamp), msg->gyroscope.x, msg->gyroscope.y, msg->gyroscope.z,
              msg->accelerator.x, msg->accelerator.y, msg->accelerator.z)) {
          return;
        }
      });
    RCLCPP_INFO(get_logger(), "subscribe %s (fog_topic)", fog_sub_->get_topic_name());
  }
  fog_state_sub_ = subscribe_input<naiad_interfaces::msg::FogState>(
    "fog_state_topic", sensor_qos, latest_fog_state_, have_fog_state_);
  depth_sub_ = subscribe_input<naiad_interfaces::msg::DepthData>(
    "depth_topic", sensor_qos, latest_depth_, have_depth_);
  ais_sub_ = subscribe_input<naiad_interfaces::msg::AisData>(
    "ais_topic", sensor_qos, latest_ais_, have_ais_);
  profile_sub_ = subscribe_input<std_msgs::msg::String>(
    "profile_topic", rclcpp::QoS(1).transient_local(), latest_profile_, have_profile_);
  task_parameter_sub_ = subscribe_input<naiad_interfaces::msg::TaskParameter>(
    "task_parameter_topic", sensor_qos, latest_task_parameter_, have_task_parameter_);
  pumpbox_sub_ = subscribe_input<naiad_interfaces::msg::PumpboxState>(
    "pumpbox_topic", sensor_qos, latest_pumpbox_, have_pumpbox_);
  pushbox_sub_ = subscribe_input<naiad_interfaces::msg::PushboxState>(
    "pushbox_topic", sensor_qos, latest_pushbox_, have_pushbox_);
  set_position_service_ = create_service<naiad_interfaces::srv::RobotSetPosition>(
    "algo/set_position", std::bind(
                           &NodeLocalizationV2::on_set_position, this, std::placeholders::_1,
                           std::placeholders::_2));
  RCLCPP_INFO(
    get_logger(), "create service %s", set_position_service_->get_service_name());
  status_pub_ = create_publisher<std_msgs::msg::String>(
    get_parameter("status_topic").as_string(), 10);
  odom_pub_ = create_publisher<naiad_interfaces::msg::RobotOdometry>(
    get_parameter("odometry_topic").as_string(), 10);
  timer_ = create_wall_timer(period, std::bind(&NodeLocalizationV2::on_tick, this));
}

void NodeLocalizationV2::load_lever(const char * name, double out[3])
{
  const auto values = get_parameter(name).as_double_array();
  if (values.size() != 3) {
    RCLCPP_ERROR(get_logger(), "lever %s needs 3 values, got %zu", name, values.size());
    out[0] = out[1] = out[2] = 0.0;
    return;
  }
  out[0] = values[0];
  out[1] = values[1];
  out[2] = values[2];
}

bool NodeLocalizationV2::wheel_speed_is_zero() const
{
  const auto & wheel = latest_motion_.wheel_velocity;
  return std::fabs(wheel[0]) <= wheel_speed_zero_ && std::fabs(wheel[1]) <= wheel_speed_zero_ &&
         std::fabs(wheel[2]) <= wheel_speed_zero_ &&
         std::fabs(latest_motion_.real_velocity) <= wheel_speed_zero_;
}

void NodeLocalizationV2::on_set_position(
  const std::shared_ptr<naiad_interfaces::srv::RobotSetPosition::Request> request,
  std::shared_ptr<naiad_interfaces::srv::RobotSetPosition::Response> response)
{
  using ApiState = naiad_interfaces::msg::ApiState;
  using FogState = naiad_interfaces::msg::FogState;
  if (rtk_available()) {
    if (!gnss_tc_.start_quasi_static_align()) {
      response->result.status = ApiState::RET_NOT_IMPLEMENTED;
      response->result.message = "Fail (RTK/IMU quasi-static align is not wired)";
      RCLCPP_INFO(get_logger(), "algo/set_position: RTK present, align interface returned false");
      return;
    }
    response->result.status = ApiState::RET_SUCCESS;
    response->result.message = "Success";
    return;
  }
  // 没有 motion/data，以及轮速 0.001 静止门槛，暂时都不拦截设位。
  if (!have_fog_) {
    response->result.status = ApiState::RET_NO_FOG_DATA;
    response->result.message = "Fail (没有fog数据)";
    RCLCPP_INFO(get_logger(), "algo/set_position rejected: no fog/data");
    return;
  }
  const uint32_t fog_state = latest_fog_state_.state;
  if (!have_fog_state_ || (fog_state & FogState::STATE_ALIGNMENT_DONE) == 0) {
    response->result.status = ApiState::RET_NO_FOG_DATA;
    response->result.message = "Fail (fog未对准)";
    RCLCPP_INFO(get_logger(), "algo/set_position rejected: fog alignment not done");
    return;
  }
  const double ship_heading_deg = latest_ais_.heading;
  if (!have_ais_ || !std::isfinite(ship_heading_deg) ||
    ship_heading_deg < 0.0 || ship_heading_deg > 360.0)
  {
    response->result.status = ApiState::RET_INVALID_DATA;
    response->result.message = "Fail (没有大船航向)";
    RCLCPP_INFO(get_logger(), "algo/set_position rejected: no ais heading");
    return;
  }
  const auto & q = latest_fog_.orientation;
  const double qn = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
  if (!std::isfinite(qn) || qn < 1e-6) {
    response->result.status = ApiState::RET_INVALID_DATA;
    response->result.message = "Fail (fog姿态无效)";
    RCLCPP_INFO(get_logger(), "algo/set_position rejected: fog quaternion invalid");
    return;
  }

  ModeAOrigin origin;
  origin.valid = true;
  const double p_rfu[3] = {request->position.x, request->position.y, request->position.z};
  double p_frd[3];
  ship_rfu_to_frd(p_rfu, p_frd);
  origin.p_s[0] = p_frd[0];
  origin.p_s[1] = p_frd[1];
  origin.p_s[2] = p_frd[2];
  origin.latitude_deg = 30.0;
  origin.longitude_deg = 114.0;
  origin.height_m = 0.0;
  origin.quat_x = q.x;
  origin.quat_y = q.y;
  origin.quat_z = q.z;
  origin.quat_w = q.w;
  origin.ship_heading_deg = ship_heading_deg;
  hull_.set_mode_a_origin(origin);
  const ModeAOrigin & set = hull_.mode_a_origin();

  response->result.status = ApiState::RET_SUCCESS;
  response->result.message = "Success";
  RCLCPP_INFO(
    get_logger(),
    "algo/set_position ship=(%.3f,%.3f,%.3f) llh=(%.1f,%.1f,%.1f) v=0 "
    "heading=%.3f rel_rpy=(%.3f,%.3f,%.3f)",
    p_rfu[0], p_rfu[1], p_rfu[2], origin.latitude_deg, origin.longitude_deg,
    origin.height_m, ship_heading_deg, set.roll_deg, set.pitch_deg, set.yaw_deg);
}

double NodeLocalizationV2::stamp_sec(const builtin_interfaces::msg::Time & stamp) const
{
  if (stamp.sec == 0 && stamp.nanosec == 0) {
    return now().seconds();
  }
  return static_cast<double>(stamp.sec) + static_cast<double>(stamp.nanosec) * 1e-9;
}

void NodeLocalizationV2::on_tick()
{
  const bool depth_valid =
    have_depth_ && latest_depth_.state == naiad_interfaces::msg::DepthData::STATE_VALID;
  const bool submerged = depth_gate_.update(
    depth_valid, latest_depth_.pressure, dive_on_bar_, dive_off_bar_);

  bool still = false;
  if (have_fog_ && have_wheel_) {
    const auto & g = latest_fog_.gyroscope;
    const double gyro_norm = std::sqrt(g.x * g.x + g.y * g.y + g.z * g.z);
    const auto & d = latest_wheel_.wheel_delta_distance;
    const double wheel_span = std::fabs(d[0]) + std::fabs(d[1]) + std::fabs(d[2]);
    still = gyro_norm <= gyro_static_ && wheel_span <= wheel_static_m_;
  }
  const double align_t = have_fog_ ? stamp_sec(latest_fog_.header.stamp) : now().seconds();
  const bool align_ready = static_gate_.update(align_t, have_fog_ && have_wheel_, still, align_hold_sec_);
  ModeInput mode_input;
  mode_input.submerged = submerged;
  mode_input.rtk_fixed = rtk_available();
  mode_input.align_ready = align_ready;
  const WorkMode mode = mode_manager_.update(mode_input);
  hull_.set_wheel_aiding(mode != WorkMode::SurfaceTight);
  const NavStatus status = hull_.step(false);

  builtin_interfaces::msg::Time stamp;
  if (have_fog_ && (latest_fog_.header.stamp.sec != 0 || latest_fog_.header.stamp.nanosec != 0)) {
    stamp = latest_fog_.header.stamp;
  } else {
    const int64_t ns = now().nanoseconds();
    stamp.sec = static_cast<int32_t>(ns / 1000000000);
    stamp.nanosec = static_cast<uint32_t>(ns % 1000000000);
  }
  char stamp_text[32];
  std::snprintf(
    stamp_text, sizeof(stamp_text), "%d.%09u", stamp.sec, stamp.nanosec);

  std_msgs::msg::String text;
  std::ostringstream os;
  os << "t=" << stamp_text
     << " backend=" << status.backend
     << " filter_wired=" << (status.filter_wired ? 1 : 0)
     << " solution_valid=" << (status.solution_valid ? 1 : 0)
     << " submerged=" << (submerged ? 1 : 0)
     << " align_ready=" << (align_ready ? 1 : 0)
     << " rx="
     << (have_motion_ ? "motion," : "")
     << (have_wheel_ ? "wheel," : "")
     << (have_motor_info_ ? "motor," : "")
     << (have_fog_ ? "fog," : "")
     << (have_fog_state_ ? "fog_state," : "")
     << (have_depth_ ? "depth," : "")
     << (have_ais_ ? "ais," : "")
     << (have_profile_ ? "profile," : "")
     << (have_task_parameter_ ? "task," : "")
     << (have_pumpbox_ ? "pumpbox," : "")
     << (have_pushbox_ ? "pushbox," : "");
  if (status.solution_valid) {
    double p_rfu[3];
    double v_rfu[3];
    ship_frd_to_rfu(status.p_s, p_rfu);
    ship_frd_to_rfu(status.v_s, v_rfu);
    os << " ps=(" << p_rfu[0] << "," << p_rfu[1] << "," << p_rfu[2] << ")"
       << " vs=(" << v_rfu[0] << "," << v_rfu[1] << "," << v_rfu[2] << ")"
       << " rpy_s=(" << status.roll_s_deg << "," << status.pitch_s_deg << ","
       << status_yaw_ccw_deg(status.yaw_s_deg) << ")"
       << " llh=(" << status.latitude_deg << "," << status.longitude_deg << "," << status.height_m
       << ")"
       << " vn=(" << status.v_n[0] << "," << status.v_n[1] << "," << status.v_n[2] << ")"
       << " rpy_n=(" << status.roll_n_deg << "," << status.pitch_n_deg << ","
       << status_yaw_ccw_deg(status.yaw_n_deg) << ")";
  }
  text.data = os.str();
  status_pub_->publish(text);

  if (!status.solution_valid) {
    return;
  }
  double p_rfu[3];
  double v_rfu[3];
  ship_frd_to_rfu(status.p_s, p_rfu);
  ship_frd_to_rfu(status.v_s, v_rfu);
  naiad_interfaces::msg::RobotOdometry odom;
  odom.header.stamp = now();
  odom.header.frame_id = "ship";
  odom.child_frame_id = "imu";
  odom.pose.position.x = p_rfu[0];
  odom.pose.position.y = p_rfu[1];
  odom.pose.position.z = p_rfu[2];
  odom.pose.euler.roll = static_cast<float>(status.roll_s_deg);
  odom.pose.euler.pitch = static_cast<float>(status.pitch_s_deg);
  odom.pose.euler.yaw = static_cast<float>(status.yaw_s_deg);
  odom.pose.orientation.w = status.quat_s_w;
  odom.pose.orientation.x = status.quat_s_x;
  odom.pose.orientation.y = status.quat_s_y;
  odom.pose.orientation.z = status.quat_s_z;
  for (int i = 0; i < 9; ++i) {
    odom.pose.dcm[i] = status.C_bs[i];
  }
  odom.twist.linear.x = v_rfu[0];
  odom.twist.linear.y = v_rfu[1];
  odom.twist.linear.z = v_rfu[2];
  odom_pub_->publish(odom);
}

}  // namespace localization_v2
}  // namespace naiad
