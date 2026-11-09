#include "localization_v2/wheel_speed.hpp"

#include <cmath>
#include <string>

#include "naiad_interfaces/msg/ais_data.hpp"
#include "naiad_interfaces/msg/fog_data.hpp"
#include "naiad_interfaces/msg/motion_odometer.hpp"
#include "naiad_interfaces/msg/robot_odometry.hpp"
#include "rclcpp/rclcpp.hpp"

namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg2Rad = kPi / 180.0;
constexpr double kRad2Deg = 180.0 / kPi;

double stamp_sec(const builtin_interfaces::msg::Time & stamp)
{
  return static_cast<double>(stamp.sec) + static_cast<double>(stamp.nanosec) * 1e-9;
}

/// Hamilton 四元数到方向余弦，列主序。与 hull_fusion::quat_to_dcm 相同。
void quat_to_dcm(double w, double x, double y, double z, double * c)
{
  const double n = std::sqrt(w * w + x * x + y * y + z * z);
  if (n < 1e-12) {
    c[0] = c[4] = c[8] = 1.0;
    c[1] = c[2] = c[3] = c[5] = c[6] = c[7] = 0.0;
    return;
  }
  w /= n;
  x /= n;
  y /= n;
  z /= n;
  c[0] = w * w + x * x - y * y - z * z;
  c[1] = 2.0 * (x * y + w * z);
  c[2] = 2.0 * (x * z - w * y);
  c[3] = 2.0 * (x * y - w * z);
  c[4] = w * w - x * x + y * y - z * z;
  c[5] = 2.0 * (y * z + w * x);
  c[6] = 2.0 * (x * z + w * y);
  c[7] = 2.0 * (y * z - w * x);
  c[8] = w * w - x * x - y * y + z * z;
}

void nte_to_ned(const double nte[3], double ned[3])
{
  ned[0] = nte[0];
  ned[1] = nte[2];
  ned[2] = -nte[1];
}

/// 列主序，机体前右下到北东地，转成 3-2-1 角，弧度。
void dcm_to_rpy(const double * c_bn, double * roll, double * pitch, double * yaw)
{
  *roll = std::atan2(c_bn[5], c_bn[8]);
  *pitch = -std::asin(std::max(-1.0, std::min(1.0, c_bn[2])));
  *yaw = std::atan2(c_bn[1], c_bn[0]);
}

void dcm_to_quat(const double * c, double * q)
{
  const double tr = c[0] + c[4] + c[8];
  if (tr > 0.0) {
    const double s = 0.5 / std::sqrt(tr + 1.0);
    q[0] = 0.25 / s;
    q[1] = (c[5] - c[7]) * s;
    q[2] = (c[6] - c[2]) * s;
    q[3] = (c[1] - c[3]) * s;
    return;
  }
  if (c[0] > c[4] && c[0] > c[8]) {
    const double s = 2.0 * std::sqrt(1.0 + c[0] - c[4] - c[8]);
    q[0] = (c[5] - c[7]) / s;
    q[1] = 0.25 * s;
    q[2] = (c[3] + c[1]) / s;
    q[3] = (c[6] + c[2]) / s;
    return;
  }
  if (c[4] > c[8]) {
    const double s = 2.0 * std::sqrt(1.0 + c[4] - c[0] - c[8]);
    q[0] = (c[6] - c[2]) / s;
    q[1] = (c[3] + c[1]) / s;
    q[2] = 0.25 * s;
    q[3] = (c[7] + c[5]) / s;
    return;
  }
  const double s = 2.0 * std::sqrt(1.0 + c[8] - c[0] - c[4]);
  q[0] = (c[1] - c[3]) / s;
  q[1] = (c[6] + c[2]) / s;
  q[2] = (c[7] + c[5]) / s;
  q[3] = 0.25 * s;
}

void ship_frd_to_rfu(const double frd[3], double rfu[3])
{
  rfu[0] = frd[1];
  rfu[1] = frd[0];
  rfu[2] = -frd[2];
}

}  // namespace

namespace naiad
{
namespace localization_v2
{

/// 后轮左右大轮平均速度，用 FOG 四元数和 AIS 船首向转到大船系后积分。
class TestOdo : public rclcpp::Node
{
public:
  explicit TestOdo(const std::string & name)
  : rclcpp::Node(name)
  {
    declare_parameter<std::string>("fog_topic", "fog/data");
    declare_parameter<std::string>("wheel_topic", "motion/odometer");
    declare_parameter<std::string>("ais_topic", "ais/data");
    declare_parameter<std::string>("odometry_topic", "localization_v2/wheel_odometry");

    const rclcpp::QoS qos(10);
    fog_sub_ = create_subscription<naiad_interfaces::msg::FogData>(
      get_parameter("fog_topic").as_string(), qos,
      [this](const naiad_interfaces::msg::FogData::SharedPtr msg) {
        const auto & q = msg->orientation;
        if (!std::isfinite(q.w) || !std::isfinite(q.x) || !std::isfinite(q.y) ||
          !std::isfinite(q.z))
        {
          return;
        }
        const double n = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
        if (n < 1e-6) {
          return;
        }
        fog_q_[0] = q.w / n;
        fog_q_[1] = q.x / n;
        fog_q_[2] = q.y / n;
        fog_q_[3] = q.z / n;
        have_fog_ = true;
      });
    ais_sub_ = create_subscription<naiad_interfaces::msg::AisData>(
      get_parameter("ais_topic").as_string(), qos,
      [this](const naiad_interfaces::msg::AisData::SharedPtr msg) {
        if (!std::isfinite(msg->heading) || msg->heading < 0.0f || msg->heading > 360.0f) {
          return;
        }
        ais_heading_deg_ = msg->heading;
        have_ais_ = true;
      });
    wheel_sub_ = create_subscription<naiad_interfaces::msg::MotionOdometer>(
      get_parameter("wheel_topic").as_string(), qos,
      [this](const naiad_interfaces::msg::MotionOdometer::SharedPtr msg) {
        on_wheel(*msg);
      });
    odom_pub_ = create_publisher<naiad_interfaces::msg::RobotOdometry>(
      get_parameter("odometry_topic").as_string(), 10);

    RCLCPP_INFO(
      get_logger(), "test_odo fog=%s wheel=%s ais=%s publish=%s",
      fog_sub_->get_topic_name(), wheel_sub_->get_topic_name(), ais_sub_->get_topic_name(),
      odom_pub_->get_topic_name());
  }

private:
  void on_wheel(const naiad_interfaces::msg::MotionOdometer & msg)
  {
    const auto & motor = msg.motor_odometer;
    const auto & small = msg.wheel_odometer;
    const WheelSpeed speed = wheels_.update(
      stamp_sec(msg.header.stamp), motor[0], motor[1], motor[2], small[0], small[1]);
    if (!speed.valid || !have_fog_ || !have_ais_) {
      return;
    }

    const double v = 0.5 * (speed.motor_right_mps + speed.motor_left_mps);
    // orientation 是前-上-右到北-天-东。轮子前向是机体 x。
    double c_fur_nte[9];
    quat_to_dcm(fog_q_[0], fog_q_[1], fog_q_[2], fog_q_[3], c_fur_nte);
    double v_nte[3];
    double v_ned[3];
    for (int row = 0; row < 3; ++row) {
      v_nte[row] = c_fur_nte[row] * v;
    }
    nte_to_ned(v_nte, v_ned);
    const double psi = ais_heading_deg_ * kDeg2Rad;
    const double cpsi = std::cos(psi);
    const double spsi = std::sin(psi);
    const double v_frd[3] = {
      cpsi * v_ned[0] + spsi * v_ned[1],
      -spsi * v_ned[0] + cpsi * v_ned[1],
      v_ned[2]};
    for (int i = 0; i < 3; ++i) {
      p_frd_[i] += v_frd[i] * speed.dt_sec;
    }

    // 前右下在北东地的三列：x 与前-上-右的 x 相同，y 取 z，z 取 -y。
    double c_frd_ned[9];
    for (int axis = 0; axis < 3; ++axis) {
      double nte[3];
      double ned[3];
      for (int row = 0; row < 3; ++row) {
        nte[row] = c_fur_nte[row + 3 * axis];
      }
      nte_to_ned(nte, ned);
      const double sign = axis == 1 ? -1.0 : 1.0;
      const int col = axis == 0 ? 0 : (axis == 1 ? 2 : 1);
      for (int row = 0; row < 3; ++row) {
        c_frd_ned[row + 3 * col] = sign * ned[row];
      }
    }
    const double c_sn[9] = {
      cpsi, spsi, 0.0,
      -spsi, cpsi, 0.0,
      0.0, 0.0, 1.0};
    double c_bs[9];
    for (int col = 0; col < 3; ++col) {
      for (int row = 0; row < 3; ++row) {
        c_bs[row + 3 * col] =
          c_sn[row] * c_frd_ned[3 * col] +
          c_sn[row + 3] * c_frd_ned[1 + 3 * col] +
          c_sn[row + 6] * c_frd_ned[2 + 3 * col];
      }
    }
    double roll = 0.0;
    double pitch = 0.0;
    double yaw = 0.0;
    dcm_to_rpy(c_bs, &roll, &pitch, &yaw);

    double p_rfu[3];
    double v_rfu[3];
    ship_frd_to_rfu(p_frd_, p_rfu);
    ship_frd_to_rfu(v_frd, v_rfu);
    double q[4];
    dcm_to_quat(c_bs, q);

    naiad_interfaces::msg::RobotOdometry odom;
    odom.header.stamp = msg.header.stamp;
    odom.header.frame_id = "ship";
    odom.child_frame_id = "rear_axle";
    odom.pose.position.x = p_rfu[0];
    odom.pose.position.y = p_rfu[1];
    odom.pose.position.z = p_rfu[2];
    odom.pose.euler.roll = static_cast<float>(roll * kRad2Deg);
    odom.pose.euler.pitch = static_cast<float>(pitch * kRad2Deg);
    odom.pose.euler.yaw = static_cast<float>(yaw * kRad2Deg);
    odom.pose.orientation.w = q[0];
    odom.pose.orientation.x = q[1];
    odom.pose.orientation.y = q[2];
    odom.pose.orientation.z = q[3];
    for (int col = 0; col < 3; ++col) {
      for (int row = 0; row < 3; ++row) {
        odom.pose.dcm[row * 3 + col] = c_bs[row + 3 * col];
      }
    }
    odom.twist.linear.x = v_rfu[0];
    odom.twist.linear.y = v_rfu[1];
    odom.twist.linear.z = v_rfu[2];
    odom_pub_->publish(odom);
  }

  rclcpp::Subscription<naiad_interfaces::msg::FogData>::SharedPtr fog_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::AisData>::SharedPtr ais_sub_;
  rclcpp::Subscription<naiad_interfaces::msg::MotionOdometer>::SharedPtr wheel_sub_;
  rclcpp::Publisher<naiad_interfaces::msg::RobotOdometry>::SharedPtr odom_pub_;
  WheelOdomSpeed wheels_;
  bool have_fog_{false};
  bool have_ais_{false};
  double fog_q_[4]{1.0, 0.0, 0.0, 0.0};
  double ais_heading_deg_{0.0};
  double p_frd_[3]{0.0, 0.0, 0.0};
};

}  // namespace localization_v2
}  // namespace naiad

int main(int argc, char ** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<naiad::localization_v2::TestOdo>("test_odo"));
  rclcpp::shutdown();
  return 0;
}
