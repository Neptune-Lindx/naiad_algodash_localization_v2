#ifndef LOCALIZATION_V2_WHEEL_SPEED_HPP_
#define LOCALIZATION_V2_WHEEL_SPEED_HPP_

namespace naiad
{
namespace localization_v2
{

/// motion/odometer 的轮速。大轮来自 motor_odometer（右、左、前），
/// 小轮来自 wheel_odometer 的右、左。单位 m/s。
struct WheelSpeed
{
  bool valid{false};
  double dt_sec{0.0};
  double motor_right_mps{0.0};
  double motor_left_mps{0.0};
  double motor_front_mps{0.0};
  double small_right_mps{0.0};
  double small_left_mps{0.0};
};

/// 配置文件里的杆臂，fog 机体：前 x、上 y、右 z，单位 m。
struct WheelLeverArms
{
  double fog_to_rear[3]{0.0, 0.0, 0.0};
  double motor_right[3]{0.0, 0.0, 0.4235};
  double motor_left[3]{0.0, 0.0, -0.4235};
  double motor_front[3]{0.68, 0.0, 0.0};
  double small_right[3]{0.0, -0.075, 0.24};
  double small_left[3]{0.0, -0.075, -0.24};
};

/// 各轮速折到 fog 原点后的速度，单位 m/s。分量是算法机体系前、右、下。
struct WheelSpeedAtFog
{
  bool valid{false};
  double motor_right[3]{0.0, 0.0, 0.0};
  double motor_left[3]{0.0, 0.0, 0.0};
  double motor_front[3]{0.0, 0.0, 0.0};
  double small_right[3]{0.0, 0.0, 0.0};
  double small_left[3]{0.0, 0.0, 0.0};
};

/// 话题与杆臂是前-上-右，先换成算法机体系前-右-下。轮速沿前向，v = (v,0,0) - ω×l。
WheelSpeedAtFog project_wheel_speed_to_fog(
  const WheelSpeed & speed, const WheelLeverArms & levers, double gyro_x, double gyro_y,
  double gyro_z);

/// 用累积里程差分得到速度。两帧时间差固定为 WHEEL_ODOM_FRAME_DT_SEC。
/// 第一帧只记里程，从第二帧起 valid。
class WheelOdomSpeed
{
public:
  WheelSpeed update(
    double t_sec, double motor_right, double motor_left, double motor_front, double small_right,
    double small_left);

  const WheelSpeed & speed() const { return speed_; }

private:
  static constexpr double kMinDt = 1.0e-4;
  static constexpr double kMaxDt = 1.0;
  static constexpr int kCount = 5;

  bool have_prev_{false};
  double t_prev_{0.0};
  double odometer_prev_[kCount]{0.0, 0.0, 0.0, 0.0, 0.0};
  WheelSpeed speed_{};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
