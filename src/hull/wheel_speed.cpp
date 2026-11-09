#include "localization_v2/wheel_speed.hpp"

namespace naiad
{
namespace localization_v2
{

namespace
{

void cross3(const double * a, const double * b, double * c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/// 与 hull_fusion 的 fur_to_frd 相同：前-上-右到前-右-下。
void fur_to_frd(const double * fur, double * frd)
{
  frd[0] = fur[0];
  frd[1] = fur[2];
  frd[2] = -fur[1];
}

void project_one(
  double v_mps, const double * fog_to_rear, const double * rear_to_wheel, const double * gyro_frd,
  double * body)
{
  const double lever_fur[3] = {
    fog_to_rear[0] + rear_to_wheel[0], fog_to_rear[1] + rear_to_wheel[1],
    fog_to_rear[2] + rear_to_wheel[2]};
  double lever[3];
  fur_to_frd(lever_fur, lever);
  double wxr[3];
  cross3(gyro_frd, lever, wxr);
  body[0] = v_mps - wxr[0];
  body[1] = -wxr[1];
  body[2] = -wxr[2];
}

void midpoint(const double * a, const double * b, double * mid)
{
  mid[0] = 0.5 * (a[0] + b[0]);
  mid[1] = 0.5 * (a[1] + b[1]);
  mid[2] = 0.5 * (a[2] + b[2]);
}

}  // namespace

// 轮速计算用的两帧时间差，单位秒。按 50 Hz 固定，后期可迁到 yaml。
#define WHEEL_ODOM_FRAME_DT_SEC (1.0 / 50.0)

WheelSpeed WheelOdomSpeed::update(
  double t_sec, double motor_right, double motor_left, double motor_front, double small_right,
  double small_left)
{
  const double odometer[kCount] = {motor_right, motor_left, motor_front, small_right, small_left};
  if (!have_prev_) {
    have_prev_ = true;
    t_prev_ = t_sec;
    for (int i = 0; i < kCount; ++i) {
      odometer_prev_[i] = odometer[i];
    }
    speed_ = WheelSpeed{};
    return speed_;
  }

  const double stamp_dt = t_sec - t_prev_;
  t_prev_ = t_sec;
  if (!(stamp_dt > kMinDt) || stamp_dt > kMaxDt) {
    for (int i = 0; i < kCount; ++i) {
      odometer_prev_[i] = odometer[i];
    }
    speed_ = WheelSpeed{};
    return speed_;
  }

  const double dt = WHEEL_ODOM_FRAME_DT_SEC;
  const double rate[kCount] = {
    (odometer[0] - odometer_prev_[0]) / dt, (odometer[1] - odometer_prev_[1]) / dt,
    (odometer[2] - odometer_prev_[2]) / dt, (odometer[3] - odometer_prev_[3]) / dt,
    (odometer[4] - odometer_prev_[4]) / dt};
  speed_.valid = true;
  speed_.dt_sec = dt;
  speed_.motor_right_mps = rate[0];
  speed_.motor_left_mps = rate[1];
  speed_.motor_front_mps = rate[2];
  speed_.small_right_mps = rate[3];
  speed_.small_left_mps = rate[4];
  for (int i = 0; i < kCount; ++i) {
    odometer_prev_[i] = odometer[i];
  }
  return speed_;
}

WheelSpeedAtFog project_wheel_speed_to_fog(
  const WheelSpeed & speed, const WheelLeverArms & levers, double gyro_x, double gyro_y,
  double gyro_z)
{
  WheelSpeedAtFog out;
  if (!speed.valid) {
    return out;
  }
  const double gyro_fur[3] = {gyro_x, gyro_y, gyro_z};
  double gyro[3];
  fur_to_frd(gyro_fur, gyro);
  double rear_mid[3];
  double small_mid[3];
  midpoint(levers.motor_right, levers.motor_left, rear_mid);
  midpoint(levers.small_right, levers.small_left, small_mid);
  const double v_rear = 0.5 * (speed.motor_right_mps + speed.motor_left_mps);
  const double v_small = 0.5 * (speed.small_right_mps + speed.small_left_mps);
  project_one(v_rear, levers.fog_to_rear, rear_mid, gyro, out.motor_right);
  project_one(v_small, levers.fog_to_rear, small_mid, gyro, out.small_right);
  out.valid = true;
  return out;
}

}  // namespace localization_v2
}  // namespace naiad
