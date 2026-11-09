#ifndef LOCALIZATION_V2_HULL_FUSION_HPP_
#define LOCALIZATION_V2_HULL_FUSION_HPP_

#include <string>

#include "localization_v2/nav_status.hpp"
#include "localization_v2/wheel_speed.hpp"

namespace naiad
{
namespace localization_v2
{

/// 无 RTK 时模式 A 的初始状态。全局锚点固定，机器人位置在大船系。
struct ModeAOrigin
{
  bool valid{false};
  double p_s[3]{0.0, 0.0, 0.0};
  double v_s[3]{0.0, 0.0, 0.0};
  double roll_deg{0.0};
  double pitch_deg{0.0};
  double yaw_deg{0.0};
  double quat_x{0.0};
  double quat_y{0.0};
  double quat_z{0.0};
  double quat_w{1.0};
  double latitude_deg{30.0};
  double longitude_deg{114.0};
  double height_m{0.0};
  double v_global[3]{0.0, 0.0, 0.0};
};

/// 大船坐标系 IMU + 轮速融合入口。
/// 公式见 man/ins-hull-rel.md。不要调用陆地 odo() / nhc() / zvu()。
enum class HullVariant
{
  Auto,          ///< 壁虎 1 Hz 新鲜则模式 B，否则模式 A
  HeadingOnly,   ///< 模式 A：仅低频艏向
  ShipNav        ///< 模式 B：壁虎位置速度姿态
};

class HullFusion
{
public:
  void set_variant(HullVariant variant) { variant_ = variant; }
  HullVariant variant() const { return variant_; }

  /// xcar_fresh 只在 Auto 下用来选择 A/B。
  NavStatus step(bool xcar_fresh) const;

  void set_mode_a_origin(const ModeAOrigin & origin);
  const ModeAOrigin & mode_a_origin() const { return origin_; }

  /// 话题轴前-上-右，IMU 到后轮中心。输出的大船位置和经纬高落在这个点上。
  void set_fog_to_rear_fur(const double lever[3]);

  /// 话题 IMU：前 x、上 y、右 z，陀螺 rad/s，加表 m/s^2 比力。
  /// 内部先换成前-右-下，再按模式 A 做时间积分。
  bool integrate_fog(
    double t_sec, double gyro_x, double gyro_y, double gyro_z, double acc_x, double acc_y,
    double acc_z);

  /// 水平静止、比力沿话题 y 轴向上时，短时间速度应留在零附近。
  static bool level_specific_force_holds();

  /// 轮速观测阵：大轮轴中心在第 0 行，小轮轴中心在第 3 行。两行姿态三列为 0，速度列形式相同。
  /// 大轮比例因子占 s_0、s_1，小轮占 s_3、s_4，同一对的两列都是 v_x/2。其余三行为 0。
  static bool wheel_observation_holds();

  /// 模式 A 非完整性约束：侧向和法向都只改速度，姿态三列为 0。
  static bool nhc_observation_holds();

  /// SurfaceTight 时关掉。从关到开时重置协方差，比例因子估计保留。
  void set_wheel_aiding(bool enabled);

  /// 轮速融合和非完整性约束的更新频率，单位 Hz。两个值都要大于 0。
  void set_update_hz(double wheel_hz, double nhc_hz);

  /// 轮速更新。频率由 set_update_hz 决定。大轮和小轮各自左右平均后补到 IMU，只用前向分量。姿态三列为 0。t_sec 用轮速消息时间。
  bool update_wheels(double t_sec, const WheelSpeedAtFog & speed);

  static HullVariant parse(const std::string & text);

private:
  void project_state();
  void propagate_error(double dt, const double * f_b);
  void reset_error_covariance();
  void feedback_error(const double * dx);
  /// 模式 A。频率由 set_update_hz 决定。侧向把后轮中心右向速度拉回 0，法向把 IMU 下向速度拉回 0。两行都不改姿态。
  /// omega_b、lever_frd 为前-右-下；omega_b 已减陀螺零偏，lever_frd 是 FOG 到后轮中心。速度门和角速度门已屏蔽。
  void update_nhc(double t_sec, const double * omega_b, const double * lever_frd);
  bool apply_scalar_update(
    const double * row, double prediction, double measurement, double variance);

  HullVariant variant_{HullVariant::Auto};
  ModeAOrigin origin_{};
  NavStatus nav_{};
  bool armed_{false};
  double t_sec_{-1.0};
  double re_[3]{0.0, 0.0, 0.0};
  double ve_[3]{0.0, 0.0, 0.0};
  double Cbe_[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double Cse_[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double rs_[3]{0.0, 0.0, 0.0};
  bool have_prev_{false};
  double prev_dt_{0.01};
  double prev_w_[3]{0.0, 0.0, 0.0};
  double prev_f_[3]{0.0, 0.0, 0.0};

  static constexpr int kErrorN = 20;
  bool wheel_aiding_{true};
  double ba_[3]{0.0, 0.0, 0.0};
  double bg_[3]{0.0, 0.0, 0.0};
  double scale_[5]{0.0, 0.0, 0.0, 0.0, 0.0};
  double f_b_[3]{0.0, 0.0, 0.0};
  double lever_fur_[3]{0.0, 0.0, 0.0};
  double P_[kErrorN * kErrorN]{};
  double wheel_period_sec_{0.02};
  double nhc_period_sec_{0.01};
  double wheel_sec_{-1.0};
  double nhc_sec_{-1.0};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
