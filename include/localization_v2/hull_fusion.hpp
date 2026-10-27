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

  /// 话题 IMU：前 x、上 y、右 z，陀螺 rad/s，加表 m/s^2 比力。
  /// 内部先换成前-右-下，再按模式 A 做时间积分。
  bool integrate_fog(
    double t_sec, double gyro_x, double gyro_y, double gyro_z, double acc_x, double acc_y,
    double acc_z);

  /// 水平静止、比力沿话题 y 轴向上时，短时间速度应留在零附近。
  static bool level_specific_force_holds();

  /// 五轮观测阵：公共列在比例因子为 0 时五行相同，比例因子块为对角。
  static bool wheel_observation_holds();

  /// 模式 A 非完整性约束的两行观测阵，右扰动姿态列与速度有限差分一致。
  static bool nhc_observation_holds();

  /// SurfaceTight 时关掉。从关到开时重置协方差，比例因子估计保留。
  void set_wheel_aiding(bool enabled);

  /// 五维轮速更新，2 Hz。speed 已在算法机体系，只用前向分量。t_sec 用轮速消息时间。
  bool update_wheels(double t_sec, const WheelSpeedAtFog & speed);

  static HullVariant parse(const std::string & text);

private:
  void project_state();
  void propagate_error(double dt, const double * f_b);
  void reset_error_covariance();
  void feedback_error(const double * dx);
  /// 模式 A：约 1 Hz，约束机体系侧向和法向速度为 0。omega_b 为前-右-下角速度。
  void update_nhc(double t_sec, const double * omega_b);
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
  double P_[kErrorN * kErrorN]{};
  double wheel_sec_{-1.0};
  double nhc_sec_{-1.0};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
