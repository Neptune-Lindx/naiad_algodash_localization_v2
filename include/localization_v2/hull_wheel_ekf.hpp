#ifndef LOCALIZATION_V2_HULL_WHEEL_EKF_HPP_
#define LOCALIZATION_V2_HULL_WHEEL_EKF_HPP_

namespace naiad
{
namespace localization_v2
{

/// 模式 A 的 20 维误差状态。顺序：姿态、速度、位置、加表零偏、陀螺零偏、五个轮速比例因子。
class HullWheelEkf
{
public:
  static constexpr int kNx = 20;
  static constexpr int kNv = 5;

  HullWheelEkf();

  void reset();

  /// c_bs、f_b 为前-右-下，列主序。ω_is^s = 0。
  void propagate(double dt, const double * c_bs, const double * f_b);

  /// z 为五个轮折到 fog 后的前向速度。成功时 dx 是要反馈的误差，随后状态视为 0。
  bool update(
    const double * c_bs, const double * v_s, const double * scale, const double * z, double * dx);

private:
  double p_[kNx * kNx];
};

}  // namespace localization_v2
}  // namespace naiad

#endif
