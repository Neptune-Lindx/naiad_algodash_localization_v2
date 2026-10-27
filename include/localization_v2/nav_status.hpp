#ifndef LOCALIZATION_V2_NAV_STATUS_HPP_
#define LOCALIZATION_V2_NAV_STATUS_HPP_

namespace naiad
{
namespace localization_v2
{

/// 算法入口的返回。filter_wired 为假时不得把 p_s / v_s 当成导航解发布。
struct NavStatus
{
  bool filter_wired{false};
  bool solution_valid{false};
  const char * backend{"none"};
  bool fast_search_armed{false};
  double p_s[3]{0.0, 0.0, 0.0};
  double v_s[3]{0.0, 0.0, 0.0};
  /// 机体相对准船体的姿态，度。机体轴为前-右-下。
  double roll_s_deg{0.0};
  double pitch_s_deg{0.0};
  double yaw_s_deg{0.0};
  /// C_b^s，行主序，把机体系向量转到准船体。
  double C_bs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  /// 四元数 wxyz，同样是机体到准船体。
  double quat_s_w{1.0};
  double quat_s_x{0.0};
  double quat_s_y{0.0};
  double quat_s_z{0.0};
  double latitude_deg{0.0};
  double longitude_deg{0.0};
  double height_m{0.0};
  /// 北东地速度，m/s。
  double v_n[3]{0.0, 0.0, 0.0};
  double roll_n_deg{0.0};
  double pitch_n_deg{0.0};
  double yaw_n_deg{0.0};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
