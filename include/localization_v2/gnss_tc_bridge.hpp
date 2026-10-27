#ifndef LOCALIZATION_V2_GNSS_TC_BRIDGE_HPP_
#define LOCALIZATION_V2_GNSS_TC_BRIDGE_HPP_

#include "localization_v2/nav_status.hpp"

namespace naiad
{
namespace localization_v2
{

/// RTK/IMU 紧组合入口。
/// 实现落在 src/ins-gnss/ins-gnss-tc.cc 与 tcpostpos.cc，本类只保留调用边界。
/// 出水搜星时 arm_fast_search() 置位，后续接入时保持星历与模糊度状态，不做冷启动。
class GnssTcBridge
{
public:
  void reset();
  void arm_fast_search();
  /// 有 RTK 时的 IMU 准静态对准入口。滤波方程尚未接入，返回 false。
  bool start_quasi_static_align();
  NavStatus step_align() const;
  NavStatus step_coast() const;
  NavStatus step_tight() const;
  bool fast_search_armed() const { return fast_search_armed_; }

private:
  bool fast_search_armed_{false};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
