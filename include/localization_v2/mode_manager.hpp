#ifndef LOCALIZATION_V2_MODE_MANAGER_HPP_
#define LOCALIZATION_V2_MODE_MANAGER_HPP_

#include "localization_v2/work_mode.hpp"

namespace naiad
{
namespace localization_v2
{

struct ModeInput
{
  bool submerged{false};
  bool rtk_fixed{false};
  bool align_ready{false};
};

/// 全自动相位切换。水下优先于 RTK：没出水面不进紧组合。
class ModeManager
{
public:
  WorkMode update(const ModeInput & in);
  WorkMode mode() const { return mode_; }

private:
  WorkMode mode_{WorkMode::Align};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
