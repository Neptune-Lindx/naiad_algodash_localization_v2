#include "localization_v2/static_gate.hpp"

namespace naiad
{
namespace localization_v2
{

bool StaticGate::update(double t_sec, bool sample_ok, bool is_still, double hold_sec)
{
  if (ready_ || !sample_ok || !is_still) {
    holding_ = false;
    return ready_;
  }
  if (!holding_) {
    holding_ = true;
    hold_start_sec_ = t_sec;
  }
  if ((t_sec - hold_start_sec_) >= hold_sec) {
    ready_ = true;
    holding_ = false;
  }
  return ready_;
}

}  // namespace localization_v2
}  // namespace naiad
