#include "localization_v2/gnss_tc_bridge.hpp"

namespace naiad
{
namespace localization_v2
{

void GnssTcBridge::reset() { fast_search_armed_ = false; }

void GnssTcBridge::arm_fast_search() { fast_search_armed_ = true; }

bool GnssTcBridge::start_quasi_static_align() { return false; }

NavStatus GnssTcBridge::step_align() const
{
  NavStatus status;
  status.backend = "align";
  status.fast_search_armed = fast_search_armed_;
  return status;
}

NavStatus GnssTcBridge::step_coast() const
{
  NavStatus status;
  status.backend = "gnss_coast";
  status.fast_search_armed = fast_search_armed_;
  return status;
}

NavStatus GnssTcBridge::step_tight() const
{
  NavStatus status;
  status.backend = "gnss_tc";
  status.fast_search_armed = fast_search_armed_;
  return status;
}

}  // namespace localization_v2
}  // namespace naiad
