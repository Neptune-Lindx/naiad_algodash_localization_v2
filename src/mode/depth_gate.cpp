#include "localization_v2/depth_gate.hpp"

namespace naiad
{
namespace localization_v2
{

bool DepthGate::update(
  bool sensor_valid, double pressure_bar, double dive_on_bar, double dive_off_bar)
{
  if (!sensor_valid) {
    return submerged_;
  }
  have_sample_ = true;
  if (!submerged_ && pressure_bar >= dive_on_bar) {
    submerged_ = true;
  } else if (submerged_ && pressure_bar <= dive_off_bar) {
    submerged_ = false;
  }
  return submerged_;
}

}  // namespace localization_v2
}  // namespace naiad
