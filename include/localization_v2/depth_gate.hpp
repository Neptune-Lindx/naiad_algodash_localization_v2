#ifndef LOCALIZATION_V2_DEPTH_GATE_HPP_
#define LOCALIZATION_V2_DEPTH_GATE_HPP_

namespace naiad
{
namespace localization_v2
{

/// 用压力计做水面 / 水下滞回。传感器无效时保持上一状态。
class DepthGate
{
public:
  bool update(bool sensor_valid, double pressure_bar, double dive_on_bar, double dive_off_bar);
  bool submerged() const { return submerged_; }
  bool have_sample() const { return have_sample_; }

private:
  bool submerged_{false};
  bool have_sample_{false};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
