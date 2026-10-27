#ifndef LOCALIZATION_V2_STATIC_GATE_HPP_
#define LOCALIZATION_V2_STATIC_GATE_HPP_

namespace naiad
{
namespace localization_v2
{

/// 准静态保持计时。连续静止达到 hold_sec 后锁存 align_ready。
class StaticGate
{
public:
  bool update(double t_sec, bool sample_ok, bool is_still, double hold_sec);
  bool ready() const { return ready_; }

private:
  bool ready_{false};
  bool holding_{false};
  double hold_start_sec_{0.0};
};

}  // namespace localization_v2
}  // namespace naiad

#endif
