#ifndef LOCALIZATION_V2_WORK_MODE_HPP_
#define LOCALIZATION_V2_WORK_MODE_HPP_

namespace naiad
{
namespace localization_v2
{

/// 水面 / 水下自动切换的四种工作相位。
/// 滤波方程尚未接入；模式机只决定下一拍该调用哪条算法入口。
enum class WorkMode
{
  Align,             ///< 水面准静态初始对准
  SurfaceTight,      ///< 水面，RTK 固定，RTK/IMU 紧组合
  UnderwaterHull,    ///< 水下无 RTK，船体坐标系 IMU + 轮速
  ResurfaceSearch    ///< 出水或水面失锁，保持导航并快速搜星
};

const char * work_mode_name(WorkMode mode);

}  // namespace localization_v2
}  // namespace naiad

#endif
