#include "localization_v2/mode_manager.hpp"

namespace naiad
{
namespace localization_v2
{

const char * work_mode_name(WorkMode mode)
{
  switch (mode) {
    case WorkMode::Align:
      return "align";
    case WorkMode::SurfaceTight:
      return "surface_tight";
    case WorkMode::UnderwaterHull:
      return "underwater_hull";
    case WorkMode::ResurfaceSearch:
      return "resurface_search";
  }
  return "unknown";
}

WorkMode ModeManager::update(const ModeInput & in)
{
  switch (mode_) {
    case WorkMode::Align:
      if (in.submerged) {
        mode_ = WorkMode::UnderwaterHull;
      } else if (in.align_ready && in.rtk_fixed) {
        mode_ = WorkMode::SurfaceTight;
      } else if (in.align_ready) {
        mode_ = WorkMode::ResurfaceSearch;
      }
      break;
    case WorkMode::SurfaceTight:
      if (in.submerged) {
        mode_ = WorkMode::UnderwaterHull;
      } else if (!in.rtk_fixed) {
        mode_ = WorkMode::ResurfaceSearch;
      }
      break;
    case WorkMode::UnderwaterHull:
      if (!in.submerged) {
        mode_ = WorkMode::ResurfaceSearch;
      }
      break;
    case WorkMode::ResurfaceSearch:
      if (in.submerged) {
        mode_ = WorkMode::UnderwaterHull;
      } else if (!in.align_ready) {
        mode_ = WorkMode::Align;
      } else if (in.rtk_fixed) {
        mode_ = WorkMode::SurfaceTight;
      }
      break;
  }
  return mode_;
}

}  // namespace localization_v2
}  // namespace naiad
