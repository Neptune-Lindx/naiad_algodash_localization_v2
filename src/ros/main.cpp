#include "localization_v2/node_localization_v2.hpp"

int main(int argc, char ** argv)
{
  rclcpp::init(argc, argv);
  auto node = std::make_shared<naiad::localization_v2::NodeLocalizationV2>("naiad_localization_v2");
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
