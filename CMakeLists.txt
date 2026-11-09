cmake_minimum_required(VERSION 3.10)

execute_process(
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/package_version.py ${CMAKE_CURRENT_SOURCE_DIR}/package.xml
  OUTPUT_VARIABLE PACKAGE_VERSION
)

project(naiad_localization_v2 VERSION ${PACKAGE_VERSION})

set(RULES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/rules)
include(${RULES_DIR}/ws_common.cmake)

if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 拷贝来的 ins-gnss 仍是文件后处理库。默认不编进 ROS 节点，避免在滤波尚未接上时链接 lapack 主程序。
option(NAIAD_BUILD_INS_GNSS "Build the copied ins-gnss static library" OFF)

if(WS_ROS)
  find_package(ament_cmake REQUIRED)
  find_package(rclcpp REQUIRED)
  find_package(std_msgs REQUIRED)
  find_package(geometry_msgs REQUIRED)
  find_package(naiad_interfaces REQUIRED PATHS ${WS_INSTALL_PREFIX})
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

if(WS_ROS)
  set(NODE_LOCALIZATION_V2 naiad_localization_v2)
  set(library_name ${NODE_LOCALIZATION_V2}_core)
  set(dependencies rclcpp std_msgs geometry_msgs naiad_interfaces)

  add_library(${library_name}
    src/mode/depth_gate.cpp
    src/mode/static_gate.cpp
    src/mode/mode_manager.cpp
    src/bridge/gnss_tc_bridge.cpp
    src/hull/hull_fusion.cpp
    src/hull/hull_ekf.cpp
    src/hull/wheel_speed.cpp
    src/ros/node_localization_v2.cpp
  )
  target_include_directories(${library_name} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  )
  ament_target_dependencies(${library_name} ${dependencies})

  add_executable(${NODE_LOCALIZATION_V2} src/ros/main.cpp)
  target_link_libraries(${NODE_LOCALIZATION_V2} ${library_name})
  ament_target_dependencies(${NODE_LOCALIZATION_V2} ${dependencies})

  add_executable(test_odo tools/test_odo.cpp)
  target_link_libraries(test_odo ${library_name})
  ament_target_dependencies(test_odo ${dependencies})

  if(NAIAD_BUILD_INS_GNSS)
    add_subdirectory(src/ins-gnss)
  endif()

  install(TARGETS ${library_name}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
  )
  install(TARGETS ${NODE_LOCALIZATION_V2} test_odo
    RUNTIME DESTINATION lib/${PROJECT_NAME}
  )
  install(DIRECTORY launch config
    DESTINATION share/${PROJECT_NAME}
  )
  install(DIRECTORY include/localization_v2
    DESTINATION include
  )

  ament_package()
else()
  message(STATUS "WS_ROS is off; naiad_localization_v2 node is not built")
endif()
