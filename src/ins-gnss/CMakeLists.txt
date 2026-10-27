# 可选静态库。由包根目录 -DNAIAD_BUILD_INS_GNSS=ON 打开。
# 不包含 app/：那些文件自带 main，是文件后处理程序，输入输出已改由 ROS2 节点承担。
file(GLOB INS_GNSS_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rcv/*.cc
)

add_library(ins_gnss STATIC ${INS_GNSS_SRCS})
target_include_directories(ins_gnss PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/rcv
  ${PROJECT_SOURCE_DIR}/include/ins-gnss
)
target_link_libraries(ins_gnss PUBLIC pthread)
