cmake_minimum_required(VERSION 2.8.3)
project(dmp)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(catkin REQUIRED COMPONENTS
eigen
roscpp
rospy
)


catkin_package(
INCLUDE_DIRS include
CATKIN_DEPENDS common_msgs image_common
)


include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
)

add_library(dmp src/dmp.cpp src/fourier_approx.cpp src/radial_approx.cpp)
add_executable(dmp_server nodes/dmp_server.cpp)
target_link_libraries(dmp_server dmp)
