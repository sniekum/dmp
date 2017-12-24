cmake_minimum_required(VERSION 2.8.3)
project(dmp)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(catkin REQUIRED COMPONENTS
cmake_modules
roscpp
rospy
message_generation
)
find_package(Eigen REQUIRED)

add_message_files(
   FILES
   DMPData.msg
   DMPPoint.msg
   DMPTraj.msg
 )

add_service_files(
   FILES
   GetDMPPlan.srv
   LearnDMPFromDemo.srv
   SetActiveDMP.srv
 )

generate_messages(
   DEPENDENCIES
)

catkin_package(
INCLUDE_DIRS include
CATKIN_DEPENDS common_msgs image_common message_runtime
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
)

add_library(dmp src/dmp.cpp src/fourier_approx.cpp src/radial_approx.cpp src/linear_approx.cpp)
add_executable(dmp_server nodes/dmp_server.cpp)
add_dependencies(dmp
                 ${catkin_EXPORTED_TARGETS}
                 dmp_gencpp
                 dmp_generate_messages
                 dmp_generate_messages_cpp
                 dmp_generate_messages_py
                 dmp_generate_messages_lisp)
target_link_libraries(dmp_server dmp ${catkin_LIBRARIES} ${Eigen_LIBRARIES})
