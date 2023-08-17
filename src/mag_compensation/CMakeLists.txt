cmake_minimum_required(VERSION 2.8.3)
project(mag_compensation)

#find_package(yaml-cpp REQUIRED)
#find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
#include_directories(${EIGEN3_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS}) 

list(APPEND LIB_MAGCMP filter)

file(GLOB_RECURSE SRC_MAGCMP ./*.cpp)
add_library(mag_compensation SHARED ${SRC_MAGCMP})
target_link_libraries(mag_compensation ${LIB_MAGCMP})

