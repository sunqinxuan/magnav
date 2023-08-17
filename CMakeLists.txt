cmake_minimum_required(VERSION 2.8.3)
project(magnav)

SET(CMAKE_BUILD_TYPE "Debug") #Debug/Release
MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -fPIC -O3  -Wno-error=deprecated-declarations -Wno-deprecated-declarations -g") #-g

add_compile_options(-std=c++14)
add_definitions(-std=c++14)
set(CMAKE_CXX_STANDARD 14) 

INCLUDE_DIRECTORIES(src)
link_directories(${PROJECT_SOURCE_DIR}/lib) # ${PROJECT_SOURCE_DIR}/3rdparty/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})
#set(HDF5LIBS hdf5 hdf5_cpp)
MESSAGE("HDF5_INCLUDE_DIR: " ${HDF5_INCLUDE_DIR})
link_directories(/usr/lib/x86_64-linux-gnu/hdf5/serial/)

add_subdirectory(src/mag_compensation)
add_subdirectory(src/filter)

add_executable(magnav src/magnav.cpp)
target_link_libraries(magnav mag_compensation filter hdf5 hdf5_cpp)
 


