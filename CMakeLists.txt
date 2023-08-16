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

add_subdirectory(src/mag_compensation)

add_executable(magnav src/magnav.cpp)
#target_link_libraries(magnav mag_compensation)

