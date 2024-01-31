cmake_minimum_required(VERSION 3.10)

project(accelerateRSPlanner)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the build type to Release by default
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto=auto -std=c++20 -D_USE_MATH_DEFINES")

# Set the default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

include_directories(include)
add_executable(accelerateRSPlanner src/main.cpp src/ConfigParser.cpp src/HelperFunctions.cpp src/Solver.cpp src/DesaulniersAlgorithm.cpp src/Plotter.cpp)

find_package(ompl)
target_link_libraries(accelerateRSPlanner PRIVATE ompl)
find_package(SFML 2.5 COMPONENTS graphics REQUIRED)
target_link_libraries(accelerateRSPlanner PRIVATE sfml-graphics)
