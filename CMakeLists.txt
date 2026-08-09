cmake_minimum_required(VERSION 3.10)

project(acceleratedRSPlanner)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the build type to Release by default
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -Wpedantic -Wextra -Wnull-dereference -fexpensive-optimizations -mfpmath=sse -fno-trapping-math -fno-math-errno -flto=auto -std=c++20 -D_USE_MATH_DEFINES")

# Set the default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

include_directories(include)
add_executable(acceleratedRSPlanner src/main.cpp src/ConfigParser.cpp src/HelperFunctions.cpp src/Solver.cpp src/DesaulniersAlgorithm.cpp src/Plotter.cpp)

find_package(ompl)
target_link_libraries(acceleratedRSPlanner PRIVATE ompl)
find_package(SFML 2.5 COMPONENTS graphics REQUIRED)
target_link_libraries(acceleratedRSPlanner PRIVATE sfml-graphics)

include(CTest)
if(BUILD_TESTING)
    add_executable(
        translationInvarianceTest
        tests/TranslationInvarianceTest.cpp
        src/ConfigParser.cpp
        src/HelperFunctions.cpp
        src/Solver.cpp
        src/DesaulniersAlgorithm.cpp
        src/Plotter.cpp
    )
    target_link_libraries(translationInvarianceTest PRIVATE ompl sfml-graphics)
    set_target_properties(
        translationInvarianceTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME translation_invariance COMMAND translationInvarianceTest)
    set_tests_properties(
        translation_invariance
        PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(
        randomStateGenerationTest
        tests/RandomStateGenerationTest.cpp
        src/ConfigParser.cpp
        src/HelperFunctions.cpp
    )
    target_compile_definitions(randomStateGenerationTest PRIVATE _GLIBCXX_ASSERTIONS)
    set_target_properties(
        randomStateGenerationTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME random_state_generation COMMAND randomStateGenerationTest)
    set_tests_properties(
        random_state_generation
        PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()
