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
find_package(Threads REQUIRED)

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

    add_executable(
        inputValidationTest
        tests/InputValidationTest.cpp
        src/ConfigParser.cpp
        src/HelperFunctions.cpp
        src/Solver.cpp
        src/DesaulniersAlgorithm.cpp
        src/Plotter.cpp
    )
    target_link_libraries(inputValidationTest PRIVATE ompl sfml-graphics)
    set_target_properties(
        inputValidationTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME input_validation COMMAND inputValidationTest)
    set_tests_properties(
        input_validation
        PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(
        analyticDomainTest
        tests/AnalyticDomainTest.cpp
        src/ConfigParser.cpp
        src/HelperFunctions.cpp
        src/Solver.cpp
        src/DesaulniersAlgorithm.cpp
        src/Plotter.cpp
    )
    target_link_libraries(analyticDomainTest PRIVATE ompl sfml-graphics)
    set_target_properties(
        analyticDomainTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME analytic_domains COMMAND analyticDomainTest)
    set_tests_properties(
        analytic_domains
        PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(
        variabilityCharacterizationTest
        tests/VariabilityCharacterizationTest.cpp
        src/ConfigParser.cpp
        src/HelperFunctions.cpp
        src/Solver.cpp
        src/DesaulniersAlgorithm.cpp
        src/Plotter.cpp
    )
    target_link_libraries(variabilityCharacterizationTest PRIVATE ompl sfml-graphics)
    target_compile_definitions(
        variabilityCharacterizationTest
        PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
    )
    set_target_properties(
        variabilityCharacterizationTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME variability_characterization COMMAND variabilityCharacterizationTest)
    set_tests_properties(
        variability_characterization
        PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_executable(
        reconstructionTest
        tests/ReconstructionTest.cpp
        src/ConfigParser.cpp
        src/HelperFunctions.cpp
        src/Solver.cpp
        src/DesaulniersAlgorithm.cpp
        src/Plotter.cpp
    )
    target_link_libraries(reconstructionTest PRIVATE ompl sfml-graphics)
    set_target_properties(
        reconstructionTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME reconstruction_idempotence COMMAND reconstructionTest)
    set_tests_properties(
        reconstruction_idempotence
        PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(
        publicApiTest
        tests/PublicApiTest.cpp
        src/ConfigParser.cpp
        src/HelperFunctions.cpp
        src/Solver.cpp
        src/DesaulniersAlgorithm.cpp
        src/Plotter.cpp
    )
    target_link_libraries(
        publicApiTest PRIVATE ompl sfml-graphics Threads::Threads
    )
    set_target_properties(
        publicApiTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME public_api COMMAND publicApiTest)
    set_tests_properties(
        public_api
        PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(
        benchmarkTest
        tests/BenchmarkTest.cpp
        src/ConfigParser.cpp
        src/HelperFunctions.cpp
        src/Solver.cpp
        src/DesaulniersAlgorithm.cpp
        src/Plotter.cpp
    )
    target_link_libraries(benchmarkTest PRIVATE ompl sfml-graphics)
    set_target_properties(
        benchmarkTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME benchmark_methodology COMMAND benchmarkTest)
    set_tests_properties(
        benchmark_methodology
        PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(
        configParserTest
        tests/ConfigParserTest.cpp
        src/ConfigParser.cpp
    )
    set_target_properties(
        configParserTest
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME config_parser COMMAND configParserTest)
    set_tests_properties(
        config_parser
        PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endif()
