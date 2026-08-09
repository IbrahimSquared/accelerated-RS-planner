#include "solver/Solver.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>

int main(const int argc, char *argv[]) {
  std::size_t sampleCount = 100000;
  if (argc == 2) {
    sampleCount = std::stoull(argv[1]);
  }
  constexpr std::uint32_t seed = 0x5eed5eedU;
  constexpr std::size_t repetitions = 3;

  accelerated::Solver solver;
  const accelerated::BenchmarkResult result =
      solver.benchmarkPlanners(sampleCount, seed, repetitions);

  const bool metadataValid = result.sample_count == sampleCount &&
                             result.repetitions == repetitions &&
                             result.seed == seed;
  const bool timingsValid = result.accelerated_microseconds > 0.0 &&
                            result.ompl_microseconds > 0.0 &&
                            result.desaulniers_microseconds > 0.0 &&
                            result.accelerated_trials_microseconds.size() ==
                                repetitions &&
                            result.ompl_trials_microseconds.size() ==
                                repetitions &&
                            result.desaulniers_trials_microseconds.size() ==
                                repetitions;
  const bool errorsValid =
      std::isfinite(result.accelerated_average_error) &&
      std::isfinite(result.desaulniers_average_error) &&
      result.accelerated_maximum_error <= 1e-10 &&
      result.desaulniers_maximum_error <= 1e-10;

  std::cout << std::boolalpha
            << "Benchmark metadata is valid: " << metadataValid << '\n'
            << "Median timings are positive: " << timingsValid << '\n'
            << "Streaming error statistics are valid: " << errorsValid
            << '\n';
  return metadataValid && timingsValid && errorsValid ? 0 : 1;
}
