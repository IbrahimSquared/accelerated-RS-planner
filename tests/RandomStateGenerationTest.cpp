#include "helper/HelperFunctions.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <numbers>
#include <vector>

namespace {

bool equal(const accelerated::State &lhs, const accelerated::State &rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.theta == rhs.theta;
}

bool validBenchmarkState(const accelerated::State &state,
                         const double radius) {
  const double circleX = state.x + radius * std::cos(state.theta);
  const double circleY = state.y + radius * std::sin(state.theta);
  return state.x >= -5.0 && state.x <= 1.0 && state.y >= 0.0 &&
         state.y <= 5.0 && state.theta >= -std::numbers::pi &&
         state.theta <= std::numbers::pi && circleX <= 1.0 && circleY >= 0.0;
}

} // namespace

int main() {
  constexpr std::size_t sampleCount = 128;
  constexpr double radius = 1.0;
  constexpr std::uint32_t seed = 0x5eedU;

  const accelerated::HelperFunctions &helper =
      accelerated::HelperFunctions::getInstance();

  std::vector<accelerated::State> first = {{999.0, 999.0, 999.0}};
  first.reserve(sampleCount * 2);
  helper.generateRandomStates(first, sampleCount, false, radius, seed);

  std::vector<accelerated::State> repeated;
  helper.generateRandomStates(repeated, sampleCount, false, radius, seed);

  std::vector<accelerated::State> differentSeed;
  helper.generateRandomStates(differentSeed, sampleCount, false, radius,
                              seed + 1);

  std::vector<accelerated::State> empty = {{1.0, 2.0, 3.0}};
  helper.generateRandomStates(empty, 0, false, radius, seed);

  bool valid = first.size() == sampleCount &&
               repeated.size() == sampleCount &&
               differentSeed.size() == sampleCount && empty.empty();
  bool deterministic = valid;
  bool seedChangesOutput = false;
  for (std::size_t i = 0; i < first.size() && valid; ++i) {
    valid = validBenchmarkState(first[i], radius);
    deterministic = deterministic && equal(first[i], repeated[i]);
    seedChangesOutput = seedChangesOutput || !equal(first[i], differentSeed[i]);
  }

  std::cout << "Requested states: " << sampleCount
            << ", generated states: " << first.size() << '\n'
            << "All states satisfy the benchmark domain: " << std::boolalpha
            << valid << '\n'
            << "Repeated seed is deterministic: " << deterministic << '\n'
            << "Different seed changes output: " << seedChangesOutput << '\n'
            << "Zero-count request is empty: " << empty.empty() << '\n';

  return valid && deterministic && seedChangesOutput && empty.empty() ? 0 : 1;
}
