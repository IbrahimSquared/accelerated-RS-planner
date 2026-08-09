#include "solver/Solver.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numbers>
#include <ompl/base/ScopedState.h>
#include <ompl/base/spaces/ReedsSheppStateSpace.h>
#include <random>
#include <thread>
#include <vector>

namespace ob = ompl::base;

namespace {

struct Query {
  accelerated::State from;
  accelerated::State to;
  double expected_distance;
};

bool nearlyEqual(const double lhs, const double rhs) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= 1e-10 * scale;
}

bool validResult(const accelerated::SolveResult &result, const Query &query) {
  if (!nearlyEqual(result.distance, query.expected_distance) ||
      result.condition < 1 || result.condition > 20 || result.quadrant < 1 ||
      result.quadrant > 4 || result.path.empty()) {
    return false;
  }

  std::size_t segmentCount = 0;
  while (segmentCount < result.motion_types.size() &&
         result.motion_types[segmentCount] != 'N') {
    ++segmentCount;
  }
  if (result.path.size() != segmentCount + 1) {
    return false;
  }

  const accelerated::State &endpoint = result.path.back();
  const double positionError =
      std::hypot(endpoint.x - query.to.x, endpoint.y - query.to.y);
  const double headingError = std::remainder(endpoint.theta - query.to.theta,
                                             2.0 * std::numbers::pi);
  return positionError <= 1e-9 && std::abs(headingError) <= 1e-9;
}

} // namespace

int main() {
  constexpr std::size_t queryCount = 256;
  constexpr std::size_t threadCount = 4;
  constexpr std::size_t repetitions = 10;
  constexpr double radius = 10.0;

  const ob::StateSpacePtr space =
      std::make_shared<ob::ReedsSheppStateSpace>(radius);
  std::mt19937_64 generator(0xa915eedULL);
  std::uniform_real_distribution<double> position(-250.0, 250.0);
  std::uniform_real_distribution<double> heading(-std::numbers::pi,
                                                 std::numbers::pi);

  std::vector<Query> queries;
  queries.reserve(queryCount);
  for (std::size_t i = 0; i < queryCount; ++i) {
    Query query = {{position(generator), position(generator), heading(generator)},
                   {position(generator), position(generator), heading(generator)},
                   0.0};
    ob::ScopedState<> fromState(space);
    ob::ScopedState<> toState(space);
    fromState[0] = query.from.x;
    fromState[1] = query.from.y;
    fromState[2] = query.from.theta;
    toState[0] = query.to.x;
    toState[1] = query.to.y;
    toState[2] = query.to.theta;
    query.expected_distance = space->distance(fromState(), toState());
    queries.push_back(query);
  }

  bool sequentialValid = true;
  for (const Query &query : queries) {
    sequentialValid = sequentialValid &&
                      validResult(accelerated::Solver::solve(
                                      query.from, query.to, radius),
                                  query);
  }

  std::atomic<bool> concurrentValid = true;
  std::vector<std::thread> workers;
  workers.reserve(threadCount);
  for (std::size_t thread = 0; thread < threadCount; ++thread) {
    workers.emplace_back([&] {
      for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
        for (const Query &query : queries) {
          if (!validResult(
                  accelerated::Solver::solve(query.from, query.to, radius),
                  query)) {
            concurrentValid = false;
            return;
          }
        }
      }
    });
  }
  for (std::thread &worker : workers) {
    worker.join();
  }

  std::cout << std::boolalpha
            << "Sequential structured results are valid: " << sequentialValid
            << '\n'
            << "Concurrent structured results are valid: "
            << concurrentValid.load() << '\n'
            << "Concurrent queries checked: "
            << threadCount * repetitions * queryCount << '\n';
  return sequentialValid && concurrentValid ? 0 : 1;
}
