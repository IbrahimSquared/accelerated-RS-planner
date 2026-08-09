#include "solver/Solver.hpp"

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <ompl/base/ScopedState.h>
#include <ompl/base/spaces/ReedsSheppStateSpace.h>
#include <random>

namespace ob = ompl::base;

namespace accelerated {

class SolverTestAccess {
public:
  static double distance(Solver &solver, const State &from, const State &to,
                         const double radius) {
    double distance = 0.0;
    int condition = 0;
    int quadrant = 0;
    solver.acceleratedRSPlanner(from, to, distance, radius, condition,
                                quadrant);
    return distance;
  }
};

} // namespace accelerated

namespace {

constexpr std::size_t randomQueryCount = 200000;
constexpr double relativeTolerance = 1e-10;

bool nearlyEqual(const double lhs, const double rhs) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= relativeTolerance * scale;
}

double omplDistance(const ob::StateSpacePtr &space,
                    const accelerated::State &from,
                    const accelerated::State &to) {
  ob::ScopedState<> fromState(space);
  ob::ScopedState<> toState(space);
  fromState[0] = from.x;
  fromState[1] = from.y;
  fromState[2] = from.theta;
  toState[0] = to.x;
  toState[1] = to.y;
  toState[2] = to.theta;
  return space->distance(fromState(), toState());
}

} // namespace

int main() {
  using accelerated::SolverTestAccess;
  using accelerated::State;

  accelerated::Solver solver;
  const State from = {0.0, 0.0, 0.0};
  constexpr double radius = 1.0;
  const ob::StateSpacePtr space =
      std::make_shared<ob::ReedsSheppStateSpace>(radius);

  std::mt19937_64 generator(0xd0a1a5eedULL);
  std::uniform_real_distribution<double> position(-5.0, 5.0);
  std::uniform_real_distribution<double> heading(-std::numbers::pi,
                                                 std::numbers::pi);

  std::size_t floatingPointFailures = 0;
  std::size_t nonfiniteOutputs = 0;
  std::size_t omplFailures = 0;
  double maximumOmplError = 0.0;

  std::size_t boundaryFloatingPointFailures = 0;
  std::size_t boundaryOmplFailures = 0;
  constexpr double boundaryHeading = 0.5;
  for (int step = -1000; step <= 1000; ++step) {
    const double rl = 2.0 * radius + step * 1e-9;
    const State to = {std::sin(boundaryHeading),
                      -radius + rl - std::cos(boundaryHeading),
                      boundaryHeading};

    std::feclearexcept(FE_ALL_EXCEPT);
    const double accelerated =
        SolverTestAccess::distance(solver, from, to, radius);
    if (std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) != 0) {
      ++boundaryFloatingPointFailures;
    }
    const double reference = omplDistance(space, from, to);
    if (!std::isfinite(accelerated) || !nearlyEqual(accelerated, reference)) {
      ++boundaryOmplFailures;
    }
  }

  for (std::size_t i = 0; i < randomQueryCount; ++i) {
    const State to = {position(generator), position(generator),
                      heading(generator)};

    std::feclearexcept(FE_ALL_EXCEPT);
    const double accelerated =
        SolverTestAccess::distance(solver, from, to, radius);
    const int exceptions =
        std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    if (exceptions != 0) {
      ++floatingPointFailures;
    }
    if (!std::isfinite(accelerated)) {
      ++nonfiniteOutputs;
    }

    const double reference = omplDistance(space, from, to);
    const double error = std::abs(accelerated - reference);
    maximumOmplError = std::max(maximumOmplError, error);
    if (!nearlyEqual(accelerated, reference)) {
      ++omplFailures;
    }
  }

  std::cout << std::setprecision(12)
            << "Random queries: " << randomQueryCount << '\n'
            << "RL=2r boundary floating-point failures: "
            << boundaryFloatingPointFailures << '\n'
            << "RL=2r boundary OMPL failures: " << boundaryOmplFailures
            << '\n'
            << "Floating-point domain failures: " << floatingPointFailures
            << '\n'
            << "Nonfinite outputs: " << nonfiniteOutputs << '\n'
            << "OMPL failures: " << omplFailures
            << ", maximum error=" << maximumOmplError << '\n';

  return boundaryFloatingPointFailures == 0 && boundaryOmplFailures == 0 &&
                 floatingPointFailures == 0 && nonfiniteOutputs == 0 &&
                 omplFailures == 0
             ? 0
             : 1;
}
