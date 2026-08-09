#include "solver/Solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <ompl/base/ScopedState.h>
#include <ompl/base/spaces/ReedsSheppStateSpace.h>
#include <random>
#include <vector>

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

  static State endpointErrors(Solver &solver, const State &from,
                              const State &to, const double radius) {
    double distance = 0.0;
    int condition = 0;
    int quadrant = 0;
    solver.acceleratedRSPlanner(from, to, distance, radius, condition,
                                quadrant);

    char motionTypes[5] = {'N', 'N', 'N', 'N', 'N'};
    solver.getMotionTypes(condition, motionTypes);
    solver.backProjectMotion(condition, quadrant, motionTypes);
    std::vector<State> path = {from};
    solver.getPath(condition, motionTypes, path, radius);

    State errors = {};
    solver.getErrors(path, to, errors);
    return errors;
  }
};

} // namespace accelerated

namespace {

constexpr std::size_t randomQueryCount = 100000;
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
  constexpr std::array<double, 4> radii = {0.5, 1.0, 10.0, 25.0};
  std::array<ob::StateSpacePtr, radii.size()> spaces;
  for (std::size_t i = 0; i < radii.size(); ++i) {
    spaces[i] = std::make_shared<ob::ReedsSheppStateSpace>(radii[i]);
  }

  std::size_t omplFailures = 0;
  std::size_t translationFailures = 0;
  double maximumOmplError = 0.0;
  double maximumTranslationError = 0.0;

  const State knownFrom = {159.776, 107.651, -0.346374};
  const State knownTo = {153.0, 51.0, -0.193722};
  const double knownAccelerated =
      SolverTestAccess::distance(solver, knownFrom, knownTo, 10.0);
  const double knownOmpl = omplDistance(spaces[2], knownFrom, knownTo);
  const State knownEndpointErrors =
      SolverTestAccess::endpointErrors(solver, knownFrom, knownTo, 10.0);
  const double knownEndpointErrorNorm =
      std::hypot(knownEndpointErrors.x, knownEndpointErrors.y);
  if (!nearlyEqual(knownAccelerated, knownOmpl)) {
    ++omplFailures;
    maximumOmplError = std::abs(knownAccelerated - knownOmpl);
  }
  const bool endpointValid = knownEndpointErrorNorm <= 1e-10 &&
                             std::abs(knownEndpointErrors.theta) <= 1e-10;

  std::mt19937_64 generator(0x5eed5eedULL);
  std::uniform_real_distribution<double> position(-250.0, 250.0);
  std::uniform_real_distribution<double> translation(-500.0, 500.0);
  std::uniform_real_distribution<double> heading(-std::numbers::pi,
                                                 std::numbers::pi);
  std::uniform_int_distribution<std::size_t> radiusIndex(0,
                                                         radii.size() - 1);

  for (std::size_t i = 0; i < randomQueryCount; ++i) {
    const State from = {position(generator), position(generator),
                        heading(generator)};
    const State to = {position(generator), position(generator),
                      heading(generator)};
    const double dx = translation(generator);
    const double dy = translation(generator);
    const State translatedFrom = {from.x + dx, from.y + dy, from.theta};
    const State translatedTo = {to.x + dx, to.y + dy, to.theta};
    const std::size_t selectedRadius = radiusIndex(generator);
    const double radius = radii[selectedRadius];

    const double accelerated =
        SolverTestAccess::distance(solver, from, to, radius);
    const double translated = SolverTestAccess::distance(
        solver, translatedFrom, translatedTo, radius);
    const double reference = omplDistance(spaces[selectedRadius], from, to);

    const double omplError = std::abs(accelerated - reference);
    const double translationError = std::abs(accelerated - translated);
    maximumOmplError = std::max(maximumOmplError, omplError);
    maximumTranslationError =
        std::max(maximumTranslationError, translationError);

    if (!nearlyEqual(accelerated, reference)) {
      ++omplFailures;
    }
    if (!nearlyEqual(accelerated, translated)) {
      ++translationFailures;
    }
  }

  std::cout << std::setprecision(12)
            << "Known query: accelerated=" << knownAccelerated
            << ", OMPL=" << knownOmpl
            << ", endpoint position error=" << knownEndpointErrorNorm
            << ", endpoint heading error=" << knownEndpointErrors.theta
            << '\n'
            << "Random queries: " << randomQueryCount << '\n'
            << "OMPL failures: " << omplFailures
            << ", maximum error=" << maximumOmplError << '\n'
            << "Translation failures: " << translationFailures
            << ", maximum error=" << maximumTranslationError << '\n';

  return omplFailures == 0 && translationFailures == 0 && endpointValid ? 0
                                                                        : 1;
}
