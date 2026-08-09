#include "solver/Solver.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

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

  static State errors(Solver &solver, const State &reconstructed,
                      const State &target) {
    const std::vector<State> path = {reconstructed};
    State result = {};
    solver.getErrors(path, target, result);
    return result;
  }
};

} // namespace accelerated

namespace {

template <typename Function> bool throwsInvalidArgument(Function &&function) {
  try {
    std::forward<Function>(function)();
  } catch (const std::invalid_argument &) {
    return true;
  } catch (...) {
    return false;
  }
  return false;
}

} // namespace

int main() {
  using accelerated::SolverTestAccess;
  using accelerated::State;

  accelerated::Solver solver;
  const State origin = {0.0, 0.0, 0.0};
  const State goal = {1.0, 2.0, 0.5};
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();

  const bool rejectsZeroRadius = throwsInvalidArgument(
      [&] { SolverTestAccess::distance(solver, origin, goal, 0.0); });
  const bool rejectsNegativeRadius = throwsInvalidArgument(
      [&] { SolverTestAccess::distance(solver, origin, goal, -1.0); });
  const bool rejectsNanRadius = throwsInvalidArgument(
      [&] { SolverTestAccess::distance(solver, origin, goal, nan); });
  const bool rejectsInfiniteRadius = throwsInvalidArgument(
      [&] { SolverTestAccess::distance(solver, origin, goal, infinity); });
  const bool publicApiRejectsZeroRadius =
      throwsInvalidArgument([&] { solver.evaluateQuery(origin, goal, 0.0); });
  const bool rejectsNanStart = throwsInvalidArgument([&] {
    SolverTestAccess::distance(solver, {nan, 0.0, 0.0}, goal, 1.0);
  });
  const bool rejectsInfiniteStart = throwsInvalidArgument([&] {
    SolverTestAccess::distance(solver, {0.0, infinity, 0.0}, goal, 1.0);
  });
  const bool rejectsNanGoal = throwsInvalidArgument([&] {
    SolverTestAccess::distance(solver, origin, {1.0, 2.0, nan}, 1.0);
  });

  const State equivalentAngleErrors = SolverTestAccess::errors(
      solver, {0.0, 0.0, -std::numbers::pi},
      {0.0, 0.0, std::numbers::pi});
  const bool equivalentAnglesMatch = equivalentAngleErrors.theta <= 1e-12;

  const bool valid = rejectsZeroRadius && rejectsNegativeRadius &&
                     rejectsNanRadius && rejectsInfiniteRadius &&
                     publicApiRejectsZeroRadius &&
                     rejectsNanStart && rejectsInfiniteStart &&
                     rejectsNanGoal && equivalentAnglesMatch;

  std::cout << std::boolalpha
            << "Rejects zero radius: " << rejectsZeroRadius << '\n'
            << "Rejects negative radius: " << rejectsNegativeRadius << '\n'
            << "Rejects NaN radius: " << rejectsNanRadius << '\n'
            << "Rejects infinite radius: " << rejectsInfiniteRadius << '\n'
            << "Public API rejects zero radius: "
            << publicApiRejectsZeroRadius << '\n'
            << "Rejects NaN start state: " << rejectsNanStart << '\n'
            << "Rejects infinite start state: " << rejectsInfiniteStart << '\n'
            << "Rejects NaN goal state: " << rejectsNanGoal << '\n'
            << "Equivalent -pi/+pi headings match: " << equivalentAnglesMatch
            << " (error=" << equivalentAngleErrors.theta << ")\n";

  return valid ? 0 : 1;
}
