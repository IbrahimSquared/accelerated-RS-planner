#include "solver/Solver.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <random>
#include <vector>

namespace accelerated {

struct TestSolution {
  int condition;
  int quadrant;
};

class SolverTestAccess {
public:
  static TestSolution solve(Solver &solver, const State &from, const State &to,
                            const double radius) {
    double distance = 0.0;
    TestSolution result = {};
    solver.acceleratedRSPlanner(from, to, distance, radius, result.condition,
                                result.quadrant);
    return result;
  }

  static std::vector<State> reconstruct(Solver &solver, const State &from,
                                        const TestSolution &solution,
                                        const double radius) {
    char motionTypes[5] = {'N', 'N', 'N', 'N', 'N'};
    double motionLengths[5] = {};
    solver.getMotionTypes(solution.condition, motionTypes);
    solver.backProjectMotion(solution.condition, solution.quadrant, motionTypes,
                             motionLengths);
    std::vector<State> path = {from};
    solver.getPath(motionTypes, motionLengths, path, radius);
    return path;
  }

  static std::array<double, 5> canonicalLengths(
      const Solver &solver, const int condition) {
    std::array<double, 5> result = {};
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = solver.lengths_[condition - 1][i];
    }
    return result;
  }

  static State errors(Solver &solver, const std::vector<State> &path,
                      const State &to) {
    State result = {};
    solver.getErrors(path, to, result);
    return result;
  }
};

} // namespace accelerated

namespace {

bool samePath(const std::vector<accelerated::State> &lhs,
              const std::vector<accelerated::State> &rhs) {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i].x != rhs[i].x || lhs[i].y != rhs[i].y ||
        lhs[i].theta != rhs[i].theta) {
      return false;
    }
  }
  return true;
}

bool reachesTarget(const accelerated::State &errors) {
  return std::hypot(errors.x, errors.y) <= 1e-10 && errors.theta <= 1e-10;
}

} // namespace

int main() {
  using accelerated::SolverTestAccess;
  using accelerated::State;

  accelerated::Solver solver;
  const State from = {0.0, 0.0, 0.0};
  const State quadrantTwoTarget = {-4.0, 4.0, 0.5};
  constexpr double radius = 1.0;

  const accelerated::TestSolution quadrantTwoSolution =
      SolverTestAccess::solve(solver, from, quadrantTwoTarget, radius);
  const std::array<double, 5> originalLengths =
      SolverTestAccess::canonicalLengths(solver,
                                         quadrantTwoSolution.condition);
  const std::vector<State> firstPath = SolverTestAccess::reconstruct(
      solver, from, quadrantTwoSolution, radius);
  const std::vector<State> secondPath = SolverTestAccess::reconstruct(
      solver, from, quadrantTwoSolution, radius);
  const std::array<double, 5> lengthsAfterReconstruction =
      SolverTestAccess::canonicalLengths(solver,
                                         quadrantTwoSolution.condition);

  const bool quadrantIsMirrored = quadrantTwoSolution.quadrant == 2;
  const bool reconstructionIsIdempotent = samePath(firstPath, secondPath);
  const bool canonicalLengthsUnchanged =
      originalLengths == lengthsAfterReconstruction;
  const bool mirroredPathIsValid = reachesTarget(
      SolverTestAccess::errors(solver, firstPath, quadrantTwoTarget));

  std::mt19937_64 generator(0x5eedU);
  std::uniform_real_distribution<double> position(-5.0, 5.0);
  std::uniform_real_distribution<double> heading(-std::numbers::pi,
                                                 std::numbers::pi);
  bool foundFiveSegmentPath = false;
  bool fiveSegmentPathIsValid = false;
  for (std::size_t attempt = 0; attempt < 100000; ++attempt) {
    const State to = {position(generator), position(generator),
                      heading(generator)};
    const accelerated::TestSolution solution =
        SolverTestAccess::solve(solver, from, to, radius);
    if (solution.condition != 9 && solution.condition != 12) {
      continue;
    }

    const std::vector<State> path =
        SolverTestAccess::reconstruct(solver, from, solution, radius);
    foundFiveSegmentPath = path.size() == 6;
    fiveSegmentPathIsValid =
        reachesTarget(SolverTestAccess::errors(solver, path, to));
    break;
  }

  const bool valid = quadrantIsMirrored && reconstructionIsIdempotent &&
                     canonicalLengthsUnchanged && mirroredPathIsValid &&
                     foundFiveSegmentPath && fiveSegmentPathIsValid;
  std::cout << std::boolalpha
            << "Mirrored quadrant exercised: " << quadrantIsMirrored << '\n'
            << "Repeated reconstruction is identical: "
            << reconstructionIsIdempotent << '\n'
            << "Canonical lengths remain unchanged: "
            << canonicalLengthsUnchanged << '\n'
            << "Mirrored path reaches target: " << mirroredPathIsValid << '\n'
            << "Five-segment path contains six states: "
            << foundFiveSegmentPath << '\n'
            << "Five-segment path reaches target: " << fiveSegmentPathIsValid
            << '\n';
  return valid ? 0 : 1;
}
