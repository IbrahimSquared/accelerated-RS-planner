#include "helper/HelperFunctions.hpp"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

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
  const accelerated::HelperFunctions &helper =
      accelerated::HelperFunctions::getInstance();
  const std::vector<double> reference = {1.0, 2.0, 3.0};
  const std::vector<double> candidate = {1.0, 2.0, 3.0};
  const std::vector<double> shortCandidate = {1.0, 2.0};
  const std::vector<double> nonfinite = {
      1.0, std::numeric_limits<double>::quiet_NaN(), 3.0};

  const bool emptyRejected = throwsInvalidArgument(
      [&] { helper.errorStats(std::vector<double>{}, std::vector<double>{}); });
  const bool mismatchRejected = throwsInvalidArgument(
      [&] { helper.errorStats(shortCandidate, reference); });
  const bool threeWayMismatchRejected = throwsInvalidArgument(
      [&] { helper.errorStats(candidate, reference, shortCandidate); });
  const bool nonfiniteRejected = throwsInvalidArgument(
      [&] { helper.errorStats(nonfinite, reference); });
  const bool zeroDurationRejected =
      throwsInvalidArgument([&] { helper.timeStats(0, 10, 10, 3); });
  const bool zeroCountRejected =
      throwsInvalidArgument([&] { helper.timeStats(10, 10, 10, 0); });

  bool validInputsAccepted = true;
  try {
    helper.errorStats(candidate, reference);
    helper.errorStats(candidate, reference, candidate);
    helper.timeStats(10, 20, 15, 3);
  } catch (...) {
    validInputsAccepted = false;
  }

  const bool valid = emptyRejected && mismatchRejected &&
                     threeWayMismatchRejected && nonfiniteRejected &&
                     zeroDurationRejected && zeroCountRejected &&
                     validInputsAccepted;
  std::cout << std::boolalpha
            << "Empty vectors rejected: " << emptyRejected << '\n'
            << "Mismatched vectors rejected: " << mismatchRejected << '\n'
            << "Three-way mismatch rejected: " << threeWayMismatchRejected
            << '\n'
            << "Nonfinite values rejected: " << nonfiniteRejected << '\n'
            << "Zero durations rejected: " << zeroDurationRejected << '\n'
            << "Zero sample counts rejected: " << zeroCountRejected << '\n'
            << "Valid statistics accepted: " << validInputsAccepted << '\n';
  return valid ? 0 : 1;
}
