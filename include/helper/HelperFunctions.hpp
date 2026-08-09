#ifndef HELPERFUNCTIONS_HPP
#define HELPERFUNCTIONS_HPP

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "parser/ConfigParser.hpp"

namespace accelerated {

struct State {
  double x, y, theta;
};

template <typename T> auto durationInMicroseconds(T start, T end) {
  return std::chrono::duration_cast<std::chrono::microseconds>(end - start)
      .count();
}

class HelperFunctions {
public:
  static HelperFunctions &getInstance();

  void generateRandomStates(
      std::vector<State> &states, std::size_t number_of_final_states,
      bool use_config, double r = 1.0,
      std::optional<std::uint32_t> seed = std::nullopt) const;
  void errorStats(const std::vector<double> &distancesAcceleratedRS,
                  const std::vector<double> &distancesOMPL,
                  const std::vector<double> &distancesDesaulniers) const;
  void errorStats(const std::vector<double> &distancesAcceleratedRS,
                  const std::vector<double> &distancesOMPL) const;

  void timeStats(const long durationAcceleratedRS, const long durationOMPL,
                 const long durationDesaulniers,
                 const int number_of_final_states) const;

  HelperFunctions(HelperFunctions const &) = delete;
  void operator=(HelperFunctions const &) = delete;

private:
  HelperFunctions() = default;
  const Config &config_ = ConfigParser::getInstance().getConfig();
};

} // namespace accelerated

#endif // HELPERFUNCTIONS_HPP
