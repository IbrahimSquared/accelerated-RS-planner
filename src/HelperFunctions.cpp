#include "helper/HelperFunctions.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>

using namespace std::numbers;

namespace accelerated {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void HelperFunctions::generateRandomStates(
    std::vector<State> &states, const std::size_t number_of_final_states,
    const bool use_config, double r,
    const std::optional<std::uint32_t> seed) const {
  double min_x, max_x, min_y, max_y, min_thetaf, max_thetaf;
  if (use_config) {
    min_x = config_.min_x;
    max_x = config_.max_x;
    min_y = config_.min_y;
    max_y = config_.max_y;
    min_thetaf = config_.min_thetaf * pi / 180;
    max_thetaf = config_.max_thetaf * pi / 180;
    r = config_.r;
  } else {
    min_x = -5.0;
    max_x = 1.0;
    min_y = 0.0;
    max_y = 5.0;
    min_thetaf = -pi;
    max_thetaf = pi;
  }

  if (use_config) {
    std::cout << "\nUsing configuration to generate " << number_of_final_states
              << " random states" << std::endl;
  } else {
    std::cout << "\nGenerating " << number_of_final_states
              << " random states for benchmarking the 3 planners" << std::endl;
  }
  std::cout << "x: [" << min_x << ", " << max_x << "]" << std::endl;
  std::cout << "y: [" << min_y << ", " << max_y << "]" << std::endl;
  std::cout << "theta: [" << min_thetaf << ", " << max_thetaf << "]"
            << std::endl;
  std::cout << "r: " << r << std::endl;

  std::uint32_t resolved_seed;
  if (seed.has_value()) {
    resolved_seed = *seed;
  } else {
    std::random_device rd;
    resolved_seed = rd();
  }
  std::cout << "seed: " << resolved_seed << std::endl;
  std::mt19937 gen(resolved_seed);
  std::uniform_real_distribution<> x_dis(min_x, max_x);
  std::uniform_real_distribution<> y_dis(min_y, max_y);
  std::uniform_real_distribution<> theta_dis(min_thetaf, max_thetaf);

  states.clear();
  states.reserve(number_of_final_states);

  const std::size_t progress_print_rate =
      std::max<std::size_t>(1, std::ceil(0.1 * number_of_final_states));

  std::size_t counter = 0;
  std::size_t counter_old = 0;
  while (counter < number_of_final_states) {
    const double x = x_dis(gen);
    const double y = y_dis(gen);
    const double theta = theta_dis(gen);
    const State state(x, y, theta);

    bool accept_state = true;
    if (!use_config) { // for benchmarking the 3 solvers
      const double cx = x + r * cos(theta);
      const double cy = y + r * sin(theta);
      accept_state = (cx <= 1 && cy >= 0);
    }

    if ((counter % progress_print_rate) == 0 && counter != counter_old) {
      counter_old = counter;
      std::cout << "State generation progress: "
                << 10.0 * counter / progress_print_rate << "%" << std::endl;
    }
    if (accept_state) {
      states.push_back(state);
      ++counter;
    }
  }
  std::cout << "State generation progress: 100% \n" << std::endl;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
HelperFunctions &HelperFunctions::getInstance() {
  static HelperFunctions instance;
  return instance;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void HelperFunctions::timeStats(const long duration_acceleratedRS,
                                const long duration_OMPL,
                                const long duration_Desaulniers,
                                const int number_of_final_states) const {
  if (duration_acceleratedRS <= 0 || duration_OMPL <= 0 ||
      duration_Desaulniers <= 0) {
    throw std::invalid_argument("Timing durations must be positive");
  }
  if (number_of_final_states <= 0) {
    throw std::invalid_argument("Timing statistics require at least one state");
  }

  std::cout << "------------------" << std::endl;
  std::cout << "Time stats: " << std::endl;
  std::cout << "******************" << std::endl;
  std::cout << "acceleratedRS (Proposed) total time: " << duration_acceleratedRS
            << " [µs] \n";
  std::cout << "DesaulniersAlgorithm total time: " << duration_Desaulniers
            << " [µs]\n";
  std::cout << "OMPL total time: " << duration_OMPL << " [µs]\n";

  std::cout << "******************" << std::endl;

  std::cout << "Ratio acceleratedRS (Proposed): "
            << 1.0 * duration_OMPL / duration_acceleratedRS << "\n";
  std::cout << "Ratio DesaulniersAlgorithm: "
            << 1.0 * duration_OMPL / duration_Desaulniers << "\n";

  std::cout << "******************" << std::endl;

  std::cout << "acceleratedRS (Proposed): "
            << 1.0 * duration_acceleratedRS / number_of_final_states
            << " [µs] per state\n";
  std::cout << "DesaulniersAlgorithm: "
            << 1.0 * duration_Desaulniers / number_of_final_states
            << " [µs] per state\n";
  std::cout << "OMPL: " << 1.0 * duration_OMPL / number_of_final_states
            << " [µs] per state\n";
  std::cout << "------------------" << std::endl;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void HelperFunctions::errorStats(
    const std::vector<double> &distances_acceleratedRS,
    const std::vector<double> &distances_OMPL,
    const std::vector<double> &distances_Desaulniers) const {

  if (distances_OMPL.empty() ||
      distances_acceleratedRS.size() != distances_OMPL.size() ||
      distances_Desaulniers.size() != distances_OMPL.size()) {
    throw std::invalid_argument(
        "Distance vectors must be nonempty and have equal sizes");
  }

  double avg_error_accelerated = 0.0;
  double avg_error_Desaulniers = 0.0;
  double max_error_accelerated = 0.0;
  double max_error_Desaulniers = 0.0;
  for (size_t i = 0; i < distances_OMPL.size(); i++) {
    if (!std::isfinite(distances_acceleratedRS[i]) ||
        !std::isfinite(distances_OMPL[i]) ||
        !std::isfinite(distances_Desaulniers[i])) {
      throw std::invalid_argument("Distance vectors must contain finite values");
    }
    double error_accelerated =
        std::abs(distances_acceleratedRS[i] - distances_OMPL[i]);
    double error_Desaulniers =
        std::abs(distances_Desaulniers[i] - distances_OMPL[i]);
    avg_error_accelerated += error_accelerated;
    avg_error_Desaulniers += error_Desaulniers;
    max_error_accelerated = std::max(max_error_accelerated, error_accelerated);
    max_error_Desaulniers = std::max(max_error_Desaulniers, error_Desaulniers);
  }
  avg_error_accelerated /= distances_OMPL.size();
  avg_error_Desaulniers /= distances_OMPL.size();
  std::cout << "\n------------------" << std::endl;
  std::cout << "Error stats: " << std::endl;
  std::cout << "Average error acceleratedRS (Proposed): "
            << avg_error_accelerated << std::endl;
  std::cout << "Average error DesaulniersAlgorithm: " << avg_error_Desaulniers
            << std::endl;
  std::cout << "******************" << std::endl;
  std::cout << "Maximum error acceleratedRS (Proposed): "
            << max_error_accelerated << std::endl;
  std::cout << "Maximum error DesauliersAgorithm: " << max_error_Desaulniers
            << std::endl;
  std::cout << "-----------------" << std::endl;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void HelperFunctions::errorStats(
    const std::vector<double> &distances_acceleratedRS,
    const std::vector<double> &distances_OMPL) const {

  if (distances_OMPL.empty() ||
      distances_acceleratedRS.size() != distances_OMPL.size()) {
    throw std::invalid_argument(
        "Distance vectors must be nonempty and have equal sizes");
  }

  double avg_error_accelerated = 0.0;
  double max_error_accelerated = 0.0;
  for (size_t i = 0; i < distances_OMPL.size(); i++) {
    if (!std::isfinite(distances_acceleratedRS[i]) ||
        !std::isfinite(distances_OMPL[i])) {
      throw std::invalid_argument("Distance vectors must contain finite values");
    }
    double error_accelerated =
        std::abs(distances_acceleratedRS[i] - distances_OMPL[i]);
    avg_error_accelerated += error_accelerated;
    max_error_accelerated = std::max(max_error_accelerated, error_accelerated);
  }
  avg_error_accelerated /= distances_OMPL.size();
  std::cout << "\n------------------" << std::endl;
  std::cout << "Error stats: " << std::endl;
  std::cout << "Average error acceleratedRS (Proposed): "
            << avg_error_accelerated << std::endl;
  std::cout << "******************" << std::endl;
  std::cout << "Maximum error acceleratedRS (Proposed): "
            << max_error_accelerated << std::endl;
  std::cout << "-----------------" << std::endl;
}

} // namespace accelerated
