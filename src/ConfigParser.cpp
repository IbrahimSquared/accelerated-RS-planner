#include "parser/ConfigParser.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

namespace accelerated {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
bool ConfigParser::parse(const std::string &filename) {
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Failed to open " << filename << '\n';
    return false;
  }
  std::cout << "Attempting to parse config file: " << filename << '\n';
  std::string line;
  while (std::getline(file, line)) {
    // Ignore comments and blank lines
    if (line.empty() || line[0] == '#') {
      continue;
    }

    // Split the line into key and value
    std::istringstream iss(line);
    std::string key;
    if (!std::getline(iss, key, '=')) {
      continue;
    }
    std::string value;
    if (!std::getline(iss, value)) {
      continue;
    }

    // Trim leading and trailing whitespace from key and value
    key.erase(0, key.find_first_not_of(" \t"));
    key.erase(key.find_last_not_of(" \t") + 1);
    value.erase(0, value.find_first_not_of(" \t"));
    value.erase(value.find_last_not_of(" \t") + 1);

    if (key == "min_x") {
      try {
        config_.min_x = std::stod(value);
      } catch (...) {
        std::cerr << "Invalid value for " << key << ": " << value << '\n';
        return false;
      }
    } else if (key == "max_x") {
      try {
        config_.max_x = std::stod(value);
      } catch (...) {
        std::cerr << "Invalid value for " << key << ": " << value << '\n';
        return false;
      }
    } else if (key == "min_y") {
      try {
        config_.min_y = std::stod(value);
      } catch (...) {
        std::cerr << "Invalid value for " << key << ": " << value << '\n';
        return false;
      }
    } else if (key == "max_y") {
      try {
        config_.max_y = std::stod(value);
      } catch (...) {
        std::cerr << "Invalid value for " << key << ": " << value << '\n';
        return false;
      }
    } else if (key == "min_thetaf") {
      try {
        config_.min_thetaf = std::stod(value);
      } catch (...) {
        std::cerr << "Invalid value for " << key << ": " << value << '\n';
        return false;
      }
    } else if (key == "max_thetaf") {
      try {
        config_.max_thetaf = std::stod(value);
      } catch (...) {
        std::cerr << "Invalid value for " << key << ": " << value << '\n';
        return false;
      }
    } else if (key == "number_of_final_states") {
      try {
        config_.number_of_final_states = std::stoi(value);
      } catch (...) {
        std::cerr << "Invalid value for " << key << ": " << value << '\n';
        return false;
      }
    } else if (key == "r") {
      try {
        config_.r = std::stod(value);
      } catch (...) {
        std::cerr << "Invalid value for " << key << ": " << value << '\n';
        return false;
      }
    } else {
      std::cerr << "Invalid key: " << key << '\n';
      return false;
    }
    if (config_.min_x > config_.max_x) {
      std::cerr << "min_x must be less than or equal to max_x\n";
      return false;
    }
    if (config_.min_y > config_.max_y) {
      std::cerr << "min_y must be less than or equal to max_y\n";
      return false;
    }
    if (config_.min_thetaf > config_.max_thetaf) {
      std::cerr << "min_thetaf must be less than or equal to max_thetaf\n";
      return false;
    }
    if (config_.number_of_final_states < 0) {
      std::cerr << "number_of_final_states must be non-negative\n";
      return false;
    }
  }
  std::cout << "Config file parsed successfully" << std::endl;

  return true;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
const Config &ConfigParser::getConfig() const { return config_; }

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ConfigParser::printConfig() const {
  std::cout << "min_x: " << config_.min_x << '\n';
  std::cout << "max_x: " << config_.max_x << '\n';
  std::cout << "min_y: " << config_.min_y << '\n';
  std::cout << "max_y: " << config_.max_y << '\n';
  std::cout << "min_thetaf: " << config_.min_thetaf << '\n';
  std::cout << "max_thetaf: " << config_.max_thetaf << '\n';
  std::cout << "r: " << config_.r << '\n';
  std::cout << "number_of_final_states: " << config_.number_of_final_states
            << '\n';
}

} // namespace accelerated
