#include "parser/ConfigParser.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace accelerated {

namespace {

void trim(std::string &text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    text.clear();
    return;
  }
  const std::size_t last = text.find_last_not_of(" \t\r\n");
  text = text.substr(first, last - first + 1);
}

bool parseFiniteDouble(const std::string &key, const std::string &value,
                       double &result) {
  try {
    std::size_t consumed = 0;
    const double parsed = std::stod(value, &consumed);
    if (consumed != value.size() || !std::isfinite(parsed)) {
      throw std::invalid_argument("nonfinite or trailing input");
    }
    result = parsed;
    return true;
  } catch (...) {
    std::cerr << "Invalid value for " << key << ": " << value << '\n';
    return false;
  }
}

bool parseInteger(const std::string &key, const std::string &value,
                  int &result) {
  try {
    std::size_t consumed = 0;
    const int parsed = std::stoi(value, &consumed);
    if (consumed != value.size()) {
      throw std::invalid_argument("trailing input");
    }
    result = parsed;
    return true;
  } catch (...) {
    std::cerr << "Invalid value for " << key << ": " << value << '\n';
    return false;
  }
}

} // namespace

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
  Config parsed_config = config_;
  std::string line;
  while (std::getline(file, line)) {
    trim(line);
    if (line.empty() || line.front() == '#') {
      continue;
    }

    const std::size_t separator = line.find('=');
    if (separator == std::string::npos) {
      std::cerr << "Expected key=value entry: " << line << '\n';
      return false;
    }
    std::string key = line.substr(0, separator);
    std::string value = line.substr(separator + 1);
    trim(key);
    trim(value);
    if (key.empty() || value.empty()) {
      std::cerr << "Expected nonempty key and value: " << line << '\n';
      return false;
    }

    if (key == "min_x") {
      if (!parseFiniteDouble(key, value, parsed_config.min_x)) {
        return false;
      }
    } else if (key == "max_x") {
      if (!parseFiniteDouble(key, value, parsed_config.max_x)) {
        return false;
      }
    } else if (key == "min_y") {
      if (!parseFiniteDouble(key, value, parsed_config.min_y)) {
        return false;
      }
    } else if (key == "max_y") {
      if (!parseFiniteDouble(key, value, parsed_config.max_y)) {
        return false;
      }
    } else if (key == "min_thetaf") {
      if (!parseFiniteDouble(key, value, parsed_config.min_thetaf)) {
        return false;
      }
    } else if (key == "max_thetaf") {
      if (!parseFiniteDouble(key, value, parsed_config.max_thetaf)) {
        return false;
      }
    } else if (key == "number_of_final_states") {
      if (!parseInteger(key, value, parsed_config.number_of_final_states)) {
        return false;
      }
    } else if (key == "r") {
      if (!parseFiniteDouble(key, value, parsed_config.r)) {
        return false;
      }
    } else {
      std::cerr << "Invalid key: " << key << '\n';
      return false;
    }
  }

  if (parsed_config.min_x >= parsed_config.max_x) {
    std::cerr << "min_x must be less than max_x\n";
    return false;
  }
  if (parsed_config.min_y >= parsed_config.max_y) {
    std::cerr << "min_y must be less than max_y\n";
    return false;
  }
  if (parsed_config.min_thetaf > parsed_config.max_thetaf) {
    std::cerr << "min_thetaf must be less than or equal to max_thetaf\n";
    return false;
  }
  if (parsed_config.r <= 0.0) {
    std::cerr << "r must be greater than zero\n";
    return false;
  }
  if (parsed_config.number_of_final_states <= 0) {
    std::cerr << "number_of_final_states must be positive\n";
    return false;
  }

  config_ = parsed_config;
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
