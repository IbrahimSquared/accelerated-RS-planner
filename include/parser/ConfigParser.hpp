#ifndef CONFIGPARSER_HPP
#define CONFIGPARSER_HPP

#include <string>

namespace accelerated {

struct Config {
  double min_x = -100;
  double max_x = 100;
  double min_y = -100;
  double max_y = 100;
  double min_thetaf = -180;
  double max_thetaf = 180;
  double r = 1;
  int number_of_final_states = 1e5;
};

class ConfigParser {
public:
  static ConfigParser &getInstance() {
    static ConfigParser instance;
    return instance;
  }

  bool parse(const std::string &filename);
  const Config &getConfig() const;

  ConfigParser(ConfigParser const &) = delete;
  void operator=(ConfigParser const &) = delete;

  void printConfig() const;

private:
  ConfigParser() = default;
  Config config_;
};

} // namespace accelerated

#endif // CONFIGPARSER_HPP
