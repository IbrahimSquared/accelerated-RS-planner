#include "parser/ConfigParser.hpp"
#include "solver/Solver.hpp"
#include <iostream>

int main() {
  accelerated::ConfigParser &config_parser =
      accelerated::ConfigParser::getInstance();
  if (!config_parser.parse("config/settings.config")) {
    std::cerr << "Error parsing config file" << std::endl;
    return 1;
  }
  config_parser.printConfig();

  accelerated::Solver solver;
  // solver.RandomPathsValidityChecks();
  solver.drawRandomPath();
  // solver.benchmarkPlanners();
  // solver.characterizeVariability();

  return 0;
}
