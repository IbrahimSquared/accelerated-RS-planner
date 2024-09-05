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
  // config_parser.printConfig();

  accelerated::Solver solver;
  // solver.RandomPathsValidityChecks();
  // solver.drawRandomPath();
  solver.benchmarkPlanners();
  // solver.characterizeVariability();

  // const double pi = 3.1415926535897932;
  // // query a single path
  // accelerated::State from = {100, 100, 0};
  // accelerated::State to = {101.24796, 101.64845, -0.9306};
  // const double radius = 1;
  // solver.evaluateQuery(from, to, radius);

  return 0;
}
