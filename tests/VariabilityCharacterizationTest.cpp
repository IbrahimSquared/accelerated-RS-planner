#include "parser/ConfigParser.hpp"
#include "solver/Solver.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

int main(const int argc, char *argv[]) {
  std::size_t sampleCount = 100000;
  if (argc == 2) {
    sampleCount = std::stoull(argv[1]);
  }
  constexpr std::uint32_t seed = 0x5eed5eedU;
  const std::string configPath =
      std::string(PROJECT_SOURCE_DIR) + "/config/settings.config";

  accelerated::ConfigParser &parser = accelerated::ConfigParser::getInstance();
  if (!parser.parse(configPath)) {
    std::cerr << "Failed to load the variability test configuration\n";
    return 1;
  }

  accelerated::Solver solver;
  solver.characterizeVariability(sampleCount, seed);

  std::ifstream report("variability.txt");
  if (!report) {
    std::cerr << "variability.txt was not generated\n";
    return 1;
  }

  std::size_t totalCount = 0;
  std::size_t rowCount = 0;
  bool validRows = true;
  std::size_t count = 0;
  double acceleratedTime = 0.0;
  double omplTime = 0.0;
  double speedup = 0.0;
  while (report >> count >> acceleratedTime >> omplTime >> speedup) {
    ++rowCount;
    totalCount += count;
    if (count == 0) {
      validRows = validRows && acceleratedTime == 0.0 && omplTime == 0.0 &&
                  speedup == 0.0;
    } else {
      validRows = validRows && acceleratedTime > 0.0 && omplTime > 0.0 &&
                  std::isfinite(speedup) && speedup > 0.0;
    }
  }

  const bool valid = rowCount == 20 && totalCount == sampleCount && validRows;
  std::cout << "Report rows: " << rowCount << '\n'
            << "Categorized samples: " << totalCount << '\n'
            << "Rows have consistent timing fields: " << std::boolalpha
            << validRows << '\n';
  return valid ? 0 : 1;
}
