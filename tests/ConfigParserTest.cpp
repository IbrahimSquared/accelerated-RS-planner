#include "parser/ConfigParser.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace {

bool writeConfig(const std::string &path, const std::string &contents) {
  std::ofstream file(path, std::ios::trunc);
  file << contents;
  return file.good();
}

bool sameConfig(const accelerated::Config &lhs,
                const accelerated::Config &rhs) {
  return lhs.min_x == rhs.min_x && lhs.max_x == rhs.max_x &&
         lhs.min_y == rhs.min_y && lhs.max_y == rhs.max_y &&
         lhs.min_thetaf == rhs.min_thetaf &&
         lhs.max_thetaf == rhs.max_thetaf && lhs.r == rhs.r &&
         lhs.number_of_final_states == rhs.number_of_final_states;
}

} // namespace

int main() {
  accelerated::ConfigParser &parser = accelerated::ConfigParser::getInstance();

  const std::string reorderedPath = "config_parser_reordered.config";
  const bool wroteReordered = writeConfig(
      reorderedPath,
      "   # A minimum above the defaults is valid once the later maximum is "
      "read\n"
      "min_x = 200\n"
      "max_x = 300\n"
      "min_y = 400\n"
      "max_y = 500\n"
      "min_thetaf = -90\n"
      "max_thetaf = 90\n"
      "r = 2\n"
      "number_of_final_states = 10\n");
  const bool reorderedAccepted = wroteReordered && parser.parse(reorderedPath);
  const accelerated::Config accepted = parser.getConfig();

  const std::string nanPath = "config_parser_nan.config";
  const bool nanRejected =
      writeConfig(nanPath, "r=nan\n") && !parser.parse(nanPath);
  const bool unchangedAfterNan = sameConfig(accepted, parser.getConfig());

  const std::string trailingPath = "config_parser_trailing.config";
  const bool trailingRejected =
      writeConfig(trailingPath, "r=2oops\n") && !parser.parse(trailingPath);
  const bool unchangedAfterTrailing = sameConfig(accepted, parser.getConfig());

  const std::string zeroRadiusPath = "config_parser_zero_radius.config";
  const bool zeroRadiusRejected =
      writeConfig(zeroRadiusPath, "r=0\n") && !parser.parse(zeroRadiusPath);
  const bool unchangedAfterZeroRadius =
      sameConfig(accepted, parser.getConfig());

  const std::string equalBoundsPath = "config_parser_equal_bounds.config";
  const bool equalBoundsRejected =
      writeConfig(equalBoundsPath, "min_x=250\nmax_x=250\n") &&
      !parser.parse(equalBoundsPath);
  const bool unchangedAfterEqualBounds =
      sameConfig(accepted, parser.getConfig());

  const bool valid = reorderedAccepted && accepted.min_x == 200.0 &&
                     accepted.max_x == 300.0 && nanRejected &&
                     unchangedAfterNan && trailingRejected &&
                     unchangedAfterTrailing && zeroRadiusRejected &&
                     unchangedAfterZeroRadius && equalBoundsRejected &&
                     unchangedAfterEqualBounds;

  std::cout << std::boolalpha
            << "Order-independent configuration accepted: "
            << reorderedAccepted << '\n'
            << "NaN rejected atomically: "
            << (nanRejected && unchangedAfterNan) << '\n'
            << "Trailing junk rejected atomically: "
            << (trailingRejected && unchangedAfterTrailing) << '\n'
            << "Zero radius rejected atomically: "
            << (zeroRadiusRejected && unchangedAfterZeroRadius) << '\n'
            << "Equal plot bounds rejected atomically: "
            << (equalBoundsRejected && unchangedAfterEqualBounds) << '\n';
  return valid ? 0 : 1;
}
