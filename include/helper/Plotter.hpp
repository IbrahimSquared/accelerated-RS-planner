#ifndef PLOTTER_HPP
#define PLOTTER_HPP

#include "helper/HelperFunctions.hpp"
#include "parser/ConfigParser.hpp"
#include <SFML/Graphics.hpp>

namespace accelerated {

class Plotter {
public:
  static Plotter &getInstance();
  void draw(sf::RenderWindow &window, const State &to,
            const std::vector<State> &path,
            const std::vector<char> &motionTypes) const;

  Plotter(Plotter const &) = delete;
  void operator=(Plotter const &) = delete;

private:
  Plotter() = default;
  const Config &config_ = ConfigParser::getInstance().getConfig();
};

} // namespace accelerated

#endif // PLOTTER_HPP
