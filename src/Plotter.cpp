#include "helper/Plotter.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace std::numbers;

namespace accelerated {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
Plotter &Plotter::getInstance() {
  static Plotter instance;
  return instance;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Plotter::draw(sf::RenderWindow &window, const State &to,
                   const std::vector<State> &path,
                   const std::vector<char> &motionTypes) const {

  const auto firstInactiveMotion =
      std::find(motionTypes.begin(), motionTypes.end(), 'N');
  const std::size_t segmentCount =
      std::distance(motionTypes.begin(), firstInactiveMotion);
  if (path.size() != segmentCount + 1) {
    throw std::invalid_argument(
        "A plotted path must contain one more state than motion segments");
  }

  const double min_x = config_.min_x;
  const double max_x = config_.max_x;
  const double min_y = config_.min_y;
  const double max_y = config_.max_y;
  if (max_x <= min_x || max_y <= min_y) {
    throw std::invalid_argument("Plot bounds must have positive width and height");
  }

  const double xScale = window.getSize().x / (max_x - min_x);
  const double yScale = window.getSize().y / (max_y - min_y);

  const double marker_size = 20.0;
  sf::CircleShape firstPathPointMarker(marker_size);
  firstPathPointMarker.setPosition((path[0].x - min_x) * xScale - marker_size,
                                   (path[0].y - min_y) * yScale - marker_size);
  firstPathPointMarker.setFillColor(sf::Color::Green);
  window.draw(firstPathPointMarker);

  sf::CircleShape fromMarker(marker_size / 2.0);
  fromMarker.setPosition((path[0].x - min_x) * xScale - marker_size / 2.0,
                         (path[0].y - min_y) * yScale - marker_size / 2.0);
  fromMarker.setFillColor(sf::Color::White);
  window.draw(fromMarker);

  sf::CircleShape lastPathPointMarker(marker_size);
  const std::size_t lastIdx = segmentCount;
  lastPathPointMarker.setPosition(
      (path[lastIdx].x - min_x) * xScale - marker_size,
      (path[lastIdx].y - min_y) * yScale - marker_size);
  lastPathPointMarker.setFillColor(sf::Color::Red);
  window.draw(lastPathPointMarker);

  sf::CircleShape toMarker(marker_size / 2.0);
  toMarker.setPosition((to.x - min_x) * xScale - marker_size / 2.0,
                       (to.y - min_y) * yScale - marker_size / 2.0);
  toMarker.setFillColor(sf::Color::White);
  window.draw(toMarker);

  const double r = config_.r;
  const double L = 2;
  const double x_ext_start = path[0].x + L * r * cos(path[0].theta);
  const double y_ext_start = path[0].y + L * r * sin(path[0].theta);
  const double x_ext_end = path[lastIdx].x + L * r * cos(path[lastIdx].theta);
  const double y_ext_end = path[lastIdx].y + L * r * sin(path[lastIdx].theta);

  sf::Vertex line_start[] = {
      sf::Vertex(sf::Vector2f((path[0].x - min_x) * xScale,
                              (path[0].y - min_y) * yScale),
                 sf::Color::Green),
      sf::Vertex(sf::Vector2f((x_ext_start - min_x) * xScale,
                              (y_ext_start - min_y) * yScale),
                 sf::Color::Green)};
  line_start[0].color = sf::Color::Green;
  line_start[1].color = sf::Color::Green;
  window.draw(line_start, 2, sf::Lines);

  sf::Vertex line_end[] = {
      sf::Vertex(sf::Vector2f((path[lastIdx].x - min_x) * xScale,
                              (path[lastIdx].y - min_y) * yScale),
                 sf::Color::Green),
      sf::Vertex(sf::Vector2f((x_ext_end - min_x) * xScale,
                              (y_ext_end - min_y) * yScale),
                 sf::Color::Green)};
  line_end[0].color = sf::Color::Green;
  line_end[1].color = sf::Color::Green;
  window.draw(line_end, 2, sf::Lines);

  std::cout << "******************" << std::endl;
  std::cout << "Drawing motion type: ";
  for (size_t i = 0; i < motionTypes.size(); ++i) {
    std::cout << motionTypes[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "******************" << std::endl;

  // Draw path
  constexpr int numArcPoints = 100;
  for (std::size_t i = 0; i < segmentCount; ++i) {
    if (motionTypes[i] == 'S') {
      sf::Vertex line[] = {
          sf::Vertex(sf::Vector2f((path[i].x - min_x) * xScale,
                                  (path[i].y - min_y) * yScale),
                     sf::Color::Yellow),
          sf::Vertex(sf::Vector2f((path[i + 1].x - min_x) * xScale,
                                  (path[i + 1].y - min_y) * yScale),
                     sf::Color::Yellow)};
      line[0].color = sf::Color::Yellow;
      line[1].color = sf::Color::Yellow;
      window.draw(line, 2, sf::Lines);
    } else if (motionTypes[i] == 'N') {
      break;
    } else if (motionTypes[i] == 'L') {
      const double startAngle = path[i].theta;
      const double cx = path[i].x - r * cos(startAngle - pi / 2);
      const double cy = path[i].y - r * sin(startAngle - pi / 2);

      sf::VertexArray arc(sf::LineStrip, numArcPoints + 1);

      arc[0].position.x = (path[i].x - min_x) * xScale;
      arc[0].position.y = (path[i].y - min_y) * yScale;
      arc[0].color = sf::Color::Green;

      const double endAngle = path[i + 1].theta;
      for (int j = 1; j <= numArcPoints; ++j) {
        const double angle =
            startAngle + j * (endAngle - startAngle) / numArcPoints;
        const double pt_x = cx + r * cos(angle - pi / 2);
        const double pt_y = cy + r * sin(angle - pi / 2);
        arc[j].position =
            sf::Vector2f((pt_x - min_x) * xScale, (pt_y - min_y) * yScale);
        arc[j].color = sf::Color::Magenta;
      }
      window.draw(arc);
    } else if (motionTypes[i] == 'R') {
      const double startAngle = path[i].theta;
      const double cx = path[i].x + r * cos(startAngle - pi / 2);
      const double cy = path[i].y + r * sin(startAngle - pi / 2);

      sf::VertexArray arc(sf::LineStrip, numArcPoints + 1);

      arc[0].position.x = (path[i].x - min_x) * xScale;
      arc[0].position.y = (path[i].y - min_y) * yScale;
      arc[0].color = sf::Color::Green;

      const double endAngle = path[i + 1].theta;
      for (int j = 1; j <= numArcPoints; ++j) {
        const double angle =
            startAngle + j * (endAngle - startAngle) / numArcPoints;
        const double pt_x = cx - r * cos(angle - pi / 2);
        const double pt_y = cy - r * sin(angle - pi / 2);
        arc[j].position =
            sf::Vector2f((pt_x - min_x) * xScale, (pt_y - min_y) * yScale);
        arc[j].color = sf::Color::Cyan;
      }
      window.draw(arc);
    }
  }
}

} // namespace accelerated
