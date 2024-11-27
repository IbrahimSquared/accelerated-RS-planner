#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "helper/HelperFunctions.hpp"
#include "helper/Plotter.hpp"
#include "parser/ConfigParser.hpp"
#include <SFML/Graphics.hpp>
#include <cmath>

namespace accelerated {

class Solver {
public:
  Solver() = default;

  /*!
   * @brief
   */
  void benchmarkPlanners();
  // Simulates random start/end config and checks validity
  void RandomPathsValidityChecks();
  // Press enter to solve and draw a random path
  void drawRandomPath();
  // Characterizes variability in speedup between types/conditions
  void characterizeVariability();
  // Evaluate a single query
  void evaluateQuery(const State &from, const State &to, const double r = 1);

private:
  const Config &config_ = ConfigParser::getInstance().getConfig();
  const HelperFunctions &helper_functions_ = HelperFunctions::getInstance();
  const Plotter &plotter = Plotter::getInstance();

  // Main partitions of the algorithm
  void setA(const double x, const double y, const double thetaf, const double r,
            const double beta0, int &cond, double &d, const double xn,
            const double x0, const double yn);
  void setB(const double x, const double y, const double thetaf, const double r,
            int &cond, double &d, const double beta0);

  // Forward simulate the path & check if it reaches the final config
  void validityCheck(const State &from, const State &to, const int cond,
                     const int Q, State &errors, const double r);

  void getMotionTypes(const int condition, char motionTypes[5]) const;
  void backProjectMotion(const int Q, const int condition, char motionTypes[5]);
  void getPath(const int condition, char motionTypes[5],
               std::vector<State> &path, const double r) const;
  void getErrors(const std::vector<State> &path, const State &to,
                 State &errors) const;

  // Proposed method
  void acceleratedRSPlanner(const State &from, const State &to, double &d,
                            const double r, int &cond, int &Q);

  void computeAndDrawRandomPath(sf::RenderWindow &window);

  inline double wrapToPi(const double angle) const;
  inline double wrapTo2Pi(const double angle) const;
  inline double euclideanDistance(const double x1, const double y1,
                                  const double x2, const double y2) const {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
  }

  /*!
   * @brief Compute paths.
   */
  inline void P1(const double x, const double y, const double thetaf,
                 double &d);
  inline void P2(const double x, const double y, const double thetaf,
                 double &d);
  inline void P3(const double x, const double y, const double thetaf,
                 double &d);
  inline void P4(const double x, const double y, const double thetaf,
                 double &d);
  inline void P5(const double x, const double y, const double thetaf,
                 double &d);
  inline void P6(const double x, const double y, const double thetaf,
                 double &d);
  inline void P7(const double x, const double y, const double thetaf,
                 double &d);
  inline void P8(const double x, const double y, const double thetaf,
                 double &d);
  inline void P9(const double x, const double y, const double thetaf,
                 double &d);
  inline void P10(const double x, const double y, const double thetaf,
                  double &d);
  inline void P11(const double x, const double y, const double thetaf,
                  double &d);
  inline void P12(const double x, const double y, const double thetaf,
                  double &d);
  inline void P13(const double x, const double y, const double thetaf,
                  double &d);
  inline void P14(const double x, const double y, const double thetaf,
                  double &d);
  inline void P15(const double x, const double y, const double thetaf,
                  double &d);
  inline void P16(const double x, const double y, const double thetaf,
                  double &d);
  inline void P17(const double x, const double y, const double thetaf,
                  double &d);
  inline void P18(const double x, const double y, const double thetaf,
                  double &d);
  inline void P19(const double x, const double y, const double thetaf,
                  double &d);
  inline void P20(const double x, const double y, const double thetaf,
                  double &d);

  // 21x5 matrix, holds segments lengths
  double lengths_[21][5] = {
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
  // 21x5 matrix, hold segment directions (+ or -)
  int directions_[21][5] = {
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
  // Centers of the circles
  double c0Lx_ = 0, c0Rx_ = 0, c0Ly_ = 0, c0Ry_ = 0;
  double cfLx_ = 0, cfRx_ = 0, cfLy_ = 0, cfRy_ = 0;
  // Distances between circles & threshold
  double RR = 0, LL = 0, RL = 0, LR = 0, limit = 0;
};

} // namespace accelerated

#endif // SOLVER_HPP
