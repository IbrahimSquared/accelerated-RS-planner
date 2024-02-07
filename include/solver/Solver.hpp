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
   * @brief Updates accessibility/visibility to a point using PDE advection.
   * @param [in] lightSourceNumber number of the lightsource whose visibility of
   * the point we are checking.
   * @param [in] lightSource_x x position of the lightsource.
   * @param [in] lightSource_y y position of the lightsource.
   * @param [in] x position of our queried pixel.
   * @param [in] y position of our queried pixel.
   */
  void benchmarkPlanners();
  // Simulates random start/end config and checks validity
  void RandomPathsValidityChecks();
  // Press enter to solve and draw a random path
  void drawRandomPath();

private:
  const Config &config_ = ConfigParser::getInstance().getConfig();
  const HelperFunctions &helper_functions_ = HelperFunctions::getInstance();
  const Plotter &plotter = Plotter::getInstance();

  /*!
   * @brief Transforms final config to local rotational frame
   * @param [in] start config
   * @param [in] final config
   * @param [in] final local config
   */
  void globalToLocal(const State &from, const State &to, State &to_local) const;

  // Main partitions of the algorithm
  void setA(const double x, const double y, const double thetaf, const double r,
            const double beta0, int &cond, double &d);
  void setB(const double x, const double y, const double thetaf, const double r,
            const bool pred_1, const bool pred_2, int &cond, double &d);

  // Forward simulate the path & check if it reaches the final config
  void validityCheck(const State &from, const State &to, const int cond,
                     const int Q, State &errors);

  void getMotionTypes(const int condition, char motionTypes[5]) const;
  void backProjectMotion(const int Q, const int condition, char motionTypes[5]);
  void getPath(const int condition, char motionTypes[5],
               std::vector<State> &path) const;
  void getErrors(const std::vector<State> &path, const State &to,
                 State &errors) const;

  // Proposed method
  void acceleratedRSPlanner(const State &from, const State &to, double &d,
                            const double r, int &cond, int &Q);

  void computeAndDrawRandomPath(sf::RenderWindow &window);

  inline const double wrapToPi(const double angle) const;
  inline const double euclideanDistance(const double x1, const double y1,
                                        const double x2,
                                        const double y2) const {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
  }

  inline void polar(const double x, const double y, double &rho,
                    double &theta) const {
    rho = sqrt(x * x + y * y);
    theta = atan2(y, x);
  }

  /*!
   * @brief Compute paths.
   */
  inline void T3(const double x, const double y, const double thetaf,
                 double &d);
  inline void T6(const double x, const double y, const double thetaf,
                 double &d);
  inline void T8(const double x, const double y, const double thetaf,
                 double &d);
  inline void T9(const double x, const double y, const double thetaf,
                 double &d);
  inline void T9(const double xi, const double eta, const double rho,
                 const double thetaf, double &d);
  inline void T10(const double x, const double y, const double thetaf,
                  double &d);
  inline void T11(const double x, const double y, const double thetaf,
                  double &d);
  inline void T12(const double x, const double y, const double thetaf,
                  double &d);
  inline void T13(const double x, const double y, const double thetaf,
                  double &xi, double &eta, double &rho, double &t, double &d);
  inline void T13(const double x, const double y, const double thetaf,
                  double &xi, double &eta, double &rho, double &t1, double &t2,
                  double &u1, double &u2, double &v1, double &v2) const;
  inline void T14(const double x, const double y, const double thetaf,
                  double &d);
  inline void T15(const double x, const double y, const double thetaf,
                  double &t, double &u, double &v) const;
  inline void T16(const double x, const double y, const double thetaf,
                  double &d);
  inline void T18(const double x, const double y, const double thetaf,
                  double &d);
  inline void T19(const double xi, const double eta, const double rho,
                  const double thetaf, double &d);
  inline void T19(const double x, const double y, const double thetaf,
                  double &d);
  inline void T19_rho(const double x, const double y, const double thetaf,
                      double &xi, double &eta, double &rho) const;
  inline void T20(const double x, const double y, const double thetaf,
                  double &d);
  inline void T21(const double x, const double y, const double thetaf,
                  double &d);
  /*!
   * @brief Compute angles.
   */
  inline void Beta1(const double x, const double y, const double thetaf,
                    double &t, double &u) const;
  inline void Beta4(const double x, const double y, const double thetaf,
                    double &t, double &u, double &v) const;
  inline void Beta5(const double x, const double y, const double thetaf,
                    double &t, double &u) const;
  inline void Beta7(const double x, const double y, const double thetaf,
                    double &t, double &u) const;
  inline void Beta9(const double x, const double y, const double thetaf,
                    double &t, double &u) const;

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
  double cxl0_ = 0, cxr0_ = 0, cyl0_ = 0, cyr0_ = 0;
  double cxlf_ = 0, cxrf_ = 0, cylf_ = 0, cyrf_ = 0;
  // Distances between circles & threshold
  double RR = 0, LL = 0, RL = 0, LR = 0, limit = 0;
};

} // namespace accelerated

#endif // SOLVER_HPP
