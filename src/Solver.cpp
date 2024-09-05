#include "solver/Solver.hpp"
#include "classic/DesaulniersAlgorithm.hpp"

#include <SFML/Graphics.hpp>
#include <chrono>
#include <complex>
#include <fstream>
#include <numbers>
#include <ompl/base/ScopedState.h>
#include <ompl/base/spaces/ReedsSheppStateSpace.h>

using namespace std::numbers;
namespace ob = ompl::base;

namespace accelerated {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double Solver::wrapToPi(double angle) const {
  angle = fmod(angle + pi, 2 * pi);
  if (angle < 0) {
    angle += 2 * pi;
  }
  return angle - pi;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double Solver::wrapTo2Pi(double angle) const {
  angle = fmod(angle + 2 * pi, 2 * pi);
  if (angle < 0) {
    angle += 2 * pi;
  }
  return angle;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P1(const double x, const double y, const double thetaf,
                       double &d) {
  const double dx = x + sin(thetaf);
  const double dy = y - 1 - cos(thetaf);
  const double phi = atan2(dy, dx);
  double rho = dx * dx + dy * dy;
  double t, u;
  if (rho < 4) {
    u = 0;
    t = wrapToPi(phi + pi / 2);
  } else {
    u = sqrt(rho - 4);
    t = atan2(2, u);
    t = wrapToPi(phi + t);
  }
  const double v = wrapToPi(t - thetaf);
  d = fabs(t) + fabs(u) + fabs(v);
  lengths_[0][0] = t;
  lengths_[0][1] = u;
  lengths_[0][2] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P2(const double x, const double y, const double thetaf,
                       double &d) {
  const double dx = x - sin(thetaf);
  const double dy = y - 1 + cos(thetaf);

  const double u = sqrt(dx * dx + dy * dy);
  const double t = atan2(dy, dx);
  const double v = wrapToPi(thetaf - t);

  d = fabs(t) + fabs(u) + fabs(v);

  lengths_[1][0] = t;
  lengths_[1][1] = u;
  lengths_[1][2] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P3(const double x, const double y, const double thetaf,
                       double &d) {
  // backwards
  const double xb = -(x * cos(thetaf) + y * sin(thetaf));
  const double yb = -(x * sin(thetaf) - y * cos(thetaf));
  const double xi = xb + sin(thetaf);
  const double eta = yb - 1 - cos(thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double phi = atan2(xi, -eta);
  const double u = 2 - rho;
  const double v = wrapToPi(phi + 0.5 * pi - thetaf);

  d = fabs(phi) + fabs(u) + fabs(v) + pi / 2;

  lengths_[2][0] = -v;
  lengths_[2][1] = -u;
  lengths_[2][2] = pi / 2;
  lengths_[2][3] = -phi;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P4(const double x, const double y, const double thetaf,
                       double &d) {
  const double xb = -(x * cos(thetaf) + y * sin(thetaf));
  const double yb = x * sin(thetaf) - y * cos(thetaf);

  const double xi = xb - sin(-thetaf);
  const double eta = yb - 1 + cos(-thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double phi = atan2(eta, xi);
  const double rr = sqrt(rho * rho - 4);

  const double u = 2 - rr;
  const double t = wrapToPi(phi + atan2(rr, -2));
  const double v = wrapToPi(-thetaf - 0.5 * pi - t);

  d = fabs(t) + fabs(u) + fabs(v) + pi / 2;

  lengths_[3][0] = -v;
  lengths_[3][1] = -u;
  lengths_[3][2] = pi / 2;
  lengths_[3][3] = -t;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P5(const double x, const double y, const double thetaf,
                       double &d) {
  const double xi = x - sin(-thetaf);
  const double eta = -y - 1 + cos(-thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double phi = atan2(eta, xi);
  const double rr = sqrt(rho * rho - 4);
  const double u = 2 - rr;
  const double t = wrapToPi(phi + atan2(rr, -2));
  const double v = wrapToPi(-thetaf - 0.5 * pi - t);

  d = fabs(t) + fabs(u) + fabs(v) + pi / 2;

  lengths_[4][0] = t;
  lengths_[4][1] = -pi / 2;
  lengths_[4][2] = u;
  lengths_[4][3] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P6(const double x, const double y, const double thetaf,
                       double &d) {
  const double xi = x + sin(-thetaf);
  const double eta = -y - 1 - cos(-thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double phi = atan2(xi, -eta);
  const double t = phi;
  const double u = 2 - rho;
  const double v = wrapToPi(phi + 0.5 * pi + thetaf);

  d = fabs(t) + fabs(u) + fabs(v) + pi / 2;

  lengths_[5][0] = t;
  lengths_[5][1] = -pi / 2;
  lengths_[5][2] = u;
  lengths_[5][3] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P7(const double x, const double y, const double thetaf,
                       double &d) {
  const double dx = x + sin(-thetaf);
  const double dy = -y - 1 - cos(-thetaf);
  const double rho = (dx * dx + dy * dy);
  const double phi = atan2(dy, dx);
  double u, t;
  if (rho < 4) {
    u = 0;
    t = wrapToPi(phi + pi / 2);
  } else {
    u = sqrt(rho - 4);
    t = wrapToPi(phi + atan2(2, u));
  }
  const double v = wrapToPi(t + thetaf);

  d = fabs(t) + fabs(u) + fabs(v);

  lengths_[6][0] = t;
  lengths_[6][1] = u;
  lengths_[6][2] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P8(const double x, const double y, const double thetaf,
                       double &d) {
  // must be backwards, check later
  const double xb = -(x * cos(thetaf) + y * sin(thetaf));
  const double yb = -(x * sin(thetaf) - y * cos(thetaf));
  const double xi = xb - sin(thetaf);
  const double eta = yb - 1 + cos(thetaf);
  const double rho = xi * xi + eta * eta;
  const double phi = atan2(eta, xi);
  const double rr = sqrt(rho - 4);
  const double u = 2 - rr;
  const double t = wrapToPi(phi + atan2(rr, -2));
  const double v = wrapToPi(thetaf - 0.5 * pi - t);

  d = fabs(t) + fabs(u) + fabs(v) + pi / 2;

  lengths_[7][0] = -v;
  lengths_[7][1] = -u;
  lengths_[7][2] = pi / 2;
  lengths_[7][3] = -t;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P9(const double x, const double y, const double thetaf,
                       double &d) {
  const double xi = -x + sin(thetaf);
  const double eta = -y - 1 - cos(thetaf);
  const double rho = xi * xi + eta * eta;
  const double phi = atan2(eta, xi);
  double t, u;
  u = 4 - sqrt(rho - 4);
  t = wrapToPi(atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta));
  const double v = wrapToPi(t - thetaf);

  d = fabs(t) + fabs(u) + fabs(v) + pi;

  lengths_[8][0] = -t;
  lengths_[8][1] = pi / 2;
  lengths_[8][2] = -u;
  lengths_[8][3] = pi / 2;
  lengths_[8][4] = -v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P10(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = -x - sin(thetaf);
  const double eta = -y - 1 + cos(thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double phi = atan2(eta, xi);
  const double rr = sqrt(rho * rho - 4);
  const double u = 2 - rr;
  const double t = wrapToPi(phi + atan2(rr, -2));
  const double v = wrapToPi(thetaf - 0.5 * pi - t);

  d = fabs(t) + fabs(u) + fabs(v) + pi / 2;

  lengths_[9][0] = -t;
  lengths_[9][1] = pi / 2;
  lengths_[9][2] = -u;
  lengths_[9][3] = -v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P11(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = -x + sin(thetaf);
  const double eta = -y - 1 - cos(thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double phi = atan2(xi, -eta);
  const double t = phi;
  const double u = 2 - rho;
  const double v = wrapToPi(phi + 0.5 * pi - thetaf);

  d = fabs(t) + fabs(u) + fabs(v) + pi / 2;

  lengths_[10][0] = -t;
  lengths_[10][1] = pi / 2;
  lengths_[10][2] = -u;
  lengths_[10][3] = -v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P12(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = x + sin(-thetaf);
  const double eta = -y - 1 - cos(-thetaf);
  const double rho = xi * xi + eta * eta;
  const double u = 4 - sqrt(rho - 4);
  const double t =
      wrapToPi(atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta));
  const double v = wrapToPi(t + thetaf);

  d = fabs(t) + fabs(u) + fabs(v) + pi;

  lengths_[11][0] = t;
  lengths_[11][1] = -pi / 2;
  lengths_[11][2] = u;
  lengths_[11][3] = -pi / 2;
  lengths_[11][4] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P13(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = -x + sin(thetaf);
  const double eta = -y - 1 - cos(thetaf);
  const double rho = xi * xi + eta * eta;
  const double rho_CCCC = (20 - rho) / 16;
  double u;
  if (rho_CCCC < 0 || rho_CCCC > 1) {
    std::complex<double> a;
    a = -acos(rho_CCCC);
    u = std::real(a);
  } else {
    u = -acos(rho_CCCC);
  }
  const double A = sin(u);
  const double B = cos(u) - 2;
  const double t = atan2(eta * A - xi * B, xi * A + eta * B);
  const double v = (t - thetaf);

  d = fabs(t) + 2 * fabs(u) + fabs(v);

  lengths_[12][0] = -t;
  lengths_[12][1] = -u;
  lengths_[12][2] = -u;
  lengths_[12][3] = -v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P14(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = -x - sin(thetaf);
  const double eta = -y - 1 + cos(thetaf);
  const double u1 = sqrt(xi * xi + eta * eta);
  const double theta = atan2(eta, xi);
  const double u = -2 * asin(0.25 * u1);
  const double t = wrapToPi(theta + 0.5 * u + pi);
  const double v = wrapToPi(thetaf - t + u);

  d = fabs(t) + fabs(u) + fabs(v);

  lengths_[13][0] = -t;
  lengths_[13][1] = -u;
  lengths_[13][2] = -v;

  // std::cout << "T14: " << (-t < 0 ? "-" : "+") << "t, " << (-u < 0 ? "-" :
  // "+")
  //           << "u, " << (-v < 0 ? "-" : "+") << "v" << std::endl;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P15(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = x - sin(thetaf);
  const double eta = y - 1 + cos(thetaf);
  const double u1 = sqrt(xi * xi + eta * eta);
  const double theta = atan2(eta, xi);
  const double u = -2 * asin(0.25 * u1);
  const double t = wrapToPi(theta + 0.5 * u + pi);
  const double v = wrapToPi(thetaf - t + u);

  d = fabs(t) + fabs(u) + fabs(v);

  lengths_[14][0] = t;
  lengths_[14][1] = u;
  lengths_[14][2] = v;

  // std::cout << "T15: " << (t < 0 ? "-" : "+") << "t, " << (u < 0 ? "-" : "+")
  //           << "u, " << (v < 0 ? "-" : "+") << "v" << std::endl;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P16(const double x, const double y, const double thetaf,
                        double &d) {
  const double xb = x * cos(thetaf) + y * sin(thetaf);
  const double yb = x * sin(thetaf) - y * cos(thetaf);
  const double xi = -xb + sin(thetaf);
  const double eta = yb - 1 + cos(thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double theta = atan2(eta, xi);
  const double u = -2 * asin(0.25 * rho);
  const double t = wrapToPi(theta + 0.5 * u + pi);
  const double v = wrapToPi(-thetaf - t + u);

  d = fabs(t) + fabs(u) + fabs(v);

  lengths_[15][0] = -v;
  lengths_[15][1] = -u;
  lengths_[15][2] = -t;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P17(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = x + sin(-thetaf);
  const double eta = -y - 1 - cos(-thetaf);
  const double rho = xi * xi + eta * eta;
  const double rho_CCCC = (20 - rho) / 16;
  double u;
  if (rho_CCCC < 0 || rho_CCCC > 1) {
    std::complex<double> a;
    a = -acos(rho_CCCC);
    u = std::real(a);
  } else {
    u = -acos(rho_CCCC);
  }
  const double A = sin(u);
  const double B = cos(u) - 2;
  const double t = atan2(eta * A - xi * B, xi * A + eta * B);
  const double v = (t + thetaf);
  d = fabs(t) + 2 * fabs(u) + fabs(v);
  lengths_[16][0] = t;
  lengths_[16][1] = u;
  lengths_[16][2] = u;
  lengths_[16][3] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P18(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = x + sin(thetaf);
  const double eta = -y - 1 + cos(-thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double theta = atan2(eta, xi);
  const double u = -2 * asin(0.25 * rho);
  const double t = wrapToPi(theta + 0.5 * u + pi);
  const double v = wrapToPi(-thetaf - t + u);

  d = fabs(t) + fabs(u) + fabs(v);

  lengths_[17][0] = t;
  lengths_[17][1] = u;
  lengths_[17][2] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P19(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = -x + sin(-thetaf);
  const double eta = y - 1 - cos(-thetaf);
  const double rho = 0.25 * (2 + sqrt(xi * xi + eta * eta));
  const double u = acos(rho);
  const double delta = wrapToPi(2 * u);
  const double A = sin(u) - sin(delta);
  const double B = cos(u) - cos(delta) - 1;
  const double t = atan2(eta * A - xi * B, xi * A + eta * B);
  const double v = (t - 2 * u + thetaf);

  d = fabs(t) + 2 * fabs(u) + fabs(v);

  lengths_[18][0] = -t;
  lengths_[18][1] = -u;
  lengths_[18][2] = u;
  lengths_[18][3] = -v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::P20(const double x, const double y, const double thetaf,
                        double &d) {
  const double xi = x + sin(thetaf);
  const double eta = y - 1 - cos(thetaf);
  const double rho = 0.25 * (2 + sqrt(xi * xi + eta * eta));
  const double u = acos(rho);
  const double delta = wrapToPi(2 * u);
  const double A = sin(u) - sin(delta);
  const double B = cos(u) - cos(delta) - 1;
  const double t1 = atan2(eta * A - xi * B, xi * A + eta * B);
  const double t = (t1);
  const double v = (t - 2 * u - thetaf);

  d = fabs(t) + 2 * fabs(u) + fabs(v);

  lengths_[19][0] = t;
  lengths_[19][1] = u;
  lengths_[19][2] = -u;
  lengths_[19][3] = v;

  // std::cout << "T20: " << (t < 0 ? "-" : "+") << "t, " << (u < 0 ? "-" : "+")
  //           << "u, " << (-u < 0 ? "-" : "+") << "u, " << (v < 0 ? "-" : "+")
  //           << "v" << std::endl;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::setA(const double x, const double y, const double thetaf,
                  const double r, const double beta0, int &cond, double &d,
                  const double x0, const double xn, const double yn) {
  // Set A
  if (thetaf >= 0) {
    // 7, 8 - singled out
    if (cfLy_ <= c0Ly_ && cfRy_ <= c0Ly_) {
      double t2 = (c0Rx_ - xn) * cos(thetaf) + (c0Ry_ - yn) * sin(thetaf);
      const double t2x = t2 * cos(thetaf) + xn;
      const double t2y = t2 * sin(thetaf) + yn;
      const double d1 = euclideanDistance(t2x, t2y, c0Rx_, c0Ry_);
      if (t2 <= -2 * r || d1 <= r) {
        // P7
        P7(x, y, thetaf, d);
        d *= r;
        cond = 7;
        return;
      }
      // P8
      P8(x, y, thetaf, d);
      d *= r;
      cond = 8;
      return;
    }

    // 1, 9, 10, 11
    const double LfL0 = atan2(cfLy_ - c0Ly_, cfLx_ - c0Lx_);
    if (thetaf < fabs(LfL0)) { // omega
      // 11 directly singled out
      if (thetaf > atan2(cfLy_ - c0Ry_, cfLx_ - c0Rx_)) {
        // CCSC | -R+L+S+L
        P11(x, y, thetaf, d);
        d *= r;
        cond = 11;
        return;
      }

      // 1 , 9, 10
      // 1
      double u;
      if (cfRx_ >= 2 * r + x0 || cfRy_ <= c0Ly_) {
        // CSC | +L+S+R
        P1(x, y, thetaf, d);
        d *= r;
        cond = 1;
        return;
      }
      double t2 = (c0Rx_ - xn) * cos(thetaf) + (c0Ry_ - yn) * sin(thetaf);
      // 9
      if (fabs(t2) <= 2 * r) {
        // CCSCC | -R+L+S+R-L
        P9(x, y, thetaf, d);
        d *= r;
        cond = 9;
        return;
      }
      // 10
      // CCSC | -R+L+S+R
      P10(x, y, thetaf, d);
      d *= r;
      cond = 10;
      return;
    }

    // thetaf >= LfL0
    if (cfLx_ < 0) {
      // CCSC | -R+L+S+L
      P11(x, y, thetaf, d);
      d *= r;
      cond = 11;
      return;
    }
    // 2, 3
    if (thetaf > LfL0 + pi / 2) {
      P3(x, y, thetaf, d);
      d *= r;
      cond = 3;
      return;
    } else {
      P2(x, y, thetaf, d);
      d *= r;
      cond = 2;
      return;
    }
  }
  // const double RfL0 = atan2(cfRy_ - c0Ly_, cfRx_ - c0Lx_);
  // const double LfR0 = atan2(cfLy_ - c0Ry_, cfLx_ - c0Rx_);
  // if (thetaf < 2 * beta0 - pi) {
  // } else {
  //   // std::cout << (LfR0 >= RfL0) << std::endl;
  // }

  // thetaf < 0
  double u;
  // 5, 6, 12
  if (thetaf < 2 * beta0 - pi) {
    // 6
    const double R0Lf = atan2(c0Ry_ - cfLy_, c0Rx_ - cfLx_);
    if (thetaf < R0Lf) {
      // CCSC | +R-L-S-L
      P6(x, y, thetaf, d);
      d *= r;
      cond = 6;
      return;
    }
    // 5
    double t2 = (c0Rx_ - xn) * cos(thetaf) + (c0Ry_ - yn) * sin(thetaf);
    if (fabs(t2) <= 2 * r) {
      // 12
      // CCSCC | +R-L-S-R+L
      P12(x, y, thetaf, d);
      d *= r;
      cond = 12;
      return;
    }
    // CCSC | +R-L-S-R
    P5(x, y, thetaf, d);
    d *= r;
    cond = 5;
    return;
  }

  // 1
  const double RfL0 = atan2(cfRy_ - c0Ly_, cfRx_ - c0Lx_);
  double t1 = (c0Lx_ - xn) * cos(thetaf) + (c0Ly_ - yn) * sin(thetaf);
  if (RfL0 <= thetaf || t1 <= -2 * r) {
    // CSC | +L+S+R
    P1(x, y, thetaf, d);
    d *= r;
    cond = 1;
    return;
  }
  if (cfLx_ >= 2 * r) {
    // 4
    P4(x, y, thetaf, d);
    d *= r;
    cond = 4;
    return;
  }
  // 9
  P9(x, y, thetaf, d);
  d *= r;
  cond = 9;
  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::setB(const double x, const double y, const double thetaf,
                  const double xn, const double yn, const double r, int &cond,
                  double &d, const double beta0) {
  // 9, 12
  if (RL >= sqrt(20) * r) {
    // const double RfL0 = atan2(cfRy_ - c0Ly_, cfRx_ - c0Lx_);
    // const double LfR0 = atan2(cfLy_ - c0Ry_, cfLx_ - c0Rx_);
    // (LfR0 >= RfL0)
    if (thetaf > 2 * beta0 - pi) {
      // CCSCC
      P9(x, y, thetaf, d);
      d *= r;
      cond = 9;
      return;
    }
    // CCSCC
    P12(x, y, thetaf, d);
    d *= r;
    cond = 12;
    return;
  }

  if (thetaf >= 0) {
    // cond 13, 14, 19
    if (thetaf < pi / 2) {
      double alpha = acos((3 * r * r + 0.25 * RL * RL) / (2 * r * RL));
      const double R0Lf = atan2(c0Ry_ - cfLy_, c0Rx_ - cfLx_);
      double beta = thetaf - pi / 2 - R0Lf;

      // 13, 19
      if (alpha >= beta) {
        const double c_13 = sqrt((RL * RL + 4 * r * r) / 2 - 4 * r * r);
        const double beta_13 =
            acos((4 * r * r + RL * RL - c_13 * c_13) / (4 * r * RL));
        const double T_13 = R0Lf + pi / 2 + beta_13;
        const double U_13 = acos((8 * r * r - c_13 * c_13) / (8 * r * r));

        const double U_19 = acos((0.5 * LR + r) / (2 * r));
        const double L0Rf = atan2(c0Ly_ - cfRy_, c0Lx_ - cfRx_);
        const double beta_19 = thetaf + pi / 2 - L0Rf;
        const double V_19 = wrapTo2Pi(U_19 - beta_19);

        if (T_13 <= V_19 || T_13 + U_13 <= 2 * U_19) {
          // U_13 is always greater than U_19 here
          d = r * (2 * T_13 - thetaf + 2 * U_13);
          cond = 13;
          const double V_13 = -T_13 + thetaf;
          lengths_[12][0] = -T_13;
          lengths_[12][1] = U_13;
          lengths_[12][2] = U_13;
          lengths_[12][3] = V_13;
          return;
        }
        d = r * (4 * U_19 - thetaf);
        cond = 19;
        const double T_19 = 2 * U_19 - V_19 - thetaf;
        lengths_[18][0] = -T_19;
        lengths_[18][1] = -U_19;
        lengths_[18][2] = U_19;
        lengths_[18][3] = V_19;
        return;
      }

      // 14, 19
      double gamma = acos((0.5 * LR + r) / (2 * r));
      double beta_3 = atan2(cfRy_ - c0Ly_, cfRx_ - c0Lx_) + pi / 2;
      if (RL <= 2 * r || beta_3 >= gamma) {
        // CCC
        P14(x, y, thetaf, d);
        d *= r;
        cond = 14;
        return;
      }
      // CCC or CCCC
      P19(x, y, thetaf, d);
      d *= r;
      cond = 19;
      return;
    }

    // thetaf > pi/2
    // 14, 15 split
    // CCC
    if (LR <= 2 * r && RL <= 2 * r) {
      P15(x, y, thetaf, d);
      d *= r;
      cond = 15;
      return;
    }
    P14(x, y, thetaf, d);
    d *= r;
    cond = 14;
    return;
  }

  // thetaf < 0, two cases
  double alpha;
  if (RL < 2.0 * r) {
    alpha = 0;
  } else {
    alpha = acos((3 * r * r + 0.25 * RL * RL) / (2 * r * RL));
  }
  const double LfR0 = atan2(cfLy_ - c0Ry_, cfLx_ - c0Rx_);
  const double beta_1 = pi / 2 - LfR0;
  const double beta_2 = -thetaf - beta_1;
  // Case 2: 13 20, 16, 20
  if (thetaf >= 2 * LfR0 - pi) {
    if (alpha > beta_1) {
      // CCCC two possibilities
      const double rl = LfR0 - pi;
      const double c_13 = sqrt((RL * RL + 4 * r * r) / 2 - 4 * r * r);
      const double beta_13 =
          acos((4 * r * r + RL * RL - c_13 * c_13) / (4 * r * RL));
      const double T_13 = rl + pi / 2 + beta_13;
      const double V_13 = T_13 - thetaf;
      const double U_13 = acos((8 * r * r - c_13 * c_13) / (8 * r * r));

      const double U_20 = acos((0.5 * LR + r) / (2 * r));
      const double L0Rf = atan2(c0Ly_ - cfRy_, c0Lx_ - cfRx_);
      const double T_20 = wrapTo2Pi(U_20 - (pi / 2 - L0Rf));

      if (V_13 <= T_20 || T_13 + U_13 <= thetaf + 2 * U_20) {
        d = r * (2 * T_13 - thetaf + 2 * U_13);
        cond = 13;
        lengths_[12][0] = -T_13;
        lengths_[12][1] = U_13;
        lengths_[12][2] = U_13;
        lengths_[12][3] = -V_13;
        return;
      }
      d = r * (4 * U_20 + thetaf);
      cond = 20;
      lengths_[19][0] = T_20;
      lengths_[19][1] = U_20;
      lengths_[19][2] = -U_20;
      const double V_20 = thetaf + 2 * U_20 - T_20;
      lengths_[19][3] = -V_20;
      return;
    }
    // CCC or CCCC
    const double gamma = acos((0.5 * LR + r) / (2 * r));
    const double O = 4 * r * sin(gamma / 2.0);
    const bool cond_1 = O <= LL;
    const bool cond_2 = O <= RR;
    const bool cond_3 = cond_1 || cond_2;
    if (cond_3) {
      P16(x, y, thetaf, d);
      d *= r;
      cond = 16;
      return;
    }
    P20(x, y, thetaf, d);
    d *= r;
    cond = 20;
    return;
  }

  // Case 1: 17, 20, 18, 20
  // beta_1 <= beta2 || d_ < 2 * r - 17, 20
  if (alpha >= beta_2) {
    // CCCC two possibilities
    const double c_17 = sqrt((RL * RL + 4 * r * r) / 2 - 4 * r * r);
    const double U_17 = acos((8 * r * r - c_17 * c_17) / (8 * r * r));
    const double beta_17 =
        acos((4 * r * r + RL * RL - c_17 * c_17) / (4 * r * RL));
    const double T_17 = pi / 2 - LfR0 + beta_17;

    const double U_20 = acos((0.5 * LR + r) / (2 * r));
    const double L0Rf = atan2(c0Ly_ - cfRy_, c0Lx_ - cfRx_);
    const double T_20 = wrapTo2Pi(U_20 - (pi / 2 - L0Rf));

    if (T_17 <= T_20 || T_17 + U_17 <= 2 * U_20) {
      d = r * (2 * T_17 + thetaf + 2 * U_17);
      const double V_17 = T_17 + thetaf;
      lengths_[16][0] = T_17;
      lengths_[16][1] = -U_17;
      lengths_[16][2] = -U_17;
      lengths_[16][3] = V_17;
      cond = 17;
      return;
    }
    d = r * (4 * U_20 + thetaf);
    cond = 20;
    lengths_[19][0] = T_20;
    lengths_[19][1] = U_20;
    lengths_[19][2] = -U_20;
    const double V_20 = thetaf + 2 * U_20 - T_20;
    lengths_[19][3] = -V_20;
    return;
  }
  // alpha < beta_2 - 18 20
  const double gamma = acos((0.5 * LR + r) / (2 * r));
  const double O = 4 * r * sin(gamma / 2.0);
  if (O <= RR || RL <= 2 * r) {
    P18(x, y, thetaf, d);
    cond = 18;
  } else {
    P20(x, y, thetaf, d);
    cond = 20;
  }
  d *= r;
  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::acceleratedRSPlanner(const State &from, const State &to, double &d,
                                  const double r, int &cond, int &Q) {
  const double x0 = from.x;
  const double y0 = from.y;
  const double theta0g = from.theta;
  const double xng = to.x;
  const double yng = to.y;
  const double thetafg = to.theta;

  // Projection into local rotational frame
  double xn = x0 - (x0 - xng) * cos(theta0g) - (y0 - yng) * sin(theta0g);
  double yn = y0 + (x0 - xng) * sin(theta0g) - (y0 - yng) * cos(theta0g);
  double thetaf = (thetafg - theta0g);

  // Project into Q1 - obtaining the mirrored configuration
  const double dx_ = xn - x0;
  const double dy_ = yn - y0;
  const double theta0 = 0;
  Q = 1;
  if (dx_ <= 0 && dy_ >= 0) {
    xn = 2 * x0 - xn;
    thetaf = 2 * pi - (thetaf);
    Q = 2;
  } else if (dx_ <= 0 && dy_ <= 0) {
    xn = 2 * x0 - xn;
    yn = 2 * y0 - yn;
    Q = 3;
  } else if (dx_ >= 0 && dy_ <= 0) {
    yn = 2 * y0 - yn;
    thetaf = 2 * pi - (thetaf);
    Q = 4;
  }
  thetaf = wrapToPi(thetaf);

  // LHC and RHC
  c0Rx_ = x0 + r * cos(theta0 - pi / 2);
  c0Ry_ = y0 + r * sin(theta0 - pi / 2);
  c0Lx_ = x0 - r * cos(theta0 - pi / 2);
  c0Ly_ = y0 - r * sin(theta0 - pi / 2);
  cfRx_ = xn + r * cos(thetaf - pi / 2);
  cfRy_ = yn + r * sin(thetaf - pi / 2);
  cfLx_ = xn - r * cos(thetaf - pi / 2);
  cfLy_ = yn - r * sin(thetaf - pi / 2);

  // Algorithm 3 (IsInSetB)
  limit = 2 * r * sqrt(2);
  LL = euclideanDistance(c0Lx_, c0Ly_, cfLx_, cfLy_);
  LR = euclideanDistance(c0Lx_, c0Ly_, cfRx_, cfRy_);
  RL = euclideanDistance(c0Rx_, c0Ry_, cfLx_, cfLy_);
  RR = euclideanDistance(c0Rx_, c0Ry_, cfRx_, cfRy_);
  const bool pred_1 = RR <= limit && LL <= limit && LR <= 2 * r;
  const bool pred_2 = RR <= limit && LL <= limit && RL <= 2 * r;
  const bool pred_3 = LR <= 2 * r && LL <= limit && RL <= 2 * r;

  const double x = (xn - x0) / r;
  const double y = (yn - y0) / r;
  const double beta0 = atan2(yn - y0, xn - x0);

  // Main Algorithm (4)
  if (pred_1 || pred_2 || pred_3) {
    setB(x, y, thetaf, xn, yn, r, cond, d, beta0);
    return;
  } else {
    setA(x, y, thetaf, r, beta0, cond, d, x0, xn, yn);
    return;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::evaluateQuery(const State &from, const State &to,
                           const double radius) {

  std::cout << "\n******************" << std::endl;
  std::cout << "State_from: {" << from.x << ", " << from.y << ", " << from.theta
            << "}" << std::endl;
  std::cout << "State_to: {" << to.x << ", " << to.y << ", " << to.theta << "}"
            << std::endl;

  const double r = radius;
  // Reeds-Shepp State Space (OMPL's implementation)
  ob::StateSpacePtr space(std::make_shared<ob::ReedsSheppStateSpace>(r));
  ob::ScopedState<> fromRSS(space), toRSS(space);

  double d_accelerated = 0;
  int condition = 0;
  int Q = 0;
  acceleratedRSPlanner(from, to, d_accelerated, r, condition, Q);

  std::cout << "******************" << std::endl;
  std::cout << "AcceleratedRSPlanner condition: " << condition << std::endl;
  std::cout << "******************" << std::endl;
  for (int i = 0; i < 5; ++i) {
    std::cout << lengths_[condition - 1][i] << " ";
  }

  fromRSS[0] = from.x;
  fromRSS[1] = from.y;
  fromRSS[2] = from.theta;
  toRSS[0] = to.x;
  toRSS[1] = to.y;
  toRSS[2] = to.theta;

  ob::ReedsSheppStateSpace::ReedsSheppPath reedsSheppPath =
      space->as<ob::ReedsSheppStateSpace>()->reedsShepp(fromRSS(), toRSS());
  std::cout << "\n"
            << reedsSheppPath.length_[0] << " " << reedsSheppPath.length_[1]
            << " " << reedsSheppPath.length_[2] << " "
            << reedsSheppPath.length_[3] << " " << reedsSheppPath.length_[4]
            << std::endl;

  std::cout << "******************" << std::endl;
  std::cout << "AcceleratedRSPlanner path length: " << d_accelerated
            << std::endl;
  std::cout << "OMPL path length: " << reedsSheppPath.length() * r << std::endl;
  State errors;
  double d = 0;
  acceleratedRSPlanner(from, to, d, r, condition, Q);
  std::cout << "Condition: " << condition << std::endl;

  char motionTypes[5] = {'N', 'N', 'N', 'N', 'N'};
  getMotionTypes(condition, motionTypes);

  std::cout << "Quadrant: " << Q << std::endl;
  backProjectMotion(condition, Q, motionTypes);

  std::vector<State> path;
  path.push_back(from);

  getPath(condition, motionTypes, path, r);
  getErrors(path, to, errors);

  std::cout << "Path: " << std::endl;
  for (int i = 0; i < path.size(); ++i) {
    std::cout << path[i].x << " " << path[i].y << " " << path[i].theta
              << std::endl;
  }

  std::cout << "******************" << std::endl;
  std::cout << "Errors: " << errors.x << " " << errors.y << " " << errors.theta
            << std::endl;
  double errors_norm = sqrt(errors.x * errors.x + errors.y * errors.y +
                            errors.theta * errors.theta);
  std::cout << "Errors norm ||end-start||: " << errors_norm << std::endl;

  if (errors_norm > 1e-10) {
    std::cout << "Error in validity check!" << std::endl;
  }

  // cast motion types to string
  std::string motionTypesStr;
  for (int i = 0; i < 5; ++i) {
    motionTypesStr += motionTypes[i];
  }
  std::cout << "******************" << std::endl;
  std::cout << "Drawing motion type: ";
  for (int i = 0; i < motionTypesStr.size(); ++i) {
    std::cout << motionTypesStr[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "******************" << std::endl;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::benchmarkPlanners() {
  const int number_of_final_states = 1000000;
  // std::random_device rd;
  // std::mt19937 gen(rd());
  // std::uniform_real_distribution<> r_dis(0.01, 10);
  // const double r = r_dis(gen);
  const double r = 1;

  std::vector<State> toStates;
  toStates.reserve(number_of_final_states);
  const bool use_config = false;
  helper_functions_.generateRandomStates(toStates, use_config, r);

  std::vector<double> distances_acceleratedRS(number_of_final_states);
  std::vector<double> distances_Desaulniers(number_of_final_states);
  std::vector<double> distances_OMPL(number_of_final_states);

  const State from = {0, 0, pi / 2};

  // AcceleratedRS (Proposed)
  double d_accelerated = 0;
  State to;
  int Q, condition;
  auto start_accelerated = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < number_of_final_states; ++i) {
    to.x = toStates[i].x;
    to.y = toStates[i].y;
    to.theta = toStates[i].theta - pi / 2;
    acceleratedRSPlanner(from, to, d_accelerated, r, condition, Q);
    distances_acceleratedRS[i] = d_accelerated;
  }
  auto end_accelerated = std::chrono::high_resolution_clock::now();
  auto duration_accelerated =
      durationInMicroseconds(start_accelerated, end_accelerated);

  // Reeds-Shepp State Space (OMPL's implementation)
  ob::StateSpacePtr space(std::make_shared<ob::ReedsSheppStateSpace>(r));
  ob::ScopedState<> fromOMPL(space), toOMPL(space);
  fromOMPL[0] = from.x;
  fromOMPL[1] = from.y;
  fromOMPL[2] = from.theta;
  double d_OMPL = 0;
  auto start_OMPL = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < number_of_final_states; ++i) {
    toOMPL[0] = toStates[i].x;
    toOMPL[1] = toStates[i].y;
    toOMPL[2] = toStates[i].theta - pi / 2;
    d_OMPL = space->distance(fromOMPL(), toOMPL());
    distances_OMPL[i] = d_OMPL;
  }
  auto end_OMPL = std::chrono::high_resolution_clock::now();
  auto duration_OMPL = durationInMicroseconds(start_OMPL, end_OMPL);

  // Our implmentation of An Efficient Algorithm to find a shortest path for a
  // car-like robot
  DesaulniersAlgorithm desaulniersAlgorithm;
  double d_Desaulniers = 0;
  auto start_Desaulniers = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < number_of_final_states; ++i) {
    d_Desaulniers = desaulniersAlgorithm.getDistance(
        toStates[i].x, toStates[i].y, toStates[i].theta);
    distances_Desaulniers[i] = d_Desaulniers;
  }
  auto end_Desaulniers = std::chrono::high_resolution_clock::now();
  auto duration_Desaulniers =
      durationInMicroseconds(start_Desaulniers, end_Desaulniers);

  helper_functions_.timeStats(duration_accelerated, duration_OMPL,
                              duration_Desaulniers, number_of_final_states);
  helper_functions_.errorStats(distances_acceleratedRS, distances_OMPL,
                               distances_Desaulniers);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::computeAndDrawRandomPath(sf::RenderWindow &window) {
  const double r = config_.r;

  const double min_x = config_.min_x;
  const double max_x = config_.max_x;
  const double min_y = config_.min_y;
  const double max_y = config_.max_y;
  const double min_thetaf = config_.min_thetaf * pi / 180;
  const double max_thetaf = config_.max_thetaf * pi / 180;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> x_dis(min_x, max_x);
  std::uniform_real_distribution<> y_dis(min_y, max_y);
  std::uniform_real_distribution<> theta_dis(min_thetaf, max_thetaf);

  std::cout << "\nGenerating random start and end states" << std::endl;
  std::cout << "x: [" << min_x << ", " << max_x << "]" << std::endl;
  std::cout << "y: [" << min_y << ", " << max_y << "]" << std::endl;
  std::cout << "theta: [" << min_thetaf << ", " << max_thetaf << "]"
            << std::endl;
  std::cout << "r: " << r << std::endl;

  State from, to;
  from.x = x_dis(gen);
  from.y = y_dis(gen);
  from.theta = theta_dis(gen);
  to.x = x_dis(gen);
  to.y = y_dis(gen);
  to.theta = theta_dis(gen);

  std::cout << "\n******************" << std::endl;
  std::cout << "State_from: {" << from.x << ", " << from.y << ", " << from.theta
            << "}" << std::endl;
  std::cout << "State_to: {" << to.x << ", " << to.y << ", " << to.theta << "}"
            << std::endl;

  // Reeds-Shepp State Space (OMPL's implementation)
  ob::StateSpacePtr space(std::make_shared<ob::ReedsSheppStateSpace>(r));
  ob::ScopedState<> fromRSS(space), toRSS(space);

  double d_accelerated = 0;
  int condition = 0;
  int Q = 0;
  acceleratedRSPlanner(from, to, d_accelerated, r, condition, Q);

  std::cout << "******************" << std::endl;
  std::cout << "AcceleratedRSPlanner condition: " << condition << std::endl;
  std::cout << "******************" << std::endl;
  for (int i = 0; i < 5; ++i) {
    std::cout << lengths_[condition - 1][i] << " ";
  }

  fromRSS[0] = from.x;
  fromRSS[1] = from.y;
  fromRSS[2] = from.theta;
  toRSS[0] = to.x;
  toRSS[1] = to.y;
  toRSS[2] = to.theta;

  ob::ReedsSheppStateSpace::ReedsSheppPath reedsSheppPath =
      space->as<ob::ReedsSheppStateSpace>()->reedsShepp(fromRSS(), toRSS());
  std::cout << "\n"
            << reedsSheppPath.length_[0] << " " << reedsSheppPath.length_[1]
            << " " << reedsSheppPath.length_[2] << " "
            << reedsSheppPath.length_[3] << " " << reedsSheppPath.length_[4]
            << std::endl;

  std::cout << "******************" << std::endl;
  std::cout << "AcceleratedRSPlanner path length: " << d_accelerated
            << std::endl;
  std::cout << "OMPL path length: " << reedsSheppPath.length() * r << std::endl;
  State errors;

  char motionTypes[5] = {'N', 'N', 'N', 'N', 'N'};
  getMotionTypes(condition, motionTypes);

  std::cout << "Quadrant: " << Q << std::endl;
  backProjectMotion(condition, Q, motionTypes);
  std::vector<State> path;
  path.push_back(from);
  getPath(condition, motionTypes, path, r);
  getErrors(path, to, errors);
  std::cout << "******************" << std::endl;
  std::cout << "Errors: " << errors.x << " " << errors.y << " " << errors.theta
            << std::endl;
  double errors_norm = sqrt(errors.x * errors.x + errors.y * errors.y +
                            errors.theta * errors.theta);
  std::cout << "Errors norm ||end-start||: " << errors_norm << std::endl;

  if (errors_norm > 1e-10) {
    std::cout << "Error in validity check!" << std::endl;
  }

  std::vector<char> motionTypesVector(motionTypes, motionTypes + 5);
  std::vector<double> motionLengths = {
      lengths_[condition - 1][0], lengths_[condition - 1][1],
      lengths_[condition - 1][2], lengths_[condition - 1][3],
      lengths_[condition - 1][4]};
  plotter.draw(window, from, to, path, motionTypesVector, motionLengths,
               condition);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::drawRandomPath() {
  std::cout << "Press enter to generate and plot new random path" << std::endl;
  sf::RenderWindow window(sf::VideoMode(1000, 1000), "Path Visualization");
  window.setFramerateLimit(60);
  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed ||
          event.key.code == sf::Keyboard::Escape) {
        window.close();
      }

      if (event.type == sf::Event::KeyPressed &&
          event.key.code == sf::Keyboard::Enter) {
        // Clear previous drawings
        window.clear();
        computeAndDrawRandomPath(window);

        // Display the window
        window.display();
      }
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::RandomPathsValidityChecks() {
  const bool use_config = true;
  const int number_of_final_states = config_.number_of_final_states;
  const double r = config_.r;
  // const double r between 0 and config_.r
  // std::random_device rd;
  // std::mt19937 gen(rd());
  // std::uniform_real_distribution<> r_dis(0.01, config_.r);
  // double r = r_dis(gen);

  std::vector<State> toStates;
  toStates.reserve(number_of_final_states);
  helper_functions_.generateRandomStates(toStates, use_config);
  State from;
  std::vector<State> fromStates;
  fromStates.reserve(number_of_final_states);
  helper_functions_.generateRandomStates(fromStates, use_config);

  double d_accelerated = 0;
  int condition = 0;
  int Q = 0;
  for (int i = 0; i < number_of_final_states; ++i) {
    std::cout << "******************************************" << std::endl;
    std::cout << "Checking state: " << i + 1 << std::endl;
    std::cout << "State_from: {" << fromStates[i].x << ", " << fromStates[i].y
              << ", " << fromStates[i].theta << "}" << std::endl;
    std::cout << "State_to: {" << toStates[i].x << ", " << toStates[i].y << ", "
              << toStates[i].theta << "}" << std::endl;
    // r = r_dis(gen);
    // std::cout << "State_radius: " << r << std::endl;
    from.x = fromStates[i].x;
    from.y = fromStates[i].y;
    from.theta = fromStates[i].theta;
    // AcceleratedRS (Proposed)
    acceleratedRSPlanner(from, toStates[i], d_accelerated, r, condition, Q);

    std::cout << "AcceleratedRSPlanner condition: " << condition << std::endl;
    for (int i = 0; i < 5; ++i) {
      std::cout << lengths_[condition - 1][i] << " ";
    }

    // Reeds-Shepp State Space (OMPL's implementation)
    ob::StateSpacePtr space(std::make_shared<ob::ReedsSheppStateSpace>(r));
    ob::ScopedState<> fromRSS(space), toRSS(space);
    toRSS[0] = toStates[i].x;
    toRSS[1] = toStates[i].y;
    toRSS[2] = toStates[i].theta;
    fromRSS[0] = fromStates[i].x;
    fromRSS[1] = fromStates[i].y;
    fromRSS[2] = fromStates[i].theta;
    ob::ReedsSheppStateSpace::ReedsSheppPath reedsSheppPath =
        space->as<ob::ReedsSheppStateSpace>()->reedsShepp(fromRSS(), toRSS());

    std::cout << "\n"
              << reedsSheppPath.length_[0] << " " << reedsSheppPath.length_[1]
              << " " << reedsSheppPath.length_[2] << " "
              << reedsSheppPath.length_[3] << " " << reedsSheppPath.length_[4]
              << std::endl;
    std::cout << "AcceleratedRSPlanner path length: " << d_accelerated
              << std::endl;
    std::cout << "OMPL path length: " << reedsSheppPath.length() * r
              << std::endl;

    State errors;
    validityCheck(from, toStates[i], condition, Q, errors, r);
    double errors_norm = sqrt(errors.x * errors.x + errors.y * errors.y +
                              errors.theta * errors.theta);
    std::cout << "Errors norm ||end-start||: " << errors_norm << std::endl;
    if (errors_norm > 1e-10) {
      std::cout << errors.x << " " << errors.y << " " << errors.theta
                << std::endl;
      std::cout << "Error in validity check" << std::endl;
      std::cout << "******************************************" << std::endl;
      break;
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::validityCheck(const State &from, const State &to, const int cond,
                           const int Q, State &errors, const double r) {
  char motionTypes[5] = {'N', 'N', 'N', 'N', 'N'};
  getMotionTypes(cond, motionTypes);
  backProjectMotion(cond, Q, motionTypes);
  std::vector<State> path;
  path.push_back(from);
  getPath(cond, motionTypes, path, r);
  getErrors(path, to, errors);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::getMotionTypes(const int condition, char motionTypes[5]) const {
  switch (condition) {
  case 1:
    motionTypes[0] = 'L';
    motionTypes[1] = 'S';
    motionTypes[2] = 'R';
    break;
  case 2:
    motionTypes[0] = 'L';
    motionTypes[1] = 'S';
    motionTypes[2] = 'L';
    break;
  case 3:
    motionTypes[0] = 'L';
    motionTypes[1] = 'S';
    motionTypes[2] = 'L';
    motionTypes[3] = 'R';
    break;
  case 4:
    motionTypes[0] = 'L';
    motionTypes[1] = 'S';
    motionTypes[2] = 'R';
    motionTypes[3] = 'L';
    break;
  case 5:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'S';
    motionTypes[3] = 'R';
    break;
  case 6:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'S';
    motionTypes[3] = 'L';
    break;
  case 7:
    motionTypes[0] = 'R';
    motionTypes[1] = 'S';
    motionTypes[2] = 'L';
    break;
  case 8:
    motionTypes[0] = 'R';
    motionTypes[1] = 'S';
    motionTypes[2] = 'L';
    motionTypes[3] = 'R';
    break;
  case 9:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'S';
    motionTypes[3] = 'R';
    motionTypes[4] = 'L';
    break;
  case 10:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'S';
    motionTypes[3] = 'R';
    break;
  case 11:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'S';
    motionTypes[3] = 'L';
    break;
  case 12:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'S';
    motionTypes[3] = 'R';
    motionTypes[4] = 'L';
    break;
  case 13:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'R';
    motionTypes[3] = 'L';
    break;
  case 14:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'R';
    break;
  case 15:
    motionTypes[0] = 'L';
    motionTypes[1] = 'R';
    motionTypes[2] = 'L';
    break;
  case 16:
    motionTypes[0] = 'L';
    motionTypes[1] = 'R';
    motionTypes[2] = 'L';
    break;
  case 17:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'R';
    motionTypes[3] = 'L';
    break;
  case 18:
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'R';
    break;
  case 19:
    motionTypes[0] = 'L';
    motionTypes[1] = 'R';
    motionTypes[2] = 'L';
    motionTypes[3] = 'R';
    break;
  case 20:
    motionTypes[0] = 'L';
    motionTypes[1] = 'R';
    motionTypes[2] = 'L';
    motionTypes[3] = 'R';
    break;
  default:
    std::cout << "Invalid condition" << std::endl;
    break;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::backProjectMotion(const int condition, const int Q,
                               char motionTypes[5]) {
  if (Q == 1) {
    return;
  };
  for (int k = 0; k < 5; ++k) {
    if (motionTypes[k] == 'N') {
      break;
    };
    if (Q == 4) {
      if (motionTypes[k] == 'L') {
        motionTypes[k] = 'R';
      } else if (motionTypes[k] == 'R') {
        motionTypes[k] = 'L';
      }
    } else if (Q == 3) {
      lengths_[condition - 1][k] = -lengths_[condition - 1][k];
      if (motionTypes[k] == 'L') {
        motionTypes[k] = 'R';
      } else if (motionTypes[k] == 'R') {
        motionTypes[k] = 'L';
      }
    } else if (Q == 2) {
      lengths_[condition - 1][k] = -lengths_[condition - 1][k];
    }
  }
};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::getPath(const int condition, char motionTypes[5],
                     std::vector<State> &path, const double r) const {
  double x = path[0].x;
  double y = path[0].y;
  double t = path[0].theta;
  for (int k = 0; k < 5; ++k) {
    char motionType = motionTypes[k];
    double motionLength = fabs(lengths_[condition - 1][k]);
    int motionDirection = (lengths_[condition - 1][k] > 0) ? 1 : -1;
    if (motionType == 'N') {
      break;
    } else if (motionType == 'L') {
      const double cx = x - r * cos(t - pi / 2);
      const double cy = y - r * sin(t - pi / 2);
      if (motionDirection == 1) {
        t = t + motionLength;
      } else {
        t = t - motionLength;
      }
      x = cx + r * cos(t - pi / 2);
      y = cy + r * sin(t - pi / 2);
    } else if (motionType == 'R') {
      const double cx = x + r * cos(t - pi / 2);
      const double cy = y + r * sin(t - pi / 2);
      if (motionDirection == 1) {
        t = t - motionLength;
      } else {
        t = t + motionLength;
      }
      x = cx - r * cos(t - pi / 2);
      y = cy - r * sin(t - pi / 2);
    } else if (motionType == 'S') {
      x = x + r * motionDirection * motionLength * cos(t);
      y = y + r * motionDirection * motionLength * sin(t);
    }
    path.push_back({x, y, t});
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::getErrors(const std::vector<State> &path, const State &to,
                       State &errors) const {
  const int lastIdx = path.size() - 1;
  double x = path[lastIdx].x;
  double y = path[lastIdx].y;
  double t = path[lastIdx].theta;

  errors.x = fabs(x - to.x);
  errors.y = fabs(y - to.y);
  t = wrapToPi(t);
  errors.theta = fabs(t - to.theta);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::characterizeVariability() {
  const int number_of_final_states = 100000000;
  std::vector<State> toStates;
  toStates.reserve(number_of_final_states);
  const bool use_config = true;
  helper_functions_.generateRandomStates(toStates, use_config);

  std::vector<std::vector<State>> stateCategory(20);
  const State from = {0, 0, 0};

  const double r = 400;
  // Catergorize the final states into 20 categories
  double d_accelerated = 0;
  State to;
  int Q, condition;
  for (int i = 0; i < number_of_final_states; ++i) {
    to.x = toStates[i].x;
    to.y = toStates[i].y;
    to.theta = toStates[i].theta - pi / 2;
    acceleratedRSPlanner(from, to, d_accelerated, r, condition, Q);
    stateCategory[condition - 1].push_back(toStates[i]);
  }

  // Display number of states in each category
  for (int i = 0; i < 20; ++i) {
    std::cout << "Category " << i + 1 << " has " << stateCategory[i].size()
              << " states" << std::endl;
  }
  std::cout << std::endl;

  // For each category, compute the time by proposed and OMPL then compare
  std::vector<double> distances_acceleratedRS(number_of_final_states);
  std::vector<double> distances_OMPL(number_of_final_states);

  std::vector<double> time_accelerated(20);

  ob::StateSpacePtr space(std::make_shared<ob::ReedsSheppStateSpace>(r));
  ob::ScopedState<> fromOMPL(space), toOMPL(space);
  fromOMPL[0] = from.x;
  fromOMPL[1] = from.y;
  fromOMPL[2] = from.theta;
  double d_OMPL = 0;

  std::vector<double> time_OMPL(20);
  for (int i = 0; i < 20; ++i) {

    auto start_accelerated = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < stateCategory[i].size(); ++j) {
      to.x = stateCategory[i][j].x;
      to.y = stateCategory[i][j].y;
      to.theta = stateCategory[i][j].theta;
      acceleratedRSPlanner(from, to, d_accelerated, r, condition, Q);
      distances_acceleratedRS[j] = d_accelerated;
    }
    auto end_accelerated = std::chrono::high_resolution_clock::now();
    auto duration_accelerated =
        durationInMicroseconds(start_accelerated, end_accelerated);
    time_accelerated[i] = duration_accelerated;

    auto start_OMPL = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < stateCategory[i].size(); ++j) {
      toOMPL[0] = stateCategory[i][j].x;
      toOMPL[1] = stateCategory[i][j].y;
      toOMPL[2] = stateCategory[i][j].theta;
      d_OMPL = space->distance(fromOMPL(), toOMPL());
      distances_OMPL[j] = d_OMPL;
    }
    auto end_OMPL = std::chrono::high_resolution_clock::now();
    auto duration_OMPL = durationInMicroseconds(start_OMPL, end_OMPL);
    time_OMPL[i] = duration_OMPL;

    // display error stats for each category
    std::cout << "Distance error for category " << i + 1 << std::endl;
    helper_functions_.errorStats(distances_acceleratedRS, distances_OMPL);
  }

  // Display the time speedup for each category
  for (int i = 0; i < 20; ++i) {
    std::cout << "Category " << i + 1 << " has a time speedup of "
              << time_OMPL[i] / time_accelerated[i] << std::endl;
  }

  // Save the number of states in each category, the time for each category for
  // each method and the speedup
  std::ofstream file;
  file.open("variability.txt");
  for (int i = 0; i < 20; ++i) {
    file << stateCategory[i].size() << " " << time_accelerated[i] << " "
         << time_OMPL[i] << " " << time_OMPL[i] / time_accelerated[i] << "\n";
  }

  // weighted speedup factor
  double total_time_accelerated = 0;
  double total_time_OMPL = 0;
  for (int i = 0; i < 20; ++i) {
    total_time_accelerated += time_accelerated[i] * stateCategory[i].size();
    total_time_OMPL += time_OMPL[i] * stateCategory[i].size();
  }
  double weighted_speedup = total_time_OMPL / total_time_accelerated;
  std::cout << "Weighted speedup: " << weighted_speedup << std::endl;

}

} // namespace accelerated
