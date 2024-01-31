#include "solver/Solver.hpp"
#include "classic/DesaulniersAlgorithm.hpp"

#include <SFML/Graphics.hpp>
#include <chrono>
#include <complex>
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
inline void Solver::Beta1(const double x, const double y, const double thetaf,
                          double &t, double &u) const {
  const double dx = x + sin(thetaf);
  const double dy = y - 1 - cos(thetaf);
  const double phi = atan2(dy, dx);
  double rho = dx * dx + dy * dy;
  if (rho < 4) {
    std::complex<double> a;
    a = sqrt(std::complex<double>(rho - 4, 0));
    u = std::real(a);
  } else {
    u = sqrt(rho - 4);
  }
  t = atan2(2, u);
  t = wrapToPi(phi + t);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::Beta4(const double x, const double y, const double thetaf,
                          double &t, double &u, double &v) const {
  const double xb = -(x * cos(thetaf) + y * sin(thetaf));
  const double yb = x * sin(thetaf) - y * cos(thetaf);
  const double xi = xb - sin(-thetaf);
  const double eta = yb - 1 + cos(-thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double phi = atan2(eta, xi);
  const double rr = sqrt(rho * rho - 4);
  u = 2 - rr;
  t = wrapToPi(phi + atan2(rr, -2));
  v = wrapToPi(-thetaf - 0.5 * pi - t);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::Beta5(const double x, const double y, const double thetaf,
                          double &t, double &u) const {
  const double xi = x - sin(-thetaf);
  const double eta = -y - 1 + cos(-thetaf);
  const double rho = sqrt(xi * xi + eta * eta);
  const double phi = atan2(eta, xi);
  const double rr = sqrt(rho * rho - 4);
  u = 2 - rr;
  t = wrapToPi(phi + atan2(rr, -2));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::Beta7(const double x, const double y, const double thetaf,
                          double &t, double &u) const {
  const double dx = x + sin(-thetaf);
  const double dy = -y - 1 - cos(-thetaf);
  const double rho = (dx * dx + dy * dy);
  const double phi = atan2(dy, dx);
  if (rho < 4) {
    std::complex<double> a;
    a = sqrt(std::complex<double>(rho - 4, 0));
    u = std::real(a);
  } else {
    u = sqrt(rho - 4);
  }
  t = wrapToPi(phi + atan2(2, u));
};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::Beta9(const double x, const double y, const double thetaf,
                          double &t, double &u) const {
  const double xi = -x + sin(thetaf);
  const double eta = -y - 1 - cos(thetaf);
  const double rho = xi * xi + eta * eta;
  const double phi = atan2(eta, xi);
  u = 4 - sqrt(rho - 4);
  t = wrapToPi(atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::T3(const double x, const double y, const double thetaf,
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
inline void Solver::T6(const double x, const double y, const double thetaf,
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
inline void Solver::T8(const double x, const double y, const double thetaf,
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
inline void Solver::T9(const double x, const double y, const double thetaf,
                       double &d) {
  double beta9, u;
  Beta9(x, y, thetaf, beta9, u);
  const double v = wrapToPi(beta9 - thetaf);

  d = fabs(beta9) + fabs(u) + fabs(v) + pi;

  lengths_[8][0] = -beta9;
  lengths_[8][1] = pi / 2;
  lengths_[8][2] = -u;
  lengths_[8][3] = pi / 2;
  lengths_[8][4] = -v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::T9(const double xi, const double eta, const double rho,
                       const double thetaf, double &d) {
  const double phi = atan2(eta, xi);
  const double u = 4 - sqrt(rho - 4);
  const double t =
      wrapToPi(atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta));
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
inline void Solver::T10(const double x, const double y, const double thetaf,
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
inline void Solver::T11(const double x, const double y, const double thetaf,
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
inline void Solver::T12(const double x, const double y, const double thetaf,
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
inline void Solver::T13(const double x, const double y, const double thetaf,
                        double &xi, double &eta, double &rho, double &t,
                        double &d) {
  xi = -x + sin(thetaf);
  eta = -y - 1 - cos(thetaf);
  rho = xi * xi + eta * eta;
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
  t = atan2(eta * A - xi * B, xi * A + eta * B);
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
inline void Solver::T13(const double x, const double y, const double thetaf,
                        double &xi, double &eta, double &rho, double &t1,
                        double &t2, double &u1, double &u2, double &v1,
                        double &v2) const {
  xi = x + sin(-thetaf);
  eta = -y - 1 - cos(-thetaf);
  rho = xi * xi + eta * eta;
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
  t1 = atan2(eta * A - xi * B, xi * A + eta * B);
  u1 = u;
  v1 = (t1 + thetaf);
  t2 = atan2(eta * A + xi * B, -xi * A + eta * B);
  u2 = u;
  v2 = (t2 - thetaf);
  xi = -xi;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::T14(const double x, const double y, const double thetaf,
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
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::T15(const double x, const double y, const double thetaf,
                        double &t, double &u, double &v) const {
  const double xi = x - sin(thetaf);
  const double eta = y - 1 + cos(thetaf);
  const double u1 = sqrt(xi * xi + eta * eta);
  const double theta = atan2(eta, xi);
  u = -2 * asin(0.25 * u1);
  t = wrapToPi(theta + 0.5 * u + pi);
  v = wrapToPi(thetaf - t + u);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::T16(const double x, const double y, const double thetaf,
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
inline void Solver::T18(const double x, const double y, const double thetaf,
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
inline void Solver::T20(const double x, const double y, const double thetaf,
                        double &d) {
  const double xb = x * cos(thetaf) + y * sin(thetaf);
  const double yb = x * sin(thetaf) - y * cos(thetaf);
  const double xi = -xb - sin(thetaf);
  const double eta = -yb - 1 + cos(thetaf);
  const double u1 = sqrt(xi * xi + eta * eta);
  const double theta = atan2(eta, xi);
  const double u = -2 * asin(0.25 * u1);
  const double t = wrapToPi(theta + 0.5 * u + pi);
  const double v = wrapToPi(thetaf - t + u);

  d = fabs(t) + fabs(u) + fabs(v);

  lengths_[19][0] = -v;
  lengths_[19][1] = -u;
  lengths_[19][2] = -t;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::T19_rho(const double x, const double y, const double thetaf,
                            double &xi, double &eta, double &rho) const {
  xi = -x + sin(-thetaf);
  eta = y - 1 - cos(-thetaf);
  rho = 0.25 * (2 + sqrt(xi * xi + eta * eta));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::T19(const double x, const double y, const double thetaf,
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
inline void Solver::T19(const double xi, const double eta, const double rho,
                        const double thetaf, double &d) {
  const double u1 = acos(rho);
  const double delta = wrapToPi(2 * u1);
  const double A = sin(u1) - sin(delta);
  const double B = cos(u1) - cos(delta) - 1;
  const double t1 = atan2(eta * A - xi * B, xi * A + eta * B);
  const double t = wrapToPi(t1);
  const double v = (t - 2 * u1 + thetaf);

  d = fabs(t) + 2 * fabs(u1) + fabs(v);

  lengths_[18][0] = -t;
  lengths_[18][1] = -u1;
  lengths_[18][2] = u1;
  lengths_[18][3] = -v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void Solver::T21(const double x, const double y, const double thetaf,
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

  lengths_[20][0] = t;
  lengths_[20][1] = u;
  lengths_[20][2] = -u;
  lengths_[20][3] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::setA(const double x, const double y, const double thetaf,
                  const double r, const double beta0, int &cond, double &d) {
  // Set A
  if (thetaf >= 0) {
    if (thetaf <= atan2(fabs(cylf_ - cyl0_), cxlf_ - cxl0_)) { // omega
      if (thetaf >= atan2(cylf_ - cyr0_, cxlf_ - cxr0_)) {
        // CCSC | -R+L+S+L
        T11(x, y, thetaf, d);
        d *= r;
        cond = 11;
        return;
      }
      double beta1, u;
      Beta1(x, y, thetaf, beta1, u);
      if (beta1 <= pi / 2) {
        // CSC | +L+S+R
        const double v = wrapToPi(beta1 - thetaf);
        d = fabs(beta1) + fabs(u) + fabs(v);
        d *= r;
        cond = 1;
        lengths_[0][0] = beta1;
        lengths_[0][1] = u;
        lengths_[0][2] = v;
        return;
      }
      double beta9;
      Beta9(x, y, thetaf, beta9, u);
      if (beta9 > thetaf) {
        // CCSCC | -R+L+S+R-L
        double v = wrapToPi(beta9 - thetaf);
        d = fabs(beta9) + fabs(u) + fabs(v) + pi;
        d *= r;
        cond = 9;
        lengths_[8][0] = -beta9;
        lengths_[8][1] = pi / 2;
        lengths_[8][2] = -u;
        lengths_[8][3] = pi / 2;
        lengths_[8][4] = -v;
        return;
      }
      // CCSC | -R+L+S+R
      T10(x, y, thetaf, d);
      d *= r;
      cond = 10;
      return;
    }

    double dx = x - sin(thetaf);
    if (dx < 0) {
      // CCSC | -R+L+S+L
      T11(x, y, thetaf, d);
      d *= r;
      cond = 11;
      return;
    }

    const double dy = y - 1 + cos(thetaf);
    double rho, beta2;
    polar(dx, dy, rho, beta2);
    if (beta2 < 0 && atan2(cyrf_ - cyl0_, cxrf_ - cxl0_) < 0) {
      double beta7, u;
      Beta7(x, y, thetaf, beta7, u);
      if (beta7 + thetaf <= pi / 2) {
        // CSC | +R+S+L
        double v = wrapToPi(beta7 + thetaf);
        d = fabs(beta7) + fabs(u) + fabs(v);
        d *= r;
        cond = 7;
        lengths_[6][0] = beta7;
        lengths_[6][1] = u;
        lengths_[6][2] = v;
        return;
      }
      // CSCC | +R+S+L-R
      T8(x, y, thetaf, d);
      d *= r;
      cond = 8;
      return;
    }
    if (thetaf > pi / 2 + beta2) {
      // CSCC | +L+S+L-R
      T3(x, y, thetaf, d);
      d *= r;
      cond = 3;
      return;
    }
    // CSC | +L+S+L
    const double v_ = wrapToPi(thetaf - beta2);
    d = fabs(beta2) + fabs(rho) + fabs(v_);
    d *= r;
    lengths_[1][0] = beta2;
    lengths_[1][1] = rho;
    lengths_[1][2] = v_;
    cond = 2;
    return;
  }

  // thetaf < 0
  double beta1, u;
  Beta1(x, y, thetaf, beta1, u);
  if (thetaf < 2 * beta0 - pi) {
    if (thetaf < atan2(cylf_ - cyr0_, cxlf_ - cxr0_) - pi) {
      // CCSC | +R-L-S-L
      T6(x, y, thetaf, d);
      d *= r;
      cond = 6;
      return;
    }
    double beta5;
    Beta5(x, y, thetaf, beta5, u);
    if (thetaf <= -beta5) {
      // CCSC | +R-L-S-R
      double v = wrapToPi(-thetaf - 0.5 * pi - beta5);
      d = fabs(beta5) + fabs(u) + fabs(v) + pi / 2;
      d *= r;
      cond = 5;
      lengths_[4][0] = beta5;
      lengths_[4][1] = -pi / 2;
      lengths_[4][2] = u;
      lengths_[4][3] = v;
      return;
    }
    // CCSCC | +R-L-S-R+L
    T12(x, y, thetaf, d);
    d *= r;
    cond = 12;
    return;
  }
  if (fabs(thetaf) + beta1 >= pi / 2) {
    double t, u, beta4;
    Beta4(x, y, thetaf, t, u, beta4);
    if (fabs(beta4) < pi / 2) {
      // CSCC | +L+S+R-L
      d = fabs(t) + fabs(u) + fabs(beta4) + pi / 2;
      d *= r;
      cond = 4;
      lengths_[3][0] = -beta4;
      lengths_[3][1] = -u;
      lengths_[3][2] = pi / 2;
      lengths_[3][3] = -t;
      return;
    }
    // CCSCC | -R+L+S+R-L
    T9(x, y, thetaf, d);
    d *= r;
    cond = 9;
    return;
  }
  // CSC | +L+S+R
  const double v = wrapToPi(beta1 - thetaf);
  d = fabs(beta1) + fabs(u) + fabs(v);
  d *= r;
  cond = 1;
  lengths_[0][0] = beta1;
  lengths_[0][1] = u;
  lengths_[0][2] = v;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::setB(const double x, const double y, const double thetaf,
                  const double r, const bool pred_1, const bool pred_2,
                  int &cond, double &d) {
  // Set B
  if (thetaf >= 0) {
    double rho, beta13, d2, xi, eta;
    T13(x, y, thetaf, xi, eta, rho, beta13, d2);
    d2 = d2 * r;
    if (beta13 >= 0 && beta13 >= thetaf) {
      if (rho >= 20) {
        // CCSCC
        T9(xi, eta, rho, thetaf, d);
        d *= r;
        cond = 9;
        return;
      }
      // CCCC two possibilities
      double d1;
      T19(x, y, thetaf, d1);
      d1 *= r;
      if (d1 < d2) {
        d = d1;
        cond = 19;
        return;
      }
      d = d2;
      cond = 13;
      return;
    }

    if (thetaf < pi / 2) {
      double d2;
      T14(x, y, thetaf, d2);
      d2 *= r;
      if (pred_2) {
        // CCC
        d = d2;
        cond = 14;
        return;
      }
      const bool pred_5 = RL <= limit && RR <= 2 * r && LL <= 2 * r;
      const bool pred_6 = LR <= limit && RR <= 2 * r && LL <= 2 * r;
      if (!(pred_5 && pred_6)) {
        cond = 14;
        d = d2;
        return;
      }
      // CCC or CCCC
      double xi, eta, rho;
      T19_rho(x, y, thetaf, xi, eta, rho);
      double d1;
      T19(xi, eta, rho, thetaf, d1);
      d1 *= r;
      if (d1 < d2) {
        d = d1;
        cond = 19;
        return;
      }
      d = d2;
      cond = 14;
      return;
    }

    // thetaf > pi/2
    // CCC
    if (pred_1 == pred_2) { // v > 0 && t > 0) {
      double t, u, v;
      T15(x, y, thetaf, t, u, v);
      d = fabs(t) + fabs(u) + fabs(v);
      d *= r;
      lengths_[14][0] = t;
      lengths_[14][1] = u;
      lengths_[14][2] = v;
      cond = 15;
      return;
    }
    T20(x, y, thetaf, d);
    d *= r;
    cond = 20;
    return;
  }

  // thetaf < 0, two cases
  double xi, eta, rho, t1, t2, u1, u2, v1, v2;
  T13(x, y, thetaf, xi, eta, rho, t1, t2, u1, u2, v1, v2);
  // Case 2
  if (-t1 + t2 < thetaf) {
    if (t2 >= 0) {
      if (rho >= 20) {
        // CCSCC
        T9(xi, eta, rho, thetaf, d);
        d *= r;
        cond = 9;
        return;
      }
      // CCCC two possibilities
      double d1 = fabs(t2) + 2 * fabs(u2) + fabs(v2);
      d1 *= r;
      double d2;
      T21(x, y, thetaf, d2);
      d2 *= r;
      if (d2 >= d1) {
        d = d1;
        cond = 13;
        lengths_[12][0] = -t2;
        lengths_[12][1] = -u2;
        lengths_[12][2] = -u2;
        lengths_[12][3] = -v2;
        return;
      }
      d = d2;
      cond = 21;
      return;
    }
    // CCC or CCCC
    double d1;
    T16(x, y, thetaf, d1);
    d1 *= r;
    double d2;
    T21(x, y, thetaf, d2);
    d2 = d2 * r;
    if (d1 <= d2) {
      d = d1;
      cond = 16;
      return;
    }
    d = d2;
    cond = 21;
    return;
  }

  // Case 1
  if (t1 >= 0 && t1 >= -thetaf) {
    if (rho >= 20) {
      // CCSCC
      T12(x, y, thetaf, d);
      d *= r;
      cond = 12;
      return;
    }
    // CCCC two possibilities
    double d1 = fabs(t1) + 2 * fabs(u1) + fabs(v1);
    lengths_[16][0] = t1;
    lengths_[16][1] = u1;
    lengths_[16][2] = u1;
    lengths_[16][3] = v1;
    d1 *= r;
    double d2;
    T21(x, y, thetaf, d2);
    d2 *= r;
    if (d1 <= d2) {
      // CCCC
      d = d1;
      cond = 17;
      return;
    }
    // CCCC
    d = d2;
    cond = 21;
    return;
  }

  double d1;
  T18(x, y, thetaf, d1);
  d1 *= r;
  double d2;
  T21(x, y, thetaf, d2);
  d2 *= r;
  if (d2 >= d1) {
    // CCC
    d = d1;
    cond = 18;
    return;
  }
  // CCCC
  d = d2;
  cond = 21;
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
  cxr0_ = x0 + r * cos(theta0 - pi / 2);
  cyr0_ = y0 + r * sin(theta0 - pi / 2);
  cxl0_ = x0 - r * cos(theta0 - pi / 2);
  cyl0_ = y0 - r * sin(theta0 - pi / 2);
  cxrf_ = xn + r * cos(thetaf - pi / 2);
  cyrf_ = yn + r * sin(thetaf - pi / 2);
  cxlf_ = xn - r * cos(thetaf - pi / 2);
  cylf_ = yn - r * sin(thetaf - pi / 2);

  // Algorithm 3 (IsInSetB)
  limit = 2 * r * sqrt(2);
  LL = euclideanDistance(cxl0_, cyl0_, cxlf_, cylf_);
  LR = euclideanDistance(cxl0_, cyl0_, cxrf_, cyrf_);
  RL = euclideanDistance(cxr0_, cyr0_, cxlf_, cylf_);
  RR = euclideanDistance(cxr0_, cyr0_, cxrf_, cyrf_);
  const bool pred_1 = RR <= limit && LL <= limit && LR <= 2 * r;
  const bool pred_2 = RR <= limit && LL <= limit && RL <= 2 * r;
  const bool pred_3 = LL <= limit && LR <= 2 * r && RL <= 2 * r;

  const std::complex<double> i(0.0, 1.0);
  const double x = (xn - x0) / r;
  const double y = (yn - y0) / r;
  const double beta0 = atan2(yn - y0, xn - x0);

  // Main Algorithm (4)
  if (pred_1 || pred_2 || pred_3) {
    setB(x, y, thetaf, r, pred_1, pred_2, cond, d);
    return;
  } else {
    setA(x, y, thetaf, r, beta0, cond, d);
    return;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Solver::benchmarkPlanners() {
  const int number_of_final_states = 10000000;
  const double r = 1;

  std::vector<State> toStates;
  toStates.reserve(number_of_final_states);
  const bool use_config = false;
  helper_functions_.generateRandomStates(toStates, use_config);

  std::vector<double> distances_acceleratedRS(number_of_final_states);
  std::vector<double> distances_Desaulniers(number_of_final_states);
  std::vector<double> distances_OMPL(number_of_final_states);

  const State from = {0, 0, pi - pi / 2};

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
  ob::StateSpacePtr space(std::make_shared<ob::ReedsSheppStateSpace>(r, false));
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
  ob::StateSpacePtr space(std::make_shared<ob::ReedsSheppStateSpace>(r, false));
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
  getPath(condition, motionTypes, path);
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
    ob::StateSpacePtr space(
        std::make_shared<ob::ReedsSheppStateSpace>(r, false));
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
    validityCheck(from, toStates[i], condition, Q, errors);
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
                           const int Q, State &errors) {
  char motionTypes[5] = {'N', 'N', 'N', 'N', 'N'};
  getMotionTypes(cond, motionTypes);
  backProjectMotion(cond, Q, motionTypes);
  std::vector<State> path;
  path.push_back(from);
  getPath(cond, motionTypes, path);
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
    motionTypes[0] = 'R';
    motionTypes[1] = 'L';
    motionTypes[2] = 'R';
    break;
  case 21:
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
                     std::vector<State> &path) const {
  double x = path[0].x;
  double y = path[0].y;
  double t = path[0].theta;
  const double r = config_.r;
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

} // namespace accelerated
