#include "classic/DesaulniersAlgorithm.hpp"
#include <cmath>
#include <numbers>

using namespace std::numbers;

namespace accelerated {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
const double DesaulniersAlgorithm::getDistance(double xf, double yf,
                                               double thetaf) const {
  // Algorithm Step 1
  thetaf = wrapTo2Pi(thetaf);
  double cx = xf + cos(thetaf);
  double cy = yf + sin(thetaf);

  // Algorithm Step 2
  const bool q1 = (cx > 1 && cy >= 0);
  const bool q2 = (cx <= 1 && cy >= 0);
  const bool q3 = (cx <= 1 && cy < 0);
  const bool q4 = (cx > 1 && cy < 0);
  if (q1) {
    xf = -xf;
    yf = yf;
    thetaf = 2 * pi - thetaf;
  } else if (q3) {
    yf = -yf;
    thetaf = 2 * pi - thetaf;
  } else if (q4) {
    xf = -xf;
    yf = -yf;
  }
  cx = xf + cos(thetaf);
  cy = yf + sin(thetaf);

  // ALgorithm step 3
  const int R = getRegion(cx, cy);
  int T1 = 0, T2 = 0, T3 = 0;
  getType(cx, cy, thetaf, R, T1, T2, T3);

  // Algorithm step 4
  double t1, u1, v1, L1;
  double t2, u2, v2, L2;
  double t3, u3, v3, L3;
  L2 = L3 = std::numeric_limits<double>::infinity();
  thetaf = wrapToPi(thetaf - pi / 2);
  const double dx = xf - 0;
  const double dy = yf - 0;
  const double theta0 = pi - pi / 2;
  const double phi = wrapToPi(thetaf - theta0);
  const double x = dx * cos(theta0) + dy * sin(theta0);
  const double y = -dx * sin(theta0) + dy * cos(theta0);
  getTUV(x, y, phi, T1, t1, u1, v1, L1);
  if (T2) {
    getTUV(x, y, phi, T2, t2, u2, v2, L2);
  }
  if (T3) {
    getTUV(x, y, phi, T3, t3, u3, v3, L3);
  }
  // If multiple paths, select the shortest one
  double L;
  if (L1 <= L2 && L1 <= L3) {
    L = L1;
  } else if (L2 <= L1 && L2 <= L3) {
    L = L2;
  } else {
    L = L3;
  }
  // Compute path of opposite type in case it's shorter
  double L_backwards;
  const double xb = x * cos(phi) + y * sin(phi);
  const double yb = x * sin(phi) - y * cos(phi);
  getTUV(-xb, -yb, phi, T1, t1, u1, v1, L_backwards);
  if (L_backwards < L) {
    L = L_backwards;
  }

  return L;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::wrapToPi(double angle) const {
  angle = fmod(angle + pi, 2 * pi);
  if (angle < 0) {
    angle += 2 * pi;
  }
  return angle - pi;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::wrapTo2Pi(double angle) const {
  angle = fmod(angle, 2 * pi);
  if (angle < 0) {
    angle += 2 * pi;
  }
  return angle;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const int DesaulniersAlgorithm::getRegion(const double x,
                                                 const double y) const {
  if (x <= -1) {
    // R1, R2, R5, R6, R7, R8, R9, R10, R11
    if (y >= 2) {
      const double d_1 = sqrt((x + 1) * (x + 1) + y * y);
      if (d_1 > sqrt(8)) {
        // R1
        return 1;
      } else {
        // R2
        return 2;
      }
    } else {
      const double d_1 = sqrt((x - 1) * (x - 1) + y * y);
      if (d_1 > sqrt(20)) {
        // R5
        return 5;
      } else {
        const double d_2 = sqrt((x + 1) * (x + 1) + (y - 2) * (y - 2));
        const double d_3 = sqrt((x + 1) * (x + 1) + y * y);
        if (d_2 > 2) {
          if (d_3 > 2) {
            // R6
            return 6;
          } else {
            // R8
            return 8;
          }
        } else {
          const double d_4 = sqrt((x + 3) * (x + 3) + y * y);
          if (d_4 >= 2) {
            if (d_3 > 2) {
              // R2
              return 2;
            } else {
              // R11
              return 11;
            }
          } else {
            if (d_3 > 2) {
              // R7
              return 7;
            } else {
              const bool c1 =
                  (64 * y) <= pow((x - 1) * (x - 1) + y * y - 4, 1.5) *
                                  sqrt(36 - ((x - 1) * (x - 1) + y * y));
              if (c1) {
                // R9
                return 9;
              } else {
                // R10
                return 10;
              }
            }
          }
        }
      }
    }
  } else {
    // R3, R4, R12, R13, R14, R15, R16, R17, R18, R19
    const double d_1 = sqrt((x + 1) * (x + 1) + y * y);
    if (d_1 > sqrt(8)) {
      // R3
      return 3;
    } else {
      const double d_2 = sqrt((x - 1) * (x - 1) + y * y);
      if (d_1 > 2) {
        if (d_2 > 2) {
          // R4
          return 4;
        } else {
          if (d_2 < 2 * (sqrt(2) - 1)) {
            // R15
            return 15;
          } else {
            const double d_3 = sqrt((x + 1) * (x + 1) + (y - 2) * (y - 2));
            if (d_3 >= 2) {
              // R14
              return 14;
            } else {
              // R13
              return 13;
            }
          }
        }
      } else {
        if (d_2 > 2) {
          // R12
          return 12;
        } else {
          const bool c1 = y >= sqrt(5 - 2 * x - x * x - 2 * sqrt(5 - 4 * x));
          if (c1) {
            const double d_3 = sqrt((x + 1) * (x + 1) + (y - 2) * (y - 2));
            if (d_3 < 2) {
              // R13
              return 13;
            } else {
              if (d_2 > 2 * (sqrt(2) - 1)) {
                // R14
                return 14;
              } else {
                // R15
                return 15;
              }
            }
          } else {
            const double d_3 = sqrt((x + 1) * (x + 1) + (y - 2) * (y - 2));
            const bool c2 =
                (8 * pow(sqrt((x + 1) * (x + 1) + y * y), 3) -
                 pow(sqrt((x + 1) * (x + 1) + y * y), 4)) *
                        pow(4 - pow(sqrt((x + 1) * (x + 1) + y * y), 2), 2) -
                    256 * y * y * pow(1 + sqrt((x + 1) * (x + 1) + y * y), 2) <=
                0;
            if (d_3 < 2) {
              if (c2) {
                // R17
                return 17;
              } else {
                // R18
                return 18;
              }
            } else {
              if (c2) {
                // R16
                return 16;
              } else {
                // R19
                return 19;
              }
            }
          }
        }
      }
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::getType(const double cx, const double cy,
                                          const double thetaf, const int R,
                                          int &T1, int &T2, int &T3) const {
  switch (R) {
  case 1:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf <= theta2(cx, cy)) {
      T1 = 8;
    } else if (thetaf > theta2(cx, cy) && thetaf < theta3(cx, cy)) {
      T1 = 18;
    } else if (thetaf >= theta3(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 2:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf <= theta7(cx, cy)) {
      T1 = 10;
    } else if (thetaf > theta7(cx, cy) && thetaf <= theta8(cx, cy)) {
      T1 = 13;
    } else if (thetaf > theta8(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 3:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta2(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta2(cx, cy) && thetaf < theta3(cx, cy)) {
      T1 = 18;
    } else if (thetaf > theta3(cx, cy) && thetaf <= theta9(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta9(cx, cy) && thetaf <= theta10(cx, cy)) {
      T1 = 21;
    } else if (thetaf > theta10(cx, cy) && thetaf <= theta11(cx, cy)) {
      T1 = 22;
    } else if (thetaf > theta11(cx, cy) && thetaf <= theta12(cx, cy)) {
      T1 = 24;
    } else if (thetaf > theta12(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 4:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf <= theta7(cx, cy)) {
      T1 = 10;
    } else if (thetaf > theta7(cx, cy) && thetaf <= theta8(cx, cy)) {
      T1 = 13;
    } else if (thetaf > theta8(cx, cy) && thetaf <= theta9(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta9(cx, cy) && thetaf <= theta10(cx, cy)) {
      T1 = 21;
    } else if (thetaf > theta10(cx, cy) && thetaf <= theta11(cx, cy)) {
      T1 = 22;
    } else if (thetaf > theta11(cx, cy) && thetaf <= theta12(cx, cy)) {
      T1 = 24;
    } else if (thetaf > theta12(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 5:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf <= theta13(cx, cy)) {
      T1 = 8;
    } else if (thetaf > theta13(cx, cy) && thetaf <= theta14(cx, cy)) {
      T1 = 9;
    } else if (thetaf > theta14(cx, cy) && thetaf < theta15(cx, cy)) {
      T1 = 14;
    } else if (thetaf >= theta15(cx, cy) && thetaf < theta16(cx, cy)) {
      T1 = 17;
    } else if (thetaf >= theta16(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 6:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf <= theta17(cx, cy)) {
      T1 = 10;
    } else if (thetaf > theta17(cx, cy) && thetaf < theta18(cx, cy)) {
      T1 = 11;
    } else if (thetaf >= theta18(cx, cy) && thetaf <= theta19(cx, cy)) {
      T1 = 15;
    } else if (thetaf > theta19(cx, cy) && thetaf < theta20(cx, cy)) {
      T1 = 16;
    } else if (thetaf >= theta20(cx, cy) && thetaf < theta16(cx, cy)) {
      T1 = 17;
    } else if (thetaf >= theta16(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 7:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf <= theta17(cx, cy)) {
      T1 = 10;
    } else if (thetaf > theta17(cx, cy) && thetaf < theta18(cx, cy)) {
      T1 = 11;
    } else if (thetaf >= theta18(cx, cy) && thetaf <= theta19(cx, cy)) {
      T1 = 15;
    } else if (thetaf > theta19(cx, cy) && thetaf <= theta21(cx, cy)) {
      T1 = 16;
    } else if (thetaf > theta21(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 8:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf <= theta17(cx, cy)) {
      T1 = 10;
    } else if (thetaf > theta17(cx, cy) && thetaf < theta18(cx, cy)) {
      T1 = 11;
      T2 = 12;
    } else if (thetaf >= theta18(cx, cy) && thetaf <= theta23(cx, cy)) {
      T1 = 12;
      T2 = 15;
    } else if (thetaf > theta23(cx, cy) && thetaf <= theta19(cx, cy)) {
      T1 = 15;
    } else if (thetaf > theta19(cx, cy) && thetaf < theta20(cx, cy)) {
      T1 = 16;
    } else if (thetaf >= theta20(cx, cy) && thetaf < theta16(cx, cy)) {
      T1 = 17;
    } else if (thetaf >= theta16(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 9:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf <= theta17(cx, cy)) {
      T1 = 10;
    } else if (thetaf > theta17(cx, cy) && thetaf < theta18(cx, cy)) {
      T1 = 11;
      T2 = 12;
    } else if (thetaf >= theta18(cx, cy) && thetaf <= theta23(cx, cy)) {
      T1 = 12;
      T2 = 15;
    } else if (thetaf > theta23(cx, cy) && thetaf <= theta19(cx, cy)) {
      T1 = 15;
    } else if (thetaf > theta19(cx, cy) && thetaf <= theta21(cx, cy)) {
      T1 = 16;
    } else if (thetaf > theta21(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 10:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf < theta24(cx, cy)) {
      T1 = 10;
      T2 = 12;
    } else if (thetaf >= theta24(cx, cy) && thetaf <= theta23(cx, cy)) {
      T1 = 12;
      T2 = 15;
    } else if (thetaf > theta23(cx, cy) && thetaf <= theta19(cx, cy)) {
      T1 = 15;
    } else if (thetaf > theta19(cx, cy) && thetaf <= theta21(cx, cy)) {
      T1 = 16;
    } else if (thetaf > theta21(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 11:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf < theta25(cx, cy)) {
      T1 = 10;
      T2 = 12;
    } else if (thetaf >= theta25(cx, cy) && thetaf <= theta21(cx, cy)) {
      T1 = 13;
    } else if (thetaf > theta21(cx, cy) && thetaf <= theta4(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta4(cx, cy) && thetaf <= theta5(cx, cy)) {
      T1 = 20;
    } else if (thetaf > theta5(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 12:
    if (thetaf >= 0 && thetaf <= theta1(cx, cy)) {
      T1 = 7;
    } else if (thetaf > theta1(cx, cy) && thetaf < theta6(cx, cy)) {
      T1 = 8;
    } else if (thetaf >= theta6(cx, cy) && thetaf < theta25(cx, cy)) {
      T1 = 10;
      T2 = 12;
    } else if (thetaf >= theta25(cx, cy) && thetaf <= theta21(cx, cy)) {
      T1 = 13;
    } else if (thetaf > theta21(cx, cy) && thetaf <= theta9(cx, cy)) {
      T1 = 19;
    } else if (thetaf > theta9(cx, cy) && thetaf <= theta10(cx, cy)) {
      T1 = 21;
    } else if (thetaf > theta10(cx, cy) && thetaf <= theta11(cx, cy)) {
      T1 = 22;
    } else if (thetaf > theta11(cx, cy) && thetaf <= theta12(cx, cy)) {
      T1 = 24;
    } else if (thetaf > theta12(cx, cy) && thetaf < 2 * pi) {
      T1 = 25;
    }
    break;
  case 13:
    if (thetaf >= 0 && thetaf <= theta26(cx, cy)) {
      T1 = 1;
    } else if (thetaf > theta26(cx, cy) && thetaf <= theta27(cx, cy)) {
      T1 = 4;
    } else if (thetaf > theta27(cx, cy) && thetaf < theta28(cx, cy)) {
      T1 = 6;
    } else if (thetaf >= theta28(cx, cy) && thetaf <= theta20(cx, cy)) {
      T1 = 23;
    } else if (thetaf > theta20(cx, cy) && thetaf <= theta12(cx, cy)) {
      T1 = 24;
    } else if (thetaf > theta12(cx, cy) && thetaf < theta29(cx, cy)) {
      T1 = 25;
    } else if (thetaf > theta29(cx, cy) && thetaf < 2 * pi) {
      T1 = 26;
    }
    break;
  case 14:
    if (thetaf >= 0 && thetaf <= theta26(cx, cy)) {
      T1 = 1;
    } else if (thetaf > theta26(cx, cy) && thetaf <= theta27(cx, cy)) {
      T1 = 4;
    } else if (thetaf > theta27(cx, cy) && thetaf < theta28(cx, cy)) {
      T1 = 6;
    } else if (thetaf >= theta28(cx, cy) && thetaf <= theta20(cx, cy)) {
      T1 = 23;
    } else if (thetaf > theta20(cx, cy) && thetaf < theta30(cx, cy)) {
      T1 = 24;
    } else if (thetaf >= theta30(cx, cy) && thetaf < theta29(cx, cy)) {
      T1 = 23;
    } else if (thetaf > theta29(cx, cy) && thetaf < 2 * pi) {
      T1 = 26;
    }
    break;
  case 15:
    if (thetaf >= 0 && thetaf <= theta26(cx, cy)) {
      T1 = 1;
    } else if (thetaf > theta26(cx, cy) && thetaf <= theta27(cx, cy)) {
      T1 = 4;
    } else if (thetaf > theta27(cx, cy) && thetaf < theta28(cx, cy)) {
      T1 = 6;
    } else if (thetaf >= theta28(cx, cy) && thetaf <= theta29(cx, cy)) {
      T1 = 23;
    } else if (thetaf > theta29(cx, cy) && thetaf < 2 * pi) {
      T1 = 26;
    }
    break;
  case 16:
    if (thetaf >= 0 && thetaf <= theta26(cx, cy)) {
      T1 = 1;
    } else if (thetaf > theta26(cx, cy) && thetaf <= theta31(cx, cy)) {
      T1 = 2;
    } else if (thetaf > theta31(cx, cy) && thetaf <= theta33(cx, cy)) {
      T1 = 3;
      T2 = 6;
    } else if (thetaf > theta33(cx, cy) && thetaf < theta28(cx, cy)) {
      T1 = 6;
    } else if (thetaf >= theta28(cx, cy) && thetaf < theta29(cx, cy)) {
      T1 = 23;
    } else if (thetaf > theta29(cx, cy) && thetaf < 2 * pi) {
      T1 = 26;
    }
    break;
  case 17:
    if (thetaf >= 0 && thetaf <= theta26(cx, cy)) {
      T1 = 1;
    } else if (thetaf > theta26(cx, cy) && thetaf <= theta31(cx, cy)) {
      T1 = 2;
    } else if (thetaf > theta31(cx, cy) && thetaf <= theta33(cx, cy)) {
      T1 = 3;
      T2 = 6;
    } else if (thetaf > theta33(cx, cy) && thetaf < theta28(cx, cy)) {
      T1 = 6;
    } else if (thetaf >= theta28(cx, cy) && thetaf <= theta20(cx, cy)) {
      T1 = 23;
    } else if (thetaf > theta20(cx, cy) && thetaf <= theta12(cx, cy)) {
      T1 = 24;
    } else if (thetaf > theta12(cx, cy) && thetaf < theta29(cx, cy)) {
      T1 = 25;
    } else if (thetaf > theta29(cx, cy) && thetaf < 2 * pi) {
      T1 = 26;
    }
    break;
  case 18:
    if (thetaf >= 0 && thetaf <= theta26(cx, cy)) {
      T1 = 1;
    } else if (thetaf > theta26(cx, cy) && thetaf <= theta31(cx, cy)) {
      T1 = 2;
    } else if (thetaf > theta31(cx, cy) && thetaf <= theta33(cx, cy)) {
      T1 = 3;
      T2 = 5;
      T3 = 6;
    } else if (thetaf > theta33(cx, cy) && thetaf < theta28(cx, cy)) {
      T1 = 6;
    } else if (thetaf >= theta28(cx, cy) && thetaf <= theta20(cx, cy)) {
      T1 = 23;
    } else if (thetaf > theta20(cx, cy) && thetaf <= theta12(cx, cy)) {
      T1 = 24;
    } else if (thetaf > theta12(cx, cy) && thetaf < theta29(cx, cy)) {
      T1 = 25;
    } else if (thetaf > theta29(cx, cy) && thetaf < 2 * pi) {
      T1 = 26;
    }
    break;
  case 19:
    if (thetaf >= 0 && thetaf <= theta26(cx, cy)) {
      T1 = 1;
    } else if (thetaf > theta26(cx, cy) && thetaf <= theta31(cx, cy)) {
      T1 = 2;
    } else if (thetaf > theta31(cx, cy) && thetaf <= theta33(cx, cy)) {
      T1 = 3;
      T2 = 5;
      T3 = 6;
    } else if (thetaf > theta33(cx, cy) && thetaf < theta28(cx, cy)) {
      T1 = 6;
    } else if (thetaf >= theta28(cx, cy) && thetaf < theta29(cx, cy)) {
      T1 = 23;
    } else if (thetaf > theta29(cx, cy) && thetaf < 2 * pi) {
      T1 = 26;
    }
    break;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::D1(const double cx,
                                             const double cy) const {
  return sqrt((cx - 1) * (cx - 1) + cy * cy);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::D2(const double cx,
                                             const double cy) const {
  return sqrt((cx + 1) * (cx + 1) + cy * cy);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::D3(const double cx,
                                             const double cy) const {
  return sqrt((cx - 3) * (cx - 3) + cy * cy);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::K(const double cx,
                                            const double cy) const {
  return D1(cx, cy) * D1(cx, cy) + 4 -
         16 * (sin(F3(D2(cx, cy), cx + 1, cy, 4) / 2) *
               sin(F3(D2(cx, cy), cx + 1, cy, 4) / 2));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::F1(const double a, const double b,
                                             const double c,
                                             const double d) const {
  return acos((a * b + c * sqrt(d * d - a * a)) / (d * d));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::F2(const double a, const double b,
                                             const double c,
                                             const double d) const {
  return acos((a * b + c * sqrt(4 * d * d - a * a)) / (d * d));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::F3(const double a, const double b,
                                             const double c,
                                             const double d) const {
  return acos((a * b + c * sqrt(d * d - a * a)) / (a * d));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::Q1(const double cx,
                                             const double cy) const {
  return (cx - 1) / 2 +
         (cy * sqrt(16 - D2(cx, cy) * D2(cx, cy))) / (2 * D2(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::Q2(const double cx,
                                             const double cy) const {
  return (cx - 1) / 2 -
         (cy * sqrt(16 - D2(cx, cy) * D2(cx, cy))) / (2 * D2(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::polar(const double x, const double y,
                                        double &rho, double &theta) const {
  rho = sqrt(x * x + y * y);
  theta = atan2(y, x);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::tauOmega(const double u, const double v,
                                           const double xi, const double eta,
                                           const double phi, double &tau,
                                           double &omega) const {
  const double delta = wrapToPi(u - v);
  const double A = sin(u) - sin(delta);
  const double B = cos(u) - cos(delta) - 1;
  const double t1 = atan2(eta * A - xi * B, xi * A + eta * B);
  const double t2 = 2 * (cos(delta) - cos(v) - cos(u)) + 3;
  if (t2 < 0) {
    tau = wrapToPi(t1 + pi);
  } else {
    tau = wrapToPi(t1);
  }
  omega = wrapToPi(tau - u + v - phi);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_1(const double x, const double y,
                                        const double phi, double &t, double &u,
                                        double &v) const {
  polar(x - sin(phi), y - 1 + cos(phi), u, t);
  v = wrapToPi(phi - t);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_2(const double x, const double y,
                                        const double phi, double &t, double &u,
                                        double &v) const {
  double u1, t1;
  polar(x + sin(phi), y - 1 - cos(phi), u1, t1);
  u1 = u1 * u1;
  u = sqrt(u1 - 4);
  const double theta = atan2(2, u);
  t = wrapToPi(t1 + theta);
  v = wrapToPi(t - phi);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_3(const double x, const double y,
                                        const double phi, double &t, double &u,
                                        double &v) const {
  const double xi = x - sin(phi);
  const double eta = y - 1 + cos(phi);
  double rho, theta;
  polar(xi, eta, rho, theta);
  u = -2 * asin(0.25 * rho);
  t = wrapToPi(theta + 0.5 * u + pi);
  v = wrapToPi(phi - t + u);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_4(const double x, const double y,
                                        const double phi, double &t, double &u,
                                        double &v) const {
  const double xi = x - sin(phi);
  const double eta = y - 1 + cos(phi);
  double rho, theta;
  polar(xi, eta, rho, theta);
  u = -2 * asin(0.25 * rho);
  t = wrapToPi(theta + 0.5 * u + pi);
  v = wrapToPi(phi - t + u);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_7(const double x, const double y,
                                        const double phi, double &t, double &u,
                                        double &v) const {
  const double xi = x + sin(phi);
  const double eta = y - 1 - cos(phi);
  const double rho = 0.25 * (2 + sqrt(xi * xi + eta * eta));
  u = acos(rho);
  tauOmega(u, -u, xi, eta, phi, t, v);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_8(const double x, const double y,
                                        const double phi, double &t, double &u,
                                        double &v) const {
  const double xi = x + sin(phi);
  const double eta = y - 1 - cos(phi);
  const double rho = (20 - xi * xi - eta * eta) / 16;
  u = -acos(rho);
  tauOmega(u, u, xi, eta, phi, t, v);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_9(const double x, const double y,
                                        const double phi, double &t, double &u,
                                        double &v) const {
  const double xi = x - sin(phi);
  const double eta = y - 1 + cos(phi);
  double rho, theta;
  polar(xi, eta, rho, theta);
  const double r = sqrt(rho * rho - 4);
  u = 2 - r;
  t = wrapToPi(theta + atan2(r, -2));
  v = wrapToPi(phi - 0.5 * pi - t);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_10(const double x, const double y,
                                         const double phi, double &t, double &u,
                                         double &v) const {
  const double xi = x + sin(phi);
  const double eta = y - 1 - cos(phi);
  double rho, theta;
  polar(-eta, xi, rho, theta);
  t = theta;
  u = 2 - rho;
  v = wrapToPi(t + 0.5 * pi - phi);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::eq8_11(const double x, const double y,
                                         const double phi, double &t, double &u,
                                         double &v) const {
  const double xi = x + sin(phi);
  const double eta = y - 1 - cos(phi);
  double rho, theta;
  polar(xi, eta, rho, theta);
  u = 4 - sqrt(rho * rho - 4);
  t = wrapToPi(atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta));
  v = wrapToPi(t - phi);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta1(const double cx,
                                                 const double cy) const {
  return pi / 2 - asin(cy / D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta2(const double cx,
                                                 const double cy) const {
  return theta3(cx, cy) -
         F1(2 * cos(F1(2, cx + 1, cy, D2(cx, cy))) - 1, cy, 1 - cx, D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta3(const double cx,
                                                 const double cy) const {
  return pi / 2 + F1(2, cx + 1, cy, D2(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta4(const double cx,
                                                 const double cy) const {
  return 3 * pi / 2 - asin(cy / D2(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta5(const double cx,
                                                 const double cy) const {
  return theta4(cx, cy) + pi / 2;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta6(const double cx,
                                                 const double cy) const {
  return F2(D1(cx, cy) * D1(cx, cy) - 4, cx - 1, cy, 2 * D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta7(const double cx,
                                                 const double cy) const {
  return pi + F3(D2(cx, cy), cx + 1, cy, 4) -
         F2(K(cx, cy), 1 - cx, -cy, 2 * D1(cx, cy)) - 2 * asin(D2(cx, cy) / 4);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta8(const double cx,
                                                 const double cy) const {
  return theta7(cx, cy) + F2(K(cx, cy), 1 - cx, -cy, 2 * D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta9(const double cx,
                                                 const double cy) const {
  return pi - acos((1 - cx) / 2);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta10(const double cx,
                                                  const double cy) const {
  return pi - F1(2, 1 - cx, cy, D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta11(const double cx,
                                                  const double cy) const {
  return theta10(cx, cy) + pi / 2;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta12(const double cx,
                                                  const double cy) const {
  return 3 * pi / 2 + asin((cx + 1) / 2);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta13(const double cx,
                                                  const double cy) const {
  return pi - F1(2, -cy, 1 - cx, D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta14(const double cx,
                                                  const double cy) const {
  return pi - 2 * asin(cy / D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta15(const double cx,
                                                  const double cy) const {
  return theta14(cx, cy) + F1(2, -cy, 1 - cx, D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta16(const double cx,
                                                  const double cy) const {
  double input = (2 - cy) / 2;
  if (input > 1.0) {
    input = 1.0; // Limit input to maximum value if it exceeds 1
  } else if (input < -1.0) {
    input = -1.0; // Limit input to minimum value if it's less than -1
  }
  return pi + asin(input);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta17(const double cx,
                                                  const double cy) const {
  return pi - acos((12 + D1(cx, cy) * D1(cx, cy)) / (8 * D1(cx, cy))) -
         asin(cy / D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta18(const double cx,
                                                  const double cy) const {
  return pi - 2 * asin(cy / D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta19(const double cx,
                                                  const double cy) const {
  return theta17(cx, cy) +
         2 * acos((12 + D1(cx, cy) * D1(cx, cy)) / (8 * D1(cx, cy)));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta20(const double cx,
                                                  const double cy) const {
  return pi + F2(D1(cx, cy) * D1(cx, cy) - 4, 1 - cx, cy, 2 * D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta21(const double cx,
                                                  const double cy) const {
  return pi + F2(D2(cx, cy) * D2(cx, cy), cx + 1, cy, 2 * D2(cx, cy)) -
         2 * asin(D2(cx, cy) / 4);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta23(const double cx,
                                                  const double cy) const {
  return pi - asin(cy / 2);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta24(const double cx,
                                                  const double cy) const {
  return theta19(cx, cy) - acos((20 - D1(cx, cy) * D1(cx, cy)) / 16);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta25(const double cx,
                                                  const double cy) const {
  return theta21(cx, cy) - 2 * asin(D2(cx, cy) / 4);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta26(const double cx,
                                                  const double cy) const {
  return acos((cx - Q1(cx, cy)) / 2);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta27(const double cx,
                                                  const double cy) const {
  return theta26(cx, cy) + 2 * asin(D2(cx, cy) / 4);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta28(const double cx,
                                                  const double cy) const {
  return pi + asin((2 * sin(F2(D3(cx, cy) * D3(cx, cy), 3 - cx, -cy,
                               2 * D3(cx, cy))) -
                    cy) /
                   2);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta29(const double cx,
                                                  const double cy) const {
  return pi + acos((Q2(cx, cy) - cx) / 2);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta30(const double cx,
                                                  const double cy) const {
  return pi + F2(D1(cx, cy) * D1(cx, cy) - 4, 1 - cx, -cy, 2 * D1(cx, cy));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta31(const double cx,
                                                  const double cy) const {
  return pi - asin(cy / D1(cx, cy)) - acos((2 + D1(cx, cy)) / 4);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline const double DesaulniersAlgorithm::theta33(const double cx,
                                                  const double cy) const {
  return theta31(cx, cy) + acos((2 + D1(cx, cy)) / 4);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
inline void DesaulniersAlgorithm::getTUV(const double x, const double y,
                                         const double phi, const int T,
                                         double &t, double &u, double &v,
                                         double &L) const {

  double xb, yb;
  switch (T) {
  case 1:
    eq8_3(-x, y, -phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 2:
    xb = x * cos(phi) + y * sin(phi);
    yb = x * sin(phi) - y * cos(phi);
    eq8_4(xb, -yb, -phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 3:
    eq8_7(-x, -y, phi, t, u, v);
    L = fabs(t) + 2 * fabs(u) + fabs(v);
    break;
  case 4:
    eq8_4(-x, y, -phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 5:
    eq8_8(x, y, phi, t, u, v);
    L = fabs(t) + 2 * fabs(u) + fabs(v);
    break;
  case 6:
    eq8_8(-x, y, -phi, t, u, v);
    L = fabs(t) + 2 * fabs(u) + fabs(v);
    break;
  case 7:
    eq8_10(x, -y, -phi, t, u, v);
    L = fabs(t) + pi / 2 + fabs(u) + fabs(v);
    break;
  case 8:
    eq8_9(x, -y, -phi, t, u, v);
    L = fabs(t) + pi / 2 + fabs(u) + fabs(v);
    break;
  case 9:
    eq8_11(x, -y, -phi, t, u, v);
    L = fabs(t) + pi / 2 + fabs(u) + pi / 2 + fabs(v);
    break;
  case 10:
    eq8_4(x, -y, -phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 11:
    eq8_8(x, -y, -phi, t, u, v);
    L = fabs(t) + 2 * fabs(u) + fabs(v);
    break;
  case 12:
    eq8_7(x, y, phi, t, u, v);
    L = fabs(t) + 2 * fabs(u) + fabs(v);
    break;
  case 13:
    xb = x * cos(phi) + y * sin(phi);
    yb = x * sin(phi) - y * cos(phi);
    eq8_4(-xb, yb, -phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 14:
    eq8_11(-x, -y, phi, t, u, v);
    L = fabs(t) + pi / 2 + fabs(u) + pi / 2 + fabs(v);
    break;
  case 15:
    eq8_8(-x, -y, phi, t, u, v);
    L = fabs(t) + 2 * fabs(u) + fabs(v);
    break;
  case 16:
    eq8_4(-x, -y, phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 17:
    eq8_9(-x, -y, phi, t, u, v);
    L = fabs(t) + pi / 2 + fabs(u) + fabs(v);
    break;
  case 18:
    xb = x * cos(phi) + y * sin(phi);
    yb = x * sin(phi) - y * cos(phi);
    eq8_9(-xb, yb, -phi, t, u, v);
    L = fabs(t) + fabs(u) + pi / 2 + fabs(v);
    break;
  case 19:
    eq8_2(x, y, phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 20:
    eq8_1(x, y, phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 21:
    eq8_1(x, -y, -phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 22:
    eq8_2(x, -y, -phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 23:
    xb = x * cos(phi) + y * sin(phi);
    yb = x * sin(phi) - y * cos(phi);
    eq8_4(-xb, -yb, phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  case 24:
    xb = x * cos(phi) + y * sin(phi);
    yb = x * sin(phi) - y * cos(phi);
    eq8_9(-xb, -yb, phi, t, u, v);
    L = fabs(t) + fabs(u) + pi / 2 + fabs(v);
    break;
  case 25:
    xb = x * cos(phi) + y * sin(phi);
    yb = x * sin(phi) - y * cos(phi);
    eq8_10(-xb, -yb, phi, t, u, v);
    L = fabs(t) + fabs(u) + pi / 2 + fabs(v);
    break;
  case 26:
    eq8_3(x, y, phi, t, u, v);
    L = fabs(t) + fabs(u) + fabs(v);
    break;
  }
}

} // namespace accelerated
