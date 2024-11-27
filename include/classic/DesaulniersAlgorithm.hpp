#ifndef DESAULNIERS_HPP
#define DESAULNIERS_HPP

namespace accelerated {

class DesaulniersAlgorithm {
public:
  DesaulniersAlgorithm() {};
  double getDistance(double xf, double yf, double thetaf) const;

private:
  inline double wrapToPi(const double angle) const;
  inline double wrapTo2Pi(const double angle) const;
  inline double D1(const double cx, const double cy) const;
  inline double D2(const double cx, const double cy) const;
  inline double D3(const double cx, const double cy) const;
  inline double K(const double cx, const double cy) const;
  inline double F1(const double a, const double b, const double c,
                   const double d) const;
  inline double F2(const double a, const double b, const double c,
                   const double d) const;
  inline double F3(const double a, const double b, const double c,
                   const double d) const;
  inline double Q1(const double cx, const double cy) const;
  inline double Q2(const double cx, const double cy) const;

  inline void eq8_1(const double x, const double y, const double phi, double &t,
                    double &u, double &v) const;
  inline void eq8_2(const double x, const double y, const double phi, double &t,
                    double &u, double &v) const;
  inline void eq8_3(const double x, const double y, const double phi, double &t,
                    double &u, double &v) const;
  inline void eq8_4(const double x, const double y, const double phi, double &t,
                    double &u, double &v) const;
  inline void eq8_7(const double x, const double y, const double phi, double &t,
                    double &u, double &v) const;
  inline void eq8_8(const double x, const double y, const double phi, double &t,
                    double &u, double &v) const;
  inline void eq8_9(const double x, const double y, const double phi, double &t,
                    double &u, double &v) const;
  inline void eq8_10(const double x, const double y, const double phi,
                     double &t, double &u, double &v) const;
  inline void eq8_11(const double x, const double y, const double phi,
                     double &t, double &u, double &v) const;

  inline void polar(const double x, const double y, double &rho,
                    double &theta) const;
  inline void tauOmega(const double u, const double v, const double xi,
                       const double eta, const double phi, double &tau,
                       double &omega) const;

  inline int getRegion(const double x, const double y) const;
  inline void getType(const double cx, const double cy, const double thetaf,
                      const int R, int &T1, int &T2, int &T3) const;

  inline void getTUV(const double x, const double y, const double phi,
                     const int T, double &t, double &u, double &v,
                     double &L) const;

  inline double theta1(const double cx, const double cy) const;
  inline double theta2(const double cx, const double cy) const;
  inline double theta3(const double cx, const double cy) const;
  inline double theta4(const double cx, const double cy) const;
  inline double theta5(const double cx, const double cy) const;
  inline double theta6(const double cx, const double cy) const;
  inline double theta7(const double cx, const double cy) const;
  inline double theta8(const double cx, const double cy) const;
  inline double theta9(const double cx) const;
  inline double theta10(const double cx, const double cy) const;
  inline double theta11(const double cx, const double cy) const;
  inline double theta12(const double cx) const;
  inline double theta13(const double cx, const double cy) const;
  inline double theta14(const double cx, const double cy) const;
  inline double theta15(const double cx, const double cy) const;
  inline double theta16(const double cy) const;
  inline double theta17(const double cx, const double cy) const;
  inline double theta18(const double cx, const double cy) const;
  inline double theta19(const double cx, const double cy) const;
  inline double theta20(const double cx, const double cy) const;
  inline double theta21(const double cx, const double cy) const;
  inline double theta22(const double cx, const double cy) const;
  inline double theta23(const double cy) const;
  inline double theta24(const double cx, const double cy) const;
  inline double theta25(const double cx, const double cy) const;
  inline double theta26(const double cx, const double cy) const;
  inline double theta27(const double cx, const double cy) const;
  inline double theta28(const double cx, const double cy) const;
  inline double theta29(const double cx, const double cy) const;
  inline double theta30(const double cx, const double cy) const;
  inline double theta31(const double cx, const double cy) const;
  inline double theta33(const double cx, const double cy) const;
};

} // namespace accelerated
#endif // DESAULNIERS_HPP
