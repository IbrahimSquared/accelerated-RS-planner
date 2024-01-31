#ifndef DESAULNIERS_HPP
#define DESAULNIERS_HPP

namespace accelerated {

class DesaulniersAlgorithm {
public:
  DesaulniersAlgorithm(){};
  const double getDistance(double xf, double yf, double thetaf) const;

private:
  inline const double wrapToPi(const double angle) const;
  inline const double wrapTo2Pi(const double angle) const;
  inline const double D1(const double cx, const double cy) const;
  inline const double D2(const double cx, const double cy) const;
  inline const double D3(const double cx, const double cy) const;
  inline const double K(const double cx, const double cy) const;
  inline const double F1(const double a, const double b, const double c,
                         const double d) const;
  inline const double F2(const double a, const double b, const double c,
                         const double d) const;
  inline const double F3(const double a, const double b, const double c,
                         const double d) const;
  inline const double Q1(const double cx, const double cy) const;
  inline const double Q2(const double cx, const double cy) const;

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

  inline const int getRegion(const double x, const double y) const;
  inline void getType(const double cx, const double cy, const double thetaf,
                      const int R, int &T1, int &T2, int &T3) const;

  inline void getTUV(const double x, const double y, const double phi,
                     const int T, double &t, double &u, double &v,
                     double &L) const;

  inline const double theta1(const double cx, const double cy) const;
  inline const double theta2(const double cx, const double cy) const;
  inline const double theta3(const double cx, const double cy) const;
  inline const double theta4(const double cx, const double cy) const;
  inline const double theta5(const double cx, const double cy) const;
  inline const double theta6(const double cx, const double cy) const;
  inline const double theta7(const double cx, const double cy) const;
  inline const double theta8(const double cx, const double cy) const;
  inline const double theta9(const double cx, const double cy) const;
  inline const double theta10(const double cx, const double cy) const;
  inline const double theta11(const double cx, const double cy) const;
  inline const double theta12(const double cx, const double cy) const;
  inline const double theta13(const double cx, const double cy) const;
  inline const double theta14(const double cx, const double cy) const;
  inline const double theta15(const double cx, const double cy) const;
  inline const double theta16(const double cx, const double cy) const;
  inline const double theta17(const double cx, const double cy) const;
  inline const double theta18(const double cx, const double cy) const;
  inline const double theta19(const double cx, const double cy) const;
  inline const double theta20(const double cx, const double cy) const;
  inline const double theta21(const double cx, const double cy) const;
  inline const double theta22(const double cx, const double cy) const;
  inline const double theta23(const double cx, const double cy) const;
  inline const double theta24(const double cx, const double cy) const;
  inline const double theta25(const double cx, const double cy) const;
  inline const double theta26(const double cx, const double cy) const;
  inline const double theta27(const double cx, const double cy) const;
  inline const double theta28(const double cx, const double cy) const;
  inline const double theta29(const double cx, const double cy) const;
  inline const double theta30(const double cx, const double cy) const;
  inline const double theta31(const double cx, const double cy) const;
  inline const double theta33(const double cx, const double cy) const;
};

} // namespace accelerated
#endif // DESAULNIERS_HPP
