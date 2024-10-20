#ifndef EXAMPLES_2D_ADVECTOR_IRL_2D_H_
#define EXAMPLES_2D_ADVECTOR_IRL_2D_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include "irl/generic_cutting/paraboloid_intersection/moment_contributions.h"
#include "irl/geometry/general/math_vector.h"
#include "irl/helpers/expression_templates.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"

namespace IRL2D {

class Vec {
 public:
  constexpr Vec(void) : vec_m{0., 0.} {};
  constexpr Vec(const double x, const double y) : vec_m{x, y} {};
  Vec& operator=(const Vec& a) {
    vec_m[0] = a.x();
    vec_m[1] = a.y();
    return (*this);
  };
  double& operator[](const IRL::UnsignedIndex_t d) { return vec_m[d]; };
  const double& operator[](const IRL::UnsignedIndex_t d) const {
    return vec_m[d];
  };
  double& x(void) { return vec_m[0]; };
  double& y(void) { return vec_m[1]; };
  const double& x(void) const { return vec_m[0]; };
  const double& y(void) const { return vec_m[1]; };
  void normalize(void) {
    const double inv_magnitude =
        1. / IRL::safelyEpsilon(
                 std::sqrt(vec_m[0] * vec_m[0] + vec_m[1] * vec_m[1]));
    vec_m[0] *= inv_magnitude;
    vec_m[1] *= inv_magnitude;
  };
  const double magnitude(void) const {
    return std::sqrt(vec_m[0] * vec_m[0] + vec_m[1] * vec_m[1]);
  };
  Vec operator-(void) const { return Vec(-vec_m[0], -vec_m[1]); };
  Vec& operator+=(const Vec& a) {
    vec_m[0] += a.x();
    vec_m[1] += a.y();
    return (*this);
  };
  Vec& operator-=(const Vec& a) {
    vec_m[0] -= a.x();
    vec_m[1] -= a.y();
    return (*this);
  };
  Vec operator/(const double a) { return Vec(vec_m[0] / a, vec_m[1] / a); };
  Vec& operator*=(const double a) {
    vec_m[0] *= a;
    vec_m[1] *= a;
    return (*this);
  };
  Vec& operator/=(const double a) {
    vec_m[0] /= a;
    vec_m[1] /= a;
    return (*this);
  };
  friend std::ostream& operator<<(std::ostream& out, const Vec& vec) {
    out << std::setprecision(3) << std::scientific << "(" << std::setw(10)
        << std::setfill(' ') << vec.x() << "," << std::setw(10)
        << std::setfill(' ') << vec.y() << ")";
    return out;
  }

 private:
  std::array<double, 2> vec_m;
};

const double operator*(const Vec& a_vec_0, const Vec& a_vec_1);
const Vec operator+(const Vec& a_vec_0, const Vec& a_vec_1);
const Vec operator-(const Vec& a_vec_0, const Vec& a_vec_1);
const Vec operator*(const double a_scalar, const Vec& a_vec);
const Vec operator*(const Vec& a_vec, const double a_scalar);
const Vec operator/(const Vec& a_vec, const double a_scalar);

class Mat {
 public:
  constexpr Mat(void) : matrix_m{Vec(1., 0.), Vec(0., 1.)} {};
  constexpr Mat(const double angular_pos)
      : matrix_m{Vec(std::cos(angular_pos), std::sin(angular_pos)),
                 Vec(-std::sin(angular_pos), std::cos(angular_pos))} {};
  constexpr Mat(const Vec& v0, const Vec& v1) : matrix_m{v0, v1} {};
  Vec& operator[](const IRL::UnsignedIndex_t d) { return matrix_m[d]; };
  const Vec& operator[](const IRL::UnsignedIndex_t d) const {
    return matrix_m[d];
  };
  Mat& operator=(const Mat& a) {
    matrix_m[0] = a[0];
    matrix_m[1] = a[1];
    return (*this);
  };
  const Mat transpose(void) const {
    Mat mT(*this);
    mT[0][1] = (*this)[1][0];
    mT[1][0] = (*this)[0][1];
    return mT;
  }
  Mat operator-(void) const { return Mat(-matrix_m[0], -matrix_m[1]); };
  Mat& operator+=(const Mat& a) {
    matrix_m[0] += a[0];
    matrix_m[1] += a[1];
    return (*this);
  };
  Mat& operator-=(const Mat& a) {
    matrix_m[0] -= a[0];
    matrix_m[1] -= a[1];
    return (*this);
  };
  Mat operator/(const double a) {
    return Mat(matrix_m[0] / a, matrix_m[1] / a);
  };
  Mat& operator*=(const double a) {
    matrix_m[0] *= a;
    matrix_m[1] *= a;
    return (*this);
  };
  Mat& operator/=(const double a) {
    matrix_m[0] /= a;
    matrix_m[1] /= a;
    return (*this);
  };

  friend std::ostream& operator<<(std::ostream& out, const Mat& mat) {
    out << std::setprecision(3) << std::scientific << "(" << mat[0] << ","
        << mat[1] << ")";
    return out;
  }

 private:
  std::array<Vec, 2> matrix_m;
};

const Mat operator*(const Mat& a_mat_0, const Mat& a_mat_1);
const Mat operator*(const double a_scalar, const Mat& a_mat);
const Mat operator*(const Mat& a_mat, const double a_scalar);
const Mat operator+(const Mat& a_mat_0, const Mat& a_mat_1);
const Mat operator-(const Mat& a_mat_0, const Mat& a_mat_1);
const Vec operator*(const Mat& a_mat, const Vec& a_vec);
const Mat outer_product(const Vec& a_vec_0, const Vec& a_vec_1);

class Moments {
 public:
  constexpr Moments(void)
      : m0_m{0},
        m1_m{Vec(0.0, 0.0)},
        m2_m{Mat(Vec(0.0, 0.0), Vec(0.0, 0.0))} {};
  constexpr Moments(const double m0, const Vec& m1, const Mat& m2)
      : m0_m{m0}, m1_m{m1}, m2_m{m2} {};
  double& m0(void) { return m0_m; };
  Vec& m1(void) { return m1_m; };
  Mat& m2(void) { return m2_m; };
  const double& m0(void) const { return m0_m; };
  const Vec& m1(void) const { return m1_m; };
  const Mat& m2(void) const { return m2_m; };
  Moments& operator=(const Moments& a) {
    m0_m = a.m0();
    m1_m = a.m1();
    m2_m = a.m2();
    return (*this);
  };
  Moments operator-(void) const { return Moments(-m0_m, -m1_m, -m2_m); };
  Moments& operator+=(const Moments& a) {
    m0_m += a.m0();
    m1_m += a.m1();
    m2_m += a.m2();
    return (*this);
  };
  Moments& operator-=(const Moments& a) {
    m0_m -= a.m0();
    m1_m -= a.m1();
    m2_m -= a.m2();
    return (*this);
  };
  Moments operator/(const double a) {
    return Moments(m0_m / a, m1_m / a, m2_m / a);
  };
  Moments& operator*=(const double a) {
    m0_m *= a;
    m1_m *= a;
    m2_m *= a;
    return (*this);
  };
  Moments& operator/=(const double a) {
    m0_m /= a;
    m1_m /= a;
    m2_m /= a;
    return (*this);
  };

  friend std::ostream& operator<<(std::ostream& out, const Moments& moments) {
    out << std::setprecision(3) << std::scientific << "m0 = " << moments.m0()
        << "; m1 = " << moments.m1() << "; m2 = " << moments.m2();
    return out;
  }

 private:
  double m0_m;
  Vec m1_m;
  Mat m2_m;
};

const Moments operator*(const double a_scalar, const Moments& a_mom);
const Moments operator*(const Moments& a_mom, const double a_scalar);
const Moments operator+(const Moments& a_mom_0, const Moments& a_mom_1);
const Moments operator-(const Moments& a_mom_0, const Moments& a_mom_1);

using PtAndControl = std::pair<Vec, Vec>;
using BezierList = std::vector<PtAndControl>;
using ReferenceFrame = Mat;

void Print(const BezierList& list);
void ToVTK(const BezierList& list, const std::string& filename);
void ToVTK(const std::vector<BezierList>& list, const std::string& filename);

BezierList RectangleFromBounds(const Vec& x0, const Vec& x1);

class Parabola {
 public:
  constexpr Parabola(void)
      : datum_m{0., 0.},
        frame_m{Vec(1., 0.), Vec(0., 1.)},
        coeff_m{0.},
        above_m{false},
        below_m{false} {};
  constexpr Parabola(const Vec& datum, const ReferenceFrame& frame,
                     const double coeff)
      : datum_m{datum},
        frame_m{frame},
        coeff_m{coeff},
        above_m{false},
        below_m{false} {};
  double& coeff(void) { return coeff_m; };
  Vec& datum(void) { return datum_m; };
  ReferenceFrame& frame(void) { return frame_m; };
  const double& coeff(void) const { return coeff_m; };
  const Vec& datum(void) const { return datum_m; };
  const ReferenceFrame& frame(void) const { return frame_m; };
  void markAsAlwaysAbove(void) {
    above_m = true;
    below_m = false;
  }
  void markAsAlwaysBelow(void) {
    above_m = false;
    below_m = true;
  }
  bool isAlwaysAbove(void) const { return above_m; }
  bool isAlwaysBelow(void) const { return below_m; }
  Parabola createAlwaysAbove(void) {
    Parabola parabola;
    parabola.markAsAlwaysAbove();
    return parabola;
  }
  Parabola createAlwaysBelow(void) {
    Parabola parabola;
    parabola.markAsAlwaysBelow();
    return parabola;
  }
  Parabola& operator=(const Parabola& a) {
    datum_m = a.datum();
    frame_m = a.frame();
    coeff_m = a.coeff();
    return (*this);
  };
  friend std::ostream& operator<<(std::ostream& out, const Parabola& parabola) {
    out << std::setprecision(3) << std::scientific
        << "Datum = " << parabola.datum() << "; Frame = " << parabola.frame()
        << "; Coeff = " << parabola.coeff() << "; Above ? "
        << (parabola.isAlwaysAbove() ? "true" : "false") << "; Below ? "
        << (parabola.isAlwaysBelow() ? "true" : "false");
    return out;
  }

 private:
  Vec datum_m;
  ReferenceFrame frame_m;
  double coeff_m;
  bool above_m, below_m;
};

const Vec BezierPoint(const Vec& p0, const Vec& p1, const Vec& p2,
                      const double t);

std::vector<double> solve_quartic(const double a, const double b,
                                  const double c, const double d);
unsigned int solveP3(double* x, const double a, const double b, const double c);
std::vector<double> solve_cubic(const double a, const double b, const double c,
                                const double d);

std::vector<double> AnalyticIntersections(const Parabola& parabola,
                                          const Vec& p0, const Vec& p1,
                                          const Vec& p2);

double DistanceToParabola(const Parabola& parabola, const Vec& pt);
bool IsBelow(const Parabola& parabola, const Vec& pt);

std::vector<BezierList> ParabolaClipWeilerAtherton(
    const BezierList& original_cell, const Parabola& parabola);
std::vector<BezierList> ParabolaClipWeilerAtherton(
    const std::vector<BezierList>& original_cell, const Parabola& parabola);

double ArcVolume(const Vec& P0, const Vec& P1, const Vec& P2);
double ComputeVolume(const BezierList& cell);
double ComputeVolume(const std::vector<BezierList>& cell);
Moments ComputeMoments(const BezierList& cell);
Moments ComputeMoments(const BezierList& cell, const Parabola& parabola);
Moments ComputeMoments(const BezierList& cell, const Vec& x0, const Vec& x1,
                       const Parabola& parabola);

Vec RK4Point(const Vec& P, const double dt, const double time,
             const Vec (*vel)(const double t, const Vec& P));

std::pair<Vec, Vec> RK4PointAndTangent(
    const Vec& P, const Vec& T, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P));

double Determinant(const Vec& T0, const Vec& T1);

std::pair<bool, double> RayIntersection(const Vec& P0, const Vec& P1,
                                        const Vec& T0, const Vec& T1);

BezierList ConstructPathline(const Vec& P00, const double dt, const double time,
                             Vec (*vel)(const double t, const Vec& P),
                             const int rec_num);

BezierList TransportEdge(const Vec& P00, const Vec& P10, const double dt,
                         const double time,
                         const Vec (*vel)(const double t, const Vec& P),
                         const Mat (*grad_vel)(const double t, const Vec& P),
                         const int recursion_num, const bool add_pathlines,
                         const bool close_flux, const bool correct_area,
                         const double exact_area);

BezierList CreateFluxCell(const Vec& P00, const Vec& P10, const double dt,
                          const double time,
                          const Vec (*vel)(const double t, const Vec& P),
                          const Mat (*grad_vel)(const double t, const Vec& P),
                          const bool correct_area = false,
                          const double exact_area = 0.);

BezierList ParabolaClip(const BezierList& original_cell,
                        const Parabola& parabola,
                        const bool return_parabola_only = false);

BezierList ClipByRectangleAndParabola(const BezierList& original_cell,
                                      const Vec& x0, const Vec& x1,
                                      const Parabola& parabola);

double IntegrateFlux(const Vec& P0, const Vec& P1, const double dt,
                     const double time,
                     const Vec (*vel)(const double t, const Vec& P));
}  // namespace IRL2D

#endif  // EXAMPLES_2D_ADVECTOR_IRL_2D_H_