#ifndef EXAMPLES_2D_ADVECTOR_IRL_2D_H_
#define EXAMPLES_2D_ADVECTOR_IRL_2D_H_

#include <quadmath.h>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <unsupported/Eigen/Polynomials>
#include <vector>
#include "irl/generic_cutting/paraboloid_intersection/moment_contributions.h"
#include "irl/geometry/general/math_vector.h"
#include "irl/helpers/expression_templates.h"
#include "irl/helpers/mymath.h"
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
  void serialize(IRL::ByteBuffer* a_buffer) const {
    a_buffer->pack(vec_m.data(), vec_m.size());
  }
  void unpackSerialized(IRL::ByteBuffer* a_buffer) {
    a_buffer->unpack(vec_m.data(), vec_m.size());
  }

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
const Vec minvec(const Vec& a_vec_0, const Vec& a_vec_1);
const Vec maxvec(const Vec& a_vec_0, const Vec& a_vec_1);

class Mat {
 public:
  constexpr Mat(void) : matrix_m{Vec(0., 0.), Vec(0., 0.)} {};
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
      : datum_m{Vec(0., 0.)},
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
    auto parabola = Parabola();
    parabola.markAsAlwaysAbove();
    return parabola;
  }
  Parabola createAlwaysBelow(void) {
    auto parabola = Parabola();
    parabola.markAsAlwaysBelow();
    return parabola;
  }
  Parabola& operator=(const Parabola& a) {
    datum_m = a.datum();
    frame_m = a.frame();
    coeff_m = a.coeff();
    above_m = a.isAlwaysAbove();
    below_m = a.isAlwaysBelow();
    return (*this);
  }
  void serialize(IRL::ByteBuffer* a_buffer) const {
    datum_m.serialize(a_buffer);
    frame_m[0].serialize(a_buffer);
    frame_m[1].serialize(a_buffer);
    a_buffer->pack(&coeff_m, 1);
    const IRL::UnsignedIndex_t bool_to_int =
        (above_m ? 1 : 0) + 2 * (below_m ? 1 : 0);
    a_buffer->pack(&bool_to_int, 1);
  }
  void unpackSerialized(IRL::ByteBuffer* a_buffer) {
    datum_m.unpackSerialized(a_buffer);
    frame_m[0].unpackSerialized(a_buffer);
    frame_m[1].unpackSerialized(a_buffer);
    a_buffer->unpack(&coeff_m, 1);
    IRL::UnsignedIndex_t int_to_bool = 0;
    a_buffer->unpack(&int_to_bool, 1);
    above_m = (int_to_bool % 2 == 1) ? true : false;
    below_m = (int_to_bool / 2 == 1) ? true : false;
  }
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

template <class ScalarType>
std::vector<ScalarType> AnalyticIntersections(const Parabola& parabola,
                                              const Vec& p0, const Vec& p1,
                                              const Vec& p2);

template <class ScalarType>
ScalarType DistanceToParabola(const Parabola& parabola, const Vec& pt);

bool IsBelow(const Parabola& parabola, const Vec& pt);

std::vector<BezierList> ParabolaClipWeilerAtherton(
    const BezierList& original_cell, const Parabola& parabola);
std::vector<BezierList> ParabolaClipWeilerAtherton(
    const std::vector<BezierList>& original_cell, const Parabola& parabola);

double ArcVolume(const Vec& P0, const Vec& P1, const Vec& P2);
double ComputeArea(const BezierList& cell);
double ComputeArea(const BezierList& cell, const Parabola& parabola);
double ComputeArea(const std::vector<BezierList>& cell);
double ComputeVFrac(const BezierList& cell, const Parabola& parabola);

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

BezierList TransportEdgeMidPoint(
    const Vec& P00, const Vec& P10, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P), const int rec_num,
    const bool add_pathlines = false, const bool close_flux = false,
    const bool correct_area = false, const double exact_area = 0.0);

BezierList TransportEdge(const Vec& P00, const Vec& P10, const double dt,
                         const double time,
                         const Vec (*vel)(const double t, const Vec& P),
                         const Mat (*grad_vel)(const double t, const Vec& P),
                         const int recursion_num,
                         const bool add_pathlines = false,
                         const bool close_flux = false,
                         const bool correct_area = false,
                         const double exact_area = 0.0);

BezierList CreateFluxCell(const Vec& P00, const Vec& P10, const double dt,
                          const double time,
                          const Vec (*vel)(const double t, const Vec& P),
                          const Mat (*grad_vel)(const double t, const Vec& P),
                          const bool correct_area = false,
                          const double exact_area = 0.);

BezierList CreatePreImage(const Vec& X0, const Vec& X1, const double dt,
                          const double time,
                          const Vec (*vel)(const double t, const Vec& P),
                          const Mat (*grad_vel)(const double t, const Vec& P),
                          const bool correct_area,
                          const std::array<double, 4>& exact_area);

BezierList CreateFluxCellMidPoint(
    const Vec& P00, const Vec& P10, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P),
    const bool correct_area = false, const double exact_area = 0.);

BezierList CreatePreImageMidPoint(
    const Vec& X0, const Vec& X1, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P),
    const bool correct_area, const std::array<double, 4>& exact_area);

BezierList CreateLinearPreImage(
    const Vec& X0, const Vec& X1, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P),
    const bool correct_area, const std::array<double, 4>& exact_area);

BezierList TransportLinearEdge(
    const Vec& P00, const Vec& P10, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P), const int rec_num,
    const bool add_pathlines = false, const bool close_flux = false,
    const bool correct_area = false, const double exact_area = 0.0);

BezierList CreateLinearFluxCell(
    const Vec& P00, const Vec& P10, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P),
    const bool correct_area = false, const double exact_area = 0.);

BezierList ParabolaClip(const BezierList& original_cell,
                        const Parabola& parabola,
                        const bool return_parabola_only = false);

BezierList ClipByRectangleAndParabola(const BezierList& original_cell,
                                      const Vec& x0, const Vec& x1,
                                      const Parabola& parabola);

double IntegrateFlux(const Vec& P0, const Vec& P1, const double dt,
                     const double time,
                     const Vec (*vel)(const double t, const Vec& P));

Parabola MatchToVolumeFraction(const BezierList& cell, const Parabola& parabola,
                               const double vfrac);

Parabola MatchToVolumeFractionBrent(const BezierList& cell,
                                    const Parabola& parabola,
                                    const double vfrac);

Parabola MatchToVolumeFractionBisection(const BezierList& cell,
                                        const Parabola& parabola,
                                        const double vfrac,
                                        const int max_bisection_iter = 50);

Parabola MatchToVolumeFractionIllinois(const BezierList& cell,
                                       const Parabola& parabola,
                                       const double vfrac);

std::pair<Vec, Vec> BoundingBox(const BezierList& cell);

std::vector<BezierList> TriangulateCell(const BezierList& cell,
                                        const bool is_preimage);

std::pair<Mat, Vec> MappingMatVec(const BezierList& triangle1,
                                  const BezierList& triangle2);

Vec MappingPoint(const Mat& A, const Vec& b, const Vec& point);

Mat MappingM2(const Mat& A, const Vec& b, const Moments& tri_liq_moment);

Moments ComputeMappedTriangleMoments(const Moments& triangle_liq_moments,
                                     const Mat& A, const Vec& b);

BezierList ComputeTransformedCell(const BezierList& cell,
                                  const bool& toUnitCell);

Mat MappingCellMat(const BezierList& cell, const bool& toUnitCell);

Parabola ComputeTransformedParabola(const BezierList& cell,
                                    const Parabola& parabola,
                                    const bool& toUnitCell);

Moments ComputeTransformedCellMoments(const BezierList& cell,
                                      const Parabola& parabola,
                                      const bool& toUnitCell);

std::vector<Vec> ComputeParticlePositions(const int& N, const Vec& p,
                                          const double& phi,
                                          const double& theta,
                                          const double& hp);

Vec ComputeParticleForce(
    const Vec& x, const std::vector<std::pair<Vec, Vec>>& line_seg_endpoints,
    const double& eta);

std::vector<Vec> InitializeParticlePositions(
    const std::pair<Vec, Vec>& target_endpoints, const double& hp,
    const int& N);

double ComputeParticleForceProjection(const int& N, const double& phi,
                                      const double& theta, const double& hp,
                                      const bool& iswrtPhi,
                                      const std::vector<Vec> particle_forces);

double getCurvature(const Parabola& target_interface,
                    const BezierList& target_cell,
                    const std::vector<Parabola>& interfaces,
                    const std::vector<BezierList>& cells, const int& N,
                    const double& Hp, const double& h, const double& eta);

struct InterfaceEndPoints {
  int xIndex, yIndex;
  double ax, ay, bx, by, tx, ty, nx, ny, vf;
  bool mixed = false;
};

void printParticleData(const std::vector<Vec>& pp, const std::vector<Vec>& pf);

Vec findCircleCenter(const std::vector<Vec>& points);

// std::vector<int> findClosestSegmentIndex(const std::vector<Vec>&
// particle_positions,
//                                          const
//                                          std::vector<std::pair<Vec,Vec>>&
//                                          line_seg_endpoints);

// std::vector<Vec> findClosestSegmentNormal(const std::vector<int>& indices,
//                                           const
//                                           std::vector<std::vector<InterfaceEndPoints>>&
//                                           plic_data);

std::vector<Vec> findSegmentNormals(
    const std::vector<Vec>& particle_positions,
    const std::vector<std::vector<InterfaceEndPoints>>& plic_data);

std::vector<double> computeEta(const std::vector<Vec>& particle_positions,
                               const std::vector<Vec>& pointedPLIC_normals);

void particle_pf(const std::vector<Vec>& pp0, const std::vector<Vec>& pf0,
                 std::vector<Vec>& pp_final, std::vector<Vec>& pf_final,
                 const std::pair<Vec, Vec>& target_endpoints,
                 const std::vector<std::pair<Vec, Vec>>& line_seg_endpoints,
                 const int& N, const double& Hp, const double& h,
                 const double& eta);

void curvature_vareta(
    const std::vector<Vec>& pp0, const std::vector<Vec>& pf0,
    std::vector<Vec>& pp_final, std::vector<Vec>& pf_final,
    const std::pair<Vec, Vec>& target_endpoints,
    const std::vector<std::pair<Vec, Vec>>& line_seg_endpoints,
    const std::vector<std::vector<InterfaceEndPoints>>& plicDataMat,
    const int& N, const double& Hp, const double& h);

std::vector<Vec> generatePoints(
    const std::vector<std::pair<Vec, Vec>>& line_seg_endpoints);

std::vector<Vec> getPoints(const std::pair<Vec, Vec>& endpoints,
                           const int& num_points);

std::vector<Vec> generateParabolaPoints(
    const std::vector<std::pair<Vec, Vec>>& end_points,
    const std::vector<Parabola>& parabola);

Vec getParabolaCenter(const std::pair<Vec, Vec>& end_points,
                      const Parabola& parabola);

double getVfracWeight(double vfrac);

double getDistanceWeight(const Vec& pref, const Vec& ploc, const double& h);

double DistanceWeight(const Vec& pref, const Vec& ploc, const double& h,
                      const double& delta);

double getNormalWeight(const Vec& nref, const Vec& nloc);

// double getNormalGradWeight(const Vec& nref, const Vec& nloc,
//                            const Vec& pref, const Vec& ploc);

Mat estimateFrame(const Vec& circle_center, const Vec& plic_center,
                  const double& r, const Vec& plic_normal, bool& flip_coeff);

Parabola getPrattParabola(const std::vector<IRL2D::Vec>& points,
                          const std::vector<double>& vfw,
                          const std::vector<double>& dw,
                          const std::vector<double>& nw, const Mat& plic_frame,
                          const Vec& plic_center);

Parabola getPrattParabola_localframe(const std::vector<IRL2D::Vec>& points,
                                     const std::vector<double>& vfw,
                                     const std::vector<double>& dw,
                                     const std::vector<double>& nw,
                                     const Mat& plic_frame,
                                     const Vec& plic_center);

Parabola getTaubinParabola(const std::vector<IRL2D::Vec>& points,
                           const std::vector<double>& vfw,
                           const std::vector<double>& dw,
                           const std::vector<double>& nw,
                           const Vec& plic_normal, const Vec& plic_center);

Parabola getTaubinParabola_localframe(const std::vector<IRL2D::Vec>& points,
                                      const std::vector<double>& vfw,
                                      const std::vector<double>& dw,
                                      const std::vector<double>& nw,
                                      const Mat& plic_frame,
                                      const Vec& plic_center);

std::vector<double> getPrattParams(const std::vector<IRL2D::Vec>& points,
                                   const std::vector<double>& vfw,
                                   const std::vector<double>& dw,
                                   const std::vector<double>& nw);

std::vector<double> getTaubinParams(const std::vector<IRL2D::Vec>& points,
                                    const std::vector<double>& vfw,
                                    const std::vector<double>& dw,
                                    const std::vector<double>& nw);

// Parabola getParabolaJibben(const Parabola& target_interface, const
// BezierList& target_cell,
//                            const std::vector<Parabola>& interfaces, const
//                            std::vector<BezierList>& cells);

// for Jibben debugging
struct NeighborInfo {
  Parabola interface;
  BezierList cell;
  int ii_global;
  int jj_global;
  double lvf;
};

Parabola getParabolaJibben(const Parabola& target_interface,
                           const BezierList& target_cell,
                           const std::vector<NeighborInfo>& neighbors,
                           const int i_target, const int j_target);

std::vector<double> getJibbenCoeffs(const Parabola& target_interface,
                                    const BezierList& target_cell,
                                    const std::vector<NeighborInfo>& neighbors,
                                    const std::vector<double>& weights);

// connectivity of interfaces
// class InterfaceConnectivity{

//  public:
//   // connectivity matrix for interfaces
//   std::vector<std::vector<int>> c_matrix;

//   // constructor for initializing connectivity matrix
//   InterfaceConnectivity(int n_mixed) {
//     c_matrix = std::vector<std::vector<int>>(n_mixed,
//     std::vector<int>(n_mixed,0));
//   }

//   // adding more nodes

// };

// partition of unity
struct Wendland {
  const double delta;
  const Vec x_eval;

  Wendland(const double& delta_, const Vec& x_eval_)
      : delta(delta_), x_eval(x_eval_) {}

  double phi(Vec xi) const {
    double r = std::sqrt((x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) +
                         (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]));
    if (r >= delta) return 0.0;
    double s = 1.0 - r / delta;
    return std::pow(s, 4) * (4.0 * r / delta + 1.0);
  }

  double dphidx(Vec xi) const {
    double r = std::sqrt((x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) +
                         (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]));
    if (r >= delta) return 0.0;
    return 1.0 / std::pow(delta, 5.0) * 20.0 * (x_eval[0] - xi[0]) *
           std::pow(r - delta, 3.0);
  }

  double dphidy(Vec xi) const {
    double r = std::sqrt((x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) +
                         (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]));
    if (r >= delta) return 0.0;
    return 1.0 / std::pow(delta, 5.0) * 20.0 * (x_eval[1] - xi[1]) *
           std::pow(r - delta, 3.0);
  }

  double ddphidxx(Vec xi) const {
    double r = std::sqrt((x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) +
                         (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]));
    if (r >= delta) return 0.0;
    if (r < 1e-12) r = 1.0;
    double num = 20.0 * std::pow(r - delta, 2.0) *
                 ((x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) * (r - delta) +
                  x_eval[0] * x_eval[0] * (4.0 * r - delta) +
                  xi[0] * xi[0] * (4.0 * r - delta) +
                  2.0 * x_eval[0] * xi[0] * (-4.0 * r + delta));
    double denom = std::pow(r, 2.0) * std::pow(delta, 5.0);
    return num / denom;
  }

  double ddphidyy(Vec xi) const {
    double r = std::sqrt((x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) +
                         (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]));
    if (r >= delta) return 0.0;
    if (r < 1e-12) r = 1.0;
    double num =
        20.0 * std::pow(r - delta, 2.0) *
        (x_eval[0] * x_eval[0] * (r - delta) + xi[0] * xi[0] * (r - delta) +
         (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) * (4.0 * r - delta) +
         2.0 * x_eval[0] * xi[0] * (-r + delta));
    double denom = std::pow(r, 2.0) * std::pow(delta, 5.0);
    return num / denom;
  }

  double ddphidxy(Vec xi) const {
    double r = std::sqrt((x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) +
                         (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]));
    if (r >= delta) return 0.0;
    if (r < 1e-12) r = 1.0;
    double num = 60.0 * (x_eval[0] - xi[0]) * (x_eval[1] - xi[1]) *
                 std::pow(r - delta, 2.0);
    double denom = r * std::pow(delta, 5.0);
    return num / denom;
  }
};

struct ImplicitSurface {
  const std::vector<Vec> centroids, normals;
  const double kernel_size;

  ImplicitSurface(const std::vector<Vec>& centroids_,
                  const std::vector<Vec>& normals_, const double& kernel_size_)
      : centroids(centroids_), normals(normals_), kernel_size(kernel_size_) {}

  // evaluating functions at x
  double F(Vec x) {
    Wendland w(kernel_size, x);
    double num = 0.0, denom = 0.0;
    for (int i = 0; i < centroids.size(); i++) {
      double wi = w.phi(centroids[i]);
      double Fi = normals[i][0] * (x[0] - centroids[i][0]) +
                  normals[i][1] * (x[1] - centroids[i][1]);
      num += wi * Fi;
      denom += wi;
    }
    if (denom < 1e-12) {
      std::cout << "Sum of weights is too small, denominator = " << denom
                << std::endl;
    }
    return num / denom;
  }

  double Fx(Vec x) {
    Wendland w(kernel_size, x);
    double sum_phi = 0.0;
    double sum_phi_F = 0.0;
    double sum_dphidx = 0.0;
    double sum_dphidx_F = 0.0;
    double sum_phi_dFdx = 0.0;

    for (int i = 0; i < centroids.size(); i++) {
      // phi derivatives
      double phi_i = w.phi(centroids[i]);
      double dphidx_i = w.dphidx(centroids[i]);
      // Fi derivatives
      Vec n = normals[i];
      double F_i =
          n[0] * (x[0] - centroids[i][0]) + n[1] * (x[1] - centroids[i][1]);
      double dFdx_i = n[0];
      // terms for Fx
      sum_phi += phi_i;
      sum_phi_F += phi_i * F_i;
      sum_dphidx += dphidx_i;
      sum_dphidx_F += dphidx_i * F_i;
      sum_phi_dFdx += phi_i * dFdx_i;
    }
    return ((sum_dphidx_F + sum_phi_dFdx) * sum_phi - sum_phi_F * sum_dphidx) /
           (sum_phi * sum_phi);
  }

  double Fy(Vec x) {
    Wendland w(kernel_size, x);
    double sum_phi = 0.0;
    double sum_phi_F = 0.0;
    double sum_dphidy = 0.0;
    double sum_dphidy_F = 0.0;
    double sum_phi_dFdy = 0.0;

    for (int i = 0; i < centroids.size(); i++) {
      // phi derivatives
      double phi_i = w.phi(centroids[i]);
      double dphidy_i = w.dphidy(centroids[i]);
      // Fi derivatives
      Vec n = normals[i];
      double F_i =
          n[0] * (x[0] - centroids[i][0]) + n[1] * (x[1] - centroids[i][1]);
      double dFdy_i = n[1];
      // terms for Fy
      sum_phi += phi_i;
      sum_phi_F += phi_i * F_i;
      sum_dphidy += dphidy_i;
      sum_dphidy_F += dphidy_i * F_i;
      sum_phi_dFdy += phi_i * dFdy_i;
    }
    return ((sum_dphidy_F + sum_phi_dFdy) * sum_phi - sum_phi_F * sum_dphidy) /
           (sum_phi * sum_phi);
  }

  std::vector<double> HessianTerms(Vec x) {
    Wendland w(kernel_size, x);
    double N = 0.0;
    double Nx = 0.0;
    double Ny = 0.0;
    double Nxx = 0.0;
    double Nyy = 0.0;
    double Nxy = 0.0;
    double D = 0.0;
    double Dx = 0.0;
    double Dy = 0.0;
    double Dxx = 0.0;
    double Dyy = 0.0;
    double Dxy = 0.0;

    for (int i = 0; i < centroids.size(); i++) {
      // phi derivatives
      double phi_i = w.phi(centroids[i]);
      double dphidx_i = w.dphidx(centroids[i]);
      double dphidy_i = w.dphidy(centroids[i]);
      double ddphidxx_i = w.ddphidxx(centroids[i]);
      double ddphidyy_i = w.ddphidyy(centroids[i]);
      double ddphidxy_i = w.ddphidxy(centroids[i]);

      // F derivatives
      Vec n = normals[i];
      double F_i =
          n[0] * (x[0] - centroids[i][0]) + n[1] * (x[1] - centroids[i][1]);
      double dFdx_i = n[0];
      double dFdy_i = n[1];
      double ddFdxx_i = 0.0;
      double ddFdyy_i = 0.0;
      double ddFdxy_i = 0.0;

      // numerator terms
      N += phi_i * F_i;
      Nx += dphidx_i * F_i + phi_i * dFdx_i;
      Ny += dphidy_i * F_i + phi_i * dFdy_i;
      Nxx += F_i * ddphidxx_i + phi_i * ddFdxx_i + 2.0 * dphidx_i * dFdx_i;
      Nyy += F_i * ddphidyy_i + phi_i * ddFdyy_i + 2.0 * dphidy_i * dFdy_i;
      Nxy += F_i * ddphidxy_i + phi_i * ddFdxy_i + dphidx_i * dFdy_i +
             dphidy_i * dFdx_i;

      // denominator terms
      D += phi_i;
      Dx += dphidx_i;
      Dy += dphidy_i;
      Dxx += ddphidxx_i;
      Dyy += ddphidyy_i;
      Dxy += ddphidxy_i;
    }
    double Fxx =
        (D * (Nxx * D - N * Dxx) - 2.0 * Dx * (Nx * D - N * Dx)) / (D * D * D);
    double Fyy =
        (D * (Nyy * D - N * Dyy) - 2.0 * Dy * (Ny * D - N * Dy)) / (D * D * D);
    double Fxy = (D * (Nxy * D + Nx * Dy - Ny * Dx - N * Dxy) -
                  2.0 * Dy * (Nx * D - N * Dx)) /
                 (D * D * D);

    return {Fxx, Fyy, Fxy};
  }

  // finite differences
  double Fx_FD(Vec x, double h) {
    return (F({x[0] + h, x[1]}) - F({x[0] - h, x[1]})) / (2.0 * h);
  }
  double Fy_FD(Vec x, double h) {
    return (F({x[0], x[1] + h}) - F({x[0], x[1] - h})) / (2.0 * h);
  }
  double Fxx_FD(Vec x, double h) {
    return (F({x[0] + h, x[1]}) - 2.0 * F(x) + F({x[0] - h, x[1]})) / (h * h);
  }
  double Fyy_FD(Vec x, double h) {
    return (F({x[0], x[1] + h}) - 2.0 * F(x) + F({x[0], x[1] - h})) / (h * h);
  }
  double Fxy_FD(Vec x, double h) {
    return (F({x[0] + h, x[1] + h}) - F({x[0] + h, x[1] - h}) -
            F({x[0] - h, x[1] + h}) + F({x[0] - h, x[1] - h})) /
           (4.0 * h * h);
  }
};

Vec projectToImplicitSurface(const Vec& x0, const std::vector<Vec>& centroids,
                             const std::vector<Vec>& normals,
                             const double& kernel_size, bool& usePlane);

Parabola getPU_interface(const Vec& x0, const std::vector<Vec>& centroids,
                         const std::vector<Vec>& normals,
                         const double& kernel_size, bool& usePlane);

}  // namespace IRL2D

#endif  // EXAMPLES_2D_ADVECTOR_IRL_2D_H_