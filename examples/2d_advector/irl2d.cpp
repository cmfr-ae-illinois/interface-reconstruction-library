#include <boost/multiprecision/cpp_bin_float.hpp>

#include "examples/2d_advector/irl2d.h"

#include <Eigen/Dense>
#include <Eigen/QR>

namespace IRL2D {

using float128 = boost::multiprecision::cpp_bin_float_quad;

const double operator*(const Vec& a_vec_0, const Vec& a_vec_1) {
  return a_vec_0[0] * a_vec_1[0] + a_vec_0[1] * a_vec_1[1];
}
const Vec operator+(const Vec& a_vec_0, const Vec& a_vec_1) {
  return Vec(a_vec_0[0] + a_vec_1[0], a_vec_0[1] + a_vec_1[1]);
}
const Vec operator-(const Vec& a_vec_0, const Vec& a_vec_1) {
  return Vec(a_vec_0[0] - a_vec_1[0], a_vec_0[1] - a_vec_1[1]);
}
const Vec operator*(const double a_scalar, const Vec& a_vec) {
  auto vec_to_return = a_vec;
  vec_to_return *= a_scalar;
  return vec_to_return;
}
const Vec operator*(const Vec& a_vec, const double a_scalar) {
  return a_scalar * a_vec;
}
const Vec operator/(const Vec& a_vec, const double a_scalar) {
  Vec vec_to_return = a_vec;
  vec_to_return /= a_scalar;
  return vec_to_return;
}
const Vec operator*(const Mat& a_mat, const Vec& a_vec) {
  return Vec(a_mat[0] * a_vec, a_mat[1] * a_vec);
}
const Vec minvec(const Vec& a_vec_0, const Vec& a_vec_1) {
  return Vec(std::fmin(a_vec_0[0], a_vec_1[0]),
             std::fmin(a_vec_0[1], a_vec_1[1]));
}
const Vec maxvec(const Vec& a_vec_0, const Vec& a_vec_1) {
  return Vec(std::fmax(a_vec_0[0], a_vec_1[0]),
             std::fmax(a_vec_0[1], a_vec_1[1]));
}
const Mat outer_product(const Vec& a_vec_0, const Vec& a_vec_1) {
  return Mat(Vec(a_vec_0[0] * a_vec_1[0], a_vec_0[0] * a_vec_1[1]),
             Vec(a_vec_0[1] * a_vec_1[0], a_vec_0[1] * a_vec_1[1]));
}
const Mat operator*(const Mat& a_mat_0, const Mat& a_mat_1) {
  const auto a_mat_1_T = a_mat_1.transpose();
  return Mat(Vec(a_mat_0[0] * a_mat_1_T[0], a_mat_0[0] * a_mat_1_T[1]),
             Vec(a_mat_0[1] * a_mat_1_T[0], a_mat_0[1] * a_mat_1_T[1]));
}
const Mat operator+(const Mat& a_mat_0, const Mat& a_mat_1) {
  return Mat(Vec(a_mat_0[0][0] + a_mat_1[0][0], a_mat_0[0][1] + a_mat_1[0][1]),
             Vec(a_mat_0[1][0] + a_mat_1[1][0], a_mat_0[1][1] + a_mat_1[1][1]));
}
const Mat operator-(const Mat& a_mat_0, const Mat& a_mat_1) {
  return Mat(Vec(a_mat_0[0][0] - a_mat_1[0][0], a_mat_0[0][1] - a_mat_1[0][1]),
             Vec(a_mat_0[1][0] - a_mat_1[1][0], a_mat_0[1][1] - a_mat_1[1][1]));
}
const Mat operator*(const double a_scalar, const Mat& a_mat) {
  return Mat(a_scalar * a_mat[0], a_scalar * a_mat[1]);
}
const Mat operator*(const Mat& a_mat, const double a_scalar) {
  return a_scalar * a_mat;
}
const Moments operator+(const Moments& a_mom_0, const Moments& a_mom_1) {
  return Moments(a_mom_0.m0() + a_mom_1.m0(), a_mom_0.m1() + a_mom_1.m1(),
                 a_mom_0.m2() + a_mom_1.m2());
}
const Moments operator-(const Moments& a_mom_0, const Moments& a_mom_1) {
  return Moments(a_mom_0.m0() - a_mom_1.m0(), a_mom_0.m1() - a_mom_1.m1(),
                 a_mom_0.m2() - a_mom_1.m2());
}
const Moments operator*(const double a_scalar, const Moments& a_mom) {
  return Moments(a_scalar * a_mom.m0(), a_scalar * a_mom.m1(),
                 a_scalar * a_mom.m2());
}
const Moments operator*(const Moments& a_mom, const double a_scalar) {
  return a_scalar * a_mom;
}

void Print(const BezierList& list) {
  std::cout << "Cell:";
  for (int i = 0; i < list.size(); i++) {
    std::cout << "\n  Vec" << i << "  = " << list[i].first;
    std::cout << "\n  Ctrl" << i << " = " << list[i].second;
  }
  std::cout << std::endl;
}

void ToVTK(const std::vector<BezierList>& list, const std::string& filename) {
  const int nsamples = 10;
  const double w = 1. / static_cast<double>(nsamples);
  int npoints = 0;
  for (int i = 0; i < list.size(); i++) {
    npoints += nsamples * list[i].size();
  }
  std::ofstream file;
  file.open(filename + std::string(".vtu"));
  file << "<VTKFile type=\"UnstructuredGrid\">\n<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << npoints << "\" NumberOfCells=\""
       << list.size() << "\">\n";
  file << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";
  for (int i = 0; i < list.size(); i++) {
    for (int n = 0; n < list[i].size(); n++) {
      const auto P0 = list[i][n].first;
      const auto P1 = list[i][n].second;
      const auto P2 = list[i][(n + 1) % list[i].size()].first;
      for (int m = 0; m < nsamples; m++) {
        const double t = static_cast<double>(m) * w;
        const auto P =
            P0 * (1. - t) * (1. - t) + P1 * 2. * (1 - t) * t + P2 * t * t;
        file << std::scientific << std::setprecision(15) << P.x() << " "
             << P.y() << " 0. \n";
      }
    }
  }
  file << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" "
          "Name=\"connectivity\" format=\"ascii\">\n";
  int offset = 0;
  for (int i = 0; i < list.size(); i++) {
    int count = 0;
    for (int n = 0; n < list[i].size(); n++) {
      for (int m = 0; m < nsamples; m++) {
        file << offset + count++ << " ";
      }
    }
    offset += count;
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  int count = 0;
  for (int i = 0; i < list.size(); i++) {
    count += nsamples * list[i].size();
    file << count << " ";
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (int i = 0; i < list.size(); i++) {
    file << "7 ";
  }
  file << "\n</DataArray>\n";
  file << "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  file.close();

  npoints = 0;
  for (int i = 0; i < list.size(); i++) {
    npoints += 2 * list[i].size();
  }
  file.open(filename + std::string("_skeleton.vtu"));
  file << "<VTKFile type=\"UnstructuredGrid\">\n<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << npoints << "\" NumberOfCells=\""
       << list.size() << "\">\n";
  file << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";
  for (int i = 0; i < list.size(); i++) {
    for (int n = 0; n < list[i].size(); n++) {
      const auto P0 = list[i][n].first;
      const auto P1 = list[i][n].second;
      file << std::scientific << std::setprecision(15) << P0.x() << " "
           << P0.y() << " 0. \n";
      file << std::scientific << std::setprecision(15) << P1.x() << " "
           << P1.y() << " 0. \n";
    }
  }
  file << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" "
          "Name=\"connectivity\" format=\"ascii\">\n";
  offset = 0;
  for (int i = 0; i < list.size(); i++) {
    int count = 0;
    for (int n = 0; n < list[i].size(); n++) {
      file << offset + count++ << " ";
      file << offset + count++ << " ";
    }
    offset += count;
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  count = 0;
  for (int i = 0; i < list.size(); i++) {
    count += 2 * list[i].size();
    file << count << " ";
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (int i = 0; i < list.size(); i++) {
    file << "7 ";
  }
  file << "\n</DataArray>\n";
  file << "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  file.close();
}

void ToVTK(const BezierList& list, const std::string& filename) {
  ToVTK(std::vector<BezierList>{list}, filename);
}

BezierList RectangleFromBounds(const Vec& x0, const Vec& x1) {
  Vec p00(x0.x(), x0.y()), p10(x1.x(), x0.y()), p11(x1.x(), x1.y()),
      p01(x0.x(), x1.y());
  PtAndControl pc0{p00, .5 * (p00 + p10)};
  PtAndControl pc1{p10, .5 * (p10 + p11)};
  PtAndControl pc2{p11, .5 * (p11 + p01)};
  PtAndControl pc3{p01, .5 * (p01 + p00)};
  return BezierList{pc0, pc1, pc2, pc3};
}

unsigned int solveP3(double* x, const double a, const double b,
                     const double c) {
  /*---------------------------------------------------------------------------
   *	solve cubic equation x^3 + a*x^2 + b*x + c
   *	x - array of size 3
   *	In case 3 real roots: => x[0], x[1], x[2], return 3
   *	        2 real roots: x[0], x[1],          return 2
   *	        1 real root : x[0], x[1] ± i*x[2], return 1
   */
  double a2 = a * a;
  double q = (a2 - 3 * b) / 9.;
  double r = (a * (2. * a2 - 9. * b) + 27. * c) / 54.;
  double r2 = r * r;
  double q3 = q * q * q;
  const double tol = std::numeric_limits<double>::epsilon();
  double A, B;
  if (r2 < q3) {
    double t = r / std::sqrt(q3);
    if (t < -1.) t = -1.;
    if (t > 1.) t = 1.;
    t = std::acos(t);
    q = -2. * std::sqrt(q);
    x[0] = q * std::cos(t / 3.) - a / 3.;
    x[1] = q * std::cos((t + 2. * M_PI) / 3.) - a / 3.;
    x[2] = q * std::cos((t - 2. * M_PI) / 3.) - a / 3.;
    return 3;
  } else {
    A = -std::pow(std::fabs(r) + std::sqrt(r2 - q3), 1. / 3.);
    if (r < 0.) A = -A;
    B = (0. == A ? 0. : q / A);
    x[0] = (A + B) - a / 3.;
    x[1] = -0.5 * (A + B) - a / 3.;
    x[2] = 0.5 * std::sqrt(3.) * (A - B);
    if (std::fabs(x[2]) < tol) {
      x[2] = x[1];
      return 2;
    }
    return 1;
  }
}

std::vector<double> solve_quartic(const double a, const double b,
                                  const double c, const double d) {
  // Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
  const double a3 = -b;
  const double b3 = a * c - 4. * d;
  const double c3 = -a * a * d - c * c + 4. * b * d;
  const double tol = std::numeric_limits<double>::epsilon();
  double x3[3];
  const unsigned int iZeroes = solveP3(x3, a3, b3, c3);
  double q1, q2, p1, p2, D, sqD, y;

  y = x3[0];
  // THE ESSENCE - choosing Y with maximal absolute value !
  if (iZeroes != 1) {
    if (std::fabs(x3[1]) > std::fabs(y)) y = x3[1];
    if (std::fabs(x3[2]) > std::fabs(y)) y = x3[2];
  }

  // h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

  D = y * y - 4. * d;
  if (std::fabs(D) < tol) {  // in other words - D==0
    q1 = q2 = y * 0.5;
    // g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
    D = a * a - 4. * (b - y);
    if (std::fabs(D) < tol)  // in other words - D==0
    {
      p1 = p2 = a * 0.5;
    } else {
      sqD = std::sqrt(D);
      p1 = (a + sqD) * 0.5;
      p2 = (a - sqD) * 0.5;
    }
  } else {
    sqD = std::sqrt(D);
    q1 = (y + sqD) * 0.5;
    q2 = (y - sqD) * 0.5;
    // g1+g2 = a && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
    p1 = (a * q1 - c) / (q1 - q2);
    p2 = (c - a * q2) / (q1 - q2);
  }

  std::vector<double> retval;

  // solving quadratic eq. - x^2 + p1*x + q1 = 0
  D = p1 * p1 - 4. * q1;
  if (D < 0.0 && std::abs(std::sqrt(-D) * 0.5) < tol) {
    retval.push_back(-p1 * 0.5);
  } else {
    sqD = sqrt(D);
    retval.push_back((-p1 + sqD) * 0.5);
    retval.push_back((-p1 - sqD) * 0.5);
  }

  // solving quadratic eq. - x^2 + p2*x + q2 = 0
  D = p2 * p2 - 4. * q2;
  if (D < 0.0 && abs(sqrt(-D) * 0.5) < tol) {
    retval.push_back(-p2 * 0.5);
  } else {
    sqD = std::sqrt(D);
    retval.push_back((-p2 + sqD) * 0.5);
    retval.push_back((-p2 - sqD) * 0.5);
  }

  return retval;
}

std::vector<double> solve_cubic(const double a, const double b, const double c,
                                const double d) {
  double delta_zero = b * b - 3. * a * c;
  double delta_one = 2. * b * b * b - 9. * a * b * c + 27. * a * a * d;
  double tol = std::numeric_limits<double>::epsilon();
  std::vector<double> retval;

  if (delta_zero <= 0 && delta_one <= 0) {
    retval.push_back(-b / (3. * a));
    retval.push_back(-b / (3. * a));
    retval.push_back(-b / (3. * a));
  } else {
    double C_part = std::pow(
        delta_one * delta_one - 4. * delta_zero * delta_zero * delta_zero,
        1. / 2.);
    double C1 = std::pow((delta_one + C_part) / 2., 1. / 3.);
    double C2 = std::pow((delta_one - C_part) / 2., 1. / 3.);
    double C;

    if (C1 != 0.) {
      C = C1;
    } else {
      C = C2;
    }

    double epsilon = (-1. + std::pow(-3., 1. / 2.)) / 2.;
    std::vector<std::complex<double>> roots;
    for (int i = 0; i < 3; i++) {
      double roots_part = C * std::pow(epsilon, i);
      roots.push_back(-(b + roots_part + delta_zero / roots_part) / (3. * a));
    }
    for (const std::complex<double>& root : roots) {
      if (root.imag() < tol) {
        retval.push_back(root.real());
      }
    }
  }

  return retval;
}

template <class ScalarType>
std::vector<ScalarType> AnalyticIntersections(const Parabola& parabola,
                                              const Vec& p0, const Vec& p1,
                                              const Vec& p2) {
  const ScalarType double_eps =
      static_cast<ScalarType>(std::numeric_limits<double>::epsilon());
  const ScalarType eps = std::numeric_limits<ScalarType>::epsilon();
  const ScalarType tol = 1.0e-0 * eps;
  const ScalarType eps_distance = 1.0e2 * eps;
  const auto ZERO = ScalarType(0), ONE = ScalarType(1), TWO = ScalarType(2),
             FOUR = ScalarType(4), SIX = ScalarType(6), EIGHT = ScalarType(8),
             TWELVE = ScalarType(12);
  // If the bezier arc is linear and coincident with the x axis
  if constexpr (std::is_same_v<ScalarType, double>) {
    if (std::abs(p0.y()) < double_eps && std::abs(p1.y()) < double_eps &&
        std::abs(p2.y()) < double_eps) {
      return std::vector<ScalarType>{};
    }
  }
  const auto a = -static_cast<ScalarType>(parabola.coeff());
  const auto x1 = static_cast<ScalarType>(p0.x()),
             y1 = static_cast<ScalarType>(p0.y()),
             x2 = static_cast<ScalarType>(p2.x()),
             y2 = static_cast<ScalarType>(p2.y());
  const auto ctrl_to_mid = p1 - 0.5 * (p0 + p2);
  if (std::fabs(ctrl_to_mid.x()) < 10.0 * double_eps &&
      std::fabs(ctrl_to_mid.y()) < 10.0 * double_eps) {
    // The bezier arc is linear
    std::vector<ScalarType> t_vals(0);
    t_vals.reserve(2);
    auto A = a * (x2 - x1) * (x2 - x1);
    auto B = TWO * a * x1 * (x2 - x1) + y1 - y2;
    auto C = a * x1 * x1 - y1;
    ScalarType discriminant = B * B - FOUR * A * C;
    if (discriminant > ZERO) {
      if (A != ZERO) {
        discriminant = sqrt(discriminant);
        const ScalarType q = -(B + copysign(discriminant, B)) / TWO;
        const ScalarType safe_q = q == ZERO ? copysign(eps, q) : q;
        const ScalarType sol1 = q / A;
        const ScalarType sol2 = C / safe_q;
        if (!isnan(sol1) && sol1 > -eps_distance && sol1 < ONE + eps_distance) {
          t_vals.push_back(sol1);
        }
        if (!isnan(sol2) && sol2 > -eps_distance && sol2 < ONE + eps_distance) {
          t_vals.push_back(sol2);
        }
      } else {
        const ScalarType sol = -C / B;
        if (!isnan(sol) && sol > -eps_distance && sol < ONE + eps_distance) {
          t_vals.push_back(sol);
        }
      }
    }
    if (t_vals.size() > 0) {
      std::sort(t_vals.begin(), t_vals.end());
    }
    for (int i = 0; i < t_vals.size(); i++) {
      auto t = t_vals[i];
      auto froot = A * t * t + B * t + C;
      if (fabs(froot) > double_eps) {
        auto tn = t;
        for (int j = 0; j < 20; j++) {
          const auto fn = A * tn * tn + B * tn + C;
          auto fpn = 2. * A * tn + B;
          if (fabs(fpn) < eps) break;
          tn -= fn / fpn;
        }
        t_vals[i] = tn;
        t = tn;
        froot = A * t * t + B * t + C;
        // if (fabs(froot) > 10. * double_eps) {
        //   std::cout << "1 -- f(" << t << ") = " << froot << std::endl;
        // }
      }
    }
    return t_vals;
  } else {
    const auto xc = static_cast<ScalarType>(p1.x()),
               yc = static_cast<ScalarType>(p1.y());
    auto A = a * (x1 * x1 + TWO * x1 * x2 - FOUR * x1 * xc + x2 * x2 -
                  FOUR * x2 * xc + FOUR * xc * xc);
    auto B = a * (-FOUR * x1 * x1 - FOUR * x1 * x2 + TWELVE * x1 * xc +
                  FOUR * x2 * xc - EIGHT * xc * xc);
    auto C = a * (SIX * x1 * x1 + TWO * x1 * x2 - TWELVE * x1 * xc +
                  FOUR * xc * xc) -
             y1 - y2 + TWO * yc;
    auto D = FOUR * a * x1 * (xc - x1) + TWO * y1 - TWO * yc;
    auto E = a * x1 * x1 - y1;
    std::array<ScalarType, 5> coeffs{A, B, C, D, E};

    auto max_coeff = ZERO;
    for (int i = 0; i < 5; i++) {
      max_coeff = fmax(max_coeff, fabs(coeffs[i]));
    }
    for (int i = 0; i < 5; i++) {
      coeffs[i] /= max_coeff;
    }

    std::vector<ScalarType> t_vals(0);
    t_vals.reserve(4);

    // int nnonzero = 5;
    // for (int i = 0; i < 5; i++) {
    //   if (fabs(coeffs[i]) < tol) {
    //     nnonzero--;
    //   } else {
    //     break;
    //   }
    // }

    // Eigen::PolynomialSolver<ScalarType, Eigen::Dynamic> solver;
    // Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> poly_coeffs;
    // poly_coeffs.resize(nnonzero, 1);
    // for (int i = 0; i < nnonzero; i++) {
    //   poly_coeffs[i] = coeffs[4 - i];
    // }
    // solver.compute(poly_coeffs);
    // for (int i = 0; i < solver.roots().size(); ++i) {
    //   if (fabs(std::imag(solver.roots()[i])) < ScalarType(1.0e4) * eps) {
    //     const ScalarType sol = std::real(solver.roots()[i]);
    //     if (!isnan(sol) && sol >= ZERO && sol <= ONE) {
    //       t_vals.push_back(sol);
    //     }
    //   }
    // }

    // if constexpr (std::is_same_v<double, ScalarType>) {
    // if (coeffs[0] != ZERO && fabs(coeffs[0]) <= eps) {
    // std::cout << "Coeff[0] = " << coeffs[0] << std::endl;
    // std::cout << "ABCDE    = " << A << ", " << B << ", " << C << ", " << D
    //           << ", " << E << std::endl;
    // std::cout << "a        = " << a << std::endl;
    // std::cout << "p0, p1, p2  = " << p0 << p1 << p2 << std::endl;
    // std::cout << "ctrltomidmag = " << ctrl_to_mid.magnitude() << std::endl;
    // return std::vector<ScalarType>{ZERO};
    // }
    //   if (coeffs[0] == ZERO && coeffs[1] != ZERO && fabs(coeffs[1]) <=
    //   eps) {
    //     std::cout << "Coeff[0] = " << coeffs[0] << std::endl;
    //     std::cout << "Coeff[1] = " << coeffs[1] << std::endl;
    //     std::cout << "ABCDE    = " << A << ", " << B << ", " << C << ", "
    //     <<
    //     D
    //               << ", " << E << std::endl;
    //     return std::vector<ScalarType>{ZERO};
    //   }
    //   // if (coeffs[0] == ZERO && coeffs[1] == ZERO && coeffs[2] != ZERO
    //   &&
    //   //     abs(coeffs[2]) <= eps) {
    //   //   return std::vector<ScalarType>{ZERO};
    //   // }
    //   // if (coeffs[0] == ZERO && coeffs[1] == ZERO && coeffs[2] == ZERO
    //   &&
    //   //     coeffs[3] != ZERO && abs(coeffs[3]) <= eps) {
    //   //   return std::vector<ScalarType>{0.0};
    //   // }
    // }

    if (abs(A) > tol) {
      // t_solutions = solve_quartic(B / A, C / A, D / A, E / A);
      Eigen::PolynomialSolver<ScalarType, Eigen::Dynamic> solver;
      Eigen::Matrix<ScalarType, 5, 1> coeff(E, D, C, B, A);
      solver.compute(coeff);
      for (int i = 0; i < solver.roots().size(); ++i) {
        if (fabs(std::imag(solver.roots()[i])) < ScalarType(1.0e4) * eps) {
          const ScalarType sol = std::real(solver.roots()[i]);
          if (!isnan(sol) && sol >= -eps_distance &&
              sol <= ONE + eps_distance) {
            t_vals.push_back(sol);
          }
        }
      }
    } else if (abs(B) > tol) {
      // t_solutions = solve_cubic(B, C, D, E);
      Eigen::PolynomialSolver<ScalarType, Eigen::Dynamic> solver;
      Eigen::Matrix<ScalarType, 4, 1> coeff(E, D, C, B);
      solver.compute(coeff);
      for (int i = 0; i < solver.roots().size(); ++i) {
        if (fabs(std::imag(solver.roots()[i])) < ScalarType(1.0e4) * eps) {
          const ScalarType sol = std::real(solver.roots()[i]);
          if (!isnan(sol) && sol >= -eps_distance &&
              sol <= ONE + eps_distance) {
            t_vals.push_back(sol);
          }
        }
      }
    } else if (abs(C) > tol) {
      ScalarType discriminant = D * D - FOUR * C * E;
      if (discriminant > ZERO) {
        discriminant = sqrt(discriminant);
        const ScalarType q = -(D + copysign(discriminant, D)) / TWO;
        const ScalarType safe_q = fabs(q) < eps ? copysign(eps, q) : q;
        const ScalarType sol1 = q / C;
        const ScalarType sol2 = E / safe_q;
        if (!isnan(sol1) && sol1 >= -eps_distance &&
            sol1 <= ONE + eps_distance) {
          t_vals.push_back(sol1);
        }
        if (!isnan(sol2) && sol2 >= -eps_distance &&
            sol2 <= ONE + eps_distance) {
          t_vals.push_back(sol2);
        }
      }
      // const auto quadratic_solution = IRL::solveQuadratic(
      //     static_cast<double>(C), static_cast<double>(D),
      //     static_cast<double>(E));
      // for (int i = 0; i < quadratic_solution.size(); i++) {
      //   t_solutions.push_back(quadratic_solution[i]);
      // }
    } else if (abs(D) > tol) {
      const ScalarType sol = -E / D;
      if (!isnan(sol) && sol >= -eps_distance && sol <= ONE + eps_distance) {
        t_vals.push_back(sol);
      }
    }

    for (int i = 0; i < t_vals.size(); i++) {
      auto t = t_vals[i];
      auto froot = A * t * t * t * t + B * t * t * t + C * t * t + D * t + E;
      if (fabs(froot) > double_eps) {
        auto tn = t;
        for (int j = 0; j < 20; j++) {
          const auto fn = A * tn * tn * tn * tn + B * tn * tn * tn +
                          C * tn * tn + D * tn + E;
          auto fpn = 4. * A * tn * tn * tn + 3. * B * tn * tn + 2. * C * tn + D;
          if (fabs(fpn) < eps) break;
          tn -= fn / fpn;
        }
        t_vals[i] = tn;
        t = tn;
        froot = A * t * t * t * t + B * t * t * t + C * t * t + D * t + E;
        // if (fabs(froot) > 10. * double_eps) {
        //   std::cout << "2 -- f(" << t << ") = " << froot << std::endl;
        // }
      }
    }

    if (t_vals.size() > 0) {
      std::sort(t_vals.begin(), t_vals.end());
    }

    if (t_vals.size() > 1) {
      for (int i = 0; i < t_vals.size() - 1; i++) {
        if (fabs(t_vals[i] - t_vals[i + 1]) < 100. * double_eps) {
          std::cout << "t0 = " << t_vals[i] << ", t1 = " << t_vals[i + 1]
                    << std::endl;
        }
      }
    }
    return t_vals;
  }
}

const Vec BezierPoint(const Vec& p0, const Vec& p1, const Vec& p2,
                      const double t) {
  return p0 * (1. - t) * (1. - t) + p1 * 2. * (1 - t) * t + p2 * t * t;
}

bool IsBelow(const Parabola& parabola, const Vec& pt) {
  return pt.y() < -parabola.coeff() * pt.x() * pt.x();
}

template <class ScalarType>
ScalarType DistanceToParabola(const Parabola& parabola, const Vec& pt) {
  return static_cast<ScalarType>(pt.y()) +
         static_cast<ScalarType>(parabola.coeff()) *
             static_cast<ScalarType>(pt.x()) * static_cast<ScalarType>(pt.x());
}

std::vector<BezierList> ParabolaClipWeilerAtherton(
    const BezierList& original_cell, const Parabola& parabola) {
  // If empty cell, return empty cell
  if (original_cell.size() == 0 || parabola.isAlwaysBelow())
    return std::vector<BezierList>(0);
  if (parabola.isAlwaysAbove()) return std::vector<BezierList>{original_cell};

  // Specify constants
  double tol = std::numeric_limits<double>::epsilon();
  bool parabola_contains_vertex = true;
  const int itmax = 100;

  // Initialize list of clipped regions
  std::vector<BezierList> clipped_cell;
  // Reserve 5 times original size (i.e., 4 intersections per arc)
  clipped_cell.reserve(5 * original_cell.size());

  // Store parabola properties
  Vec datum = parabola.datum();
  ReferenceFrame frame = parabola.frame();
  double coeff = parabola.coeff();

  // This while loop is needed in case we might want to nudge the parabola
  int it = 0;
  while (parabola_contains_vertex && it < itmax) {
    parabola_contains_vertex = false;
    it++;
    // Copy cell in frame of reference of parabola
    const int nvert_init = original_cell.size();
    BezierList cell, cell_with_intersections, intersections;
    cell.assign(original_cell.begin(), original_cell.end());
    for (int i = 0; i < nvert_init; i++) {
      cell[i].first = frame * (cell[i].first - datum);
      cell[i].second = frame * (cell[i].second - datum);
    }

    // Compute all intersections and insert them in cell
    int nintersections = 0;
    std::vector<int>
        vertex_nature;  // 1 = below, 0 = above, 2 = entry, 3 = exit
    std::vector<std::pair<double, int>> sorted_list_for_mapping;
    for (int i = 0; i < nvert_init; i++) {
      // Get points of current bezier arc
      const auto p0 = cell[i].first;
      const auto p1 = cell[i].second;
      const auto p2 = cell[(i + 1) % cell.size()].first;
      // Add p0 and p1 to new list
      cell_with_intersections.push_back(PtAndControl{p0, p1});
      // Is p0 below or above parabola?
      vertex_nature.push_back(IsBelow(parabola, p0) ? 1 : 0);
      // Compute t-values of potential intersections
      auto t_vals = AnalyticIntersections<double>(parabola, p0, p1, p2);
      // Tmp variables needs for arc splitting
      double t0 = 0.;
      auto p0new = p0, p1new = p1;
      // Loop over all detected intersections
      for (int j = 0; j < t_vals.size(); j++) {
        const double t = t_vals[j];
        // If t = 0 or t = 1, we must nudge! Go back to beginning of function
        if (std::fabs(t) < tol || std::fabs(1. - t) < tol) {
          parabola_contains_vertex = true;
          break;
        }
        //// Spilling of the arc (using de Casteljau algorithm)
        // First: Update previous control point
        cell_with_intersections.back().second =
            p0new + (t - t0) / (1. - t0) * (p1new - p0new);
        // Second: Add intersection and new control point
        p0new = BezierPoint(p0, p1, p2, t);
        p1new = p1new + (t - t0) / (1. - t0) * (p2 - p1new);
        const int inter_id = cell_with_intersections.size();
        sorted_list_for_mapping.push_back(std::make_pair(-p0new.x(), inter_id));
        cell_with_intersections.push_back(PtAndControl{p0new, p1new});
        nintersections++;
        // Check if intersection is entry of exit
        // If: previous vertex is below or exit, then current vertex is entry
        // Otherwise: the current vertex is exit
        vertex_nature.push_back(
            (vertex_nature.back() == 1 || vertex_nature.back() == 2) ? 3 : 2);
        t0 = t;
      }
      if (parabola_contains_vertex) break;
    }

    // If non-even # of intersections or parabola intersects with vertex,
    // nudge!
    if (nintersections % 2 != 0 || parabola_contains_vertex) {
      datum += 10. * tol * frame[1];
      continue;
    }

    // If no intersections at all; leave loop
    if (nintersections == 0) {
      if (vertex_nature[0] == 1) {
        clipped_cell.push_back(cell_with_intersections);
      }
      break;
    }

    // Now create mappings to/from intersections
    const int nvertices = cell_with_intersections.size();
    std::sort(sorted_list_for_mapping.begin(), sorted_list_for_mapping.end());
    std::vector<int> mapping_inter_to_cell(nintersections);
    std::vector<int> mapping_cell_to_inter(nvertices);
    for (int i = 0; i < nvertices; i++) {
      mapping_cell_to_inter[i] = -1;
    }
    for (int i = 0; i < nintersections; i++) {
      mapping_inter_to_cell[i] = sorted_list_for_mapping[i].second;
      mapping_cell_to_inter[mapping_inter_to_cell[i]] = i;
    }

    ///////// Weiler-Atherton clipping algorithm
    int count_entries = 0;
    while (count_entries < nintersections / 2) {
      // Create new bezier list
      BezierList closed_clipped_region;
      // First, find entry
      int start = -1;
      for (int i = 0; i < nvertices; i++) {
        if (vertex_nature[i] == 2) {
          start = i;
          break;
        }
      }
      assert(start >= 0);
      // Add entry to new bezier list and set as "above"
      closed_clipped_region.push_back(cell_with_intersections[start]);
      vertex_nature[start] = 0;
      count_entries++;
      // Store start id
      const int check_start = start;
      // Loop over cell until next exit or start point itself;
      for (int i = 0; i < nvertices; i++) {
        // Fint next vertex on cell
        const int next_id = (start + 1 + i) % nvertices;
        // If next vertex is start, we are done
        if (next_id == check_start) {
          break;
        }
        // Add next vertex to list
        closed_clipped_region.push_back(cell_with_intersections[next_id]);
        // If next vertex is exit, we switch to parabola and find next entry
        if (vertex_nature[next_id] == 3) {
          // Find next entry in cell
          const int next_entry_id =
              mapping_inter_to_cell[mapping_cell_to_inter[next_id] + 1];
          // Set entry flag to zero
          vertex_nature[next_entry_id] = 0;
          if (next_entry_id != check_start) count_entries++;
          // Move i until before the next entry
          while ((start + 2 + i) % nvertices != next_entry_id) i++;
          // Create new control point for bezier arc extracted from parabola
          const auto p0 = cell_with_intersections[next_id].first;
          const auto p2 = cell_with_intersections[next_entry_id].first;
          closed_clipped_region.back().second =
              Vec(0.5 * (p0.x() + p2.x()),
                  p0.y() + coeff * (p0.x() - p2.x()) * p0.x());
        }
        vertex_nature[next_id] = 0;
      }
      clipped_cell.push_back(closed_clipped_region);
    }
  }

  if (parabola_contains_vertex) {
    std::cout << "WARNING: Parabola contains vertex!" << itmax
              << " nudges were not enough." << std::endl;
  }

  // Move clipped cell back to canonical frame of reference
  for (int i = 0; i < clipped_cell.size(); i++) {
    for (int j = 0; j < clipped_cell[i].size(); j++) {
      clipped_cell[i][j].first =
          parabola.datum() + frame.transpose() * clipped_cell[i][j].first;
      clipped_cell[i][j].second =
          parabola.datum() + frame.transpose() * clipped_cell[i][j].second;
    }
  }

  return clipped_cell;
}

std::vector<BezierList> ParabolaClipWeilerAtherton(
    const std::vector<BezierList>& original_cell, const Parabola& parabola) {
  // If empty cell, return empty cell
  if (original_cell.size() == 0 || parabola.isAlwaysBelow())
    return std::vector<BezierList>(0);
  if (parabola.isAlwaysAbove()) return original_cell;

  auto clipped_cell = ParabolaClipWeilerAtherton(original_cell[0], parabola);
  for (int i = 1; i < original_cell.size(); i++) {
    const auto tmp_clipped_cell =
        ParabolaClipWeilerAtherton(original_cell[i], parabola);
    clipped_cell.insert(clipped_cell.end(), tmp_clipped_cell.begin(),
                        tmp_clipped_cell.end());
  }
  return clipped_cell;
}

// std::vector<BezierList> ClipByRectangleAndParabola(
//     const BezierList& original_cell, const Vec& x0, const Vec& x1,
//     const Parabola& parabola) {
//   std::array<Parabola, 4> localizers;
//   localizers[0] = Parabola(x1, ReferenceFrame(0.), 0.);
//   localizers[1] = Parabola(x0, ReferenceFrame(M_PI / 2.), 0.);
//   localizers[2] = Parabola(x0, ReferenceFrame(M_PI), 0.);
//   localizers[3] = Parabola(x1, ReferenceFrame(3. * M_PI / 2.), 0.);
//   auto clipped_cell = ParabolaClip(original_cell, parabola);
//   for (int i = 0; i < 4; i++) {
//     clipped_cell = ParabolaClip(clipped_cell, localizers[i]);
//   }
//   return clipped_cell;
// }

double ArcVolume(const Vec& P0, const Vec& P1, const Vec& P2) {
  return -(P2.x() * (2. * P0.y() + P1.y()) + P1.x() * (P0.y() - P2.y()) -
           P0.x() * (P1.y() + 2. * P2.y())) /
         3.;
}

double ComputeArea(const BezierList& cell) {
  if (cell.size() == 0) return 0.0;

  double sum = 0.0, c = 0.0;
  double area = 0.0;
  for (int i = 0; i < cell.size(); i++) {
    const auto p0 = cell[i].first;
    const auto p1 = cell[i].second;
    const auto p2 = cell[(i + 1) % cell.size()].first;
    const double m0_contrib =
        -(2. * p1.x() * p0.y() + p2.x() * p0.y() - 2. * p0.x() * p1.y() +
          2. * p2.x() * p1.y() - p0.x() * p2.y() - 2. * p1.x() * p2.y()) /
        6.;
    double y = m0_contrib - c;
    double t = sum + y;
    c = (t - sum) - y;
    sum = t;
    area = sum;
    // area -= (2. * p1.x() * p0.y() + p2.x() * p0.y() - 2. * p0.x() * p1.y()
    // +
    //          2. * p2.x() * p1.y() - p0.x() * p2.y() - 2. * p1.x() * p2.y())
    //          /
    //         6.;
  }
  return area;
}

double ComputeArea(const BezierList& cell, const Parabola& parabola) {
  if (cell.size() == 0 || parabola.isAlwaysBelow()) return 0.0;

  const auto clipped_cell = ParabolaClip(cell, parabola);
  double area = 0.0;
  for (int i = 0; i < clipped_cell.size(); i++) {
    const auto p0 = clipped_cell[i].first;
    const auto p1 = clipped_cell[i].second;
    const auto p2 = clipped_cell[(i + 1) % clipped_cell.size()].first;
    area += -(2. * p1.x() * p0.y() + p2.x() * p0.y() - 2. * p0.x() * p1.y() +
              2. * p2.x() * p1.y() - p0.x() * p2.y() - 2. * p1.x() * p2.y()) /
            6.;
  }
  return area;
}

double ComputeVFrac(const BezierList& cell, const Parabola& parabola) {
  return ComputeArea(cell, parabola) / IRL::safelyEpsilon(ComputeArea(cell));
}

double ComputeArea(const std::vector<BezierList>& cell) {
  if (cell.size() == 0) return 0.0;

  double area = 0.0;
  for (int i = 0; i < cell.size(); i++) {
    area += ComputeArea(cell[i]);
  }

  return area;
}

Moments ComputeMoments(const BezierList& cell) {
  if (cell.size() == 0) return Moments();

  double sum = 0.0, c = 0.0;
  auto moments = Moments();
  for (int i = 0; i < cell.size(); i++) {
    const auto p0 = cell[i].first;
    const auto p1 = cell[i].second;
    const auto p2 = cell[(i + 1) % cell.size()].first;

    const double m0_contrib =
        -(2. * p1.x() * p0.y() + p2.x() * p0.y() - 2. * p0.x() * p1.y() +
          2. * p2.x() * p1.y() - p0.x() * p2.y() - 2. * p1.x() * p2.y()) /
        6.;
    double y = m0_contrib - c;
    double t = sum + y;
    c = (t - sum) - y;
    sum = t;
    moments.m0() = sum;
    // moments.m0() +=
    //     (-(p2.x() * (p0.y() + 2. * p1.y())) - 2. * p1.x() * (p0.y() -
    //     p2.y())
    //     +
    //      p0.x() * (2. * p1.y() + p2.y())) /
    //     6.;
    moments.m1() += IRL2D::Vec(
        (-4. * p1.x() * p2.x() * (p0.y() + p1.y() - 2. * p2.y()) -
         4. * (p1.x() * p1.x()) * (p0.y() - p2.y()) +
         2. * p0.x() * p2.x() * (-p0.y() + p2.y()) +
         4. * p0.x() * p1.x() * (-2. * p0.y() + p1.y() + p2.y()) +
         p0.x() * p0.x() * (5. * p0.y() + 8. * p1.y() + 2. * p2.y()) -
         p2.x() * p2.x() * (2. * p0.y() + 8. * p1.y() + 5. * p2.y())) /
            60.,
        (-4. * p1.x() * (p0.y() - p2.y()) *
             (2. * p0.y() + p1.y() + 2. * p2.y()) -
         p2.x() * (2. * (p0.y() * p0.y()) + 4. * p0.y() * p1.y() +
                   4. * (p1.y() * p1.y()) + 2. * p0.y() * p2.y() +
                   8. * p1.y() * p2.y() - 5. * (p2.y() * p2.y())) +
         p0.x() * (-5. * (p0.y() * p0.y()) + 8. * p0.y() * p1.y() +
                   4. * (p1.y() * p1.y()) + 2. * p0.y() * p2.y() +
                   4. * p1.y() * p2.y() + 2. * (p2.y() * p2.y()))) /
            60.);
    moments.m2() += IRL2D::Mat(
        IRL2D::Vec(
            (-4. * (p1.x() * p1.x()) * p2.x() *
                 (3. * p0.y() + 2. * p1.y() - 5. * p2.y()) -
             10. * p1.x() * (p2.x() * p2.x()) *
                 (p0.y() + 2. * p1.y() - 3. * p2.y()) -
             8. * (p1.x() * p1.x() * p1.x()) * (p0.y() - p2.y()) +
             5. * (p0.x() * p0.x() * p0.x()) *
                 (21. * p0.y() + 6. * p1.y() + p2.y()) -
             5. * (p2.x() * p2.x() * p2.x()) *
                 (p0.y() + 6. * p1.y() + 21. * p2.y()) +
             p0.x() * p0.x() *
                 (10. * p1.x() * (-3. * p0.y() + 2. * p1.y() + p2.y()) +
                  p2.x() * (-5. * p0.y() + 2. * p1.y() + 3. * p2.y())) +
             p0.x() *
                 (-12. * p1.x() * p2.x() * (p0.y() - p2.y()) +
                  p2.x() * p2.x() * (-3. * p0.y() - 2. * p1.y() + 5. * p2.y()) +
                  p1.x() * p1.x() *
                      (-20. * p0.y() + 8. * p1.y() + 12. * p2.y()))) /
                420.,
            (-4. * (p1.x() * p1.x()) * (p0.y() - p2.y()) *
                 (5. * p0.y() + 4. * p1.y() + 5. * p2.y()) -
             2. * p0.x() * p2.x() * (p0.y() - p2.y()) *
                 (5. * p0.y() + 4. * p1.y() + 5. * p2.y()) -
             p2.x() * p2.x() *
                 (3. * (p0.y() * p0.y()) + 12. * p0.y() * p1.y() +
                  20. * (p1.y() * p1.y()) + 10. * p0.y() * p2.y() +
                  60. * p1.y() * p2.y() - 105. * (p2.y() * p2.y())) +
             4. * p0.x() * p1.x() *
                 (-15. * (p0.y() * p0.y()) + 4. * (p1.y() * p1.y()) +
                  2. * p0.y() * p2.y() + 6. * p1.y() * p2.y() +
                  3. * (p2.y() * p2.y())) -
             4. * p1.x() * p2.x() *
                 (3. * (p0.y() * p0.y()) + 4. * (p1.y() * p1.y()) -
                  15. * (p2.y() * p2.y()) +
                  2. * p0.y() * (3. * p1.y() + p2.y())) +
             p0.x() * p0.x() *
                 (-105. * (p0.y() * p0.y()) + 20. * (p1.y() * p1.y()) +
                  12. * p1.y() * p2.y() + 3. * (p2.y() * p2.y()) +
                  10. * p0.y() * (6. * p1.y() + p2.y()))) /
                840.),
        IRL2D::Vec((-4. * (p1.x() * p1.x()) * (p0.y() - p2.y()) *
                        (5. * p0.y() + 4. * p1.y() + 5. * p2.y()) -
                    2. * p0.x() * p2.x() * (p0.y() - p2.y()) *
                        (5. * p0.y() + 4. * p1.y() + 5. * p2.y()) -
                    p2.x() * p2.x() *
                        (3. * (p0.y() * p0.y()) + 12. * p0.y() * p1.y() +
                         20. * (p1.y() * p1.y()) + 10. * p0.y() * p2.y() +
                         60. * p1.y() * p2.y() - 105. * (p2.y() * p2.y())) +
                    4. * p0.x() * p1.x() *
                        (-15. * (p0.y() * p0.y()) + 4. * (p1.y() * p1.y()) +
                         2. * p0.y() * p2.y() + 6. * p1.y() * p2.y() +
                         3. * (p2.y() * p2.y())) -
                    4. * p1.x() * p2.x() *
                        (3. * (p0.y() * p0.y()) + 4. * (p1.y() * p1.y()) -
                         15. * (p2.y() * p2.y()) +
                         2. * p0.y() * (3. * p1.y() + p2.y())) +
                    p0.x() * p0.x() *
                        (-105. * (p0.y() * p0.y()) + 20. * (p1.y() * p1.y()) +
                         12. * p1.y() * p2.y() + 3. * (p2.y() * p2.y()) +
                         10. * p0.y() * (6. * p1.y() + p2.y()))) /
                       840.,
                   (-2. * p1.x() * (p0.y() - p2.y()) *
                        (15. * (p0.y() * p0.y()) + 4. * (p1.y() * p1.y()) +
                         10. * p1.y() * p2.y() + 15. * (p2.y() * p2.y()) +
                         2. * p0.y() * (5. * p1.y() + 8. * p2.y())) +
                    p0.x() * (-105. * (p0.y() * p0.y() * p0.y()) +
                              8. * (p1.y() * p1.y() * p1.y()) +
                              12. * (p1.y() * p1.y()) * p2.y() +
                              10. * p1.y() * (p2.y() * p2.y()) +
                              5. * (p2.y() * p2.y() * p2.y()) +
                              5. * (p0.y() * p0.y()) * (6. * p1.y() + p2.y()) +
                              p0.y() * (20. * (p1.y() * p1.y()) +
                                        12. * p1.y() * p2.y() +
                                        3. * (p2.y() * p2.y()))) -
                    p2.x() * (5. * (p0.y() * p0.y() * p0.y()) +
                              8. * (p1.y() * p1.y() * p1.y()) +
                              20. * (p1.y() * p1.y()) * p2.y() +
                              30. * p1.y() * (p2.y() * p2.y()) -
                              105. * (p2.y() * p2.y() * p2.y()) +
                              p0.y() * p0.y() * (10. * p1.y() + 3. * p2.y()) +
                              p0.y() * (12. * (p1.y() * p1.y()) +
                                        12. * p1.y() * p2.y() +
                                        5. * (p2.y() * p2.y())))) /
                       420.));
  }
  return moments;
}

Moments ComputeMoments(const BezierList& cell, const Parabola& parabola) {
  return ComputeMoments(ParabolaClip(cell, parabola));
}

Moments ComputeMoments(const BezierList& cell, const Vec& x0, const Vec& x1,
                       const Parabola& parabola) {
  return ComputeMoments(ClipByRectangleAndParabola(cell, x0, x1, parabola));
}

Vec RK4Point(const Vec& P, const double dt, const double time,
             const Vec (*vel)(const double t, const Vec& P)) {
  Vec Pnew = P;
  const int nsteps = 10;
  const double ddt = dt / static_cast<double>(nsteps);
  for (int i = 0; i < nsteps; i++) {
    const double t0 = time + static_cast<double>(i) * ddt;
    const auto k1 = vel(t0, Pnew);
    const auto P1 = Pnew + 0.5 * ddt * k1;
    const auto k2 = vel(t0 + 0.5 * ddt, P1);
    const auto P2 = Pnew + 0.5 * ddt * k2;
    const auto k3 = vel(t0 + 0.5 * ddt, P2);
    const auto P3 = Pnew + ddt * k3;
    const auto k4 = vel(t0 + ddt, P3);
    Pnew += ddt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
  }
  return Pnew;
}

std::pair<Vec, Vec> RK4PointAndTangent(
    const Vec& P, const Vec& T, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P)) {
  Vec Pnew = P;
  Vec Tnew = T;
  const int nsteps = 10;
  const double ddt = dt / static_cast<double>(nsteps);
  for (int i = 0; i < nsteps; i++) {
    const double t0 = time + static_cast<double>(i) * ddt;
    const auto k1 = vel(t0, Pnew);
    const auto gradu1 = grad_vel(t0, Pnew);
    const auto vortu1 = 0.0;  // gradu1[1][0] - gradu1[0][1];
    const auto q1 = gradu1 * Tnew + Vec(Tnew[1], -Tnew[0]) * vortu1;
    const auto P1 = Pnew + 0.5 * ddt * k1;
    const auto T1 = Tnew + 0.5 * ddt * q1;
    const auto k2 = vel(t0 + 0.5 * ddt, P1);
    const auto gradu2 = grad_vel(t0 + 0.5 * ddt, P1);
    const auto vortu2 = 0.0;  // gradu2[1][0] - gradu2[0][1];
    const auto q2 = gradu2 * T1 + Vec(T1[1], -T1[0]) * vortu2;
    const auto P2 = Pnew + 0.5 * ddt * k2;
    const auto T2 = Tnew + 0.5 * ddt * q2;
    const auto k3 = vel(t0 + 0.5 * ddt, P2);
    const auto gradu3 = grad_vel(t0 + 0.5 * ddt, P2);
    const auto vortu3 = 0.0;  // gradu3[1][0] - gradu3[0][1];
    const auto q3 = gradu3 * T2 + Vec(T2[1], -T2[0]) * vortu3;
    const auto P3 = Pnew + ddt * k3;
    const auto T3 = Tnew + ddt * q3;
    const auto k4 = vel(t0 + ddt, P3);
    const auto gradu4 = grad_vel(t0 + ddt, P3);
    const auto vortu4 = 0.0;  // gradu4[1][0] - gradu4[0][1];
    const auto q4 = gradu4 * T3 + Vec(T3[1], -T3[0]) * vortu4;
    Pnew += ddt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    Tnew += ddt * (q1 + 2.0 * q2 + 2.0 * q3 + q4) / 6.0;
  }
  return std::make_pair(Pnew, Tnew);
}

double Determinant(const Vec& T0, const Vec& T1) {
  return T0[0] * T1[1] - T0[1] * T1[0];
}

std::pair<bool, double> RayIntersection(const Vec& P0, const Vec& P1,
                                        const Vec& T0, const Vec& T1) {
  double tol = 1.0e2 * std::numeric_limits<double>::epsilon();
  auto edge = P1 - P0;
  const auto dP = edge;
  const double edge_mag = edge.magnitude();
  edge.normalize();
  const double det = Determinant(T1, T0);
  // If T0 and T1 are parallel
  if (std::abs(det) < tol) {
    const double det_2 = Determinant(T0, edge);
    // If T0 and P1-P0 are parallel
    if (std::abs(det_2) < tol) {
      return std::pair<bool, double>(true, 0.5 * dP.magnitude());
    }
    return std::pair<bool, double>(false, 0.0);
  }
  double t0 = Determinant(T1, dP) / det;
  double t1 = Determinant(T0, dP) / det;
  if (t0 > edge_mag || t1 > edge_mag) {
    return std::pair<bool, double>(false, 0.0);
  }
  return std::pair<bool, double>(t0 > 0.0 && t1 > 0.0, t0);
}

BezierList ConstructPathline(const Vec& P00, const double dt, const double t,
                             const Vec (*vel)(const double t, const Vec& P),
                             const int rec_num) {
  const int max_recursion = 3;
  auto P01 = RK4Point(P00, dt, t, vel);
  auto T00 = std::copysign(1.0, dt) * vel(t, P00);
  T00.normalize();
  auto T01 = -std::copysign(1.0, dt) * vel(t + dt, P01);
  T01.normalize();

  auto intersection = RayIntersection(P00, P01, T00, T01);
  if (intersection.first == false) {
    const double ddt = 0.5 * dt;
    auto Pmid = RK4Point(P00, ddt, t, vel);
    // If max level not reached, split arc
    if (rec_num < max_recursion) {
      auto list1 = ConstructPathline(P00, ddt, t, vel, rec_num + 1);
      auto list2 = ConstructPathline(Pmid, ddt, t + ddt, vel, rec_num + 1);
      // Append list2 to list1
      list1.insert(list1.end(), list2.begin(), list2.end());
      return list1;
      // Else return straight line
    } else {
      return std::vector{std::make_pair(P00, Pmid)};
    }
  } else {
    return std::vector{std::make_pair(P00, P00 + intersection.second * T00)};
  }
}

BezierList TransportEdgeMidPoint(
    const Vec& P00, const Vec& P10, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P), const int rec_num,
    const bool add_pathlines, const bool close_flux, const bool correct_area,
    const double exact_area) {
  double tol = 10. * std::numeric_limits<double>::epsilon();
  const int max_recursion = 1;
  BezierList list(0);
  list.reserve(4);
  const auto P01 = RK4Point(P00, dt, time, vel);
  const auto P11 = RK4Point(P10, dt, time, vel);
  const auto Phalf0 = 0.5 * (P00 + P10);
  const auto Phalf1 = RK4Point(Phalf0, dt, time, vel);

  // Add first pathline
  if (add_pathlines == true || correct_area == true) {
    const auto P0half = RK4Point(P00, 0.5 * dt, time, vel);
    const auto P0ctrl = 2.0 * P0half - 0.5 * (P00 + P01);
    list.push_back(std::make_pair(P00, P0ctrl));
  }

  const auto Pctrl1 = 2.0 * Phalf1 - 0.5 * (P01 + P11);
  list.push_back(std::make_pair(P01, Pctrl1));

  // Add returning pathline to list
  if (add_pathlines == true || correct_area == true) {
    const auto P1half = RK4Point(P10, 0.5 * dt, time, vel);
    const auto P1ctrl = 2.0 * P1half - 0.5 * (P10 + P11);
    list.push_back(std::make_pair(P11, P1ctrl));
  }

  // Add original edge to list
  if (close_flux == true || correct_area == true) {
    list.push_back(std::make_pair(P10, 0.5 * (P00 + P10)));
  }

  // Correct control points of egde to match exact area
  if (correct_area == true) {
    const double uncorrected_area = ComputeArea(list);
    const double area_correction = exact_area - uncorrected_area;

    const auto P0 = P01;
    const auto P1 = Pctrl1;
    const auto P2 = P11;
    auto segment = P2 - P0;
    const double length = segment.magnitude();
    segment.normalize();
    const auto normal = IRL2D::Vec(-segment.y(), segment.x());
    const double desired_arc_area = ArcVolume(P0, P1, P2) + area_correction;
    // Find correction dir
    const auto Pmid = P1 - (normal * (P1 - P0)) * normal;
    // Calculate new distance correction
    const double s =
        (3. * desired_arc_area + 2. * P2.x() * P0.y() + Pmid.x() * P0.y() -
         2. * P0.x() * P2.y() - Pmid.x() * P2.y() - P0.x() * Pmid.y() +
         P2.x() * Pmid.y()) /
        IRL::safelyEpsilon(normal.y() * P0.x() - normal.y() * P2.x() -
                           normal.x() * P0.y() + normal.x() * P2.y());
    // Correct control point
    list[1].second = Pmid + s * normal;

    if (close_flux == false) {
      list.pop_back();
    }
    if (add_pathlines == false) {
      list.erase(list.begin() + 2, list.end());
      list.erase(list.begin(), list.begin() + 1);
    }
  }

  return list;
}

BezierList TransportEdge(const Vec& P00, const Vec& P10, const double dt,
                         const double time,
                         const Vec (*vel)(const double t, const Vec& P),
                         const Mat (*grad_vel)(const double t, const Vec& P),
                         const int rec_num, const bool add_pathlines,
                         const bool close_flux, const bool correct_area,
                         const double exact_area) {
  double tol = 10. * std::numeric_limits<double>::epsilon();
  const int max_recursion = 3;
  BezierList list(0);
  list.reserve(5);
  // Add first pathline
  if (add_pathlines == true || correct_area == true) {
    auto pathline = ConstructPathline(P00, dt, time, vel, 0);
    list.insert(list.end(), pathline.begin(), pathline.end());
  }

  // Transport both points and tangents
  auto T00 = P10 - P00;
  T00.normalize();
  const auto T10 = -T00;
  auto PT01back = RK4PointAndTangent(P00, T00, dt, time, vel, grad_vel);
  auto PT11back = RK4PointAndTangent(P10, T10, dt, time, vel, grad_vel);
  const auto P01 = PT01back.first;
  const auto P11 = PT11back.first;
  auto T01 = PT01back.second;
  T01.normalize();
  auto T11 = PT11back.second;
  T11.normalize();

  // Check for existance of control point; recursive slip if needed
  auto intersection = RayIntersection(P01, P11, T01, T11);
  const int start_edge = list.size();
  int end_edge = -1;
  if (intersection.first == false) {
    auto Pmid = 0.5 * (P00 + P10);
    // If max level not reached, split arc
    if (rec_num < max_recursion) {
      const auto list1 =
          TransportEdge(P00, Pmid, dt, time, vel, grad_vel, rec_num + 1);
      list.insert(list.end(), list1.begin(), list1.end());
      const auto list2 =
          TransportEdge(Pmid, P10, dt, time, vel, grad_vel, rec_num + 1);
      list.insert(list.end(), list2.begin(), list2.end());
    } else {
      list.push_back(std::make_pair(P01, 0.5 * (P01 + P11)));
    }
  } else {
    list.push_back(std::make_pair(P01, P01 + intersection.second * T01));
  }

  // If area correction required, store last id of edge
  if (correct_area == true) {
    end_edge = list.size();
  }

  // Add returning pathline to list
  if (add_pathlines == true || correct_area == true) {
    auto pathline = ConstructPathline(P10, dt, time, vel, 0);
    list.push_back(std::make_pair(P11, pathline.back().second));
    for (int i = pathline.size() - 1; i >= 1; i--) {
      list.push_back(std::make_pair(pathline[i].first, pathline[i - 1].second));
    }
  }

  // Add original edge to list
  if (close_flux == true || correct_area == true) {
    list.push_back(std::make_pair(P10, 0.5 * (P00 + P10)));
  }

  // Correct control points of egde to match exact area
  if (correct_area == true) {
    const int narcs = end_edge - start_edge;
    const double uncorrected_area = ComputeArea(list);
    const double area_correction = exact_area - uncorrected_area;

    // Compute length of polygon arc
    double total_length = 0.0;
    for (int i = start_edge; i < end_edge; i++) {
      const auto segment = list[i + 1].first - list[i].first;
      total_length += segment.magnitude();
    }

    // Correct each arc
    for (int i = start_edge; i < end_edge; i++) {
      const auto P0 = list[i].first;
      const auto P1 = list[i].second;
      const auto P2 = list[i + 1].first;
      auto segment = P2 - P0;
      const double length = segment.magnitude();
      segment.normalize();
      const auto normal = IRL2D::Vec(-segment.y(), segment.x());
      const double weight = length / total_length;
      const double desired_arc_area =
          ArcVolume(P0, P1, P2) + weight * area_correction;
      // Find correction dir
      const auto Pmid = P1 - (normal * (P1 - P0)) * normal;
      // Calculate new distance correction
      const double s =
          (3. * desired_arc_area + 2. * P2.x() * P0.y() + Pmid.x() * P0.y() -
           2. * P0.x() * P2.y() - Pmid.x() * P2.y() - P0.x() * Pmid.y() +
           P2.x() * Pmid.y()) /
          IRL::safelyEpsilon(normal.y() * P0.x() - normal.y() * P2.x() -
                             normal.x() * P0.y() + normal.x() * P2.y());
      // Correct control point
      list[i].second = Pmid + s * normal;
    }

    if (close_flux == false) {
      list.pop_back();
    }
    if (add_pathlines == false) {
      list.erase(list.begin() + end_edge, list.end());
      list.erase(list.begin(), list.begin() + start_edge);
    }
  }

  return list;
}

BezierList TransportLinearEdge(
    const Vec& P00, const Vec& P10, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P), const int rec_num,
    const bool add_pathlines, const bool close_flux, const bool correct_area,
    const double exact_area) {
  BezierList list(0);
  list.reserve(5);

  const auto P01 = RK4Point(P00, dt, time, vel);
  const auto P11 = RK4Point(P10, dt, time, vel);
  const auto PM0 = 0.5 * (P00 + P10);
  const auto PM1 = RK4Point(PM0, dt, time, vel);

  // Add first pathline
  if (add_pathlines == true || correct_area == true) {
    list.push_back(std::make_pair(P00, 0.5 * (P01 + P00)));
  }

  list.push_back(std::make_pair(P01, 0.5 * (P01 + PM1)));
  list.push_back(std::make_pair(PM1, 0.5 * (PM1 + P11)));

  // Add returning pathline to list
  if (add_pathlines == true || correct_area == true) {
    list.push_back(std::make_pair(P11, 0.5 * (P10 + P11)));
  }

  // Add original edge to list
  if (close_flux == true || correct_area == true) {
    list.push_back(std::make_pair(P10, 0.5 * (P00 + P10)));
  }

  // Correct control points of egde to match exact area
  if (correct_area == true) {
    const double uncorrected_area = ComputeArea(list);
    const double area_correction = exact_area - uncorrected_area;

    // Correct each arc
    const auto P0 = P01;
    const auto P1 = PM1;
    const auto P2 = P11;
    const double desired_arc_area =
        0.5 * (-(P1.x() * P0.y()) + P2.x() * P0.y() + P0.x() * P1.y() -
               P2.x() * P1.y() - P0.x() * P2.y() + P1.x() * P2.y()) +
        area_correction;
    auto segment = P2 - P0;
    segment.normalize();
    const auto normal = IRL2D::Vec(segment.y(), -segment.x());
    // Calculate new distance correction
    const double s =
        (2. * desired_arc_area + P1.x() * P0.y() - P2.x() * P0.y() -
         P0.x() * P1.y() + P2.x() * P1.y() + P0.x() * P2.y() -
         P1.x() * P2.y()) /
        IRL::safelyEpsilon(normal.y() * P0.x() - normal.y() * P2.x() -
                           normal.x() * P0.y() + normal.x() * P2.y());
    // Correct control point
    list[2].first = P1 + s * normal;
    list[1].second = 0.5 * (list[1].first + list[2].first);
    list[2].second = 0.5 * (list[2].first + list[3].first);

    if (close_flux == false) {
      list.pop_back();
    }
    if (add_pathlines == false) {
      list.erase(list.begin() + 3, list.end());
      list.erase(list.begin(), list.begin() + 1);
    }
  }

  return list;
}

BezierList CreateFluxCell(const Vec& P00, const Vec& P10, const double dt,
                          const double time,
                          const Vec (*vel)(const double t, const Vec& P),
                          const Mat (*grad_vel)(const double t, const Vec& P),
                          const bool correct_area, const double exact_area) {
  return TransportEdge(P00, P10, dt, time, vel, grad_vel, 0, true, true,
                       correct_area, exact_area);
}

BezierList CreateFluxCellMidPoint(
    const Vec& P00, const Vec& P10, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P),
    const bool correct_area, const double exact_area) {
  return TransportEdgeMidPoint(P00, P10, dt, time, vel, grad_vel, 0, true, true,
                               correct_area, exact_area);
}

BezierList CreatePreImage(const Vec& X0, const Vec& X1, const double dt,
                          const double time,
                          const Vec (*vel)(const double t, const Vec& P),
                          const Mat (*grad_vel)(const double t, const Vec& P),
                          const bool correct_area,
                          const std::array<double, 4>& exact_area) {
  const auto PList =
      std::array<Vec, 4>{Vec(X0.x(), X0.y()), Vec(X1.x(), X0.y()),
                         Vec(X1.x(), X1.y()), Vec(X0.x(), X1.y())};
  BezierList preimage(0);
  preimage.reserve(4);
  for (int i = 0; i < 4; i++) {
    auto tmp_edge =
        TransportEdge(PList[i], PList[(i + 1) % 4], dt, time, vel, grad_vel, 0,
                      false, false, correct_area, exact_area[i]);
    preimage.insert(preimage.end(), tmp_edge.begin(), tmp_edge.end());
  }
  return preimage;
}

BezierList CreatePreImageMidPoint(
    const Vec& X0, const Vec& X1, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P),
    const bool correct_area, const std::array<double, 4>& exact_area) {
  const auto PList =
      std::array<Vec, 4>{Vec(X0.x(), X0.y()), Vec(X1.x(), X0.y()),
                         Vec(X1.x(), X1.y()), Vec(X0.x(), X1.y())};
  BezierList preimage(0);
  preimage.reserve(4);
  for (int i = 0; i < 4; i++) {
    auto tmp_edge = TransportEdgeMidPoint(PList[i], PList[(i + 1) % 4], dt,
                                          time, vel, grad_vel, 0, false, false,
                                          correct_area, exact_area[i]);
    preimage.insert(preimage.end(), tmp_edge.begin(), tmp_edge.end());
  }
  return preimage;
}

BezierList CreateLinearPreImage(
    const Vec& X0, const Vec& X1, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P),
    const bool correct_area, const std::array<double, 4>& exact_area) {
  const auto PList =
      std::array<Vec, 4>{Vec(X0.x(), X0.y()), Vec(X1.x(), X0.y()),
                         Vec(X1.x(), X1.y()), Vec(X0.x(), X1.y())};
  BezierList preimage(0);
  preimage.reserve(4);
  for (int i = 0; i < 4; i++) {
    auto tmp_edge = TransportLinearEdge(PList[i], PList[(i + 1) % 4], dt, time,
                                        vel, grad_vel, 0, false, false,
                                        correct_area, exact_area[i]);
    preimage.insert(preimage.end(), tmp_edge.begin(), tmp_edge.end());
  }
  return preimage;
}

BezierList CreateLinearFluxCell(
    const Vec& P00, const Vec& P10, const double dt, const double time,
    const Vec (*vel)(const double t, const Vec& P),
    const Mat (*grad_vel)(const double t, const Vec& P),
    const bool correct_area, const double exact_area) {
  return TransportLinearEdge(P00, P10, dt, time, vel, grad_vel, 0, true, true,
                             correct_area, exact_area);
}

template <class ScalarType>
const bool AddIntersectionsToCell(const BezierList& original_cell,
                                  const Parabola& parabola, const int it,
                                  BezierList* cell_with_intersections,
                                  int* nintersections,
                                  std::vector<int>* vertex_nature) {
  const int nvert_init = original_cell.size();
  auto tol = std::numeric_limits<ScalarType>::epsilon();
  if constexpr (std::is_same_v<ScalarType, double>) {
    tol *= 1.0e4;
  } else {
    tol *= 1.0e8;
  }
  bool parabola_contains_vertex = false;

  // Resize list to be empty
  (*cell_with_intersections).resize(0);
  (*vertex_nature).resize(0);
  (*nintersections) = 0;

  ///// If necessary, nudge cell
  auto nudge = Vec(0.0, 0.0);
  if (it > 1) {
    // Create a random number generator and seed it with
    // the number of prior nudge iterations (for reproductibility)
    std::random_device rd;
    std::mt19937 gen(it);
    std::uniform_real_distribution distr(-1.0, 1.0);

    // This is a-hoc but works well
    double nudge_epsilon = 1.0e-1 * pow(static_cast<double>(it), 1.05) *
                           std::numeric_limits<double>::epsilon();

    // Compute polytope center to rotate about it
    const auto cell_moments = ComputeMoments(original_cell);
    const auto center =
        cell_moments.m1() / IRL::safelyEpsilon(cell_moments.m0());
    nudge = nudge_epsilon * Vec(distr(gen), distr(gen));
  }

  // Compute all intersections and insert them in cell
  for (int i = 0; i < nvert_init; i++) {
    // Get points of current bezier arc
    const auto p0 = original_cell[i].first;
    const auto p1 = original_cell[i].second;
    const auto p2 = original_cell[(i + 1) % original_cell.size()].first;
    // Add p0 and p1 to new list
    (*cell_with_intersections).push_back(PtAndControl{p0, p1});
    // Is p0 below or above parabola?
    const auto dist = DistanceToParabola<ScalarType>(parabola, p0 + nudge);
    if (fabs(dist) < tol) {
      parabola_contains_vertex = true;
      break;
    }
    (*vertex_nature).push_back(dist < 0. ? 1 : 0);
    // Compute t-values of potential intersections
    const auto t_vals = AnalyticIntersections<ScalarType>(
        parabola, p0 + nudge, p1 + nudge, p2 + nudge);
    // Tmp variables needs for arc splitting
    double t0 = 0.;
    Vec p0new = p0, p1new = p1;
    // Loop over all detected intersections
    for (int j = 0; j < t_vals.size(); j++) {
      // If t = 0 or t = 1, we must nudge! Go back to beginning of function
      if (fabs(t_vals[j]) < tol || fabs(1. - t_vals[j]) < tol) {
        parabola_contains_vertex = true;
        break;
      }
      const double t = static_cast<double>(t_vals[j]);
      //// Spilling of the arc (using de Casteljau algorithm)
      // First: Update previous control point
      // if (std::fabs(1. - t0) >
      //     100. * std::numeric_limits<ScalarType>::epsilon()) {
      (*cell_with_intersections).back().second =
          p0new + (t - t0) / (1. - t0) * (p1new - p0new) - nudge;
      // Second: Add intersection and new control point
      p0new = BezierPoint(p0, p1, p2, t) - nudge;
      p1new = p1new + (t - t0) / (1. - t0) * (p2 - p1new) - nudge;
      // } else {
      //   (*cell_with_intersections).back().second =
      //       p0new + 0.5 * (p1new - p0new) - nudge;
      //   p0new = BezierPoint(p0, p1, p2, t) - nudge;
      //   p1new = p1new + 0.5 * (p2 - p1new) - nudge;
      // }
      (*cell_with_intersections).push_back(PtAndControl{p0new, p1new});
      (*nintersections)++;
      // Check if intersection is entry of exit
      // If: previous vertex is below or exit, then current vertex is entry
      // Otherwise: the current vertex is exit
      if ((*vertex_nature).back() == 1 || (*vertex_nature).back() == 2) {
        (*vertex_nature).push_back(3);
      } else {
        (*vertex_nature).push_back(2);
      }
      t0 = t;
    }
    if (parabola_contains_vertex) break;
  }

  return parabola_contains_vertex;
}

BezierList ParabolaClip(const BezierList& original_cell,
                        const Parabola& parabola,
                        const bool return_parabola_only) {
  // If empty cell, return empty cell
  if (original_cell.size() == 0 || parabola.isAlwaysBelow()) {
    return BezierList(0);
  }
  if (parabola.isAlwaysAbove()) {
    if (return_parabola_only) {
      return BezierList(0);
    } else {
      return original_cell;
    }
  }

  // Specify constants
  const int nvert_init = original_cell.size();
  const double tol = 100. * std::numeric_limits<double>::epsilon();
  bool parabola_contains_vertex = false;
  const int itmax = 500;
  double max_dist = 0.0;
  for (int i = 0; i < nvert_init; i++) {
    const auto p0 = original_cell[i].first;
    const auto p2 = original_cell[(i + 1) % nvert_init].first;
    const auto segment = p2 - p0;
    max_dist = std::fmax(max_dist, segment.magnitude());
  }

  // Define scale so that the polygon's bounding box is O(1)
  const double length_scale = 1.0;  // std::fmax(1.0e6 * tol, max_dist);
  const double inv_ls = 1.0 / length_scale;

  // Store parabola properties
  const auto datum = parabola.datum();
  auto scaled_datum = inv_ls * datum;
  ReferenceFrame frame = parabola.frame();
  frame[0].normalize();
  frame[1][0] = -frame[0][1];
  frame[1][1] = frame[0][0];
  const double scaled_coeff = length_scale * parabola.coeff();
  const auto scaled_parabola = Parabola(scaled_datum, frame, scaled_coeff);

  // Copy cell in frame of reference of parabola
  BezierList cell(nvert_init);
  for (int i = 0; i < nvert_init; i++) {
    cell[i].first = frame * (inv_ls * (original_cell[i].first - datum));
    cell[i].second = frame * (inv_ls * (original_cell[i].second - datum));
  }

  // Compute bounding box and do simple cheap tests
  const auto bbox = BoundingBox(cell);
  const auto xmin = bbox.first.x(), ymin = bbox.first.y();
  const auto xmax = bbox.second.x(), ymax = bbox.second.y();
  const auto parabola_at_xmin = -scaled_coeff * xmin * xmin;
  const auto parabola_at_xmax = -scaled_coeff * xmax * xmax;
  bool is_entirely_below = false, is_entirely_above = false;
  if (scaled_coeff >= 0.0) {
    if (ymax <= parabola_at_xmin && ymax <= parabola_at_xmax) {
      is_entirely_below = true;
    } else if ((xmax <= 0.0 && ymin >= parabola_at_xmax) ||
               (xmin >= 0.0 && ymin >= parabola_at_xmin)) {
      is_entirely_above = true;
    }
  } else {
    if (ymin >= parabola_at_xmin && ymin >= parabola_at_xmax) {
      is_entirely_above = true;
    } else if ((xmax <= 0.0 && ymax <= parabola_at_xmax) ||
               (xmin >= 0.0 && ymax <= parabola_at_xmin)) {
      is_entirely_below = true;
    }
  }
  if (is_entirely_above) {
    return BezierList(0);
  }
  if (is_entirely_below) {
    if (return_parabola_only) {
      return BezierList(0);
    } else {
      return original_cell;
    }
  }

  // Initialize clipped cell
  BezierList clipped_cell(0), clipped_parabola(0);
  if (return_parabola_only) {
    clipped_parabola.reserve(5 * nvert_init);
  } else {
    clipped_cell.reserve(5 * nvert_init);
  }

  // Reserve 5 times original size (i.e., 4 intersections per arc)
  BezierList cell_with_intersections(0);
  cell_with_intersections.reserve(5 * nvert_init);
  std::vector<int> vertex_nature(0);
  int nintersections = 0;

  // This while loop is needed in case we might want to nudge the parabola
  int it = 0;
  while (it < itmax) {
    it++;
    if (it <= 1) {
      parabola_contains_vertex = AddIntersectionsToCell<double>(
          cell, scaled_parabola, it, &cell_with_intersections, &nintersections,
          &vertex_nature);
    } else {
      parabola_contains_vertex = AddIntersectionsToCell<float128>(
          cell, scaled_parabola, it, &cell_with_intersections, &nintersections,
          &vertex_nature);
    }

    // If non-even # of intersections or parabola X with vertex, go back to
    // start of loop and nudge!
    if (nintersections % 2 != 0 || parabola_contains_vertex) {
      continue;
    }

    // If no intersections at all; leave loop
    if (nintersections == 0) {
      if (return_parabola_only) {
        return BezierList(0);
      } else {
        if (vertex_nature[0] == 1) {
          return original_cell;
        }
        return BezierList(0);
      }
    }

    ///////// In-house clipping algorithm
    // Create new bezier list
    // First, find entry
    int start = -1;
    const int nvertices = cell_with_intersections.size();
    for (int i = 0; i < nvertices; i++) {
      if (vertex_nature[i] == 2) {
        start = i;
        break;
      }
    }
    assert(start >= 0);
    // Add entry to new bezier list
    clipped_cell.push_back(cell_with_intersections[start]);
    // Store start id
    const int check_start = start;
    // Loop over cell until next exit or start point itself;
    for (int i = 0; i < nvertices; i++) {
      // Fint next vertex on cell
      const int next_id = (start + 1 + i) % nvertices;
      // If next vertex is entry, modify previous control point
      if (vertex_nature[next_id] == 2) {
        // Find next entry in cell
        // Create new control point for bezier arc extracted from parabola
        const auto p0 = clipped_cell.back().first;
        const auto p2 = cell_with_intersections[next_id].first;
        const auto p1 = Vec(0.5 * (p0.x() + p2.x()),
                            p0.y() + scaled_coeff * (p0.x() - p2.x()) * p0.x());
        if (return_parabola_only) {
          clipped_parabola.push_back(std::make_pair(p0, p1));
          clipped_parabola.push_back(std::make_pair(p2, p1));
        } else {
          clipped_cell.back().second = p1;
        }
      }
      // If next vertex is start, we are done
      if (next_id == check_start) {
        break;
      }
      // If below or intersection, add next vertex to list
      if (vertex_nature[next_id] != 0) {
        clipped_cell.push_back(cell_with_intersections[next_id]);
      }
    }
    break;
  }

  if (parabola_contains_vertex) {
    std::cout << "WARNING: Parabola contains vertex!" << itmax
              << " nudges were not enough." << std::endl;
    std::cout << "  -> Parabola original = " << parabola << std::endl;
    std::cout << "  -> Parabola scaled   = " << scaled_parabola << std::endl;
    std::cout << "  -> Length scale      = " << length_scale << std::endl;
    ToVTK(original_cell, "cell_with_vertex_contained_in_parabola");
  } else if (it == itmax) {
    std::cout << "WARNING: Could not find an even number of intersections"
              << std::endl;
    ToVTK(original_cell, "cell_with_odd_intersections");
    std::cout << "  -> Parabola original = " << parabola << std::endl;
    std::cout << "  -> Parabola scaled   = " << scaled_parabola << std::endl;
    std::cout << "  -> Length scale      = " << length_scale << std::endl;
  }

  // Move back to canonical frame of reference and return
  if (return_parabola_only) {
    for (int i = 0; i < clipped_parabola.size(); i++) {
      clipped_parabola[i].first =
          datum +
          frame.transpose() * (length_scale * clipped_parabola[i].first);
      clipped_parabola[i].second =
          datum +
          frame.transpose() * (length_scale * clipped_parabola[i].second);
    }
    return clipped_parabola;
  } else {
    for (int i = 0; i < clipped_cell.size(); i++) {
      clipped_cell[i].first =
          datum + frame.transpose() * (length_scale * clipped_cell[i].first);
      clipped_cell[i].second =
          datum + frame.transpose() * (length_scale * clipped_cell[i].second);
    }
    return clipped_cell;
  }
}

BezierList ClipByRectangleAndParabola(const BezierList& original_cell,
                                      const Vec& x0, const Vec& x1,
                                      const Parabola& parabola) {
  const auto datum = parabola.datum();
  auto clipped_cell = ParabolaClip(original_cell, parabola);

  std::array<Parabola, 4> localizers{
      Parabola(Vec(x1.x(), x1.y()),
               ReferenceFrame(Vec(1.0, 0.0), Vec(0.0, 1.0)), 0.),
      Parabola(Vec(x0.x(), x0.y()),
               ReferenceFrame(Vec(0.0, 1.0), Vec(-1.0, 0.0)), 0.),
      Parabola(Vec(x0.x(), x0.y()),
               ReferenceFrame(Vec(-1.0, 0.0), Vec(0.0, -1.0)), 0.),
      Parabola(Vec(x1.x(), x1.y()),
               ReferenceFrame(Vec(0.0, -1.0), Vec(1.0, 0.0)), 0.)};
  for (int i = 0; i < 4; i++) {
    clipped_cell = ParabolaClip(clipped_cell, localizers[i]);
  }
  return clipped_cell;
}

double IntegrateFlux(const Vec& P0, const Vec& P1, const double dt,
                     const double time,
                     const Vec (*vel)(const double t, const Vec& P)) {
  const auto points = IRL::AbscissaeGauss<double, 10>();
  const auto weights = IRL::WeightsGauss<double, 10>();
  const auto edge = P1 - P0;
  auto normal = Vec(-edge.y(), edge.x());
  normal.normalize();
  double udotn = 0.0;
  for (int j = 0; j < 10; j++) {
    for (int i = 0; i < 10; i++) {
      udotn += weights[i] * weights[j] *
               (vel(time + 0.5 * (1.0 + points[j]) * dt,
                    P0 + 0.5 * (1.0 + points[i]) * edge) *
                normal);
    }
  }
  return 0.25 * dt * edge.magnitude() * udotn;
}

Parabola MatchToVolumeFraction(const BezierList& cell, const Parabola& parabola,
                               const double vfrac) {
  return MatchToVolumeFractionBisection(cell, parabola, vfrac);
}

Parabola MatchToVolumeFractionBrent(const BezierList& cell,
                                    const Parabola& parabola,
                                    const double vfrac) {
  const IRL::UnsignedIndex_t max_brent_iter = 50;
  const double tol = std::numeric_limits<double>::epsilon();
  const double vfrac_tolerance = 1.0e-14;

  // Calculate volume of cell
  const auto cell_area = ComputeArea(cell);
  const double cell_area_inv = 1.0 / cell_area;
  const double length_scale = std::sqrt(cell_area);

  // Move cell to local frame of reference of the paraboloid
  const auto datum = parabola.datum();
  const auto frame = parabola.frame();
  const auto coeff = parabola.coeff();

  double interval_min = -0.5;
  double interval_max = 0.5;

  Parabola moving_parabola(datum, frame, coeff);
  moving_parabola.datum() = datum + frame[1] * interval_max * length_scale;
  double vfrac_max = ComputeArea(cell, moving_parabola) * cell_area_inv;
  IRL::UnsignedIndex_t iter = 0;
  while (iter < 40 && vfrac_max < vfrac) {
    interval_max *= 2.0;
    moving_parabola.datum() = datum + frame[1] * interval_max * length_scale;
    vfrac_max = ComputeArea(cell, moving_parabola) * cell_area_inv;
    iter++;
  }

  moving_parabola.datum() = datum + frame[1] * interval_min * length_scale;
  double vfrac_min = ComputeArea(cell, moving_parabola) * cell_area_inv;
  iter = 0;
  while (iter < 40 && vfrac_min > vfrac) {
    interval_min *= 2.0;
    moving_parabola.datum() = datum + frame[1] * interval_min * length_scale;
    vfrac_min = ComputeArea(cell, moving_parabola) * cell_area_inv;
    iter++;
  }

  // Brent algorithm
  double a = interval_min, b = interval_max;
  double fa = 0.0 - vfrac, fb = 1.0 - vfrac;
  if (fa * fb > 0.0) {
    std::cout << "ERROR: Brent root-finding failure: f(a)*f(b) > 0"
              << std::endl;
    exit(-1);
  }
  if (std::fabs(fa) < std::fabs(fb)) {
    std::swap(a, b);
    std::swap(fa, fb);
  }
  double c = a, fc = fa, s = DBL_MAX, d = DBL_MAX;
  double fs = DBL_MAX;
  bool mflag = true;
  for (IRL::UnsignedIndex_t iter = 0; iter < max_brent_iter; ++iter) {
    if (std::fabs(fs) < vfrac_tolerance) {
      return moving_parabola;
    }
    if (fa != fc && fb != fc) {
      s = a * fb * fc / ((fa - fb) * (fa - fc)) +
          b * fa * fc / ((fb - fa) * (fb - fc)) +
          c * fa * fb / ((fc - fa) * (fc - fb));
    } else {
      s = b - fb * (b - a) / (fb - fa);
    }
    if (s < 0.25 * (3.0 * a + b) || s > b ||
        (mflag && std::fabs(s - b) >= 0.5 * std::fabs(b - c)) ||
        (!mflag && std::fabs(s - b) >= 0.5 * std::fabs(c - d)) ||
        (mflag && std::fabs(b - c) < tol) ||
        (!mflag && std::fabs(c - d) < tol)) {
      s = 0.5 * (a + b);
      mflag = true;
    } else {
      mflag = false;
    }
    moving_parabola.datum() = datum + frame[1] * s * length_scale;
    fs = ComputeArea(cell, moving_parabola) * cell_area_inv - vfrac;
    d = c;
    c = b;
    fc = fb;
    if (fa * fs < 0.0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if (std::fabs(fa) < std::fabs(fb)) {
      std::swap(a, b);
      std::swap(fa, fb);
    }
  }

  return moving_parabola;
}

Parabola MatchToVolumeFractionIllinois(const BezierList& cell,
                                       const Parabola& parabola,
                                       const double vfrac) {
  const IRL::UnsignedIndex_t max_illinois_iter = 100;
  const double vfrac_tolerance = 1.0e-14;

  // Calculate volume of cell
  const auto cell_area = ComputeArea(cell);
  const double cell_area_inv = 1.0 / cell_area;
  const double length_scale = std::sqrt(cell_area);

  // Move cell to local frame of reference of the paraboloid
  const auto datum = parabola.datum();
  const auto frame = parabola.frame();
  const auto coeff = parabola.coeff();

  double interval_min = -1.0;
  double interval_max = 1.0;

  Parabola moving_parabola(datum, frame, coeff);
  moving_parabola.datum() = datum + frame[1] * interval_max * length_scale;
  double vfrac_max = ComputeArea(cell, moving_parabola) * cell_area_inv;
  IRL::UnsignedIndex_t iter = 0;
  while (iter < 40 && vfrac_max < vfrac) {
    interval_max *= 2.0;
    moving_parabola.datum() = datum + frame[1] * interval_max * length_scale;
    vfrac_max = ComputeArea(cell, moving_parabola) * cell_area_inv;
    iter++;
  }

  moving_parabola.datum() = datum + frame[1] * interval_min * length_scale;
  double vfrac_min = ComputeArea(cell, moving_parabola) * cell_area_inv;
  iter = 0;
  while (iter < 40 && vfrac_min > vfrac) {
    interval_min *= 2.0;
    moving_parabola.datum() = datum + frame[1] * interval_min * length_scale;
    vfrac_min = ComputeArea(cell, moving_parabola) * cell_area_inv;
    iter++;
  }

  // Illinois algorithm
  double a = interval_min, b = interval_max;
  double fa = 0.0 - vfrac, fb = 1.0 - vfrac;
  if (fa * fb > 0.0) {
    std::cout << "ERROR: Illinois root-finding failure: f(a)*f(b) > 0"
              << std::endl;
    exit(-1);
  }
  for (IRL::UnsignedIndex_t iter = 0; iter < max_illinois_iter; ++iter) {
    double c = b - fb * (b - a) / (fb - fa);
    moving_parabola.datum() = datum + frame[1] * c * length_scale;
    const double fc =
        ComputeArea(cell, moving_parabola) * cell_area_inv - vfrac;
    if (std::fabs(fc) < vfrac_tolerance) return moving_parabola;
    if (fb * fc < 0.0) {
      a = b;
      fa = fb;
    } else {
      fa = 0.5 * fa;
    }
    b = c;
    fb = fc;
  }

  return moving_parabola;
}

Parabola MatchToVolumeFractionBisection(const BezierList& cell,
                                        const Parabola& parabola,
                                        const double vfrac,
                                        const int max_bisection_iter) {
  const double vfrac_tolerance = 1.0e-15;

  // Calculate volume of cell
  const auto cell_area = ComputeArea(cell);
  const double cell_area_inv = 1.0 / cell_area;
  const double length_scale = std::sqrt(cell_area);

  // Move cell to local frame of reference of the paraboloid
  const auto datum = parabola.datum();
  const auto frame = parabola.frame();
  const auto coeff = parabola.coeff();

  double interval_min = -0.5;
  double interval_max = 0.5;

  Parabola moving_parabola(datum, frame, coeff);
  moving_parabola.datum() = datum + frame[1] * interval_max * length_scale;
  double vfrac_max = ComputeArea(cell, moving_parabola) * cell_area_inv;
  IRL::UnsignedIndex_t iter = 0;
  while (iter < 40 && vfrac_max < vfrac) {
    interval_max *= 2.0;
    moving_parabola.datum() = datum + frame[1] * interval_max * length_scale;
    vfrac_max = ComputeArea(cell, moving_parabola) * cell_area_inv;
    iter++;
  }

  moving_parabola.datum() = datum + frame[1] * interval_min * length_scale;
  double vfrac_min = ComputeArea(cell, moving_parabola) * cell_area_inv;
  iter = 0;
  while (iter < 40 && vfrac_min > vfrac) {
    interval_min *= 2.0;
    moving_parabola.datum() = datum + frame[1] * interval_min * length_scale;
    vfrac_min = ComputeArea(cell, moving_parabola) * cell_area_inv;
    iter++;
  }

  // Bisection algorithm
  double a = interval_min;
  double b = interval_max;
  for (IRL::UnsignedIndex_t iter = 0; iter < max_bisection_iter; ++iter) {
    double c = 0.5 * (a + b);
    moving_parabola.datum() = datum + frame[1] * c * length_scale;
    const double vfrac_cut = ComputeArea(cell, moving_parabola) * cell_area_inv;
    if (vfrac_cut > vfrac + vfrac_tolerance) {
      b = c;
    } else if (vfrac_cut < vfrac - vfrac_tolerance) {
      a = c;
    } else {
      return moving_parabola;
    }
    if (a == b) return moving_parabola;
  }

  return moving_parabola;
}

std::pair<Vec, Vec> BoundingBox(const BezierList& cell) {
  if (cell.size() == 0)
    return std::make_pair(Vec(-DBL_MAX, -DBL_MAX), Vec(DBL_MAX, DBL_MAX));
  auto mi = cell[0].first;
  auto ma = cell[0].first;
  for (int i = 0; i < cell.size(); i++) {
    const auto p0 = cell[i].first;
    const auto p1 = cell[i].second;
    const auto p2 = cell[(i + 1) % cell.size()].first;
    mi = minvec(mi, p0);
    mi = minvec(mi, p2);
    ma = maxvec(ma, p0);
    ma = maxvec(ma, p2);
    if (p1.x() < mi.x() || p1.x() > ma.x() || p1.y() < mi.y() ||
        p1.y() > ma.y()) {
      const double tx =
          std::clamp((p0.x() - p1.x()) /
                         IRL::safelyEpsilon(p0.x() - 2.0 * p1.x() + p2.x()),
                     0.0, 1.0);
      const double ty =
          std::clamp((p0.y() - p1.y()) /
                         IRL::safelyEpsilon(p0.y() - 2.0 * p1.y() + p2.y()),
                     0.0, 1.0);
      const double sx = 1.0 - tx;
      const double sy = 1.0 - ty;
      const Vec q{sx * sx * p0.x() + 2.0 * sx * tx * p1.x() + tx * tx * p2.x(),
                  sy * sy * p0.y() + 2.0 * sy * ty * p1.y() + ty * ty * p2.y()};
      mi = minvec(mi, q);
      ma = maxvec(ma, q);
    }
  }

  return std::make_pair(mi, ma);
}

// Functions for mapping moments -----------------------------------------------

std::vector<BezierList> TriangulateCell(const BezierList& cell,
                                        const bool is_preimage) {
  std::vector<BezierList> triangles;
  std::vector<Vec> cell_points;
  for (int i = 0; i < cell.size(); ++i) {
    cell_points.push_back(cell[i].first);
    if (is_preimage == false) {
      cell_points.push_back(cell[i].second);
    }
  }
  Moments cell_moment = ComputeMoments(cell);
  auto centroid = cell_moment.m1() / cell_moment.m0();
  for (int i = 0; i < cell_points.size(); ++i) {
    Vec v1 = cell_points[i];
    Vec v2 = cell_points[(i + 1) % cell_points.size()];
    Vec v3 = centroid;
    Vec c1 = (v1 + v2) / 2.0;
    Vec c2 = (v2 + centroid) / 2.0;
    Vec c3 = (centroid + v1) / 2.0;
    BezierList triangle_list = {{v1, c1}, {v2, c2}, {v3, c3}};
    triangles.push_back(triangle_list);
  }
  return triangles;  // 8 bezier lists with 3 point & control pt pairs each
}

std::pair<Mat, Vec> MappingMatVec(const BezierList& triangle1,
                                  const BezierList& triangle2) {
  Mat A = Mat();
  Vec b = Vec();

  double x1, x2, x3, y1, y2, y3, x1p, x2p, x3p, y1p, y2p, y3p, denominator;
  x1 = triangle1[0].first[0];
  x2 = triangle1[1].first[0];
  x3 = triangle1[2].first[0];
  y1 = triangle1[0].first[1];
  y2 = triangle1[1].first[1];
  y3 = triangle1[2].first[1];

  x1p = triangle2[0].first[0];
  x2p = triangle2[1].first[0];
  x3p = triangle2[2].first[0];
  y1p = triangle2[0].first[1];
  y2p = triangle2[1].first[1];
  y3p = triangle2[2].first[1];

  denominator = (x2 * y1 - x3 * y1 - x1 * y2 + x3 * y2 + x1 * y3 - x2 * y3);

  if (std::abs(denominator) > 1.0e-14) {
    A[0][0] =
        -(-x2p * y1 + x3p * y1 + x1p * y2 - x3p * y2 - x1p * y3 + x2p * y3) /
        denominator;
    A[0][1] =
        -(x1p * x2 - x1 * x2p - x1p * x3 + x2p * x3 + x1 * x3p - x2 * x3p) /
        (-denominator);
    A[1][0] =
        -(y1p * y2 - y1 * y2p - y1p * y3 + y2p * y3 + y1 * y3p - y2 * y3p) /
        denominator;
    A[1][1] =
        -(-x2 * y1p + x3 * y1p + x1 * y2p - x3 * y2p - x1 * y3p + x2 * y3p) /
        denominator;
    b[0] = -(-x2p * x3 * y1 + x2 * x3p * y1 + x1p * x3 * y2 - x1 * x3p * y2 -
             x1p * x2 * y3 + x1 * x2p * y3) /
           (-denominator);
    b[1] = -(-x3 * y1p * y2 + x3 * y1 * y2p + x2 * y1p * y3 - x1 * y2p * y3 -
             x2 * y1 * y3p + x1 * y2 * y3p) /
           denominator;
  }

  return {A, b};
}

Vec MappingPoint(const Mat& A, const Vec& b, const Vec& point) {
  return A * point + b;
}

Mat MappingM2(const Mat& A, const Vec& b, const Moments& tri_liq_moment) {
  double M0 = tri_liq_moment.m0();
  Vec M1 = tri_liq_moment.m1();
  Mat M2 = tri_liq_moment.m2();

  Mat term1 = A * M2 * A.transpose();
  Mat term2 =
      A * Mat(Vec(M1[0] * b[0], M1[0] * b[1]), Vec(M1[1] * b[0], M1[1] * b[1]));
  Mat term3 =
      Mat(Vec(b[0] * M1[0], b[0] * M1[1]), Vec(b[1] * M1[0], b[1] * M1[1])) *
      A.transpose();
  Mat term4 =
      M0 * Mat(Vec(b[0] * b[0], b[0] * b[1]), Vec(b[1] * b[0], b[1] * b[1]));

  return (term1 + term2 + term3 + term4);
}

Moments ComputeMappedTriangleMoments(const Moments& triangle_liq_moments,
                                     const Mat& A, const Vec& b) {
  Moments MappedTriangleMoments;
  double detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];
  // std::cout << detA << std::endl;
  MappedTriangleMoments.m0() = triangle_liq_moments.m0() * detA;
  MappedTriangleMoments.m1() = MappingPoint(A, b * triangle_liq_moments.m0(),
                                            triangle_liq_moments.m1()) *
                               detA;
  MappedTriangleMoments.m2() = MappingM2(A, b, triangle_liq_moments) * detA;

  return MappedTriangleMoments;
}

// for MOF reconstruction on Unit cell

BezierList ComputeTransformedCell(const BezierList& cell,
                                  const bool& toUnitCell) {
  double dx = std::abs(cell[1].first[0] - cell[0].first[0]);
  Vec x0;
  Vec x1;

  if (toUnitCell == true) {
    x0 = cell[0].first - Vec((1 - dx) / 2.0, (1 - dx) / 2.0);
    x1 = cell[2].first + Vec((1 - dx) / 2.0, (1 - dx) / 2.0);
  } else {
    x0 = cell[0].first + Vec((1 - dx) / 2.0, (1 - dx) / 2.0);
    x1 = cell[2].first - Vec((1 - dx) / 2.0, (1 - dx) / 2.0);
  }

  return RectangleFromBounds(x0, x1);
}

Mat MappingCellMat(const BezierList& cell, const bool& toUnitCell) {
  Mat A = Mat();
  const auto TransformedCell = ComputeTransformedCell(cell, toUnitCell);

  double x1, x1p, y1, y1p, x2, x2p, y2, y2p;

  x1 = cell[0].first[0];
  y1 = cell[0].first[1];
  x2 = cell[1].first[0];
  y2 = cell[1].first[1];

  x1p = TransformedCell[0].first[0];
  y1p = TransformedCell[0].first[1];
  x2p = TransformedCell[1].first[0];
  y2p = TransformedCell[1].first[1];

  // std::cout << x1 << " " << x2 << " " << y1 <<  " " << y2 << std::endl;

  // mapping matrix components
  double denominator = x2 * y1 - x1 * y2;
  A[0][0] = (x2p * y1 - x1p * y2) / denominator;
  A[0][1] = -(x1 * x2p - x2 * x1p) / denominator;
  A[1][0] = (y1 * y2p - y2 * y1p) / denominator;
  A[1][1] = (x2 * y1p - x1 * y2p) / denominator;

  return A;
}

Parabola ComputeTransformedParabola(const BezierList& cell,
                                    const Parabola& parabola,
                                    const bool& toUnitCell) {
  Mat A = MappingCellMat(cell, toUnitCell);
  Vec new_datum = A * parabola.datum();
  double new_coeff = parabola.coeff() * A[1][1] / std::pow(A[0][0], 2.0);

  return Parabola(new_datum, parabola.frame(), new_coeff);
}

Moments ComputeTransformedCellMoments(const BezierList& cell,
                                      const Parabola& parabola,
                                      const bool& toUnitCell) {
  Moments TransformedCellMoments;
  Moments moments = ComputeMoments(cell, parabola);
  Mat A = MappingCellMat(cell, toUnitCell);

  // scaling 0th moment
  TransformedCellMoments.m0() = A[0][0] * A[1][1] * moments.m0();

  // scaling 1st moment
  TransformedCellMoments.m1()[0] =
      std::pow(A[0][0], 2.0) * A[1][1] * moments.m1()[0];
  TransformedCellMoments.m1()[1] =
      std::pow(A[0][0], 2.0) * A[1][1] * moments.m1()[1];

  // scaling 2nd moment
  TransformedCellMoments.m2()[0][0] =
      std::pow(A[0][0], 3.0) * A[1][1] * moments.m2()[0][0];
  TransformedCellMoments.m2()[1][1] =
      std::pow(A[1][1], 3.0) * A[0][0] * moments.m2()[1][1];
  TransformedCellMoments.m2()[0][1] =
      std::pow(A[0][0], 2.0) * std::pow(A[1][1], 2.0) * moments.m2()[0][1];
  TransformedCellMoments.m2()[1][0] = TransformedCellMoments.m2()[0][1];

  return TransformedCellMoments;
}

// for curvature estimation -------------------------------------------------

std::vector<Vec> ComputeParticlePositions(const int& N, const Vec& p,
                                          const double& phi,
                                          const double& theta,
                                          const double& hp) {
  // N: Number of particles (odd)
  // p: coordinate of central particle (origin)
  // phi: orientation angle (angle between tangent of p and x-axis)
  // theta: bending angle (turning angle between chords of the circle)
  // hp: distance between particles (or chord length)

  std::vector<Vec> particle_positions(N);

  int c = (N - 1) / 2;  // central particle index

  // computing coordinates of all other particles on the arc
  for (int i = 0; i < N; i++) {
    if (i > c) {
      for (int j = 1; j <= i - c; j++) {
        particle_positions[i] +=
            hp * Vec(std::cos(phi + (static_cast<double>(j) - 0.5) * theta),
                     std::sin(phi + (static_cast<double>(j) - 0.5) * theta));
      }
      particle_positions[i] = p + particle_positions[i];
    } else if (i < c) {
      for (int j = 1; j <= c - i; j++) {
        particle_positions[i] +=
            hp * Vec(std::cos(phi - (static_cast<double>(j) - 0.5) * theta),
                     std::sin(phi - (static_cast<double>(j) - 0.5) * theta));
      }
      particle_positions[i] = p - particle_positions[i];
    } else {
      particle_positions[i] = p;
    }
  }

  return particle_positions;
}

Vec ComputeParticleForce(
    const Vec& x, const std::vector<std::pair<Vec, Vec>>& line_seg_endpoints,
    const double& eta) {
  // x: position of particle
  // line_seg_endpoints: endpoints {a,b} of line segments that are cloest to the
  // particle

  Vec particle_force;

  // Computing closest distance to all line segments in the vicinity of the
  // particle
  for (int i = 0; i < line_seg_endpoints.size(); i++) {
    Vec a = line_seg_endpoints[i].first;
    Vec b = line_seg_endpoints[i].second;

    // finding t using projection
    Vec ab = b - a;
    Vec ax = x - a;
    double t = ax * ab / std::pow(ab.magnitude(), 2.0);
    double t_clamped = std::max(0.0, std::min(1.0, t));

    // finding closet point on the line segment to the point
    Vec y = a + t_clamped * ab;

    // Finding xy distance and keeping minimum value of "force"
    Vec xy = y - x;

    if (i == 0) {
      particle_force = xy;
    } else {
      if (xy.magnitude() < particle_force.magnitude()) {
        particle_force = xy;
      }
    }
  }
  return (eta * particle_force);
}

std::vector<Vec> InitializeParticlePositions(
    const std::pair<Vec, Vec>& target_endpoints, const double& hp,
    const int& N) {
  // target_endpoints: end points of the target interface where curvature is to
  // be estimated hp: spacing between particles along line segment N: number of
  // particles (odd)

  std::vector<Vec> initial_particle_positions(N);

  // line segment end points
  Vec a = target_endpoints.first;
  Vec b = target_endpoints.second;

  // unit vector along the line segment
  Vec unit_ab = (b - a) / (b - a).magnitude();

  // central particle at midpoint of line segment
  initial_particle_positions[(N - 1) / 2] = (a + b) / 2.0;

  // other particles are spaced by hp on either side of the central particle
  // along the line segment
  for (int i = 1; i <= (N - 1) / 2; i++) {
    initial_particle_positions[(N - 1) / 2 + i] =
        initial_particle_positions[(N - 1) / 2] + hp * unit_ab * i;
    initial_particle_positions[(N - 1) / 2 - i] =
        initial_particle_positions[(N - 1) / 2] - hp * unit_ab * i;
  }

  return initial_particle_positions;
}

double ComputeParticleForceProjection(const int& N, const double& phi,
                                      const double& theta, const double& hp,
                                      const bool& iswrtPhi,
                                      const std::vector<Vec> particle_forces) {
  // iswrtPhi: "true" will compute derivative wrt phi else wrt theta

  int c = (N - 1) / 2;  // central particle index

  std::vector<Vec> position_derivative(N);

  position_derivative[c] = Vec(0.0, 0.0);  // central particle

  if (iswrtPhi == true) {
    for (int i = 1; i <= c; i++) {
      // i > c
      position_derivative[c + i] =
          position_derivative[c + (i - 1)] +
          hp * Vec(std::cos(
                       phi +
                       (static_cast<double>(i) - static_cast<double>(c) - 0.5) *
                           theta +
                       M_PI / 2.0),
                   std::sin(
                       phi +
                       (static_cast<double>(i) - static_cast<double>(c) - 0.5) *
                           theta +
                       M_PI / 2.0));
      // i < c
      position_derivative[c - i] =
          position_derivative[c - (i - 1)] -
          hp * Vec(std::cos(
                       phi -
                       (static_cast<double>(c) - static_cast<double>(i) - 0.5) *
                           theta +
                       M_PI / 2.0),
                   std::sin(
                       phi -
                       (static_cast<double>(c) - static_cast<double>(i) - 0.5) *
                           theta +
                       M_PI / 2.0));
    }
  } else {
    for (int i = 1; i <= c; i++) {
      // i > c
      position_derivative[c + i] =
          position_derivative[c + (i - 1)] +
          hp * (static_cast<double>(i) - static_cast<double>(c) - 0.5) *
              Vec(std::cos(
                      phi +
                      (static_cast<double>(i) - static_cast<double>(c) - 0.5) *
                          theta +
                      M_PI / 2.0),
                  std::sin(
                      phi +
                      (static_cast<double>(i) - static_cast<double>(c) - 0.5) *
                          theta +
                      M_PI / 2.0));
      // i < c
      position_derivative[c - i] =
          position_derivative[c - (i - 1)] -
          hp * (static_cast<double>(c) - static_cast<double>(i) - 0.5) *
              Vec(std::cos(
                      phi -
                      (static_cast<double>(c) - static_cast<double>(i) - 0.5) *
                          theta +
                      M_PI / 2.0),
                  std::sin(
                      phi -
                      (static_cast<double>(c) - static_cast<double>(i) - 0.5) *
                          theta +
                      M_PI / 2.0));
    }
  }

  // projecting force on derivative
  double num = 0.0, denom = 0.0;
  for (int i = 0; i < N; i++) {
    num += particle_forces[i] * position_derivative[i];
    denom += position_derivative[i] * position_derivative[i];
  }

  return (num / denom);
}

double getCurvature(const Parabola& target_interface,
                    const BezierList& target_cell,
                    const std::vector<Parabola>& interfaces,
                    const std::vector<BezierList>& cells, const int& N,
                    const double& Hp, const double& h, const double& eta) {
  double theta, phi;
  std::vector<Vec> particle_positions(N), particle_forces(N),
      particle_positions_prev(N), particle_forces_prev(N),
      particle_positions_s(N), particle_positions_ss(N), particle_forces_s(N),
      particle_forces_ss(N);

  // target interface end points
  BezierList clipped_target_interface =
      ParabolaClip(target_cell, target_interface, true);
  std::pair<Vec, Vec> target_endpoints = {clipped_target_interface[0].first,
                                          clipped_target_interface[1].first};

  // end points of all linear interfaces
  std::vector<std::pair<Vec, Vec>> line_seg_endpoints(interfaces.size());
  BezierList clipped_interface;
  for (int i = 0; i < interfaces.size(); i++) {
    clipped_interface = ParabolaClip(cells[i], interfaces[i], true);
    // IRL2D::Print(clipped_interface);
    line_seg_endpoints[i] = {clipped_interface[0].first,
                             clipped_interface[1].first};
  }

  // particle spacing
  double hp = Hp * h / (static_cast<double>(N) - 1.0);

  // initializing particle positions
  particle_positions = InitializeParticlePositions(target_endpoints, hp, N);

  // initialize forces
  for (int i = 0; i < N; i++) {
    particle_forces[i] =
        ComputeParticleForce(particle_positions[i], line_seg_endpoints, eta);
  }

  // initialize orientation and bending angle
  theta = 0.0;
  Vec ab_star = target_endpoints.second - target_endpoints.first;
  phi = std::atan2(ab_star[1], ab_star[0]);
  if (phi < 0.0) {
    phi += 2.0 * M_PI;
  }

  // iteration parameters
  int max_iter = 100;
  double tol = 1e-6;
  int iter = 0;
  double residual = 1.0;

  // index of central particle
  int c = (N - 1) / 2;

  // update positions and forces
  while (std::abs(residual) > tol) {
    iter++;

    // prev iter
    particle_positions_prev = particle_positions;
    particle_forces_prev = particle_forces;

    // step 1: correct central particle position using force
    particle_positions[c] += particle_forces[c];

    // step 1: change in position for other particles
    particle_positions_s =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 1: subtracting change of position from forces
    for (int i = 0; i < N; i++) {
      particle_forces_s[i] =
          particle_forces_prev[i] -
          (particle_positions_s[i] - particle_positions_prev[i]);
    }

    // step 2: correct phi by projection of force
    phi += ComputeParticleForceProjection(N, phi, theta, hp, true,
                                          particle_forces_s);

    // step 2: change in position for other particles
    particle_positions_ss =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 2: subtracting change of position from forces
    for (int i = 0; i < N; i++) {
      particle_forces_ss[i] = particle_forces_s[i] - (particle_positions_ss[i] -
                                                      particle_positions_s[i]);
    }

    // step 3: correct theta by projection of force
    theta -= ComputeParticleForceProjection(N, phi, theta, hp, false,
                                            particle_forces_ss);

    // step 3: update particle positions
    particle_positions =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 3: update particle forces
    for (int i = 0; i < N; i++) {
      particle_forces[i] =
          ComputeParticleForce(particle_positions[i], line_seg_endpoints, eta);
    }

    // residual: change in position for all particles (max value among all
    // particles)
    residual = 0.0;
    for (int i = 0; i < N; i++) {
      residual = std::max(residual,
                          (particle_positions[i] - particle_positions_prev[i])
                              .magnitude()) /
                 (eta * h);
    }

    if (iter == max_iter) {
      break;
    }
  }

  // std::cout << "phi: " << phi * 180.0/M_PI << " theta: " << theta* 180.0/M_PI
  // << std::endl; std::cout << "residual: " << residual << " iter: " << iter <<
  // std::endl;

  return (2 * std::sin(theta / 2.0) / hp);  // curvature
}

// printing particle position and force data (in MATLAB format)
void printParticleData(const std::vector<Vec>& pp, const std::vector<Vec>& pf) {
  // printing x position
  std::cout << "x_pos = [";
  for (auto& p : pp) {
    std::cout << p.x() << " ";
  }
  std::cout << "] ;" << std::endl;

  // printing y position
  std::cout << "y_pos = [";
  for (auto& p : pp) {
    std::cout << p.y() << " ";
  }
  std::cout << "] ;" << std::endl;

  // printing x force
  std::cout << "x_force = [";
  for (auto& f : pf) {
    std::cout << f.x() << " ";
  }
  std::cout << "] ;" << std::endl;

  // printing y force
  std::cout << "y_force = [";
  for (auto& f : pf) {
    std::cout << f.y() << " ";
  }
  std::cout << "] ;" << std::endl;
}

Vec findCircleCenter(const std::vector<Vec>& points) {
  // input is 3 non-collinear points on cirlce
  Vec A = points[0], B = points[1], C = points[2];

  // midpoint of AB and BC
  Vec midAB = Vec((A.x() + B.x()) / 2.0, (A.y() + B.y()) / 2.0);
  Vec midBC = Vec((B.x() + C.x()) / 2.0, (B.y() + C.y()) / 2.0);

  // slopes
  double dxAB = B.x() - A.x();
  double dyAB = B.y() - A.y();
  double dxBC = C.x() - B.x();
  double dyBC = C.y() - B.y();

  // collinearity check
  double area = A.x() * (B.y() - C.y()) + B.x() * (C.y() - A.y()) +
                C.x() * (A.y() - B.y());
  if (std::abs(area) < 1e-10) {
    throw std::runtime_error("Points are collinear. Circle is undefined");
  }

  // slopes for perpendicular bisectors
  double slopePerpAB, slopePerpBC;

  // handling vertical line segments
  bool verticalAB = (std::abs(dxAB) < 1e-10);
  bool verticalBC = (std::abs(dxBC) < 1e-10);
  if (!verticalAB) slopePerpAB = -dxAB / dyAB;
  if (!verticalBC) slopePerpBC = -dxBC / dyBC;

  // solving for center of circle
  double cx, cy;
  if (verticalAB) {
    // AB is vertical, so its perpendicular bisector is horizontal (slope = 0)
    cy = midAB.y();
    cx = slopePerpBC * (cy - midBC.y()) + midBC.x();
  } else if (verticalBC) {
    // BC is vertical, so its perpendicular bisector is horizontal
    cy = midBC.y();
    cx = slopePerpAB * (cy - midAB.y()) + midAB.x();
  } else {
    // Solve intersection of two lines
    cx = (slopePerpAB * midAB.x() - slopePerpBC * midBC.x() + midBC.y() -
          midAB.y()) /
         (slopePerpAB - slopePerpBC);
    cy = slopePerpAB * (cx - midAB.x()) + midAB.y();
  }

  return Vec(cx, cy);
}

std::vector<Vec> findSegmentNormals(
    const std::vector<Vec>& particle_positions,
    const std::vector<std::vector<InterfaceEndPoints>>& plic_data) {
  std::vector<Vec> normals(particle_positions.size());

  // Flatten the plic segments and keep track of (i,j) indices
  std::vector<std::pair<int, int>> mixed_indices;
  std::vector<std::pair<Vec, Vec>> segments;

  for (int i = 0; i < plic_data.size(); i++) {
    for (int j = 0; j < plic_data[i].size(); j++) {
      if (plic_data[i][j].mixed) {
        Vec a(plic_data[i][j].ax, plic_data[i][j].ay);
        Vec b(plic_data[i][j].bx, plic_data[i][j].by);
        segments.emplace_back(a, b);
        mixed_indices.emplace_back(i, j);
      }
    }
  }

  // For each particle, find the closest segment and assign the corresponding
  // normal
  for (int ip = 0; ip < particle_positions.size(); ++ip) {
    const Vec& x = particle_positions[ip];
    double min_dist = std::numeric_limits<double>::infinity();
    int closest = -1;

    for (int s = 0; s < segments.size(); ++s) {
      const Vec& a = segments[s].first;
      const Vec& b = segments[s].second;
      Vec ab = b - a;
      Vec ax = x - a;
      double t = ax * ab / std::pow(ab.magnitude(), 2.0);
      double t_clamped = std::max(0.0, std::min(1.0, t));
      Vec y = a + ab * t_clamped;
      Vec xy = y - x;
      double dist = xy.magnitude();

      if (dist < min_dist) {
        min_dist = dist;
        closest = s;
      }
    }

    // Add the normal of the closest segment
    const auto& [ii, jj] = mixed_indices[closest];
    normals[ip] = Vec(plic_data[ii][jj].nx, plic_data[ii][jj].ny);
  }

  return normals;
}

std::vector<double> computeEta(const std::vector<Vec>& particle_positions,
                               const std::vector<Vec>& pointedPLIC_normals) {
  // cirlce center
  Vec center = findCircleCenter(particle_positions);
  // finding eta for each particle
  std::vector<double> etas(particle_positions.size());
  for (int i = 0; i < particle_positions.size(); ++i) {
    // outward particle normal
    Vec particle_normal = (particle_positions[i] - center);
    particle_normal.normalize();
    // normal of interface
    Vec plic_normal = pointedPLIC_normals[i];
    plic_normal.normalize();
    // eta
    etas[i] = (1.0 + (particle_normal * plic_normal)) / 2.0;
  }
  return etas;
}

void particle_pf(const std::vector<Vec>& pp0, const std::vector<Vec>& pf0,
                 std::vector<Vec>& pp_final, std::vector<Vec>& pf_final,
                 const std::pair<Vec, Vec>& target_endpoints,
                 const std::vector<std::pair<Vec, Vec>>& line_seg_endpoints,
                 const int& N, const double& Hp, const double& h,
                 const double& eta) {
  double theta, phi;
  std::vector<Vec> particle_positions(N), particle_forces(N),
      particle_positions_prev(N), particle_forces_prev(N),
      particle_positions_s(N), particle_positions_ss(N), particle_forces_s(N),
      particle_forces_ss(N);

  // particle spacing
  double hp = Hp * h / (static_cast<double>(N) - 1.0);

  // Initial particle positions and forces
  particle_positions = pp0;
  particle_forces = pf0;

  // initialize orientation and bending angle
  theta = 0.0;
  Vec ab_star = target_endpoints.second - target_endpoints.first;
  phi = std::atan2(ab_star[1], ab_star[0]);
  if (phi < 0.0) {
    phi += 2.0 * M_PI;
  }

  // iteration parameters
  int max_iter = 200;
  double tol = 1e-6;
  int iter = 0;
  double residual = 1.0;

  // index of central particle
  int c = (N - 1) / 2;

  // update positions and forces
  while (std::abs(residual) > tol) {
    iter++;

    // prev iter
    particle_positions_prev = particle_positions;
    particle_forces_prev = particle_forces;

    // step 1: correct central particle position using force
    particle_positions[c] += particle_forces[c];

    // step 1: change in position for other particles
    particle_positions_s =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 1: subtracting change of position from forces
    for (int i = 0; i < N; i++) {
      particle_forces_s[i] =
          particle_forces_prev[i] -
          (particle_positions_s[i] - particle_positions_prev[i]);
    }

    // step 2: correct phi by projection of force
    phi += ComputeParticleForceProjection(N, phi, theta, hp, true,
                                          particle_forces_s);

    // step 2: change in position for other particles
    particle_positions_ss =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 2: subtracting change of position from forces
    for (int i = 0; i < N; i++) {
      particle_forces_ss[i] = particle_forces_s[i] - (particle_positions_ss[i] -
                                                      particle_positions_s[i]);
    }

    // step 3: correct theta by projection of force
    theta -= ComputeParticleForceProjection(N, phi, theta, hp, false,
                                            particle_forces_ss);

    // step 3: update particle positions
    particle_positions =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 3: update particle forces
    for (int i = 0; i < N; i++) {
      particle_forces[i] =
          ComputeParticleForce(particle_positions[i], line_seg_endpoints, eta);
    }

    // residual: change in position for all particles (max value among all
    // particles)
    residual = 0.0;
    for (int i = 0; i < N; i++) {
      residual = std::max(residual,
                          (particle_positions[i] - particle_positions_prev[i])
                              .magnitude()) /
                 (eta * h);
    }

    if (iter == max_iter) {
      break;
    }
  }

  // final particle positions and forces
  pp_final = particle_positions;
  pf_final = particle_forces;

  // printing curvature
  std::cout << "Curvature: " << (2 * std::sin(theta / 2.0) / hp) << std::endl;
  std::cout << "residual: " << residual << std::endl;
  std::cout << "Iteration: " << iter << std::endl;
}

void curvature_vareta(
    const std::vector<Vec>& pp0, const std::vector<Vec>& pf0,
    std::vector<Vec>& pp_final, std::vector<Vec>& pf_final,
    const std::pair<Vec, Vec>& target_endpoints,
    const std::vector<std::pair<Vec, Vec>>& line_seg_endpoints,
    const std::vector<std::vector<InterfaceEndPoints>>& plicDataMat,
    const int& N, const double& Hp, const double& h) {
  double theta, phi;
  std::vector<Vec> particle_positions(N), particle_forces(N),
      particle_positions_prev(N), particle_forces_prev(N),
      particle_positions_s(N), particle_positions_ss(N), particle_forces_s(N),
      particle_forces_ss(N);

  // particle spacing
  double hp = Hp * h / (static_cast<double>(N) - 1.0);

  // Initial particle positions and forces
  particle_positions = pp0;
  particle_forces = pf0;

  // initialize orientation and bending angle
  theta = 0.0;
  Vec ab_star = target_endpoints.second - target_endpoints.first;
  phi = std::atan2(ab_star[1], ab_star[0]);
  if (phi < 0.0) {
    phi += 2.0 * M_PI;
  }

  // iteration parameters
  int max_iter = 200;
  double tol = 1e-6;
  int iter = 0;
  double residual = 1.0;

  // index of central particle
  int c = (N - 1) / 2;

  // update positions and forces
  while (std::abs(residual) > tol) {
    iter++;

    // prev iter
    particle_positions_prev = particle_positions;
    particle_forces_prev = particle_forces;

    // step 1: correct central particle position using force
    particle_positions[c] += particle_forces[c];

    // step 1: change in position for other particles
    particle_positions_s =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 1: subtracting change of position from forces
    for (int i = 0; i < N; i++) {
      particle_forces_s[i] =
          particle_forces_prev[i] -
          (particle_positions_s[i] - particle_positions_prev[i]);
    }

    // step 2: correct phi by projection of force
    phi += ComputeParticleForceProjection(N, phi, theta, hp, true,
                                          particle_forces_s);

    // step 2: change in position for other particles
    particle_positions_ss =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 2: subtracting change of position from forces
    for (int i = 0; i < N; i++) {
      particle_forces_ss[i] = particle_forces_s[i] - (particle_positions_ss[i] -
                                                      particle_positions_s[i]);
    }

    // step 3: correct theta by projection of force
    theta -= ComputeParticleForceProjection(N, phi, theta, hp, false,
                                            particle_forces_ss);

    // step 3: update particle positions
    particle_positions =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // compute eta for each particle
    // std::vector<int> pointedSegmentIndices =
    // findClosestSegmentIndex(particle_positions, line_seg_endpoints);
    // std::vector<Vec> pointedPLIC_normals =
    // findClosestSegmentNormal(pointedSegmentIndices, plicDataMat);
    std::vector<Vec> pointedPLIC_normals =
        findSegmentNormals(particle_positions, plicDataMat);
    std::vector<double> eta =
        computeEta(particle_positions, pointedPLIC_normals);

    // step 3: update particle forces
    for (int i = 0; i < N; i++) {
      particle_forces[i] = ComputeParticleForce(particle_positions[i],
                                                line_seg_endpoints, eta[i]);
    }

    // residual: change in position for all particles (max value among all
    // particles)
    residual = 0.0;
    for (int i = 0; i < N; i++) {
      residual = std::max(residual,
                          (particle_positions[i] - particle_positions_prev[i])
                              .magnitude()) /
                 (eta[i] * h);
    }

    if (iter == max_iter) {
      break;
    }
  }

  // final particle positions and forces
  pp_final = particle_positions;
  pf_final = particle_forces;

  // printing curvature
  std::cout << "Curvature: " << (2 * std::sin(theta / 2.0) / hp) << std::endl;
  std::cout << "residual: " << residual << std::endl;
  std::cout << "Iteration: " << iter << std::endl;
}

// least squares fit of a circle
// --------------------------------------------------------------------
std::vector<Vec> generatePoints(
    const std::vector<std::pair<Vec, Vec>>& line_seg_endpoints) {
  std::vector<Vec> points;
  const int num_points = 10;  // points per interface (including end points)
  for (const auto& seg : line_seg_endpoints) {
    const Vec& start = seg.first;
    const Vec& end = seg.second;
    for (int i = 0; i < num_points; i++) {
      double t = static_cast<double>(i) / (num_points - 1);
      Vec pt = start * (1.0 - t) + end * t;
      points.push_back(pt);
    }
  }
  return points;
}

std::vector<Vec> getPoints(const std::pair<Vec, Vec>& endpoints,
                           const int& num_points) {
  std::vector<Vec> points;
  const Vec& start = endpoints.first;
  const Vec& end = endpoints.second;
  for (int i = 1; i < (num_points - 1); i++) {
    double t = static_cast<double>(i) / (num_points - 1);
    Vec pt = start * (1.0 - t) + end * t;
    points.push_back(pt);
  }
  return points;
}

std::vector<Vec> generateParabolaPoints(
    const std::vector<std::pair<Vec, Vec>>& end_points,
    const std::vector<Parabola>& interfaces) {
  std::vector<Vec> points;
  const int num_points = 10;  // points to be sampled per parabola

  auto globalToLocal = [&](const Vec& p, const Parabola& parabola) -> Vec {
    Vec t = parabola.frame()[0], n = parabola.frame()[1];
    Vec ploc = {(p.x() - parabola.datum().x()) * t.x() +
                    (p.y() - parabola.datum().y()) * t.y(),
                (p.x() - parabola.datum().x()) * n.x() +
                    (p.y() - parabola.datum().y()) * n.y()};
    return ploc;
  };
  auto localToGlobal = [&](const Vec& ploc, const Parabola& parabola) -> Vec {
    Vec t = parabola.frame()[0], n = parabola.frame()[1];
    Vec p = {parabola.datum().x() + t.x() * ploc.x() + n.x() * ploc.y(),
             parabola.datum().y() + t.y() * ploc.x() + n.y() * ploc.y()};
    return p;
  };

  for (int i = 0; i < end_points.size(); i++) {
    Vec p1 = end_points[i].first;
    Vec p2 = end_points[i].second;
    Parabola parabola = interfaces[i];

    Vec p1_loc = globalToLocal(p1, parabola);
    Vec p2_loc = globalToLocal(p2, parabola);

    double x1 = p1_loc.x();
    double x2 = p2_loc.x();

    for (int j = 0; j < num_points; ++j) {
      double t = static_cast<double>(j) / (num_points - 1);
      double x = (1.0 - t) * x1 + t * x2;
      double y = -parabola.coeff() * x * x;
      Vec ploc = {x, y};
      Vec pglob = localToGlobal(ploc, parabola);
      points.push_back(pglob);
    }
  }
  return points;
}

Vec getParabolaCenter(const std::pair<Vec, Vec>& end_points,
                      const Parabola& parabola) {
  auto globalToLocal = [&](const Vec& p) -> Vec {
    Vec t = parabola.frame()[0], n = parabola.frame()[1];
    Vec ploc = {(p.x() - parabola.datum().x()) * t.x() +
                    (p.y() - parabola.datum().y()) * t.y(),
                (p.x() - parabola.datum().x()) * n.x() +
                    (p.y() - parabola.datum().y()) * n.y()};
    return ploc;
  };
  auto localToGlobal = [&](const Vec& ploc) -> Vec {
    Vec t = parabola.frame()[0], n = parabola.frame()[1];
    Vec p = {parabola.datum().x() + t.x() * ploc.x() + n.x() * ploc.y(),
             parabola.datum().y() + t.y() * ploc.x() + n.y() * ploc.y()};
    return p;
  };

  Vec p1_loc = globalToLocal(end_points.first);
  Vec p2_loc = globalToLocal(end_points.second);
  double xloc_center = 0.5 * (p1_loc.x() + p2_loc.x());
  double yloc_center = -parabola.coeff() * xloc_center * xloc_center;
  Vec center_loc = {xloc_center, yloc_center};
  Vec center = localToGlobal(center_loc);
  return center;
}

double getVfracWeight(double vfrac) {
  const double limit_vfrac = 0.1;
  if (vfrac < limit_vfrac) {
    return 0.5 - 0.5 * std::cos(M_PI * vfrac / limit_vfrac);
  } else if (vfrac > 1.0 - limit_vfrac) {
    return 0.5 - 0.5 * std::cos(M_PI * (1.0 - vfrac) / limit_vfrac);
  } else {
    return 1.0;
  }
}

double getDistanceWeight(const Vec& pref, const Vec& ploc, const double& h) {
  double distance = (ploc - pref).magnitude() / h;
  if (distance < 2.5) {
    return (1.0 + 4.0 * distance / 2.5) * pow(1.0 - distance / 2.5, 4.0);
  } else {
    return 0.0;
  }
}

double DistanceWeight(const Vec& pref, const Vec& ploc, const double& h,
                      const double& delta) {
  double distance = (ploc - pref).magnitude() / h;
  if (distance < delta) {
    return (1.0 + 4.0 * distance / delta) * pow(1.0 - distance / delta, 4.0);
  } else {
    return 0.0;
  }
}

// based on dot product
double getNormalWeight(const Vec& nref, const Vec& nloc) {
  return std::max(0.0, (nloc * nref));
}

// double getNormalWeight(const Vec& nref, const Vec& nloc) {
//   double dot = nloc * nref;
//   if (dot <= -0.5) return 0.0;
//   if (dot >= 1.0) return 1.0;
//   return (dot + 0.5) / 1.5;  // linearly map from [-0.5, 1] to [0, 1]
// }

// double getNormalWeight(const Vec& nref, const Vec& nloc) {
//   double dot = nloc * nref;
//   if (dot <= -0.5) {
//     return 0.0;
//   } else {
//     return 1.0;
//   }
// }

// double getNormalGradWeight(const Vec& nref, const Vec& nloc,
//                            const Vec& pref, const Vec& ploc) {
//   double G = 100; // maxgradient magnitude
//   Vec dn = nloc - nref;
//   Vec dp = ploc - pref;

//   double dn_mag = dn.magnitude();
//   double dp_mag = dp.magnitude();

//   if (dp_mag < 1e-10) return 1.0;

//   double grad_mag = dn_mag / dp_mag;
//   // std::cout << grad_mag << std::endl;

//   if (grad_mag < G) {
//     double q = grad_mag / G;
//     return (1.0 + 4.0 * q) * pow(1.0 - q, 4.0);
//   } else {
//     return 0.0;
//   }
// }

Mat estimateFrame(const Vec& circle_center, const Vec& plic_center,
                  const double& r, const Vec& plic_normal, bool& flip_coeff) {
  Vec dir = plic_center - circle_center;
  dir.normalize();
  Vec datum =
      circle_center + r * dir;  // point on circle closest to cell centroid
  Vec normal = datum - circle_center;
  normal.normalize();

  // comparing to plic normal
  if ((normal * plic_normal) < 0) {
    normal *= -1.0;
    flip_coeff = true;
  }

  Vec tangent = {normal.y(), -normal.x()};
  return {tangent, normal};
}

// Pratt parabola
Parabola getPrattParabola(const std::vector<IRL2D::Vec>& points,
                          const std::vector<double>& vfw,
                          const std::vector<double>& dw,
                          const std::vector<double>& nw, const Mat& plic_frame,
                          const Vec& plic_center) {
  Parabola PrattParabola;

  // moment matrix
  int n = points.size();
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    double w = vfw[i] * dw[i] * nw[i];
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();
  }
  Eigen::Matrix4d B;  // constraint matrix
  B << 0, 0, 0, -2, 0, 1, 0, 0, 0, 0, 1, 0, -2, 0, 0, 0;

  // solving the generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, B);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // extracting smallest positive eigenvalue and its eigenvector
  std::vector<std::pair<double, int>> positive_eigs;
  for (int i = 0; i < eigenvalues.size(); i++) {
    double real_part = eigenvalues[i].real();
    if (real_part > 0) {
      positive_eigs.emplace_back(real_part, i);
    }
  }
  std::sort(positive_eigs.begin(), positive_eigs.end());
  double A = eigenvectors.col(positive_eigs[0].second)[0].real();
  double B_ = eigenvectors.col(positive_eigs[0].second)[1].real();
  double C = eigenvectors.col(positive_eigs[0].second)[2].real();
  double D = eigenvectors.col(positive_eigs[0].second)[3].real();

  // scaling parameters
  double constraint = B_ * B_ + C * C - 4.0 * A * D;
  double scale_factor = 1.0 / std::sqrt(constraint);
  A *= scale_factor;
  B_ *= scale_factor;
  C *= scale_factor;
  D *= scale_factor;

  // circle radius and center
  double radius = std::sqrt((B_ * B_ + C * C - 4.0 * A * D) / (4.0 * A * A));
  Vec circle_center = IRL2D::Vec(-B_ / (2.0 * A), -C / (2.0 * A));

  // parabola coefficient
  PrattParabola.coeff() = 0.5 / radius;

  // reference frame
  bool flip_coeff = false;
  PrattParabola.frame() = estimateFrame(circle_center, plic_center, radius,
                                        plic_frame[1], flip_coeff);
  if (flip_coeff == true) {
    PrattParabola.coeff() = -PrattParabola.coeff();
  }

  // datum
  Vec direction = plic_center - circle_center;
  direction.normalize();
  PrattParabola.datum() = circle_center + radius * direction;

  return PrattParabola;
}

Parabola getPrattParabola_localframe(const std::vector<IRL2D::Vec>& points,
                                     const std::vector<double>& vfw,
                                     const std::vector<double>& dw,
                                     const std::vector<double>& nw,
                                     const Mat& plic_frame,
                                     const Vec& plic_center) {
  Parabola PrattParabola;

  // points in local frame
  std::vector<Vec> pts_local(points.size());
  double tx = plic_frame[0][0], ty = plic_frame[0][1];
  double nx = plic_frame[1][0], ny = plic_frame[1][1];
  for (int i = 0; i < points.size(); i++) {
    Vec p = points[i];
    pts_local[i] = {
        (p.x() - plic_center.x()) * tx + (p.y() - plic_center.y()) * ty,
        (p.x() - plic_center.x()) * nx + (p.y() - plic_center.y()) * ny};
  }

  // moment matrix
  int n = points.size();
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = pts_local[i].x(), yi = pts_local[i].y();
    double zi = xi * xi + yi * yi;
    double w = vfw[i] * dw[i] * nw[i];
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();
  }
  Eigen::Matrix4d B;  // constraint matrix
  B << 0, 0, 0, -2, 0, 1, 0, 0, 0, 0, 1, 0, -2, 0, 0, 0;

  // solving the generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, B);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // extracting smallest positive eigenvalue and its eigenvector
  std::vector<std::pair<double, int>> positive_eigs;
  for (int i = 0; i < eigenvalues.size(); i++) {
    double real_part = eigenvalues[i].real();
    if (real_part > 0) {
      positive_eigs.emplace_back(real_part, i);
    }
  }
  std::sort(positive_eigs.begin(), positive_eigs.end());
  double A = eigenvectors.col(positive_eigs[0].second)[0].real();
  double B_ = eigenvectors.col(positive_eigs[0].second)[1].real();
  double C = eigenvectors.col(positive_eigs[0].second)[2].real();
  double D = eigenvectors.col(positive_eigs[0].second)[3].real();

  // scaling parameters
  double constraint = B_ * B_ + C * C - 4.0 * A * D;
  double scale_factor = 1.0 / std::sqrt(constraint);
  A *= scale_factor;
  B_ *= scale_factor;
  C *= scale_factor;
  D *= scale_factor;

  // circle radius and center in local frame
  double radius = std::sqrt((B_ * B_ + C * C - 4.0 * A * D) / (4.0 * A * A));
  Vec cc_local = IRL2D::Vec(-B_ / (2.0 * A), -C / (2.0 * A));

  // circle center in global frame
  Vec circle_center = {plic_center.x() + tx * cc_local.x() + nx * cc_local.y(),
                       plic_center.y() + ty * cc_local.x() + ny * cc_local.y()};

  // parabola coefficient
  PrattParabola.coeff() = 0.5 / radius;

  // reference frame
  bool flip_coeff = false;
  PrattParabola.frame() = estimateFrame(circle_center, plic_center, radius,
                                        plic_frame[1], flip_coeff);
  if (flip_coeff == true) {
    PrattParabola.coeff() = -PrattParabola.coeff();
  }

  // datum
  Vec direction = plic_center - circle_center;
  direction.normalize();
  PrattParabola.datum() = circle_center + radius * direction;

  return PrattParabola;
}

// Taubin parabola
Parabola getTaubinParabola(const std::vector<IRL2D::Vec>& points,
                           const std::vector<double>& vfw,
                           const std::vector<double>& dw,
                           const std::vector<double>& nw,
                           const Vec& plic_normal, const Vec& plic_center) {
  Parabola TaubinParabola;

  // moment matrix
  int n = points.size();
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    double w = vfw[i] * dw[i] * nw[i];  // volume fraction and distance weight
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();  // weighted outer product
  }

  // constraint matrix
  Eigen::Matrix4d C;
  C.setZero();
  C(0, 0) = 4.0 * M(0, 3);
  C(0, 1) = 2.0 * M(1, 3);
  C(0, 2) = 2.0 * M(2, 3);
  C(1, 0) = C(0, 1);
  C(1, 1) = n;
  C(2, 0) = C(0, 2);
  C(2, 2) = n;

  // solving the generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // extracting smallest positive eigenvalue and its eigenvector
  std::vector<std::pair<double, int>> positive_eigs;
  for (int i = 0; i < eigenvalues.size(); i++) {
    double real_part = eigenvalues[i].real();
    if (real_part > 0) {
      positive_eigs.emplace_back(real_part, i);
    }
  }
  std::sort(positive_eigs.begin(), positive_eigs.end());
  double A = eigenvectors.col(positive_eigs[0].second)[0].real();
  double B = eigenvectors.col(positive_eigs[0].second)[1].real();
  double C_ = eigenvectors.col(positive_eigs[0].second)[2].real();
  double D = eigenvectors.col(positive_eigs[0].second)[3].real();

  // scaling parameters
  double constraint = 4.0 * A * A * M(0, 3) + 4.0 * A * B * M(1, 3) +
                      4.0 * A * C_ * M(2, 3) + B * B * static_cast<double>(n) +
                      C_ * C_ * static_cast<double>(n);
  double scale_factor = 1.0 / std::sqrt(constraint) * std::sqrt(80.0);
  A *= scale_factor;
  B *= scale_factor;
  C_ *= scale_factor;
  D *= scale_factor;

  // parabola coefficient
  double radius = std::sqrt((B * B + C_ * C_ - 4.0 * A * D) / (4.0 * A * A));
  TaubinParabola.coeff() = 0.5 / radius;

  // reference frame
  Vec circle_center = IRL2D::Vec(-B / (2.0 * A), -C_ / (2.0 * A));
  bool flip_coeff = false;
  TaubinParabola.frame() = estimateFrame(circle_center, plic_center, radius,
                                         plic_normal, flip_coeff);
  if (flip_coeff == true) {
    TaubinParabola.coeff() = -TaubinParabola.coeff();
  }

  // datum
  Vec direction = plic_center - circle_center;
  direction.normalize();
  TaubinParabola.datum() = circle_center + radius * direction;

  return TaubinParabola;
}

// Taubin fit in local frame of reference
Parabola getTaubinParabola_localframe(const std::vector<IRL2D::Vec>& points,
                                      const std::vector<double>& vfw,
                                      const std::vector<double>& dw,
                                      const std::vector<double>& nw,
                                      const Mat& plic_frame,
                                      const Vec& plic_center) {
  Parabola TaubinParabola;

  // points in local frame
  std::vector<Vec> pts_local(points.size());
  double tx = plic_frame[0][0], ty = plic_frame[0][1];
  double nx = plic_frame[1][0], ny = plic_frame[1][1];
  for (int i = 0; i < points.size(); i++) {
    Vec p = points[i];
    pts_local[i] = {
        (p.x() - plic_center.x()) * tx + (p.y() - plic_center.y()) * ty,
        (p.x() - plic_center.x()) * nx + (p.y() - plic_center.y()) * ny};
  }

  // moment matrix
  int n = pts_local.size();
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = pts_local[i].x(), yi = pts_local[i].y();
    double zi = xi * xi + yi * yi;
    double w = vfw[i] * dw[i] * nw[i];
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();
  }

  // constraint matrix
  Eigen::Matrix4d C;
  C.setZero();
  C(0, 0) = 4.0 * M(0, 3);
  C(0, 1) = 2.0 * M(1, 3);
  C(0, 2) = 2.0 * M(2, 3);
  C(1, 0) = C(0, 1);
  C(1, 1) = n;
  C(2, 0) = C(0, 2);
  C(2, 2) = n;

  // solving the generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // extracting smallest positive eigenvalue and its eigenvector
  std::vector<std::pair<double, int>> positive_eigs;
  for (int i = 0; i < eigenvalues.size(); i++) {
    double real_part = eigenvalues[i].real();
    if (real_part > 0) {
      positive_eigs.emplace_back(real_part, i);
    }
  }
  std::sort(positive_eigs.begin(), positive_eigs.end());
  double A = eigenvectors.col(positive_eigs[0].second)[0].real();
  double B = eigenvectors.col(positive_eigs[0].second)[1].real();
  double C_ = eigenvectors.col(positive_eigs[0].second)[2].real();
  double D = eigenvectors.col(positive_eigs[0].second)[3].real();

  // scaling parameters
  double constraint = 4.0 * A * A * M(0, 3) + 4.0 * A * B * M(1, 3) +
                      4.0 * A * C_ * M(2, 3) + B * B * static_cast<double>(n) +
                      C_ * C_ * static_cast<double>(n);
  double scale_factor = 1.0 / std::sqrt(constraint) * std::sqrt(80.0);
  A *= scale_factor;
  B *= scale_factor;
  C_ *= scale_factor;
  D *= scale_factor;

  // circle radius and center in local frame
  double radius = std::sqrt((B * B + C_ * C_ - 4.0 * A * D) / (4.0 * A * A));
  Vec cc_local = IRL2D::Vec(-B / (2.0 * A), -C_ / (2.0 * A));

  // circle center in global frame
  Vec circle_center = {plic_center.x() + tx * cc_local.x() + nx * cc_local.y(),
                       plic_center.y() + ty * cc_local.x() + ny * cc_local.y()};

  // parabola coefficient
  TaubinParabola.coeff() = 0.5 / radius;

  // reference frame
  bool flip_coeff = false;
  TaubinParabola.frame() = estimateFrame(circle_center, plic_center, radius,
                                         plic_frame[1], flip_coeff);
  if (flip_coeff == true) {
    TaubinParabola.coeff() = -TaubinParabola.coeff();
  }

  // datum
  Vec direction = plic_center - circle_center;
  direction.normalize();
  TaubinParabola.datum() = circle_center + radius * direction;

  return TaubinParabola;
}

std::vector<double> getPrattParams(const std::vector<IRL2D::Vec>& points,
                                   const std::vector<double>& vfw,
                                   const std::vector<double>& dw,
                                   const std::vector<double>& nw) {
  // Generating matrices for eigenvalue problem
  int n = points.size();
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    double w = vfw[i] * dw[i] * nw[i];  // volume fraction and distance weight
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();  // weighted outer product
  }
  Eigen::Matrix4d B;  // constraint matrix
  B << 0, 0, 0, -2, 0, 1, 0, 0, 0, 0, 1, 0, -2, 0, 0, 0;

  // solving the generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, B);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // extracting smallest positive eigenvalue and its eigenvector
  std::vector<std::pair<double, int>> positive_eigs;
  for (int i = 0; i < eigenvalues.size(); i++) {
    double real_part = eigenvalues[i].real();
    if (real_part > 0) {
      positive_eigs.emplace_back(real_part, i);
    }
  }
  std::sort(positive_eigs.begin(), positive_eigs.end());
  double A = eigenvectors.col(positive_eigs[0].second)[0].real();
  double B_ = eigenvectors.col(positive_eigs[0].second)[1].real();
  double C = eigenvectors.col(positive_eigs[0].second)[2].real();
  double D = eigenvectors.col(positive_eigs[0].second)[3].real();

  // scaling parameters
  double constraint = B_ * B_ + C * C - 4.0 * A * D;
  double scale_factor = 1.0 / std::sqrt(constraint);
  A *= scale_factor;
  B_ *= scale_factor;
  C *= scale_factor;
  D *= scale_factor;
  double theta = std::atan2(C, B_);

  return {A, D, theta};
}

std::vector<double> getTaubinParams(const std::vector<IRL2D::Vec>& points,
                                    const std::vector<double>& vfw,
                                    const std::vector<double>& dw,
                                    const std::vector<double>& nw) {
  // moment matrix
  int n = points.size();
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    double w = vfw[i] * dw[i] * nw[i];  // volume fraction and distance weight
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();  // weighted outer product
  }

  // constraint matrix
  Eigen::Matrix4d C;
  C.setZero();
  C(0, 0) = 4.0 * M(0, 3);
  C(0, 1) = 2.0 * M(1, 3);
  C(0, 2) = 2.0 * M(2, 3);
  C(1, 0) = C(0, 1);
  C(1, 1) = n;
  C(2, 0) = C(0, 2);
  C(2, 2) = n;

  // solving the generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // extracting smallest positive eigenvalue and its eigenvector
  std::vector<std::pair<double, int>> positive_eigs;
  for (int i = 0; i < eigenvalues.size(); i++) {
    double real_part = eigenvalues[i].real();
    if (real_part > 0) {
      positive_eigs.emplace_back(real_part, i);
    }
  }
  std::sort(positive_eigs.begin(), positive_eigs.end());
  double A = eigenvectors.col(positive_eigs[0].second)[0].real();
  double B = eigenvectors.col(positive_eigs[0].second)[1].real();
  double C_ = eigenvectors.col(positive_eigs[0].second)[2].real();
  double D = eigenvectors.col(positive_eigs[0].second)[3].real();

  // scaling parameters
  double constraint = 4.0 * A * A * M(0, 3) + 4.0 * A * B * M(1, 3) +
                      4.0 * A * C_ * M(2, 3) + B * B * static_cast<double>(n) +
                      C_ * C_ * static_cast<double>(n);
  double scale_factor = 1.0 / std::sqrt(constraint) * std::sqrt(80.0);
  A *= scale_factor;
  B *= scale_factor;
  C_ *= scale_factor;
  D *= scale_factor;
  double theta = std::atan2(C_, B);

  return {A, D, theta};
}

// 2D Jibben
// -------------------------------------------------------------------------

// Parabola getParabolaJibben(const Parabola& target_interface, const
// BezierList& target_cell,
//                            const std::vector<Parabola>& interfaces, const
//                            std::vector<BezierList>& cells){

//   Parabola parabolaJibben;

//   using segment = std::pair<Vec, Vec>;

//   // interface end points
//   // std::vector<segment> line_seg_endpoints(interfaces.size());
//   // BezierList clipped_interface;
//   // for (int i = 0; i < interfaces.size(); i++){
//   //   clipped_interface = ParabolaClip(cells[i], interfaces[i], true);
//   //   line_seg_endpoints[i] = {clipped_interface[0].first,
//   clipped_interface[1].first};
//   // }

//   std::vector<segment> line_seg_endpoints;
//   for (int i = 0; i < interfaces.size(); i++){
//     BezierList clipped_interface = ParabolaClip(cells[i], interfaces[i],
//     true); if (clipped_interface.size() < 2){
//       IRL2D::Print(clipped_interface);
//       continue;
//     }
//     line_seg_endpoints.push_back( {clipped_interface[0].first,
//     clipped_interface[1].first} );
//   }

//   if (line_seg_endpoints.empty()){
//     std::cout << "Empty line seg endpoints" << std::endl;
//   }

//   // local reference frame and origin
//   Vec datum = target_interface.datum();
//   Vec t = target_interface.frame()[0]; Vec n = target_interface.frame()[1];

//   // least squares problem
//   Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
//   Eigen::Vector3d b = Eigen::Vector3d::Zero();

//   for (const auto& [p1, p2]: line_seg_endpoints){
//     // translation
//     Vec r1 = p1 - datum; Vec r2 = p2 - datum;

//     // projection,  (xi,eta): local frame coordinates
//     double xi1 = r1 * t; double eta1 = r1 * n;
//     double xi2 = r2 * t; double eta2 = r2 * n;

//     // Linear form: eta = a * xi + b
//     double a_i = (eta2 - eta1) / (xi2 - xi1);
//     double b_i = eta1 - a_i * xi1;

//     // terms from integral equation
//     double s0 = xi2 - xi1;
//     double s1 = 0.5 * (std::pow(xi2, 2.0) - std::pow(xi1, 2.0));
//     double s2 = (1.0 / 3.0) * (std::pow(xi2, 3.0) - std::pow(xi1, 3.0));
//     Eigen::Vector3d S(s0, s1, s2);
//     double d_i = b_i * s0 + a_i * s1;

//     A += S * S.transpose();
//     b += d_i * S;
//   }

//   // solving Ac=b
//   Eigen::Vector3d coeffs = A.ldlt().solve(b);

//   // parabola vertex in local frame
//   double xi_v = - coeffs(1) / (2.0 * coeffs(2));
//   double eta_v = coeffs(0) + coeffs(1) * xi_v + coeffs(2) *
//   std::pow(xi_v, 2.0);

//   // fitted parabola parameters
//   parabolaJibben.coeff() = - coeffs(2);
//   parabolaJibben.frame() = target_interface.frame();
//   parabolaJibben.datum() = datum + xi_v * t + eta_v * n;

//   return parabolaJibben;
// }

Parabola getParabolaJibben(const Parabola& target_interface,
                           const BezierList& target_cell,
                           const std::vector<NeighborInfo>& neighbors,
                           const int i_target, const int j_target) {
  Parabola parabolaJibben;

  using segment = std::pair<Vec, Vec>;
  std::vector<segment> line_seg_endpoints;

  for (const auto& neighbor : neighbors) {
    const IRL2D::Parabola interface = neighbor.interface;
    const IRL2D::BezierList cell = neighbor.cell;

    BezierList clipped_interface = ParabolaClip(cell, interface, true);
    BezierList clipped_cell = ParabolaClip(cell, interface, false);

    if (clipped_interface.size() < 2) {
      std::cout << "Warning: Clipped interface has size < 2\n";
      std::cout << "Target cell global ID: (" << i_target << ", " << j_target
                << ")\n";
      std::cout << "Neighbor global ID: (" << neighbor.ii_global << ", "
                << neighbor.jj_global << ")\n";
      IRL2D::Print(clipped_interface);
      std::cout << "coefficient = " << interface.coeff() << std::endl;
      std::cout << "datum = " << interface.datum() << std::endl;
      std::cout << "frame = " << interface.frame() << std::endl;
      std::cout << "cell liquid volume fraction = " << neighbor.lvf
                << std::endl;
      std::cout << "volume fraction of clipped cell = "
                << IRL2D::ComputeArea(clipped_cell) / IRL2D::ComputeArea(cell)
                << std::endl;
      IRL2D::Print(clipped_cell);
      continue;
    }

    line_seg_endpoints.push_back(
        {clipped_interface[0].first, clipped_interface[1].first});
  }

  if (line_seg_endpoints.empty()) {
    std::cout << "Warning: No valid line segments found for Jibben fitting.\n";
    return target_interface;  // Return PLIC
  }

  // local reference frame and origin
  Vec datum = target_interface.datum();
  Vec t = target_interface.frame()[0];
  Vec n = target_interface.frame()[1];

  // least squares problem
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
  Eigen::Vector3d b = Eigen::Vector3d::Zero();

  for (const auto& [p1, p2] : line_seg_endpoints) {
    // translation
    Vec r1 = p1 - datum;
    Vec r2 = p2 - datum;

    // projection,  (xi,eta): local frame coordinates
    double xi1 = r1 * t;
    double eta1 = r1 * n;
    double xi2 = r2 * t;
    double eta2 = r2 * n;

    // Linear form: eta = a * xi + b
    double a_i = (eta2 - eta1) / (xi2 - xi1);
    double b_i = eta1 - a_i * xi1;

    // terms from integral equation
    double s0 = xi2 - xi1;
    double s1 = 0.5 * (std::pow(xi2, 2.0) - std::pow(xi1, 2.0));
    double s2 = (1.0 / 3.0) * (std::pow(xi2, 3.0) - std::pow(xi1, 3.0));
    Eigen::Vector3d S(s0, s1, s2);
    double d_i = b_i * s0 + a_i * s1;

    A += S * S.transpose();
    b += d_i * S;
  }

  // solving Ac=b
  // Eigen::Vector3d coeffs = A.ldlt().solve(b);
  Eigen::Vector3d coeffs = A.colPivHouseholderQr().solve(b);

  // parabola vertex in local frame
  double xi_v = -coeffs(1) / (2.0 * coeffs(2));
  double eta_v = coeffs(0) + coeffs(1) * xi_v + coeffs(2) * std::pow(xi_v, 2.0);

  // fitted parabola parameters
  parabolaJibben.coeff() = -coeffs(2);
  parabolaJibben.frame() = target_interface.frame();
  parabolaJibben.datum() = datum + xi_v * t + eta_v * n;

  // curvature check
  const double maxkdx = 4.0;
  const double length_scale = std::sqrt(ComputeArea(target_cell));
  const double kdx = 2.0 * parabolaJibben.coeff() * length_scale;

  if (std::abs(kdx) > maxkdx) {
    // std::cout << "Warning: Curvature too large, reverting to planar
    // interface.\n";
    parabolaJibben.coeff() = 0.0;
  }

  return parabolaJibben;
}

// std::vector<double> getJibbenCoeffs(const Parabola& target_interface, const
// BezierList& target_cell,
//                                     const std::vector<NeighborInfo>&
//                                     neighbors, const std::vector<double>&
//                                     weights){

//   using segment = std::pair<Vec, Vec>;
//   std::vector<segment> line_seg_endpoints;

//   for (const auto& neighbor : neighbors) {
//     const IRL2D::Parabola interface = neighbor.interface;
//     const IRL2D::BezierList cell = neighbor.cell;
//     BezierList clipped_interface = ParabolaClip(cell, interface, true);
//     BezierList clipped_cell = ParabolaClip(cell, interface, false);
//     line_seg_endpoints.push_back({clipped_interface[0].first,
//     clipped_interface[1].first});
//   }

//   // local reference frame and origin
//   Vec datum = target_interface.datum();
//   Vec t = target_interface.frame()[0];
//   Vec n = target_interface.frame()[1];

//   // least squares problem
//   Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
//   Eigen::Vector3d b = Eigen::Vector3d::Zero();

//   for (const auto& [p1, p2]: line_seg_endpoints){
//     // translation
//     Vec r1 = p1 - datum; Vec r2 = p2 - datum;

//     // projection,  (xi,eta): local frame coordinates
//     double xi1 = r1 * t; double eta1 = r1 * n;
//     double xi2 = r2 * t; double eta2 = r2 * n;

//     // Linear form: eta = a * xi + b
//     double a_i = (eta2 - eta1) / (xi2 - xi1);
//     double b_i = eta1 - a_i * xi1;

//     // terms from integral equation
//     double s0 = xi2 - xi1;
//     double s1 = 0.5 * (std::pow(xi2, 2.0) - std::pow(xi1, 2.0));
//     double s2 = (1.0 / 3.0) * (std::pow(xi2, 3.0) - std::pow(xi1, 3.0));
//     Eigen::Vector3d S(s0, s1, s2);
//     double d_i = b_i * s0 + a_i * s1;

//     A += S * S.transpose();
//     b += d_i * S;
//   }

//   // parabola coefficients in local frame
//   Eigen::Vector3d coeffs = A.colPivHouseholderQr().solve(b);
//   return {coeffs(2), coeffs(1), coeffs(0)};
// }

std::vector<double> getJibbenCoeffs(const Parabola& target_interface,
                                    const BezierList& target_cell,
                                    const std::vector<NeighborInfo>& neighbors,
                                    const std::vector<double>& weights) {
  using segment = std::pair<Vec, Vec>;
  std::vector<segment> line_seg_endpoints;

  for (const auto& neighbor : neighbors) {
    const IRL2D::Parabola interface = neighbor.interface;
    const IRL2D::BezierList cell = neighbor.cell;
    BezierList clipped_interface = ParabolaClip(cell, interface, true);
    BezierList clipped_cell = ParabolaClip(cell, interface, false);
    line_seg_endpoints.push_back(
        {clipped_interface[0].first, clipped_interface[1].first});
  }

  if (line_seg_endpoints.size() != weights.size()) {
    throw std::runtime_error(
        "Mismatch between number of segments and weights.");
  }

  // local reference frame and origin
  Vec datum = target_interface.datum();
  Vec t = target_interface.frame()[0];
  Vec n = target_interface.frame()[1];

  // least squares problem
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
  Eigen::Vector3d b = Eigen::Vector3d::Zero();

  for (size_t i = 0; i < line_seg_endpoints.size(); ++i) {
    const auto& [p1, p2] = line_seg_endpoints[i];
    double w = weights[i];
    if (w <= 0.0) continue;  // Skip zero or negative weights

    // translation
    Vec r1 = p1 - datum;
    Vec r2 = p2 - datum;

    // projection, (xi, eta): local frame coordinates
    double xi1 = r1 * t, eta1 = r1 * n;
    double xi2 = r2 * t, eta2 = r2 * n;

    if (std::abs(xi2 - xi1) < 1e-10) continue;  // skip degenerate segment

    // Linear form: eta = a * xi + b
    double a_i = (eta2 - eta1) / (xi2 - xi1);
    double b_i = eta1 - a_i * xi1;

    // terms from integral equation
    double s0 = xi2 - xi1;
    double s1 = 0.5 * (std::pow(xi2, 2.0) - std::pow(xi1, 2.0));
    double s2 = (1.0 / 3.0) * (std::pow(xi2, 3.0) - std::pow(xi1, 3.0));
    Eigen::Vector3d S(s0, s1, s2);
    double d_i = b_i * s0 + a_i * s1;

    // Weighted contribution
    A += w * (S * S.transpose());
    b += w * d_i * S;
  }

  Eigen::Vector3d coeffs = A.colPivHouseholderQr().solve(b);
  return {coeffs(2), coeffs(1), coeffs(0)};
}

// Functions for partition of unity
// ----------------------------------------------------------------------

Vec projectToImplicitSurface(const Vec& x0, const std::vector<Vec>& centroids,
                             const std::vector<Vec>& normals,
                             const double& kernel_size, bool& usePlane) {
  int max_iter = 200;
  double tol = 1e-12;
  Vec x = x0;
  ImplicitSurface IS(centroids, normals, kernel_size);

  // for NaN debugging
  std::ostringstream diagnostics;
  diagnostics << "x0 = " << x0 << std::endl;
  bool nan_detected = false;

  for (int i = 0; i < max_iter; i++) {
    double F = IS.F(x);
    double Fx = IS.Fx(x);
    double Fy = IS.Fy(x);
    Vec gradF = {Fx, Fy};
    Vec change = F * gradF / (Fx * Fx + Fy * Fy);
    x = x - change;

    // diagnostics for NaN
    diagnostics << "Iteration: " << i + 1 << "\n";
    diagnostics << "F(x) = " << F << "\n";
    diagnostics << "Fx = " << Fx << ", Fy = " << Fy << "\n";
    diagnostics << "gradF = " << gradF << "\n";
    diagnostics << "change = " << change
                << ", |change| = " << change.magnitude() << "\n";
    diagnostics << "x = " << x << "\n";
    diagnostics << "distance to x0 in terms of dx = "
                << (x - x0).magnitude() / (kernel_size / 2.5) << "\n";
    diagnostics << "--------------------------------------------\n";

    // NaN check
    if (std::isnan(x.magnitude()) && !nan_detected) {
      nan_detected = true;
      usePlane = true;
      std::cout << "========== NaN detected at iteration " << i + 1
                << " ==========\n";
      std::cout << diagnostics.str();
      // break;
    }

    // iteration parameter checks
    if (std::fabs(change.magnitude()) < tol && !nan_detected)
      break;  // converged solution
    if (i == (max_iter - 1)) {
      std::cout << "Max iterations reached: Projection to implicit surface is "
                   "incomplete ";
      if (std::isnan(x.magnitude())) {
        std::cout << "because NaN is detected";
      }
      std::cout << std::endl;
    }
  }
  return x;
}

Parabola getPU_interface(const Vec& x0, const std::vector<Vec>& centroids,
                         const std::vector<Vec>& normals,
                         const double& kernel_size, bool& usePlane) {
  Parabola PU_parabola;

  // projecting x0 onto implicit surface
  Vec x_proj =
      projectToImplicitSurface(x0, centroids, normals, kernel_size, usePlane);
  PU_parabola.datum() = x_proj;

  // reference frame
  ImplicitSurface IS(centroids, normals, kernel_size);
  double Fx = IS.Fx(x_proj);
  double Fy = IS.Fy(x_proj);
  IRL2D::Vec normal = {Fx, Fy};
  normal.normalize();
  IRL2D::Vec tangent = {normal.y(), -normal.x()};
  PU_parabola.frame() = {tangent, normal};

  // coefficient
  std::vector<double> HessianTerms = IS.HessianTerms(x_proj);
  double Fxx = HessianTerms[0];
  double Fxy = HessianTerms[2];
  double Fyy = HessianTerms[1];
  double curvature = (Fx * Fx * Fyy - 2.0 * Fx * Fy * Fxy + Fy * Fy * Fxx) /
                     (std::pow(Fx * Fx + Fy * Fy, 1.5));
  PU_parabola.coeff() = 0.5 * curvature;

  return PU_parabola;
}

}  // namespace IRL2D
