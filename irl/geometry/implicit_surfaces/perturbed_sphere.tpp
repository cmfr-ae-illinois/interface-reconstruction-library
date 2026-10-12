// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_IMPLICIT_SURFACES_PERTURBED_SPHERE_TPP_
#define IRL_GEOMETRY_IMPLICIT_SURFACES_PERTURBED_SPHERE_TPP_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"

namespace IRL {

template <class ScalarType, size_t MAX_REFINE_LEVEL>
class PerturbedSphere
    : public GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL> {
 public:
  using Vec3 = Eigen::Matrix<ScalarType, 3, 1>;
  using Mat3 = Eigen::Matrix<ScalarType, 3, 3>;

  PerturbedSphere()
      : xc(ScalarType(0.0)),
        yc(ScalarType(0.0)),
        zc(ScalarType(0.0)),
        R(ScalarType(0.8)),
        L(6),
        sigma0(ScalarType(5.0e-4)),
        fd_step(ScalarType(1.0e-8)) {
    generateCoefficients(0);
  }

  PerturbedSphere(const ScalarType& x_center, const ScalarType& y_center,
                  const ScalarType& z_center, const ScalarType& radius,
                  const int& order,
                  const ScalarType& variance = ScalarType(5.0e-4),
                  const unsigned int seed = 0)
      : xc(x_center),
        yc(y_center),
        zc(z_center),
        R(radius),
        L(order),
        sigma0(variance),
        fd_step(ScalarType(1.0e-5) * radius) {
    generateCoefficients(seed);
  }

  ScalarType F(const ScalarType& x, const ScalarType& y,
               const ScalarType& z) const override {
    const ScalarType X = x - xc;
    const ScalarType Y = y - yc;
    const ScalarType Z = z - zc;

    const ScalarType r2 = X * X + Y * Y + Z * Z;
    const ScalarType r = std::sqrt(r2);

    // Spherical coordinates are undefined exactly at the center
    if (r <= std::numeric_limits<ScalarType>::epsilon()) {
      return -R * R;
    }

    // converting to spherical coordiantes
    ScalarType cos_theta = Z / r;
    cos_theta =
        std::max(ScalarType(-1.0), std::min(ScalarType(1.0), cos_theta));
    const ScalarType theta = std::acos(cos_theta);  // polar angle
    ScalarType phi = std::atan2(Y, X);              // azimuthal angle
    const ScalarType pi = std::acos(ScalarType(-1.0));
    if (phi < ScalarType(0.0)) {
      phi += ScalarType(2.0) * pi;
    }

    // argument for Legendre polynomials in spherical harmonics
    const ScalarType legendre_x = std::cos(theta);

    // term with l=0 and m=0
    ScalarType R3 = R * R * R;

    size_t coeff_index = 0;

    // looping over all perturbed modes
    for (int l = 1; l <= L; ++l) {
      for (int m = -l; m <= l; ++m) {
        // to ensure continuity at poles
        if (std::abs(m) == 1) {
          continue;
        }

        const int abs_m = std::abs(m);

        // -----------------------------------------------------------
        // Associated Legendre polynomial P_l^{|m|}(cos(theta))
        // This includes Condon-Shortley phase
        ScalarType Pmm = ScalarType(1.0);

        if (abs_m > 0) {
          const ScalarType sin_theta = std::sqrt(std::max(
              ScalarType(0.0), ScalarType(1.0) - legendre_x * legendre_x));
          ScalarType factor = ScalarType(1.0);
          // P_m^m(x) = (-1)^m (2m-1)!!(1-x^2)^(m/2)
          for (int i = 1; i <= abs_m; ++i) {
            Pmm *= -factor * sin_theta;
            factor += ScalarType(2.0);
          }
        }

        ScalarType Plm = ScalarType(0.0);

        if (l == abs_m) {
          Plm = Pmm;
        } else {
          // P_{m+1}^m
          ScalarType Pmmp1 = legendre_x * ScalarType(2 * abs_m + 1) * Pmm;
          if (l == abs_m + 1) {
            Plm = Pmmp1;
          } else {
            ScalarType P_prev2 = Pmm;
            ScalarType P_prev1 = Pmmp1;
            for (int ll = abs_m + 2; ll <= l; ++ll) {
              const ScalarType Pll =
                  (ScalarType(2 * ll - 1) * legendre_x * P_prev1 -
                   ScalarType(ll + abs_m - 1) * P_prev2) /
                  ScalarType(ll - abs_m);

              P_prev2 = P_prev1;
              P_prev1 = Pll;
            }
            Plm = P_prev1;
          }
        }
        // -----------------------------------------------------------

        // normalization factors
        const ScalarType factorial_ratio =
            std::tgamma(ScalarType(l - abs_m + 1)) /
            std::tgamma(ScalarType(l + abs_m + 1));

        const ScalarType normalization = std::sqrt(
            ScalarType(2 * l + 1) / (ScalarType(4.0) * pi) * factorial_ratio);

        // Real tesseral spherical harmonic
        ScalarType Ylm;
        if (m == 0) {
          Ylm = normalization * Plm;
        } else if (m > 0) {
          Ylm = std::sqrt(ScalarType(2.0)) * normalization * Plm *
                std::cos(ScalarType(m) * phi);
        } else {
          Ylm = std::sqrt(ScalarType(2.0)) * normalization * Plm *
                std::sin(ScalarType(abs_m) * phi);
        }

        // accumulating mode contributions
        R3 += coefficients[coeff_index] * Ylm;
        ++coeff_index;
      }
    }

    const ScalarType Rlocal = std::cbrt(R3);

    return r2 - Rlocal * Rlocal;
  }

  // finite differences for now
  Vec3 gradF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    const ScalarType Fx = (F(x + fd_step, y, z) - F(x - fd_step, y, z)) /
                          (ScalarType(2.0) * fd_step);

    const ScalarType Fy = (F(x, y + fd_step, z) - F(x, y - fd_step, z)) /
                          (ScalarType(2.0) * fd_step);

    const ScalarType Fz = (F(x, y, z + fd_step) - F(x, y, z - fd_step)) /
                          (ScalarType(2.0) * fd_step);

    return Vec3(Fx, Fy, Fz);
  }

  Mat3 hessF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    const ScalarType Fxx =
        (F(x + fd_step, y, z) - ScalarType(2.0) * F(x, y, z) +
         F(x - fd_step, y, z)) /
        (fd_step * fd_step);

    const ScalarType Fyy =
        (F(x, y + fd_step, z) - ScalarType(2.0) * F(x, y, z) +
         F(x, y - fd_step, z)) /
        (fd_step * fd_step);

    const ScalarType Fzz =
        (F(x, y, z + fd_step) - ScalarType(2.0) * F(x, y, z) +
         F(x, y, z - fd_step)) /
        (fd_step * fd_step);

    const ScalarType Fxy =
        (F(x + fd_step, y + fd_step, z) - F(x + fd_step, y - fd_step, z) -
         F(x - fd_step, y + fd_step, z) + F(x - fd_step, y - fd_step, z)) /
        (ScalarType(4.0) * fd_step * fd_step);

    const ScalarType Fxz =
        (F(x + fd_step, y, z + fd_step) - F(x + fd_step, y, z - fd_step) -
         F(x - fd_step, y, z + fd_step) + F(x - fd_step, y, z - fd_step)) /
        (ScalarType(4.0) * fd_step * fd_step);

    const ScalarType Fyz =
        (F(x, y + fd_step, z + fd_step) - F(x, y + fd_step, z - fd_step) -
         F(x, y - fd_step, z + fd_step) + F(x, y - fd_step, z - fd_step)) /
        (ScalarType(4.0) * fd_step * fd_step);

    Mat3 H;

    H << Fxx, Fxy, Fxz, Fxy, Fyy, Fyz, Fxz, Fyz, Fzz;

    return H;
  }

  bool hasVolume() const override { return true; }

  ScalarType volume() const override {
    const ScalarType pi = std::acos(ScalarType(-1.0));

    return ScalarType(4.0 / 3.0) * pi * R * R * R;
  }

 private:
  void generateCoefficients(const unsigned int seed) {
    coefficients.clear();

    std::mt19937 generator(seed);  // 32-bit random number generator

    std::uniform_real_distribution<double> uniform(
        0.0, 1.0);  // uniform distribution on [0,1]

    const ScalarType pi = std::acos(ScalarType(-1.0));

    for (int l = 1; l <= L; ++l) {
      for (int m = -l; m <= l; ++m) {
        // c_l^{+/-1} = 0
        if (std::abs(m) == 1) continue;

        double gamma1 = uniform(generator);
        const double gamma2 = uniform(generator);

        // Avoid log(0)
        gamma1 = std::max(gamma1, std::numeric_limits<double>::min());

        const double gaussian =
            std::sqrt(-2.0 * std::log(gamma1)) *
            std::cos(2.0 * static_cast<double>(pi) * gamma2);

        const ScalarType c_lm = std::sqrt(sigma0) * ScalarType(gaussian);

        coefficients.push_back(c_lm);
      }
    }
  }

 private:
  ScalarType xc;
  ScalarType yc;
  ScalarType zc;
  ScalarType R;       // Base radius
  int L;              // Maximum spherical-harmonic degree
  ScalarType sigma0;  // Variance of random coefficients

  // perturbation coefficients
  std::vector<ScalarType> coefficients;

  // Finite-difference step for gradient/Hessian.
  ScalarType fd_step;
};

}  // namespace IRL

#endif  // IRL_GEOMETRY_IMPLICIT_SURFACES_PERTURBED_SPHERE_TPP_