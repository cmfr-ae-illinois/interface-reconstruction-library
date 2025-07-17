// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "gtest/gtest.h"

#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "external/NumericalIntegration/NumericalIntegration.h"

namespace {

using namespace IRL;

TEST(Integrator, Sine1D) {
  const std::size_t max_nsubdivisions = 128;
  const double epsabs = 10.0 * DBL_EPSILON;
  const double epsrel = 10.0 * DBL_EPSILON;
  const double a = 0.0;
  const double b = 10.0;

  // Define the functor.
  auto sine1d = [](const double t) { return std::sin(t); };
  auto sine1d_int = [](const double a, const double b) {
    return -std::cos(b) + std::cos(a);
  };

  // Define the integrator.
  Eigen::Integrator<double> integrator(max_nsubdivisions);

  // Integrate.
  const double integral = integrator.quadratureAdaptive(
      sine1d, a, b, epsabs, epsrel, Eigen::Integrator<double>::GaussKronrod15);
  const double exact_integral = sine1d_int(a, b);

  std::cout << "int(sin(t)) (numerical) = " << std::setprecision(16) << integral
            << std::endl;
  std::cout << "int(sin(t)) (exact)     = " << std::setprecision(16)
            << exact_integral << std::endl;

  EXPECT_NEAR(integral, exact_integral,
              std::max(epsabs, epsrel * exact_integral));
}

TEST(Integrator, Sine2D) {
  const std::size_t max_nsubdivisions = 128;
  const double epsabs = 10.0 * DBL_EPSILON;
  const double epsrel = 10.0 * DBL_EPSILON;
  const double a = 0.0;
  const double b = 10.0;

  // Define the functor.
  auto sine2d = [](const double x, const double y) {
    return std::sin(x) * std::sin(y);
  };
  auto sine2d_int = [](const double a, const double b) {
    return (std::cos(b) - std::cos(a)) * (std::cos(b) - std::cos(a));
  };
  const double exact_integral = sine2d_int(a, b);
  double integral = 0.0;

  for (std::size_t i = 1; i <= max_nsubdivisions; i *= 4) {
    // Define the integrator.
    Eigen::Integrator<double, 2> integrator(i);

    // Integrate.
    integral = integrator.quadratureAdaptive(
        sine2d, a, b, epsabs, epsrel,
        Eigen::Integrator<double, 2>::GaussKronrod15);

    std::cout << "int(sin(x)*sin(y)) (numerical, sub=" << i
              << ") = " << std::setprecision(16) << integral
              << " / error = " << std::setprecision(16)
              << std::abs(integral - exact_integral) << std::endl;
  }
  std::cout << "int(sin(x)*sin(y)) (exact)     = " << std::setprecision(16)
            << exact_integral << std::endl;

  EXPECT_NEAR(integral, exact_integral,
              std::max(epsabs, epsrel * exact_integral));
}

TEST(Integrator, Gaussian2D) {
  const std::size_t max_nsubdivisions = std::pow(2, 8);
  const double epsabs = 10.0 * DBL_EPSILON;
  const double epsrel = 10.0 * DBL_EPSILON;
  const double a = 0.0;
  const double b = 5.0;

  // Define the functor.
  auto exp2d = [](const double x, const double y) {
    return std::exp(-x * x) * std::exp(-y * y);
  };
  auto exp2d_int = [](const double a, const double b) {
    return 0.25 * M_PI * (std::erf(b) - std::erf(a)) *
           (std::erf(b) - std::erf(a));
  };
  const double exact_integral = exp2d_int(a, b);
  double integral = 0.0;

  for (std::size_t i = 1; i <= max_nsubdivisions; i *= 4) {
    // Define the integrator.
    Eigen::Integrator<double, 2> integrator(i);

    // Integrate.
    integral = integrator.quadratureAdaptive(
        exp2d, a, b, epsabs, epsrel,
        Eigen::Integrator<double, 2>::GaussKronrod15);

    std::cout << "int(exp(-x^2-y^2)) (numerical, sub=" << i
              << ") = " << std::setprecision(16) << integral
              << " / error = " << std::setprecision(16)
              << std::abs(integral - exact_integral) << std::endl;
  }
  std::cout << "int(exp(-x^2-y^2)) (exact)     = " << std::setprecision(16)
            << exact_integral << std::endl;

  EXPECT_NEAR(integral, exact_integral,
              std::max(epsabs, epsrel * exact_integral));
}

}  // namespace
