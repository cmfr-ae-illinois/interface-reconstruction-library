// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>
#include <iomanip>

#include "gtest/gtest.h"

#include "irl/quadratic_reconstruction/gauss_legendre_integrator.h"

namespace {

using namespace IRL;

TEST(GaussLegendreIntegrator, Sine1D) {
  const double a = 0.0;
  const double b = 10.0;

  // Define the functor.
  auto sine1d = [](const double t) { return std::sin(t); };
  auto sine1d_int = [](const double a, const double b) {
    return -std::cos(b) + std::cos(a);
  };

  // Define the integrator.
  IRL::GaussLegendreIntegrator<double, 1> integrator(50);

  // Integrate.
  const double integral = integrator.integrate(sine1d, a, b);
  const double exact_integral = sine1d_int(a, b);

  std::cout << "int(sin(t)) (numerical) = " << std::setprecision(16) << integral
            << std::endl;
  std::cout << "int(sin(t)) (exact)     = " << std::setprecision(16)
            << exact_integral << std::endl;

  EXPECT_NEAR(integral, exact_integral, 100.0 * DBL_EPSILON);
}

TEST(GaussLegendreIntegrator, Sine2D) {
  const double a = 0.0;
  const double b = 10.0;

  // Define the functor.
  auto sine2d = [](const double x, const double y) {
    return std::sin(x) * std::sin(y);
  };
  auto sine2d_int = [](const double a, const double b) {
    return (std::cos(b) - std::cos(a)) * (std::cos(b) - std::cos(a));
  };

  // Define the integrator
  IRL::GaussLegendreIntegrator<double, 2> integrator(50);

  // Integrate
  const double integral = integrator.integrate(sine2d, a, a, b, b);
  const double exact_integral = sine2d_int(a, b);

  std::cout << "int(sin(x)*sin(y)) (numerical)     = " << std::setprecision(16)
            << integral << std::endl;
  std::cout << "int(sin(x)*sin(y)) (exact)     = " << std::setprecision(16)
            << exact_integral << std::endl;

  EXPECT_NEAR(integral, exact_integral, 100.0 * DBL_EPSILON);
}

TEST(GaussLegendreIntegrator, Gaussian2D) {
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

  // Define the integrator
  IRL::GaussLegendreIntegrator<double, 2> integrator(50);

  // Integrate
  const double integral = integrator.integrate(exp2d, a, a, b, b);
  const double exact_integral = exp2d_int(a, b);

  std::cout << "int(exp(-x^2-y^2)) (numerical)     = " << std::setprecision(16)
            << integral << std::endl;
  std::cout << "int(exp(-x^2-y^2)) (exact)     = " << std::setprecision(16)
            << exact_integral << std::endl;

  EXPECT_NEAR(integral, exact_integral, 100.0 * DBL_EPSILON);
}

}  // namespace
