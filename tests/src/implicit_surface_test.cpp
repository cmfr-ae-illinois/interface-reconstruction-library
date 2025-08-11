// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <Eigen/Dense>

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/implicit_surfaces.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

// comparing gradient and hessian expressions against finite differences
TEST(ImplicitSurface, ImplicitSurface) {
  using BaseSurface = GeneralImplicitSurface<double, 0>;
  double x = 0.1, y = 0.2, z = 0.3;

  // vector of surfaces
  std::vector<std::unique_ptr<BaseSurface>> surface_list;
  surface_list.emplace_back(std::make_unique<Sphere<double, 0>>());
  surface_list.emplace_back(std::make_unique<Ellipsoid<double, 0>>());
  surface_list.emplace_back(std::make_unique<Genus<double, 0>>());
  surface_list.emplace_back(std::make_unique<Orthocircle<double, 0>>());

  const double h = 1e-4;
  const double tol = 1e-6;
  for (const auto& surface : surface_list) {
    // gradient terms
    auto Fx_fd = (surface->F(x + h, y, z) - surface->F(x - h, y, z)) / (2 * h);
    auto Fy_fd = (surface->F(x, y + h, z) - surface->F(x, y - h, z)) / (2 * h);
    auto Fz_fd = (surface->F(x, y, z + h) - surface->F(x, y, z - h)) / (2 * h);
    auto grad = surface->gradF(x, y, z);
    EXPECT_NEAR(grad(0), Fx_fd, tol);
    EXPECT_NEAR(grad(1), Fy_fd, tol);
    EXPECT_NEAR(grad(2), Fz_fd, tol);

    // hessian terms
    auto Fxx_fd = (surface->F(x + h, y, z) - 2 * surface->F(x, y, z) +
                   surface->F(x - h, y, z)) /
                  (h * h);
    auto Fyy_fd = (surface->F(x, y + h, z) - 2 * surface->F(x, y, z) +
                   surface->F(x, y - h, z)) /
                  (h * h);
    auto Fzz_fd = (surface->F(x, y, z + h) - 2 * surface->F(x, y, z) +
                   surface->F(x, y, z - h)) /
                  (h * h);

    auto Fxy_fd = (surface->F(x + h, y + h, z) - surface->F(x + h, y - h, z) -
                   surface->F(x - h, y + h, z) + surface->F(x - h, y - h, z)) /
                  (4 * h * h);

    auto Fxz_fd = (surface->F(x + h, y, z + h) - surface->F(x + h, y, z - h) -
                   surface->F(x - h, y, z + h) + surface->F(x - h, y, z - h)) /
                  (4 * h * h);

    auto Fyz_fd = (surface->F(x, y + h, z + h) - surface->F(x, y + h, z - h) -
                   surface->F(x, y - h, z + h) + surface->F(x, y - h, z - h)) /
                  (4 * h * h);

    auto hess = surface->hessF(x, y, z);

    EXPECT_NEAR(hess(0, 0), Fxx_fd, tol);
    EXPECT_NEAR(hess(1, 1), Fyy_fd, tol);
    EXPECT_NEAR(hess(2, 2), Fzz_fd, tol);
    EXPECT_NEAR(hess(0, 1), Fxy_fd, tol);
    EXPECT_NEAR(hess(0, 2), Fxz_fd, tol);
    EXPECT_NEAR(hess(1, 2), Fyz_fd, tol);
  }
}

}  // namespace
