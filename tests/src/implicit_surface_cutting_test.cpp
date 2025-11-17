// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <Eigen/Dense>

#include "examples/variant_advector/basic_mesh.h"
#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/implicit_surfaces.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(ImplicitSurfaceCutting, EllipsoidVolume) {
  constexpr size_t max_refine = 5;
  Ellipsoid<double, max_refine> ellipsoid;
  double volume_exact = ellipsoid.volume();

  constexpr int nx = 32;
  constexpr int ny = 32;
  constexpr int nz = 32;
  const double xmin = -0.5, ymin = -0.5, zmin = -0.5;
  const double xmax = 0.5, ymax = 0.5, zmax = 0.5;

  double dx = (xmax - xmin) / nx;
  double dy = (ymax - ymin) / ny;
  double dz = (zmax - zmin) / nz;

  double volume_cutting = 0.0;

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        Pt x0(xmin + i * dx, ymin + j * dy, zmin + k * dz);
        Pt x1(xmin + (i + 1) * dx, ymin + (j + 1) * dy, zmin + (k + 1) * dz);
        RectangularCuboid cell = RectangularCuboid::fromBoundingPts(x0, x1);
        ImplicitSurfaceCutter<Ellipsoid<double, max_refine>, Volume> cutter(
            ellipsoid, cell);
        double vol = cutter.computeVolumeMoments();
        volume_cutting += vol;
      }
    }
  }

  EXPECT_NEAR(volume_exact, volume_cutting, 1e-12);
}

TEST(ImplicitSurfaceCutting, SphereMoments) {
  constexpr size_t max_refine = 5;
  Sphere<double, max_refine> sphere(0.1, 0.2, 0.3, 0.15);
  double volume_exact = sphere.volume();

  constexpr int nx = 32;
  constexpr int ny = 32;
  constexpr int nz = 32;
  const double xmin = -0.5, ymin = -0.5, zmin = -0.5;
  const double xmax = 0.5, ymax = 0.5, zmax = 0.5;
  double dx = (xmax - xmin) / nx;
  double dy = (ymax - ymin) / ny;
  double dz = (zmax - zmin) / nz;

  VolumeMoments moments_cutting;

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        Pt x0(xmin + i * dx, ymin + j * dy, zmin + k * dz);
        Pt x1(xmin + (i + 1) * dx, ymin + (j + 1) * dy, zmin + (k + 1) * dz);
        RectangularCuboid cell = RectangularCuboid::fromBoundingPts(x0, x1);
        ImplicitSurfaceCutter<Sphere<double, max_refine>, VolumeMoments> cutter(
            sphere, cell);
        moments_cutting += cutter.computeVolumeMoments();
      }
    }
  }
  moments_cutting.normalizeByVolume();
  EXPECT_NEAR(volume_exact, moments_cutting.volume(), 1e-12);
  EXPECT_NEAR(0.1, moments_cutting.centroid()[0], 1e-12);
  EXPECT_NEAR(0.2, moments_cutting.centroid()[1], 1e-12);
  EXPECT_NEAR(0.3, moments_cutting.centroid()[2], 1e-12);
}

TEST(ImplicitSurfaceCutting, SphereSurfaceMoments) {
  constexpr size_t max_refine = 1;
  Sphere<double, max_refine> sphere(0.1, 0.2, 0.3, 0.15);
  double area_exact = sphere.surfaceArea();

  constexpr int nx = 32;
  constexpr int ny = 32;
  constexpr int nz = 32;
  const double xmin = -0.5, ymin = -0.5, zmin = -0.5;
  const double xmax = 0.5, ymax = 0.5, zmax = 0.5;
  double dx = (xmax - xmin) / nx;
  double dy = (ymax - ymin) / ny;
  double dz = (zmax - zmin) / nz;

  GeneralSurfaceMoments3D<1> moments_cutting;

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        Pt x0(xmin + i * dx, ymin + j * dy, zmin + k * dz);
        Pt x1(xmin + (i + 1) * dx, ymin + (j + 1) * dy, zmin + (k + 1) * dz);
        RectangularCuboid cell = RectangularCuboid::fromBoundingPts(x0, x1);
        ImplicitSurfaceCutter<Sphere<double, max_refine>, VolumeMoments> cutter(
            sphere, cell);
        moments_cutting += cutter.computeSurfaceMoments<1>();
      }
    }
  }
  EXPECT_NEAR(area_exact, moments_cutting[0], 1e-6);
  EXPECT_NEAR(0.1, moments_cutting[1] / moments_cutting[0], 1e-6);
  EXPECT_NEAR(0.2, moments_cutting[2] / moments_cutting[0], 1e-6);
  EXPECT_NEAR(0.3, moments_cutting[3] / moments_cutting[0], 1e-6);
}

}  // namespace