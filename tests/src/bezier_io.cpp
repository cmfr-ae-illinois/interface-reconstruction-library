// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid_parametrized_surface.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include "gtest/gtest.h"
#include "tests/src/basic_mesh.h"

namespace {

using namespace IRL;

TEST(BezierIO, ParaboloidInCube) {
  using VolumeMomentsAndSuface =
      AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

  AlignedParaboloid aligned_paraboloid({1.0, 1.0});
  Pt datum(0.25, 0.0, 0.0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Paraboloid paraboloid(datum, frame, aligned_paraboloid.a(),
                        aligned_paraboloid.b());

  RectangularCuboid cube = RectangularCuboid::fromBoundingPts(
      Pt(-0.5, -0.5, -0.5), Pt(0.5, 0.5, 0.5));

  // Compute moments and return parametrized surface
  const auto surface_and_moments =
      getVolumeMoments<VolumeMomentsAndSuface>(cube, paraboloid);
  auto surface = surface_and_moments.getSurface();
  //   surface.setLengthScale(0.01);
  //   auto triangulated_surface = surface.triangulate();
  //   triangulated_surface.write("bezier_io_test1_triangulated_surface");

  auto bezier_quadratic_surface = surface.getQuadraticBezierTriangleApprox();
  bezier_quadratic_surface.write("bezier_io_test1_bezier_quadratic_surface");

  auto bezier_cubic_surface = surface.getCubicBezierTriangleApprox();
  bezier_cubic_surface.write("bezier_io_test1_bezier_cubic_surface");

  SUCCEED();
}

TEST(BezierIO, ParaboloidAndPlaneInCubes) {
  using VolumeMomentsAndSuface =
      AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

  AlignedParaboloid aligned_paraboloid({1.0, 1.0});
  Pt datum(0.25, 0.0, 0.0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Paraboloid paraboloid(datum, frame, aligned_paraboloid.a(),
                        aligned_paraboloid.b());

  RectangularCuboid cube0 = RectangularCuboid::fromBoundingPts(
      Pt(-0.5, -0.5, -0.5), Pt(0.5, 0.5, 0.5));
  RectangularCuboid cube1 = RectangularCuboid::fromBoundingPts(
      Pt(0.5, -0.5, -0.5), Pt(1.5, 0.5, 0.5));
  RectangularCuboid cube2 = RectangularCuboid::fromBoundingPts(
      Pt(1.5, -0.5, -0.5), Pt(2.5, 0.5, 0.5));

  // Compute moments and return parametrized surface
  auto surface_and_moments1 =
      getVolumeMoments<VolumeMomentsAndSuface>(cube0, paraboloid);
  auto& surface1 = surface_and_moments1.getSurface();
  auto surface_output = surface1.getQuadraticBezierTriangleApprox();
  auto surface_and_moments2 =
      getVolumeMoments<VolumeMomentsAndSuface>(cube1, paraboloid);
  surface_output.addSurface(
      surface_and_moments2.getSurface().getQuadraticBezierTriangleApprox());

  const auto planar_separator =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 1.0), 0.0));
  const auto polygon1 = getPlanePolygonFromReconstruction<Polygon>(
      cube1, planar_separator, planar_separator[0]);
  const auto polygon2 = getPlanePolygonFromReconstruction<Polygon>(
      cube2, planar_separator, planar_separator[0]);
  surface_output.addPolygons(std::vector<Polygon>({polygon1, polygon2}));
  surface_output.write("bezier_io_test2", true);

  SUCCEED();
}

}  // namespace