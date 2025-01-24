// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// #define VALDEBUG

#include "irl/geometry/general/new_pt_calculation_functors.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/rotations.h"
// #include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
// #include "irl/paraboloid_reconstruction/hessian_paraboloid.h"

#include <sys/time.h>
#include <cmath>
#include <random>

#include "gtest/gtest.h"

#include "irl/data_structures/small_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"

#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection.h"
#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection_amr.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron_quadratic.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_quadratic.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
// #include
// "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/moments/general_moments.h"

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "irl/planar_reconstruction/planar_separator.h"
#include "tests/src/basic_mesh.h"
#include "tests/src/vtk.h"

namespace {

using namespace IRL;

TEST(CylinderIntersection, SISCPaperFig5) {
  using VolumeAndSuface =
      AddSurfaceOutput<VolumeMoments, CylinderParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({2.0, 1.2});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  auto cubes = std::array<RectangularCuboid, 5>(
      {RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(0.8, 0.8, 0.8)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(1.0, 1.0, 1.0)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(1.2, 1.2, 1.2)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, -0.5, 0.0),
                                          Pt(1.0, 0.5, 1.0)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, -1.0, 0.0),
                                          Pt(1.0, 1.0, 1.0))});
  std::array<std::string, 5> surface_filenames(
      {"surface_a", "surface_b", "surface_c", "surface_d", "surface_e"});
  std::array<std::string, 5> clipped_faces_filenames(
      {"_cube_a", "_cube_b", "_cube_c", "_cube_d", "_cube_e"});
  // Compute moments and return parametrized surface
  for (UnsignedIndex_t i = 0; i < 5; i++) {
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(cubes[i], cylinder);
    std::cout << "the " << i << "th volume is :" << std::setprecision(20)
              << temp_surface_and_moments.getMoments().volume().volume()
              << std::endl;
    std::cout << "the " << i << "th centroid is :"
              << temp_surface_and_moments.getMoments().centroid() << std::endl;
    auto& centroid = temp_surface_and_moments.getMoments().centroid();
    centroid /= temp_surface_and_moments.getMoments().volume().volume();
    std::cout << "the " << i << "th center of mass is :" << centroid
              << std::endl;
    auto temp_moments = getVolumeMoments<Volume>(cubes[i], cylinder);
    EXPECT_EQ(temp_surface_and_moments.getMoments().volume().volume(),
              temp_moments.volume());
    auto temp_param_surface = temp_surface_and_moments.getSurface();
    auto temp_tri_surface = temp_param_surface.triangulate(0.1);
    temp_tri_surface.write(surface_filenames[i]);
  }
}

TEST(HyperCylinderIntersection, SISCPaperFig5) {
  using VolumeAndSuface =
      AddSurfaceOutput<Volume, CylinderParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({-2.0, 0.5});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  auto cubes = std::array<RectangularCuboid, 5>(
      {RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(0.8, 0.8, 0.8)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(1.0, 1.0, 1.0)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(1.2, 1.2, 1.2)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, -0.5, 0.0),
                                          Pt(1.0, 0.5, 1.0)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, -1.0, 0.0),
                                          Pt(1.0, 1.0, 1.0))});
  std::array<std::string, 5> surface_filenames(
      {"surface_a", "surface_b", "surface_c", "surface_d", "surface_e"});
  std::array<std::string, 5> clipped_faces_filenames(
      {"_cube_a", "_cube_b", "_cube_c", "_cube_d", "_cube_e"});
  // Compute moments and return parametrized surface
  for (UnsignedIndex_t i = 0; i < 5; i++) {
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(cubes[i], cylinder);
    std::cout << "the " << i << "th volum is :" << std::setprecision(20)
              << temp_surface_and_moments.getMoments().volume() << std::endl;
    auto temp_moments = getVolumeMoments<Volume>(cubes[i], cylinder);
    EXPECT_EQ(temp_surface_and_moments.getMoments().volume(),
              temp_moments.volume());
    auto temp_param_surface = temp_surface_and_moments.getSurface();
    auto temp_tri_surface = temp_param_surface.triangulate(0.1);
    temp_tri_surface.write(surface_filenames[i]);
  }
}

TEST(CylinderIntersection, Debug) {
  using VolumeMomentsAndSuface =
      AddSurfaceOutput<VolumeMoments, CylinderParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({1.0, 1.0});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  std::array<Pt, 8> vertex_list{{Pt(-1.0, -0.2, 0.0), Pt(1.0, 2.2, 0.0),
                                 Pt(1.0, 0.2, 2.0), Pt(3.0, 2.6, 2.0),
                                 Pt(-2.0, -0.2, 0.0), Pt(0.0, 2.2, 0.0),
                                 Pt(0.0, 0.2, 2.0), Pt(2.0, 2.6, 2.0)}};

  std::array<std::array<UnsignedIndex_t, 4>, 6> face_mapping{{{0, 1, 3, 2},
                                                              {2, 3, 7, 6},
                                                              {3, 1, 5, 7},
                                                              {1, 0, 4, 5},
                                                              {0, 2, 6, 4},
                                                              {4, 6, 7, 5}}};

  PolyhedronConnectivity connectivity(face_mapping);
  GeneralPolyhedron prism(vertex_list, &connectivity);
  GeneralPolyhedron prism_local_frame(vertex_list, &connectivity);
  std::string surface_filename = "surface_debug";

  double th_volume_cylinder = M_PI / 4.0;
  double th_volume_triangle = 0.1;
  double th_volume_total = th_volume_cylinder + th_volume_triangle;
  std::array<double, 2> theoretical_centroid = {
      ((-0.2 / 3.0) * th_volume_triangle +
       (4.0 / (3.0 * M_PI)) * th_volume_cylinder) /
          th_volume_total,
      ((1.0 / 3.0) * th_volume_triangle +
       (4.0 / (3.0 * M_PI)) * th_volume_cylinder) /
          th_volume_total};  // the x component of the centroid is not trivial

  // Compute moments and return parametrized surface
  auto temp_surface_and_moments =
      getVolumeMoments<VolumeMomentsAndSuface>(prism, cylinder);
  std::cout << "the volume is   :" << std::setprecision(20)
            << temp_surface_and_moments.getMoments().volume() << std::endl;
  std::cout << "expected volume :" << th_volume_total << std::endl;
  auto& centroid = temp_surface_and_moments.getMoments().centroid();
  centroid /= temp_surface_and_moments.getMoments().volume();
  std::cout << "the normalize centroid is :" << centroid << std::endl;
  std::cout << "expected centroid         :( ?, " << theoretical_centroid[0]
            << ", " << theoretical_centroid[1] << ")\n";
  auto temp_param_surface = temp_surface_and_moments.getSurface();
  auto temp_tri_surface = temp_param_surface.triangulate(0.1);
  temp_tri_surface.write(surface_filename);

  EXPECT_EQ(temp_surface_and_moments.getMoments().volume(), th_volume_total);
}

TEST(HyperCylinderIntersection, Debug) {
  using VolumeAndSuface =
      AddSurfaceOutput<Volume, CylinderParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({-1.0, 0.0});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  std::array<Pt, 8> vertex_list{{Pt(0.0, 0.0, 0.0), Pt(0.0, 2.0, 0.0),
                                 Pt(0.0, 0.0, 2.0), Pt(0.0, 2.0, 2.0),
                                 Pt(-2.0, 0.0, 0.0), Pt(-2.0, 2.0, 0.0),
                                 Pt(-2.0, 0.0, 2.0), Pt(-2.0, 2.0, 2.0)}};

  std::array<std::array<UnsignedIndex_t, 4>, 6> face_mapping{{{0, 1, 3, 2},
                                                              {2, 3, 7, 6},
                                                              {3, 1, 5, 7},
                                                              {1, 0, 4, 5},
                                                              {0, 2, 6, 4},
                                                              {4, 6, 7, 5}}};

  PolyhedronConnectivity connectivity(face_mapping);
  GeneralPolyhedron prism(vertex_list, &connectivity);
  GeneralPolyhedron prism_local_frame(vertex_list, &connectivity);
  std::string surface_filename = "surface_debug";

  // Compute moments and return parametrized surface
  auto temp_surface_and_moments =
      getVolumeMoments<VolumeAndSuface>(prism, cylinder);
  std::cout << "expected volume          : " << std::setprecision(16)
            << M_PI / 4.0 + 0.1 << std::endl;
  std::cout << "irl volume               : " << std::setprecision(16)
            << temp_surface_and_moments.getMoments().volume()
            << std::setprecision(3) << " -- error: "
            << std::abs(temp_surface_and_moments.getMoments().volume() -
                        (M_PI / 4.0 + 0.1))
            << std::endl;
  auto temp_param_surface = temp_surface_and_moments.getSurface();
  auto temp_tri_surface = temp_param_surface.triangulate(0.005);
  temp_tri_surface.write(surface_filename);

  HalfEdgePolyhedronQuadratic<Pt> half_edge;
  prism_local_frame.setHalfEdgeVersion(&half_edge);
  auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

  // Calculate volume of unclipped dodecahedron using AMR
  const int nlevels = 8;
  auto amr_volume = intersectPolyhedronWithCylinderAMR<Volume>(
      &seg_half_edge, &half_edge, aligned_cylinder, nlevels, "test_cylinder");

  std::cout << "amr volume (" << nlevels
            << " levels)   : " << std::setprecision(16) << amr_volume
            << std::setprecision(3)
            << " -- error: " << std::abs(amr_volume - (M_PI / 4.0 + 0.1))
            << std::endl;
}

TEST(CylinderIntersection, Debug2) {
  using VolumeAndSuface =
      AddSurfaceOutput<Volume, CylinderParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({2.0, 1.2});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  std::array<Pt, 8> vertex_list{{Pt(1.5, -0.5, 0.0), Pt(1.5, 0.5, 0.0),
                                 Pt(1.0, -0.5, 1.0), Pt(1.0, 0.5, 1.0),
                                 Pt(0.5, -0.5, 0.0), Pt(0.5, 0.5, 0.0),
                                 Pt(0.0, -0.5, 1.0), Pt(0.0, 0.5, 1.0)}};

  std::array<std::array<UnsignedIndex_t, 4>, 6> face_mapping{{{0, 1, 3, 2},
                                                              {2, 3, 7, 6},
                                                              {3, 1, 5, 7},
                                                              {1, 0, 4, 5},
                                                              {0, 2, 6, 4},
                                                              {4, 6, 7, 5}}};

  PolyhedronConnectivity connectivity(face_mapping);
  GeneralPolyhedron prism(vertex_list, &connectivity);
  GeneralPolyhedron prism_local_frame(vertex_list, &connectivity);
  std::string surface_filename = "surface_debug";

  // Compute moments and return parametrized surface
  auto temp_surface_and_moments =
      getVolumeMoments<VolumeAndSuface>(prism, cylinder);
  auto just_volum = getVolumeMoments<Volume>(prism, cylinder);
  std::cout << "the volume is (with surface) :" << std::setprecision(20)
            << temp_surface_and_moments.getMoments().volume() << std::endl;
  std::cout << "the volume is (no   surface) :" << std::setprecision(20)
            << just_volum.volume() << std::endl;
  auto temp_param_surface = temp_surface_and_moments.getSurface();
  auto temp_tri_surface = temp_param_surface.triangulate(0.1);
  temp_tri_surface.write(surface_filename);

  EXPECT_EQ(temp_surface_and_moments.getMoments().volume(),
            just_volum.volume());
}

TEST(HyperCylinderIntersection, Debug2) {
  using VolumeAndSuface =
      AddSurfaceOutput<Volume, CylinderParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({-2.0, 0.5});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  std::array<Pt, 8> vertex_list{{Pt(1.5, -0.5, 0.0), Pt(1.5, 0.5, 0.0),
                                 Pt(1.0, -0.5, 1.0), Pt(1.0, 0.5, 1.0),
                                 Pt(0.5, -0.5, 0.0), Pt(0.5, 0.5, 0.0),
                                 Pt(0.0, -0.5, 1.0), Pt(0.0, 0.5, 1.0)}};

  std::array<std::array<UnsignedIndex_t, 4>, 6> face_mapping{{{0, 1, 3, 2},
                                                              {2, 3, 7, 6},
                                                              {3, 1, 5, 7},
                                                              {1, 0, 4, 5},
                                                              {0, 2, 6, 4},
                                                              {4, 6, 7, 5}}};

  PolyhedronConnectivity connectivity(face_mapping);
  GeneralPolyhedron prism(vertex_list, &connectivity);
  GeneralPolyhedron prism_local_frame(vertex_list, &connectivity);
  std::string surface_filename = "surface_debug";

  // Compute moments and return parametrized surface
  auto temp_surface_and_moments =
      getVolumeMoments<VolumeAndSuface>(prism, cylinder);
  auto just_volum = getVolumeMoments<Volume>(prism, cylinder);
  std::cout << "the volume is (with surface) :" << std::setprecision(20)
            << temp_surface_and_moments.getMoments().volume() << std::endl;
  std::cout << "the volume is (no   surface) :" << std::setprecision(20)
            << just_volum.volume() << std::endl;
  auto temp_param_surface = temp_surface_and_moments.getSurface();
  auto temp_tri_surface = temp_param_surface.triangulate(0.1);
  temp_tri_surface.write(surface_filename);

  EXPECT_EQ(temp_surface_and_moments.getMoments().volume(),
            just_volum.volume());
}

TEST(CylinderIntersection, DebugAMR) {
  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({4.0, 1.0});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  auto cube = RectangularCuboid::fromBoundingPts(Pt(-1.0, -1.0, -1.0),
                                                 Pt(1.0, 1.0, 1.0));

  // Compute moments and return parametrized surface
  auto temp_moments = getVolumeMoments<Volume>(cube, cylinder);
  std::cout << "expected volume          : " << std::setprecision(16) << M_PI
            << std::endl;
  std::cout << "irl volume               : " << std::setprecision(16)
            << temp_moments.volume() << std::setprecision(3)
            << " -- error: " << std::abs(temp_moments.volume() - M_PI)
            << std::endl;

  HalfEdgePolyhedronQuadratic<Pt> half_edge, dummy_half_edge;
  cube.setHalfEdgeVersion(&half_edge);
  cube.setHalfEdgeVersion(&dummy_half_edge);
  auto dummy_seg_half_edge = dummy_half_edge.generateSegmentedPolyhedron();
  auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

  // Calculate volume of unclipped dodecahedron using AMR
  int nlevels = 8;
  // This first call prints the clipped/unclipped triangles in
  // "test_cylinder_amr" files
  auto dummy_amr_volume = intersectPolyhedronWithCylinderAMR<Volume>(
      &dummy_seg_half_edge, &dummy_half_edge, aligned_cylinder, nlevels,
      "test_cylinder_amr");

  // This calculates the volume within machine zero
  nlevels = 18;
  auto amr_volume = intersectPolyhedronWithCylinderAMR<Volume>(
      &seg_half_edge, &half_edge, aligned_cylinder, nlevels);

  std::cout << "amr volume (" << nlevels
            << " levels)   : " << std::setprecision(16) << amr_volume
            << std::setprecision(3)
            << " -- error: " << std::abs(amr_volume - M_PI) << std::endl;

  EXPECT_EQ(temp_moments.volume(), M_PI);
}

}  // namespace
