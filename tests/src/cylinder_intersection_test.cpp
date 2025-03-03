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

#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection_amr.h"
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

double ExactM0TranslatingCube_cylinder(const double u) {
  const double TWO = 2.0;
  const double HALF = 1.0 / 2.0;
  if (u < 2) {
    double sqrt_1 = sqrt(u * (4.0 - u));
    double acos_1 = acos(1.0 - HALF * u);
    return (1.0/8.0) * ((u - 2.0) * sqrt_1 + 4.0 * acos_1);
  } else {
    double sqrt_1 = sqrt((6.0 - u) * (u - 2.0));
    double sqrt_2 = sqrt(u * (4.0 - u));
    double acos_1 = acos(2.0 - HALF * u);
    double acos_2 = acos(HALF * u - 1);
    return (1.0 / 8.0) * (4.0 * M_PI - 
      (u - 4.0) * sqrt_1 + (u - 2.0) * sqrt_2 - 4.0 * acos_1 - 4.0 * acos_2);
  }
}

double ExactM1yTranslatingCube_cylinder(const double u) {
  if (u < 2.0) {
    return (6.0 - u) * u * u / 48.0;
  } else {
    return (6.0 - u) * u / 8.0 - 2.0/3.0;
  }
}

double ExactM1zTranslatingCube_cylinder(const double u) {
  if (u < 2.0) {
    return pow((4.0 - u) * u, 3.0/2.0) / 24.0;
  } else {
    return (pow((4.0 - u) * u, 3.0 / 2.0) - 
            pow((2.0 - u) * (u - 6.0), 3.0 / 2.0)) / 24.0;
  }
}

using namespace IRL;

TEST(CylinderIntersection, SISCPaperFig5) {
  using VolumeAndSuface =
      AddSurfaceOutput<VolumeMoments, CylinderParametrizedSurfaceOutput>;

  const double INVSQRTWO = 1.0 / sqrt(2.0);
  const double normalization = 1.0 / (INVSQRTWO + 0.5);

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({1.0, 1.0});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  auto cubes = std::array<RectangularCuboid, 3>(
      {RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.5), Pt(1.0, 1.0, 1.5)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(1.0, 1.0, 1.0)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, -0.5), Pt(1.0, 1.0, 0.5))});
  std::array<std::string, 3> surface_filenames(
      {"surface_1", "surface_2", "surface_3"});
  std::array<std::string, 3> clipped_faces_filenames(
      {"_cube_1", "_cube_2", "_cube_3"});

  std::array<double, 3> offsets = {1., 2., 3.};


  std::array<HalfEdgePolyhedronQuadratic<Pt>, 3> half_edges;
  std::array<IRL::SegmentedHalfEdgePolyhedronQuadratic<IRL::FaceQuadratic<IRL::HalfEdgeQuadratic<IRL::VertexQuadratic<IRL::Pt>>>, IRL::VertexQuadratic<IRL::Pt>>, 3> seg_half_edges;
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    cubes[i].setHalfEdgeVersion(&(half_edges[i]));
    seg_half_edges[i] = half_edges[i].generateSegmentedPolyhedron();
  }
  const UnsignedIndex_t nlevels = 10;

  // Compute moments and return parametrized surface
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(cubes[i], cylinder);
    auto amr_moments =
        intersectPolyhedronWithCylinderAMR<VolumeMoments>(
            &(seg_half_edges[i]), &(half_edges[i]), aligned_cylinder, nlevels, clipped_faces_filenames[i]);
    auto centroid = temp_surface_and_moments.getMoments().centroid();
    double exact_volume = ExactM0TranslatingCube_cylinder(offsets[i]);
    double exact_My = ExactM1yTranslatingCube_cylinder(offsets[i]);
    double exact_Mz = ExactM1zTranslatingCube_cylinder(offsets[i]);

    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
                exact_volume, 1e-13);
    EXPECT_NEAR(centroid[0] / temp_surface_and_moments.getMoments().volume().volume(),
                0.5, 1e-13);
    // EXPECT_NEAR(centroid[1] / temp_surface_and_moments.getMoments().volume().volume(),
    //           amr_centroid[1], 1e-13);
    EXPECT_NEAR(centroid[1],
                exact_My, 1e-13);
    EXPECT_NEAR(centroid[2],
                exact_Mz, 1e-13);
    auto temp_param_surface = temp_surface_and_moments.getSurface();
    auto temp_tri_surface = temp_param_surface.triangulate(0.1);
    temp_tri_surface.write(surface_filenames[i]);
  }
}

TEST(CylinderIntersection, SISCPaperFig6) {

  AlignedCylinder aligned_cylinder({1.0, 1.0});  // DO NOT CHANGE
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(),
                        aligned_cylinder.r());

  const double VOLUME_MAX = sqrt(3.0) / 4.0 + M_PI / 6.0;
  const double M1Z_MAX = 1.0 / 3.0;
  const double M1Y_MAX = 11.0 / 24.0;

  //////////////////////////////// YOU CAN CHANGE THESE PARAMETERS
  int Ntests = 3001;  // Number of tests
  ////////////////////////////////

  double max_volume_error = 0.0, rms_volume_error = 0.0;
  // double max_surface_error = 0.0, rms_surface_error = 0.0;

  std::ofstream myfile;
  myfile.open("fig_cylinder.csv");
  myfile << "k,m0p,m0p_exact,m0p_error,m1y,m1y_exact,m1y_error,m1z,m1z_exact,m1z_error" << std::endl; // ,m1yp,m1yp_exact,m1yp_error
  myfile.close();

  for (UnsignedIndex_t i = 0; i < Ntests; i++) {
    // Create and translate unit cube
    double k = 1.5 - 1.5 * (static_cast<double>(i) /
               static_cast<double>(Ntests - 1));
    RectangularCuboid cube =
        RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, k - 0.5), Pt(1.0, 1.0, k + 0.5));

    // Compute moments and return parametrized surface
    auto our_moments =
        getVolumeMoments<VolumeMoments>(cube, cylinder);

    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;
    std::cout << "Test " << i + 1 << "/" << Ntests << std::endl;

    // Compute exact value for verification purposes
    auto exact_volume = ExactM0TranslatingCube_cylinder(3.0 * static_cast<double>(i) /
    static_cast<double>(Ntests - 1)) / VOLUME_MAX;
    auto exact_m1y = ExactM1yTranslatingCube_cylinder(3.0 * static_cast<double>(i) /
    static_cast<double>(Ntests - 1)) / M1Y_MAX;
    auto exact_m1z = ExactM1zTranslatingCube_cylinder(3.0 * static_cast<double>(i) /
    static_cast<double>(Ntests - 1)) / M1Z_MAX;
    auto exact_centroid =
        Pt(0.5*safelyEpsilon(exact_volume*VOLUME_MAX), exact_m1y, exact_m1z);
    auto our_centroid =
        Pt(our_moments.centroid()[0], our_moments.centroid()[1] / M1Y_MAX, our_moments.centroid()[2] / M1Z_MAX);
    // std::cout << std::setprecision(20)
    //           << "Surface EXACT  = " << exact_surface_area << std::endl;
    // std::cout << std::setprecision(20)
    //           << "Surface IRL    = " << our_surface_area << std::endl;
    std::cout << std::setprecision(20)
              << "Vfrac unclipped EX  = " << exact_volume
              << std::endl;
    std::cout << std::setprecision(20)
              << "Vfrac unclipped IRL = " << our_moments.volume() / VOLUME_MAX
              << std::endl;
    std::cout << std::setprecision(20)
              << "Centroid unclipped EX  = " << exact_centroid << std::endl;
    std::cout << std::setprecision(20)
              << "Centroid unclipped IRL = " << our_centroid << std::endl;
    // std::cout << "Diff Surface EX/IRL = "
    //           << std::fabs(our_surface_area - exact_surface_area) /
    //                  std::pow(poly_vol, 2.0 / 3.0)
    //           << std::endl;
    std::cout << "Diff Vfrac EX/IRL   = "
              << std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume)
              << std::endl;
    std::cout << "Diff centroid EX/IRL   = "
              << Pt(exact_centroid - our_centroid)
              << std::endl;
    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;

    myfile.open("fig_cylinder.csv", std::ios::app);
    // myfile << "k,m0p,m0p_exact,m0p_error" << std::endl; // ,m1yp,m1yp_exact,m1yp_error
    myfile << std::scientific << std::setprecision(20) << (3.0 * static_cast<double>(i) /
    static_cast<double>(Ntests - 1)) << ","
           << our_moments.volume() / VOLUME_MAX << "," << exact_volume << ","
           << std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume) << ","
           << our_centroid[1] << "," << exact_m1y << ","
           << std::fabs(our_centroid[1] - exact_m1y) << ","
           << our_centroid[2] << "," << exact_m1z << ","
           << std::fabs(our_centroid[2] - exact_m1z) << std::endl; // " "
          //  << our_surface_area << " " << exact_surface_area << " "
          //  << std::fabs(our_surface_area - exact_surface_area) << "\n";
    myfile.close();

    max_volume_error =
        max_volume_error >
                std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume)
            ? max_volume_error
            : std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume);
    // max_surface_error =
    //     max_surface_error > std::fabs(our_surface_area - exact_surface_area) /
    //                             std::pow(poly_vol, 2.0 / 3.0)
    //         ? max_surface_error
    //         : std::fabs(our_surface_area - exact_surface_area) /
    //               std::pow(poly_vol, 2.0 / 3.0);
    rms_volume_error += std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume) *
                        std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume);
    // rms_surface_error += std::fabs(our_surface_area - exact_surface_area) *
    //                      std::fabs(our_surface_area - exact_surface_area) /
    //                      std::pow(poly_vol, 4.0 / 3.0);
  }
  rms_volume_error = sqrt(rms_volume_error / static_cast<double>(Ntests));
  // rms_surface_error = sqrt(rms_surface_error / static_cast<double>(Ntests));

  // std::cout << "Max surface error = " << max_surface_error << std::endl;
  // std::cout << "RMS surface error = " << rms_surface_error << std::endl;
  std::cout << "Max volume error  = " << max_volume_error << std::endl;
  std::cout << "RMS volume error  = " << rms_volume_error << std::endl;
  std::cout << "-------------------------------------------------------------"
               "---------------------------------------------------------"
            << std::endl;

  EXPECT_NEAR(max_volume_error, 0.0, 1.0e-14);
}

TEST(CylinderIntersection, Debug1) {
  using VolumeAndSuface =
      AddSurfaceOutput<VolumeMoments, CylinderParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({2.0, 1.2});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  auto cubes = std::array<RectangularCuboid, 6>(
      {RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(0.8, 0.8, 0.8)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(1.0, 1.0, 1.0)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(1.2, 1.2, 1.2)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, -0.5, 0.0),
                                          Pt(1.0, 0.5, 1.0)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, -1.0, 0.2),
                                          Pt(1.0, 1.0, 1.2)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, -1.0, -1.0),
                                          Pt(1.0, 1.0, 1.0))});
  std::array<std::string, 6> surface_filenames(
      {"surface_a", "surface_b", "surface_c", "surface_d", "surface_e", "surface_f"});
  std::array<std::string, 6> clipped_faces_filenames(
      {"_cube_a", "_cube_b", "_cube_c", "_cube_d", "_cube_e", "_cube_f"});


  std::array<HalfEdgePolyhedronQuadratic<Pt>, 6> half_edges;
  std::array<IRL::SegmentedHalfEdgePolyhedronQuadratic<IRL::FaceQuadratic<IRL::HalfEdgeQuadratic<IRL::VertexQuadratic<IRL::Pt>>>, IRL::VertexQuadratic<IRL::Pt>>, 6> seg_half_edges;
  for (UnsignedIndex_t i = 0; i < 6; i++) {
    cubes[i].setHalfEdgeVersion(&(half_edges[i]));
    seg_half_edges[i] = half_edges[i].generateSegmentedPolyhedron();
  }
  const UnsignedIndex_t nlevels = 10;

  // Compute moments and return parametrized surface
  for (UnsignedIndex_t i = 0; i < 6; i++) {
    auto temp_moments = getVolumeMoments<VolumeMoments>(cubes[i], cylinder);
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(cubes[i], cylinder);
    auto amr_moments =
        intersectPolyhedronWithCylinderAMR<VolumeMoments>(
            &(seg_half_edges[i]), &(half_edges[i]), aligned_cylinder, nlevels, clipped_faces_filenames[i]);
    std::cout << "the " << i << "th computed volume is :" << std::setprecision(20)
              << temp_surface_and_moments.getMoments().volume().volume()
              << std::endl;
    std::cout << "the " << i << "th amr      volume is :"
              << amr_moments.volume().volume()
              << std::endl;
    std::cout << "the " << i << "th computed centroid is :"
              << temp_surface_and_moments.getMoments().centroid() << std::endl;
    std::cout << "the " << i << "th amr      centroid is :"
              << amr_moments.centroid() << std::endl;
    auto& centroid = temp_surface_and_moments.getMoments().centroid();
    centroid /= temp_surface_and_moments.getMoments().volume().volume();
    auto& centroid_No_S = temp_moments.centroid();
    centroid_No_S /= temp_moments.volume().volume();
    auto& amr_centroid = amr_moments.centroid();
    amr_centroid /= amr_moments.volume().volume();
    std::cout << "the " << i << "th computed center of mass is :" << centroid
              << std::endl;
    std::cout << "the " << i << "th amr      center of mass is :" << amr_centroid
              << std::endl << std::endl;


    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
              temp_moments.volume().volume(), 1e-13);
    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
              amr_moments.volume().volume(), 1e-13);
    EXPECT_NEAR(centroid[0],
              amr_centroid[0], 1e-13);
    EXPECT_NEAR(centroid[0],
              centroid_No_S[0], 1e-13);
    EXPECT_NEAR(centroid[1],
              amr_centroid[1], 1e-13);
    EXPECT_NEAR(centroid[1],
              centroid_No_S[1], 1e-13);
    EXPECT_NEAR(centroid[2],
              amr_centroid[2], 1e-13);
    EXPECT_NEAR(centroid[2],
              centroid_No_S[2], 1e-13);
    auto temp_param_surface = temp_surface_and_moments.getSurface();
    auto temp_tri_surface = temp_param_surface.triangulate(0.1);
    temp_tri_surface.write(surface_filenames[i]);
  }
}

TEST(CylinderIntersection, Debug2) {
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
  std::string surface_filename = "surface_debug";

  HalfEdgePolyhedronQuadratic<Pt> half_edge;
  prism.setHalfEdgeVersion(&half_edge);
  auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

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
  const int nlevels = 8;
  auto amr_volumeMoments = intersectPolyhedronWithCylinderAMR<VolumeMoments>(
      &seg_half_edge, &half_edge, aligned_cylinder, nlevels, "_cube_debug");
  std::cout << "the amr volume  :"
            << amr_volumeMoments.volume().volume() << std::endl;
  auto& centroid = temp_surface_and_moments.getMoments().centroid();
  centroid /= temp_surface_and_moments.getMoments().volume();
  auto& amr_centroid = amr_volumeMoments.centroid();
  amr_centroid /= amr_volumeMoments.volume();
  std::cout << "the normalize centroid is :" << centroid << std::endl;
  std::cout << "expected centroid         :( ??????????????????, " << theoretical_centroid[0]
            << ", " << theoretical_centroid[1] << ")\n";
  std::cout << "the amr centroid is       :" << amr_centroid << std::endl;


  auto temp_param_surface = temp_surface_and_moments.getSurface();
  auto temp_tri_surface = temp_param_surface.triangulate(0.1);
  temp_tri_surface.write(surface_filename);

  EXPECT_NEAR(temp_surface_and_moments.getMoments().volume(), th_volume_total, 1e-13);
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

  EXPECT_NEAR(temp_moments.volume(), M_PI, 1e-13);
}

TEST(HyperCylinderIntersection, SISCPaperFig5) {
  using VolumeAndSuface =
      AddSurfaceOutput<VolumeMoments, CylinderParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedCylinder aligned_cylinder({-2.0, 0.1});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // Constructing cells for each subfigure
  auto cubes = std::array<RectangularCuboid, 3>(
      {RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(0.8, 0.8, 0.8)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, -0.5, 0.0),
                                          Pt(1.0, 0.5, 1.0)),
       RectangularCuboid::fromBoundingPts(Pt(-1.0, -1.0,-1.0),
                                          Pt( 1.0,  1.0, 1.0))});
  std::array<std::string, 3> surface_filenames(
      {"hyper_surface_a", "hyper_surface_b", "hyper_surface_c"});
  std::array<std::string, 3> clipped_faces_filenames(
      {"_hyper_cube_a", "_hyper_cube_b", "_hyper_cube_c"});

  std::array<HalfEdgePolyhedronQuadratic<Pt>, 3> half_edges;
  std::array<IRL::SegmentedHalfEdgePolyhedronQuadratic<IRL::FaceQuadratic<IRL::HalfEdgeQuadratic<IRL::VertexQuadratic<IRL::Pt>>>, IRL::VertexQuadratic<IRL::Pt>>, 3> seg_half_edges;
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    cubes[i].setHalfEdgeVersion(&(half_edges[i]));
    seg_half_edges[i] = half_edges[i].generateSegmentedPolyhedron();
  }
  const UnsignedIndex_t nlevels = 10;

  // Compute moments and return parametrized surface
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(cubes[i], cylinder);
    auto temp_moments = getVolumeMoments<VolumeMoments>(cubes[i], cylinder);
    auto amr_moments =
        intersectPolyhedronWithCylinderAMR<VolumeMoments>(
            &(seg_half_edges[i]), &(half_edges[i]), aligned_cylinder, nlevels, clipped_faces_filenames[i]);
    std::cout << "the " << i << "th computed volume is :" << std::setprecision(20)
              << temp_surface_and_moments.getMoments().volume().volume()
              << std::endl;
    std::cout << "the " << i << "th amr      volume is :"
              << amr_moments.volume().volume()
              << std::endl;
    std::cout << "the " << i << "th computed centroid is :"
              << temp_surface_and_moments.getMoments().centroid() << std::endl;
    std::cout << "the " << i << "th amr      centroid is :"
              << amr_moments.centroid() << std::endl;
    auto& centroid = temp_surface_and_moments.getMoments().centroid();
    centroid /= temp_surface_and_moments.getMoments().volume().volume();
    auto& centroid_No_S = temp_moments.centroid();
    centroid_No_S /= temp_moments.volume().volume();
    auto& amr_centroid = amr_moments.centroid();
    amr_centroid /= amr_moments.volume().volume();
    std::cout << "the " << i << "th computed center of mass is :" << centroid
              << std::endl;
    std::cout << "the " << i << "th amr      center of mass is :" << amr_centroid
              << std::endl << std::endl;
    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
              temp_moments.volume().volume(), 1e-13);
    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
              amr_moments.volume().volume(), 1e-13);
    EXPECT_NEAR(centroid[0],
              amr_centroid[0], 1e-13);
    EXPECT_NEAR(centroid[0],
              centroid_No_S[0], 1e-13);
    EXPECT_NEAR(centroid[1],
              amr_centroid[1], 1e-13);
    EXPECT_NEAR(centroid[1],
              centroid_No_S[1], 1e-13);
    EXPECT_NEAR(centroid[2],
              amr_centroid[2], 1e-13);
    EXPECT_NEAR(centroid[2],
              centroid_No_S[2], 1e-13);
    auto temp_param_surface = temp_surface_and_moments.getSurface();
    auto temp_tri_surface = temp_param_surface.triangulate(0.1);
    temp_tri_surface.write(surface_filenames[i]);
  }
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

  EXPECT_NEAR(temp_surface_and_moments.getMoments().volume(),
            just_volum.volume(), 1e-13);
}

}  // namespace
