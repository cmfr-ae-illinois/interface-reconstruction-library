// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Valentin Wasquel <wasquel.valentin@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// #define DEBUG_CYL_IRL
// #define NUDGE_REGION

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
#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection_amr.h"
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

#include "irl/interface_reconstruction_methods/progressive_radius_solver_cylinder.h"

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
    return (1.0 / 8.0) * ((u - 2.0) * sqrt_1 + 4.0 * acos_1);
  } else {
    double sqrt_1 = sqrt((6.0 - u) * (u - 2.0));
    double sqrt_2 = sqrt(u * (4.0 - u));
    double acos_1 = acos(2.0 - HALF * u);
    double acos_2 = acos(HALF * u - 1);
    return (1.0 / 8.0) * (4.0 * M_PI - (u - 4.0) * sqrt_1 + (u - 2.0) * sqrt_2 -
                          4.0 * acos_1 - 4.0 * acos_2);
  }
}

double ExactM1yTranslatingCube_cylinder(const double u) {
  if (u < 2.0) {
    return (6.0 - u) * u * u / 48.0;
  } else {
    return (6.0 - u) * u / 8.0 - 2.0 / 3.0;
  }
}

double ExactM1zTranslatingCube_cylinder(const double u) {
  if (u < 2.0) {
    return pow((4.0 - u) * u, 3.0 / 2.0) / 24.0;
  } else {
    return (pow((4.0 - u) * u, 3.0 / 2.0) -
            pow((2.0 - u) * (u - 6.0), 3.0 / 2.0)) /
           24.0;
  }
}

double ExactM2yyTranslatingCube_cylinder(const double u) {
  double sqrt_1 = sqrt(u * (4.0 - u));
  if (u < 2.0) {
    return ((-u * u * u + 6.0 * u * u - 2.0 * u - 12.0) * sqrt_1 +
            24.0 * atan(sqrt_1 / (2.0 - u))) /
           192.0;
  } else {
    double sqrt_2 = sqrt(8.0 * u - 12.0 - u * u);
    return (-(u - 2.0) * (u * (u - 4.0) - 6.0) * sqrt_1 +
            (u - 4.0) * (u * (u - 8.0) + 6) * sqrt_2 +
            24.0 * atan(sqrt_1 / (u - 2.0)) + 24.0 * atan(sqrt_2 / (u - 4.0)) +
            48.0 * acos(sqrt_1 / 2.0)) /
           192.0;
  }
}

double ExactM2yzTranslatingCube_cylinder(const double u) {
  if (u < 2.0) {
    return (16.0 - 8.0 * u + u * u) * u * u / 128.0;
  } else {
    return (u * u * u - 9.0 * u * u + 24.0 * u - 18.0) / 16.0;
  }
}

double ExactM2zzTranslatingCube_cylinder(const double u) {
  double sqrt_1 = sqrt(u * (4.0 - u));
  if (u < 2.0) {
    return (3.0 * sqrt_1 * (u * u * u - 6.0 * u * u + 10.0 * u - 4.0) -
            8.0 * atan(sqrt_1 / (2.0 - u)) + 32.0 * asin(sqrt_1 / 2.0)) /
           192.0;
  } else {
    double sqrt_2 = sqrt((u - 2.0) * (6.0 - u));
    return (sqrt_1 * (3.0 * u * u * u - 18.0 * u * u + 30.0 * u - 12.0) +
            sqrt_2 * (168.0 - 138.0 * u + 36.0 * u * u - 3.0 * u * u * u) +
            48.0 * acos(sqrt_1 / 2.0) + 32.0 * asin(sqrt_1 / 2.0) -
            32.0 * asin(sqrt_2 / 2.0) - 8.0 * atan(sqrt_1 / (u - 2.0)) -
            8.0 * atan(sqrt_2 / (u - 4.0))) /
           192.0;
  }
}

using namespace IRL;

TEST(CylinderIntersection, SISCPaperFig1) {
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
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, -0.5),
                                          Pt(1.0, 1.0, 0.5))});
  std::array<std::string, 3> surface_filenames(
      {"surface_1", "surface_2", "surface_3"});
  std::array<std::string, 3> clipped_faces_filenames(
      {"_cube_1", "_cube_2", "_cube_3"});

  std::array<double, 3> offsets = {1., 2., 3.};

  std::array<HalfEdgePolyhedronQuadratic<Pt>, 3> half_edges;
  std::array<IRL::SegmentedHalfEdgePolyhedronQuadratic<
                 IRL::FaceQuadratic<
                     IRL::HalfEdgeQuadratic<IRL::VertexQuadratic<IRL::Pt>>>,
                 IRL::VertexQuadratic<IRL::Pt>>,
             3>
      seg_half_edges;
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    cubes[i].setHalfEdgeVersion(&(half_edges[i]));
    seg_half_edges[i] = half_edges[i].generateSegmentedPolyhedron();
  }
  const UnsignedIndex_t nlevels = 10;

  // Compute moments and return parametrized surface
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(cubes[i], cylinder);
    auto amr_moments = intersectPolyhedronWithCylinderAMR<VolumeMoments>(
        &(seg_half_edges[i]), &(half_edges[i]), aligned_cylinder, nlevels,
        clipped_faces_filenames[i]);
    auto centroid = temp_surface_and_moments.getMoments().centroid();
    double exact_volume = ExactM0TranslatingCube_cylinder(offsets[i]);
    double exact_My = ExactM1yTranslatingCube_cylinder(offsets[i]);
    double exact_Mz = ExactM1zTranslatingCube_cylinder(offsets[i]);

    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
                exact_volume, 1e-13);
    EXPECT_NEAR(
        centroid[0] / temp_surface_and_moments.getMoments().volume().volume(),
        0.5, 1e-13);
    // EXPECT_NEAR(centroid[1] /
    // temp_surface_and_moments.getMoments().volume().volume(),
    //           amr_centroid[1], 1e-13);
    EXPECT_NEAR(centroid[1], exact_My, 1e-13);
    EXPECT_NEAR(centroid[2], exact_Mz, 1e-13);
    auto temp_param_surface = temp_surface_and_moments.getSurface();
    auto temp_tri_surface = temp_param_surface.triangulate(0.1);
    temp_tri_surface.write(surface_filenames[i]);
  }
}

TEST(CylinderIntersection, SISCPaperFig2) {
  AlignedCylinder aligned_cylinder({1.0, 1.0});  // DO NOT CHANGE
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

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
  myfile << "k,m0p,m0p_exact,m0p_error,m1y,m1y_exact,m1y_error,m1z,m1z_exact,"
            "m1z_error"
         << std::endl;  // ,m1yp,m1yp_exact,m1yp_error
  myfile.close();

  for (UnsignedIndex_t i = 0; i < Ntests; i++) {
    // Create and translate unit cube
    double k =
        1.5 - 1.5 * (static_cast<double>(i) / static_cast<double>(Ntests - 1));
    RectangularCuboid cube = RectangularCuboid::fromBoundingPts(
        Pt(0.0, 0.0, k - 0.5), Pt(1.0, 1.0, k + 0.5));

    // Compute moments and return parametrized surface
    auto our_moments = getVolumeMoments<VolumeMoments>(cube, cylinder);

    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;
    std::cout << "Test " << i + 1 << "/" << Ntests << std::endl;

    // Compute exact value for verification purposes
    auto exact_volume =
        ExactM0TranslatingCube_cylinder(3.0 * static_cast<double>(i) /
                                        static_cast<double>(Ntests - 1)) /
        VOLUME_MAX;
    auto exact_m1y =
        ExactM1yTranslatingCube_cylinder(3.0 * static_cast<double>(i) /
                                         static_cast<double>(Ntests - 1)) /
        M1Y_MAX;
    auto exact_m1z =
        ExactM1zTranslatingCube_cylinder(3.0 * static_cast<double>(i) /
                                         static_cast<double>(Ntests - 1)) /
        M1Z_MAX;
    auto exact_centroid = Pt(0.5 * safelyEpsilon(exact_volume * VOLUME_MAX),
                             exact_m1y, exact_m1z);
    auto our_centroid =
        Pt(our_moments.centroid()[0], our_moments.centroid()[1] / M1Y_MAX,
           our_moments.centroid()[2] / M1Z_MAX);
    // std::cout << std::setprecision(20)
    //           << "Surface EXACT  = " << exact_surface_area << std::endl;
    // std::cout << std::setprecision(20)
    //           << "Surface IRL    = " << our_surface_area << std::endl;
    std::cout << std::setprecision(20)
              << "Vfrac unclipped EX  = " << exact_volume << std::endl;
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
              << Pt(exact_centroid - our_centroid) << std::endl;
    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;

    myfile.open("fig_cylinder.csv", std::ios::app);
    // myfile << "k,m0p,m0p_exact,m0p_error" << std::endl; //
    // ,m1yp,m1yp_exact,m1yp_error
    myfile << std::scientific << std::setprecision(20)
           << (3.0 * static_cast<double>(i) / static_cast<double>(Ntests - 1))
           << "," << our_moments.volume() / VOLUME_MAX << "," << exact_volume
           << "," << std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume)
           << "," << our_centroid[1] << "," << exact_m1y << ","
           << std::fabs(our_centroid[1] - exact_m1y) << "," << our_centroid[2]
           << "," << exact_m1z << "," << std::fabs(our_centroid[2] - exact_m1z)
           << std::endl;  // " "
    //  << our_surface_area << " " << exact_surface_area << " "
    //  << std::fabs(our_surface_area - exact_surface_area) << "\n";
    myfile.close();

    max_volume_error =
        max_volume_error >
                std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume)
            ? max_volume_error
            : std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume);
    // max_surface_error =
    //     max_surface_error > std::fabs(our_surface_area - exact_surface_area)
    //     /
    //                             std::pow(poly_vol, 2.0 / 3.0)
    //         ? max_surface_error
    //         : std::fabs(our_surface_area - exact_surface_area) /
    //               std::pow(poly_vol, 2.0 / 3.0);
    rms_volume_error +=
        std::fabs(our_moments.volume() / VOLUME_MAX - exact_volume) *
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

TEST(CylinderIntersection, SISCPaperFig2_M2) {
  AlignedCylinder aligned_cylinder({1.0, 1.0});  // DO NOT CHANGE
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Cylinder cylinder(datum, frame, aligned_cylinder.b(), aligned_cylinder.r());

  // const double VOLUME_MAX = sqrt(3.0) / 4.0 + M_PI / 6.0;
  // const double M1Z_MAX = 1.0 / 3.0;
  // const double M1Y_MAX = 11.0 / 24.0;

  //////////////////////////////// YOU CAN CHANGE THESE PARAMETERS
  int Ntests = 3001;  // Number of tests
  ////////////////////////////////

  double max_volume_error = 0.0, rms_volume_error = 0.0;
  // double max_surface_error = 0.0, rms_surface_error = 0.0;

  std::string file_name = "fig_cylinderM2.csv";

  std::ofstream myfile;
  myfile.open(file_name);
  myfile << "k,M0,M1x,M1y,M1z,M2xx,M2xy,M2xz,M2yy,M2yz,M2zz,M0_error,M1x_error,"
            "M1y_error,M1z_error,M2xx_error,M2xy_error,M2xz_error,M2yy_error,"
            "M2yz_error,M2zz_error"
         << std::endl;  // ,m1yp,m1yp_exact,m1yp_error
  myfile.close();

  for (UnsignedIndex_t i = 0; i < Ntests; i++) {
    // Create and translate unit cube
    const double u =
        3.0 * static_cast<double>(i) / static_cast<double>(Ntests - 1);
    double k = 1.5 - u / 2.0;
    RectangularCuboid cube = RectangularCuboid::fromBoundingPts(
        Pt(0.0, 0.0, k - 0.5), Pt(1.0, 1.0, k + 0.5));

    // Compute moments and return parametrized surface
    auto our_moments = getVolumeMoments<GeneralMoments3D<2>>(cube, cylinder);

    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;
    std::cout << "Test " << i + 1 << "/" << Ntests << std::endl;

    // Compute exact value for verification purposes
    std::array<double, 10> exact_moment = {
        ExactM0TranslatingCube_cylinder(u),
        0.5 * ExactM0TranslatingCube_cylinder(u),
        ExactM1yTranslatingCube_cylinder(u),
        ExactM1zTranslatingCube_cylinder(u),
        ExactM0TranslatingCube_cylinder(u) / 3.0,
        0.5 * ExactM1yTranslatingCube_cylinder(u),
        0.5 * ExactM1zTranslatingCube_cylinder(u),
        ExactM2yyTranslatingCube_cylinder(u),
        ExactM2yzTranslatingCube_cylinder(u),
        ExactM2zzTranslatingCube_cylinder(u)};

    std::array<double, 10> errors;
    double max_error = 0.0;
    for (int j = 0; j < 10; j++) {
      errors[j] = exact_moment[j] - our_moments[j];
      max_error = maximum(max_error, errors[j]);
    }
    std::cout << "max error = " << max_error << std::endl;
    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;

    myfile.open(file_name, std::ios::app);
    // myfile << "k,m0p,m0p_exact,m0p_error" << std::endl; //
    // ,m1yp,m1yp_exact,m1yp_error
    myfile << std::scientific << std::setprecision(20)
           << (3.0 * static_cast<double>(i) / static_cast<double>(Ntests - 1))
           << ",";
    for (int j = 0; j < 10; j++) {
      myfile << exact_moment[j] << ",";
    }
    for (int j = 0; j < 10; j++) {
      myfile << errors[j];
      if (j < 9) {
        myfile << ",";
      } else {
        myfile << std::endl;
      }
    }
    myfile.close();

    max_volume_error = max_volume_error > std::fabs(errors[0])
                           ? max_volume_error
                           : std::fabs(errors[0]);
    // max_surface_error =
    //     max_surface_error > std::fabs(our_surface_area - exact_surface_area)
    //     /
    //                             std::pow(poly_vol, 2.0 / 3.0)
    //         ? max_surface_error
    //         : std::fabs(our_surface_area - exact_surface_area) /
    //               std::pow(poly_vol, 2.0 / 3.0);
    rms_volume_error += std::fabs(errors[0]) * std::fabs(errors[0]);
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

  EXPECT_NEAR(max_volume_error, 0.0, 1.0e-12);
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
  std::array<std::string, 6> surface_filenames({"surface_a", "surface_b",
                                                "surface_c", "surface_d",
                                                "surface_e", "surface_f"});
  std::array<std::string, 6> clipped_faces_filenames(
      {"_cube_a", "_cube_b", "_cube_c", "_cube_d", "_cube_e", "_cube_f"});

  std::array<std::string, 10> moment_name({"M0", "M1x", "M1y", "M1z", "M2xx",
                                           "M2xy", "M2xz", "M2yy", "M2yz",
                                           "M2zz"});

  std::array<HalfEdgePolyhedronQuadratic<Pt>, 6> half_edges;
  std::array<IRL::SegmentedHalfEdgePolyhedronQuadratic<
                 IRL::FaceQuadratic<
                     IRL::HalfEdgeQuadratic<IRL::VertexQuadratic<IRL::Pt>>>,
                 IRL::VertexQuadratic<IRL::Pt>>,
             6>
      seg_half_edges;
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
    auto temp_moment_with_integral =
        getVolumeMoments<GeneralMoments3D<2>>(cubes[i], cylinder);
    auto amr_moments = intersectPolyhedronWithCylinderAMR<VolumeMoments>(
        &(seg_half_edges[i]), &(half_edges[i]), aligned_cylinder, nlevels,
        clipped_faces_filenames[i]);
    std::cout << "the " << i
              << "th computed volume is :" << std::setprecision(20)
              << temp_surface_and_moments.getMoments().volume().volume()
              << std::endl;
    std::cout << "the " << i
              << "th amr      volume is :" << amr_moments.volume().volume()
              << std::endl;
    std::cout << "the " << i << "th computed centroid is :"
              << temp_surface_and_moments.getMoments().centroid() << std::endl;
    std::cout << "the " << i
              << "th amr      centroid is :" << amr_moments.centroid()
              << std::endl;
    auto& centroid = temp_surface_and_moments.getMoments().centroid();
    centroid /= temp_surface_and_moments.getMoments().volume().volume();
    auto& centroid_No_S = temp_moments.centroid();
    centroid_No_S /= temp_moments.volume().volume();
    auto& amr_centroid = amr_moments.centroid();
    amr_centroid /= amr_moments.volume().volume();
    std::cout << "the " << i << "th computed center of mass is :" << centroid
              << std::endl;
    std::cout << "the " << i
              << "th amr      center of mass is :" << amr_centroid << std::endl
              << std::endl;

    std::cout << "\nmoment 2 is : " << std::endl;
    for (int i = 0; i < 10; i++) {
      std::cout << moment_name[i] << " : " << temp_moment_with_integral[i]
                << std::endl;
    }
    std::cout << std::endl;

    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
                temp_moments.volume().volume(), 1e-13);
    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
                amr_moments.volume().volume(), 1e-13);
    EXPECT_NEAR(centroid[0], amr_centroid[0], 1e-13);
    EXPECT_NEAR(centroid[0], centroid_No_S[0], 1e-13);
    EXPECT_NEAR(centroid[1], amr_centroid[1], 1e-13);
    EXPECT_NEAR(centroid[1], centroid_No_S[1], 1e-13);
    EXPECT_NEAR(centroid[2], amr_centroid[2], 1e-13);
    EXPECT_NEAR(centroid[2], centroid_No_S[2], 1e-13);
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
  std::cout << "the amr volume  :" << amr_volumeMoments.volume().volume()
            << std::endl;
  auto& centroid = temp_surface_and_moments.getMoments().centroid();
  centroid /= temp_surface_and_moments.getMoments().volume();
  auto& amr_centroid = amr_volumeMoments.centroid();
  amr_centroid /= amr_volumeMoments.volume();
  std::cout << "the normalize centroid is :" << centroid << std::endl;
  std::cout << "expected centroid         :( ??????????????????, "
            << theoretical_centroid[0] << ", " << theoretical_centroid[1]
            << ")\n";
  std::cout << "the amr centroid is       :" << amr_centroid << std::endl;

  auto temp_param_surface = temp_surface_and_moments.getSurface();
  auto temp_tri_surface = temp_param_surface.triangulate(0.1);
  temp_tri_surface.write(surface_filename);

  EXPECT_NEAR(temp_surface_and_moments.getMoments().volume(), th_volume_total,
              1e-13);
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
  nlevels = 15;
  auto amr_volume = intersectPolyhedronWithCylinderAMR<Volume>(
      &seg_half_edge, &half_edge, aligned_cylinder, nlevels);

  std::cout << "amr volume (" << nlevels
            << " levels)   : " << std::setprecision(16) << amr_volume
            << std::setprecision(3)
            << " -- error: " << std::abs(amr_volume - M_PI) << std::endl;

  EXPECT_NEAR(temp_moments.volume(), M_PI, 1e-13);
}

TEST(HyperCylinderIntersection, SISCPaperFig1) {
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
       RectangularCuboid::fromBoundingPts(Pt(-1.0, -1.0, -1.0),
                                          Pt(1.0, 1.0, 1.0))});
  std::array<std::string, 3> surface_filenames(
      {"hyper_surface_a", "hyper_surface_b", "hyper_surface_c"});
  std::array<std::string, 3> clipped_faces_filenames(
      {"_hyper_cube_a", "_hyper_cube_b", "_hyper_cube_c"});

  std::array<HalfEdgePolyhedronQuadratic<Pt>, 3> half_edges;
  std::array<IRL::SegmentedHalfEdgePolyhedronQuadratic<
                 IRL::FaceQuadratic<
                     IRL::HalfEdgeQuadratic<IRL::VertexQuadratic<IRL::Pt>>>,
                 IRL::VertexQuadratic<IRL::Pt>>,
             3>
      seg_half_edges;
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
    auto amr_moments = intersectPolyhedronWithCylinderAMR<VolumeMoments>(
        &(seg_half_edges[i]), &(half_edges[i]), aligned_cylinder, nlevels,
        clipped_faces_filenames[i]);
    std::cout << "the " << i
              << "th computed volume is :" << std::setprecision(20)
              << temp_surface_and_moments.getMoments().volume().volume()
              << std::endl;
    std::cout << "the " << i
              << "th amr      volume is :" << amr_moments.volume().volume()
              << std::endl;
    std::cout << "the " << i << "th computed centroid is :"
              << temp_surface_and_moments.getMoments().centroid() << std::endl;
    std::cout << "the " << i
              << "th amr      centroid is :" << amr_moments.centroid()
              << std::endl;
    auto& centroid = temp_surface_and_moments.getMoments().centroid();
    centroid /= temp_surface_and_moments.getMoments().volume().volume();
    auto& centroid_No_S = temp_moments.centroid();
    centroid_No_S /= temp_moments.volume().volume();
    auto& amr_centroid = amr_moments.centroid();
    amr_centroid /= amr_moments.volume().volume();
    std::cout << "the " << i << "th computed center of mass is :" << centroid
              << std::endl;
    std::cout << "the " << i
              << "th amr      center of mass is :" << amr_centroid << std::endl
              << std::endl;
    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
                temp_moments.volume().volume(), 1e-13);
    EXPECT_NEAR(temp_surface_and_moments.getMoments().volume().volume(),
                amr_moments.volume().volume(), 1e-13);
    EXPECT_NEAR(centroid[0], amr_centroid[0], 1e-13);
    EXPECT_NEAR(centroid[0], centroid_No_S[0], 1e-13);
    EXPECT_NEAR(centroid[1], amr_centroid[1], 1e-13);
    EXPECT_NEAR(centroid[1], centroid_No_S[1], 1e-13);
    EXPECT_NEAR(centroid[2], amr_centroid[2], 1e-13);
    EXPECT_NEAR(centroid[2], centroid_No_S[2], 1e-13);
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

TEST(CylinderIntersection, VFracMatching) {
  // Create random number generator and seed it with entropy
  std::random_device rd;
  std::mt19937_64 eng(rd());

  // Bounds of cylinder parameters
  std::uniform_real_distribution<double> random_rotation(-0.5 * M_PI,
                                                         0.5 * M_PI);
  std::uniform_real_distribution<double> random_translation(-0.5, 0.5);
  std::uniform_real_distribution<double> random_coeffs_b(-5.0, 5.0);
  std::uniform_real_distribution<double> random_coeffs_r(0.0, 0.25);
  std::uniform_real_distribution<double> random_vfrac(0.0, 1.0e-2);

  // Defining random cylinder
  const AlignedCylinder aligned_cylinder(
      {random_coeffs_b(eng), random_coeffs_r(eng)});
  const Pt datum(random_translation(eng), random_translation(eng),
                 random_translation(eng));
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const double angles[3] = {random_rotation(eng), random_rotation(eng), 0.0};
  const UnitQuaternion x_rotation(angles[0], frame[0]);
  const UnitQuaternion y_rotation(angles[1], frame[1]);
  const UnitQuaternion z_rotation(angles[2], frame[2]);
  frame = x_rotation * y_rotation * z_rotation * frame;

  const Cylinder cylinder(datum, frame, aligned_cylinder.b(),
                          aligned_cylinder.r());

  // Constructing cells for each subfigure
  const auto cell = RectangularCuboid::fromBoundingPts(Pt(-0.5, -0.5, -0.5),
                                                       Pt(0.5, 0.5, 0.5));

  // Compute moments and return parametrized surface
  auto vfrac_before = getVolumeMoments<Volume>(cell, cylinder);
  vfrac_before /= cell.calculateVolume();
  std::cout << "                           cylinder = " << std::setprecision(20)
            << cylinder << std::endl;
  std::cout << "the volume fraction before matching = " << std::setprecision(20)
            << vfrac_before << std::endl;

  const double vfrac_tol = 1.0e-14;
  const double liquid_volume_fraction = random_vfrac(eng);
  std::cout << "           trying to match to vfrac = " << std::setprecision(20)
            << liquid_volume_fraction << std::endl;
  std::cout << "                     with tolerance = " << std::setprecision(20)
            << vfrac_tol << std::endl;

  ProgressiveRadiusSolverCylinder<RectangularCuboid> solver_radius(
      cell, liquid_volume_fraction, vfrac_tol, cylinder);

  const auto corrected_cylinder = solver_radius.getCylinder();

  auto vfrac_after = getVolumeMoments<Volume>(cell, corrected_cylinder);
  vfrac_after /= cell.calculateVolume();
  std::cout << " the volume fraction after matching = " << std::setprecision(20)
            << vfrac_after << std::endl;

  // Print out cell
  IRL::HalfEdgePolyhedronQuadratic<IRL::Pt> half_edge;
  cell.setHalfEdgeVersion(&half_edge);
  auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
  std::ofstream myfile;
  myfile.open("test-cylinder-vfrac-matching-cell.vtu");
  myfile << seg_half_edge;
  myfile.close();

  // Generate and print initial and final cylinder surfaces
  using VolumeAndSurface =
      AddSurfaceOutput<Volume, CylinderParametrizedSurfaceOutput>;
  auto init_mom_surf = getVolumeMoments<VolumeAndSurface>(cell, cylinder);
  auto init_surface = init_mom_surf.getSurface();
  auto init_tri_surf = init_surface.triangulate(1.0e-2);
  init_tri_surf.write("test-cylinder-vfrac-matching-initial-surface");
  auto final_mom_surf =
      getVolumeMoments<VolumeAndSurface>(cell, corrected_cylinder);
  auto final_surf = final_mom_surf.getSurface();
  auto final_tri_surf = final_surf.triangulate(1.0e-2);
  final_tri_surf.write("test-cylinder-vfrac-matching-final-surface");

  EXPECT_NEAR(liquid_volume_fraction, vfrac_after, vfrac_tol);
}

TEST(CylinderIntersection, SurfaceIntegrals) {
  // Defining circular cylinder
  std::random_device rd;
  std::mt19937_64 eng(rd());

  // Bounds of cylinder parameters
  std::uniform_real_distribution<double> random_radius(0.0, 1.0);
  const double radius = random_radius(eng);
  std::cout << "Radius = " << radius << std::endl;
  const AlignedCylinder aligned_cylinder({1.0, radius * radius});
  const Pt datum(0.0, 0.0, 0.0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const Cylinder cylinder(datum, frame, aligned_cylinder.b(),
                          aligned_cylinder.r());

  // Constructing cell
  const auto cell =
      RectangularCuboid::fromBoundingPts(Pt(0, 0, 0), Pt(2, 2, 2));

  // Print out cell
  IRL::HalfEdgePolyhedronQuadratic<IRL::Pt> half_edge;
  cell.setHalfEdgeVersion(&half_edge);
  auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
  std::ofstream myfile;
  myfile.open("test-cylinder-surface-integrals-cell.vtu");
  myfile << seg_half_edge;
  myfile.close();

  // Compute moments and return parametrized surface
  using VolumeAndSurface =
      AddSurfaceOutput<Volume, CylinderParametrizedSurfaceOutput>;
  auto mom_surf = getVolumeMoments<VolumeAndSurface>(cell, cylinder);
  auto surface = mom_surf.getSurface();
  auto tri_surfare = surface.triangulate(5.0e-2);
  tri_surfare.write("test-cylinder-surface-integrals-surface");

  const auto area = surface.getSurfaceArea();
  std::cout << "Area = " << std::setprecision(20) << area
            << " (Exact = " << radius * M_PI << ")" << std::endl;
  EXPECT_NEAR(area, radius * M_PI, 1.0e-13);

  const auto mean_curvature = surface.getMeanCurvatureIntegral() / area;
  std::cout << "Mean curvature = " << std::setprecision(20) << mean_curvature
            << " (Exact = " << 1.0 / radius << ")" << std::endl;
  EXPECT_NEAR(mean_curvature, 1.0 / radius, 1.0e-13);

  const auto normal = surface.getAverageNormalNonAligned();
  std::cout << "Normal = " << std::setprecision(20) << normal
            << " (Exact = [0, " << 1.0 / std::sqrt(2.0) << ", "
            << 1.0 / std::sqrt(2.0) << "])" << std::endl;
  EXPECT_NEAR(normal[0], 0.0, 1.0e-14);
  EXPECT_NEAR(normal[1], 1.0 / std::sqrt(2.0), 1.0e-13);
  EXPECT_NEAR(normal[2], 1.0 / std::sqrt(2.0), 1.0e-13);
}

}  // namespace
