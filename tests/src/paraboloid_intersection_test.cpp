// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/geometry/general/new_pt_calculation_functors.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/rotations.h"
#include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
#include "irl/paraboloid_reconstruction/hessian_paraboloid.h"

#include <sys/time.h>
#include <cmath>
#include <random>

// quadratic_intersection_amr have some definition
// if paraboloid_intersection_test and cylinder_intersection_test
// are compile separatly, we get an error during linking.
#include "tests/src/cylinder_intersection_test.cpp"

#include "gtest/gtest.h"

#include "irl/data_structures/small_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection_amr.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron_quadratic.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_quadratic.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/moments/general_moments.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_separator.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"
#include "tests/src/basic_mesh.h"
#include "tests/src/vtk.h"

namespace {

double ExactM0TranslatingCube(const double k, const double h) {
  double exact_volume = (std::pow(k, 2.) * M_PI) / 8.;
  if (k > h) {
    exact_volume -= (std::pow(h - k, 2.) * M_PI) / 8.;
  }
  if (k > h * h) {
    exact_volume -=
        (8. * std::pow(h, 3.) * std::sqrt(-std::pow(h, 2.) + k) -
         20. * h * k * std::sqrt(-std::pow(h, 2.) + k) +
         3. * std::pow(k, 2.) * M_PI -
         6. * std::pow(k, 2.) * std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) +
         6. * std::pow(k, 2.) *
             std::atan(std::sqrt(-1. + k / std::pow(h, 2.)))) /
        24.;
  }
  if (k > 2.0 * h * h) {
    exact_volume +=
        (2. * h *
             (-4. * std::pow(h, 3.) + 6. * h * k +
              2. * std::pow(h, 2.) * std::sqrt(-std::pow(h, 2.) + k) -
              5. * k * std::sqrt(-std::pow(h, 2.) + k)) -
         3. * std::pow(k, 2.) * std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) +
         3. * std::pow(k, 2.) *
             std::atan(std::sqrt(-1. + k / std::pow(h, 2.)))) /
        12.;
  }
  if ((k - h) > h * h) {
    exact_volume +=
        (20. * std::pow(h, 2.) * std::sqrt(-h - std::pow(h, 2.) + k) +
         8. * std::pow(h, 3.) * std::sqrt(-h - std::pow(h, 2.) + k) -
         20. * h * k * std::sqrt(-h - std::pow(h, 2.) + k) +
         3. * std::pow(h, 2.) * M_PI - 6. * h * k * M_PI +
         3. * std::pow(k, 2.) * M_PI +
         6. * std::pow(h - k, 2.) *
             std::atan(
                 std::sqrt(-((h + std::pow(h, 2.) - k) / std::pow(h, 2.)))) -
         6. * std::pow(h - k, 2.) *
             std::atan(h / std::sqrt(-h - std::pow(h, 2.) + k))) /
        24.;
  }
  return exact_volume;
}

double ExactM1xTranslatingCube(const double k, const double h) {
  double exact_m1x = (2. * std::pow(k, 2.5)) / 15.;
  if (k > h) {
    exact_m1x -= (2. * std::pow(-h + k, 2.5)) / 15.;
  }
  if (k > h * h) {
    exact_m1x -= (-3. * std::pow(h, 5.) + 10. * std::pow(h, 3.) * k -
                  15. * h * std::pow(k, 2.) + 8. * std::pow(k, 2.5) +
                  8. * std::pow(-std::pow(h, 2.) + k, 2.5)) /
                 60.;
  }
  if (k > 2.0 * h * h) {
    exact_m1x +=
        (-28. * std::pow(h, 5.) + 40. * std::pow(h, 3.) * k -
         15. * h * std::pow(k, 2.) + 8. * std::pow(-std::pow(h, 2.) + k, 2.5)) /
        60.;
  }
  if ((k - h) > h * h) {
    exact_m1x += (8. * std::pow(-h - std::pow(h, 2.) + k, 2.5) +
                  5. * h * (h + std::pow(h, 2.) - k) *
                      (-3. * h + std::pow(h, 2.) + 3. * k) +
                  8. * (-std::pow(h, 5.) + std::pow(-h + k, 2.5))) /
                 60.;
  }
  return exact_m1x;
}

double ExactM1zTranslatingCube(const double k, const double h) {
  double exact_m1z = -0.08333333333333333 * (std::pow(k, 3.) * M_PI);
  if (k > h) {
    exact_m1z -= (std::pow(h - k, 3.) * M_PI) / 12.;
  }
  if (k > h * h) {
    exact_m1z -= (-16. * std::pow(h, 5.) * std::sqrt(-std::pow(h, 2.) + k) -
                  8. * std::pow(h, 3.) * k * std::sqrt(-std::pow(h, 2.) + k) +
                  84. * h * std::pow(k, 2.) * std::sqrt(-std::pow(h, 2.) + k) -
                  15. * std::pow(k, 3.) * M_PI +
                  30. * std::pow(k, 3.) *
                      std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) -
                  30. * std::pow(k, 3.) *
                      std::atan(std::sqrt(-1. + k / std::pow(h, 2.)))) /
                 180.;
  }
  if (k > 2.0 * h * h) {
    exact_m1z +=
        (28. * std::pow(h, 6.) - 45. * std::pow(h, 2.) * std::pow(k, 2.) -
         8. * std::pow(h, 5.) * std::sqrt(-std::pow(h, 2.) + k) -
         4. * std::pow(h, 3.) * k * std::sqrt(-std::pow(h, 2.) + k) +
         42. * h * std::pow(k, 2.) * std::sqrt(-std::pow(h, 2.) + k) +
         15. * std::pow(k, 3.) *
             std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) -
         15. * std::pow(k, 3.) *
             std::atan(std::sqrt(-std::pow(h, 2.) + k) / h)) /
        90.;
  }
  if ((k - h) > h * h) {
    exact_m1z +=
        (84. * std::pow(h, 3.) * std::sqrt(-h - std::pow(h, 2.) + k) +
         8. * std::pow(h, 4.) * std::sqrt(-h - std::pow(h, 2.) + k) -
         16. * std::pow(h, 5.) * std::sqrt(-h - std::pow(h, 2.) + k) -
         168. * std::pow(h, 2.) * k * std::sqrt(-h - std::pow(h, 2.) + k) -
         8. * std::pow(h, 3.) * k * std::sqrt(-h - std::pow(h, 2.) + k) +
         84. * h * std::pow(k, 2.) * std::sqrt(-h - std::pow(h, 2.) + k) +
         15. * std::pow(h, 3.) * M_PI - 45. * std::pow(h, 2.) * k * M_PI +
         45. * h * std::pow(k, 2.) * M_PI - 15. * std::pow(k, 3.) * M_PI -
         30. * std::pow(h - k, 3.) *
             std::atan(h / std::sqrt(-h - std::pow(h, 2.) + k)) +
         30. * std::pow(h - k, 3.) *
             std::atan(std::sqrt(-h - std::pow(h, 2.) + k) / h)) /
        180.;
  }
  return exact_m1z;
}

double ExactS0TranslatingCube(const double k, const double h) {
  double exact_surface_area =
      ((-1. + std::sqrt(1. + 4. * k) + 4. * k * std::sqrt(1. + 4. * k)) *
       M_PI) /
      24.;
  if (k > h) {
    exact_surface_area -= ((-1. + std::sqrt(1. - 4. * h + 4. * k) -
                            4. * h * std::sqrt(1. - 4. * h + 4. * k) +
                            4. * k * std::sqrt(1. - 4. * h + 4. * k)) *
                           M_PI) /
                          24.;
  }
  if (k > h * h) {
    exact_surface_area -=
        (-8. * h * std::sqrt(-((std::pow(h, 2.) - k) * (1. + 4. * k))) - M_PI +
         std::sqrt(1. + 4. * k) * M_PI +
         4. * k * std::sqrt(1. + 4. * k) * M_PI -
         2. * std::pow(1. + 4. * k, 1.5) *
             std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) +
         2. * std::sqrt(1. + 4. * k) *
             std::atan(std::sqrt(-std::pow(h, 2.) + k) / h) +
         8. * k * std::sqrt(1. + 4. * k) *
             std::atan(std::sqrt(-std::pow(h, 2.) + k) / h) +
         2. * std::atan(4. * h *
                        std::sqrt((-std::pow(h, 2.) + k) / (1. + 4. * k))) +
         2. * std::atan((h * std::sqrt(1. + 4. * k)) /
                        std::sqrt(-std::pow(h, 2.) + k)) -
         2. * std::atan(std::sqrt((-std::pow(h, 2.) + k) * (1. + 4. * k)) / h) -
         2. * h * (3. + 4. * std::pow(h, 2.)) *
             std::atanh(2. *
                        std::sqrt((-std::pow(h, 2.) + k) / (1. + 4. * k))) +
         3. * h * std::log(1. + 4. * std::pow(h, 2.)) +
         4. * std::pow(h, 3.) * std::log(1. + 4. * std::pow(h, 2.)) -
         6. * h *
             std::log(2. * std::sqrt(-std::pow(h, 2.) + k) +
                      std::sqrt(1. + 4. * k)) -
         8. * std::pow(h, 3.) *
             std::log(2. * std::sqrt(-std::pow(h, 2.) + k) +
                      std::sqrt(1. + 4. * k))) /
        24.;
  }
  if (k > 2.0 * h * h) {
    exact_surface_area +=
        (4. * std::pow(h, 2.) * std::sqrt(1. + 8. * std::pow(h, 2.)) -
         4. * h * std::sqrt(-((std::pow(h, 2.) - k) * (1. + 4. * k))) -
         std::atan((4. * std::pow(h, 2.)) /
                   std::sqrt(1. + 8. * std::pow(h, 2.))) -
         std::sqrt(1. + 4. * k) *
             std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) -
         4. * k * std::sqrt(1. + 4. * k) *
             std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) +
         std::sqrt(1. + 4. * k) *
             std::atan(std::sqrt(-std::pow(h, 2.) + k) / h) +
         4. * k * std::sqrt(1. + 4. * k) *
             std::atan(std::sqrt(-std::pow(h, 2.) + k) / h) +
         std::atan(4. * h * std::sqrt((-std::pow(h, 2.) + k) / (1. + 4. * k))) +
         std::atan((h * std::sqrt(1. + 4. * k)) /
                   std::sqrt(-std::pow(h, 2.) + k)) -
         std::atan(std::sqrt((-std::pow(h, 2.) + k) * (1. + 4. * k)) / h) +
         h * (3. + 4. * std::pow(h, 2.)) *
             std::atanh((2. * h) / std::sqrt(1. + 8. * std::pow(h, 2.))) -
         h * (3. + 4. * std::pow(h, 2.)) *
             std::atanh(2. *
                        std::sqrt((-std::pow(h, 2.) + k) / (1. + 4. * k))) +
         3. * h * std::log(2. * h + std::sqrt(1. + 8. * std::pow(h, 2.))) +
         4. * std::pow(h, 3.) *
             std::log(2. * h + std::sqrt(1. + 8. * std::pow(h, 2.))) -
         3. * h *
             std::log(2. * std::sqrt(-std::pow(h, 2.) + k) +
                      std::sqrt(1. + 4. * k)) -
         4. * std::pow(h, 3.) *
             std::log(2. * std::sqrt(-std::pow(h, 2.) + k) +
                      std::sqrt(1. + 4. * k))) /
        12.;
  }
  if ((k - h) > h * h) {
    exact_surface_area +=
        (-8. * h *
             std::sqrt(4. * std::pow(h, 3.) + std::pow(h, 2.) * (3. - 4. * k) +
                       k * (1. + 4. * k) - h * (1. + 8. * k)) -
         M_PI + std::sqrt(1. - 4. * h + 4. * k) * M_PI -
         4. * h * std::sqrt(1. - 4. * h + 4. * k) * M_PI +
         4. * k * std::sqrt(1. - 4. * h + 4. * k) * M_PI +
         2. * std::atan(h / std::sqrt((h + std::pow(h, 2.) - k) /
                                      (-1. + 4. * h - 4. * k))) +
         2. * std::atan(4. * h *
                        std::sqrt((h + std::pow(h, 2.) - k) /
                                  (-1. + 4. * h - 4. * k))) -
         2. * std::sqrt(1. - 4. * h + 4. * k) *
             std::atan(h / std::sqrt(-h - std::pow(h, 2.) + k)) +
         8. * h * std::sqrt(1. - 4. * h + 4. * k) *
             std::atan(h / std::sqrt(-h - std::pow(h, 2.) + k)) -
         8. * k * std::sqrt(1. - 4. * h + 4. * k) *
             std::atan(h / std::sqrt(-h - std::pow(h, 2.) + k)) +
         2. * std::sqrt(1. - 4. * h + 4. * k) *
             std::atan(std::sqrt(-h - std::pow(h, 2.) + k) / h) -
         8. * h * std::sqrt(1. - 4. * h + 4. * k) *
             std::atan(std::sqrt(-h - std::pow(h, 2.) + k) / h) +
         8. * k * std::sqrt(1. - 4. * h + 4. * k) *
             std::atan(std::sqrt(-h - std::pow(h, 2.) + k) / h) -
         2. * std::atan(std::sqrt(4. * std::pow(h, 3.) +
                                  std::pow(h, 2.) * (3. - 4. * k) +
                                  k * (1. + 4. * k) - h * (1. + 8. * k)) /
                        h) -
         2. * h * (3. + 4. * std::pow(h, 2.)) *
             std::atanh(2. * std::sqrt((h + std::pow(h, 2.) - k) /
                                       (-1. + 4. * h - 4. * k))) +
         3. * h * std::log(1. + 4. * std::pow(h, 2.)) +
         4. * std::pow(h, 3.) * std::log(1. + 4. * std::pow(h, 2.)) -
         6. * h *
             std::log(2. * std::sqrt(-h - std::pow(h, 2.) + k) +
                      std::sqrt(1. - 4. * h + 4. * k)) -
         8. * std::pow(h, 3.) *
             std::log(2. * std::sqrt(-h - std::pow(h, 2.) + k) +
                      std::sqrt(1. - 4. * h + 4. * k))) /
        24.;
  }
  return exact_surface_area;
}

using namespace IRL;

TEST(ParaboloidIntersection, Dodecahedron) {
  double tau = 0.25 * (sqrt(5.0) + 1.0);
  Pt M[12] = {Pt(0.0, tau, 0.5),   Pt(0.0, -tau, 0.5),  Pt(0.0, tau, -0.5),
              Pt(0.0, -tau, -0.5), Pt(0.5, 0.0, tau),   Pt(-0.5, 0.0, tau),
              Pt(0.5, 0.0, -tau),  Pt(-0.5, 0.0, -tau), Pt(tau, 0.5, 0.0),
              Pt(-tau, 0.5, 0.0),  Pt(tau, -0.5, 0.0),  Pt(-tau, -0.5, 0.0)};
  std::array<Pt, 20> vertex_list{{(1.0 / 3.0) * (M[0] + M[8] + M[2]),
                                  (1.0 / 3.0) * (M[0] + M[4] + M[8]),
                                  (1.0 / 3.0) * (M[0] + M[5] + M[4]),
                                  (1.0 / 3.0) * (M[0] + M[9] + M[5]),
                                  (1.0 / 3.0) * (M[0] + M[2] + M[9]),
                                  (1.0 / 3.0) * (M[2] + M[8] + M[6]),
                                  (1.0 / 3.0) * (M[8] + M[10] + M[6]),
                                  (1.0 / 3.0) * (M[8] + M[4] + M[10]),
                                  (1.0 / 3.0) * (M[4] + M[1] + M[10]),
                                  (1.0 / 3.0) * (M[4] + M[5] + M[1]),
                                  (1.0 / 3.0) * (M[5] + M[11] + M[1]),
                                  (1.0 / 3.0) * (M[5] + M[9] + M[11]),
                                  (1.0 / 3.0) * (M[9] + M[7] + M[11]),
                                  (1.0 / 3.0) * (M[9] + M[2] + M[7]),
                                  (1.0 / 3.0) * (M[2] + M[6] + M[7]),
                                  (1.0 / 3.0) * (M[3] + M[10] + M[1]),
                                  (1.0 / 3.0) * (M[3] + M[1] + M[11]),
                                  (1.0 / 3.0) * (M[3] + M[11] + M[7]),
                                  (1.0 / 3.0) * (M[3] + M[7] + M[6]),
                                  (1.0 / 3.0) * (M[3] + M[6] + M[10])}};

  std::array<std::array<UnsignedIndex_t, 5>, 12> face_mapping{
      {{4, 3, 2, 1, 0},
       {0, 1, 7, 6, 5},
       {1, 2, 9, 8, 7},
       {2, 3, 11, 10, 9},
       {3, 4, 13, 12, 11},
       {4, 0, 5, 14, 13},
       {15, 16, 17, 18, 19},
       {15, 19, 6, 7, 8},
       {8, 9, 10, 16, 15},
       {10, 11, 12, 17, 16},
       {5, 6, 19, 18, 14},
       {12, 13, 14, 18, 17}}};

  const UnsignedIndex_t Ntests = 10;
  const UnsignedIndex_t AMR_levels = 17;
  double max_error = 0.0, rms_error = 0.0;
  bool first_vertex_on_surface = false;
  HalfEdgePolyhedronQuadratic<Pt> half_edge;
  // Create random number generator and seed it with entropy
  std::random_device rd;
  std::mt19937_64 eng(0);  // rd());

  // Bounds of paraboloid parameters
  std::uniform_real_distribution<double> random_rotation(-0.5 * M_PI,
                                                         0.5 * M_PI);
  std::uniform_real_distribution<double> random_coeffs_a(-5.0, 5.0);
  std::uniform_real_distribution<double> random_coeffs_b(-5.0, 5.0);
  std::uniform_real_distribution<double> random_translation(-0.5, 0.5);

  for (UnsignedIndex_t i = 0; i < Ntests; i++) {
    PolyhedronConnectivity connectivity(face_mapping);
    GeneralPolyhedron dodeca(vertex_list, &connectivity);
    GeneralPolyhedron dodeca_local_frame(vertex_list, &connectivity);
    auto aligned_paraboloid =
        AlignedParaboloid({random_coeffs_a(eng), random_coeffs_b(eng)});
    ReferenceFrame frame(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                         Normal(0.0, 0.0, 1.0));
    double angles[3] = {random_rotation(eng), random_rotation(eng), 0.0};
    Pt datum(random_translation(eng), random_translation(eng),
             random_translation(eng));

    if (first_vertex_on_surface) datum[2] = 0.0;

    UnitQuaternion x_rotation(angles[0], frame[0]);
    UnitQuaternion y_rotation(angles[1], frame[1]);
    UnitQuaternion z_rotation(angles[2], frame[2]);
    frame = x_rotation * y_rotation * z_rotation * frame;
    for (auto& vertex : dodeca_local_frame) {
      Pt tmp_pt = vertex - datum;
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] = frame[d] * tmp_pt;
      }
    }

    double local_space_translation = 0.0;
    if (first_vertex_on_surface) {
      for (auto& vertex : dodeca_local_frame) {
        Pt tmp_pt = vertex;
        local_space_translation =
            -aligned_paraboloid.a() * vertex[0] * vertex[0] -
            aligned_paraboloid.b() * vertex[1] * vertex[1] - vertex[2];
        break;
      }
      for (auto& vertex : dodeca_local_frame) {
        vertex[2] += local_space_translation;
      }
    }
    datum -= local_space_translation * frame[2];

    dodeca_local_frame.setHalfEdgeVersion(&half_edge);
    auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

    // The AMR moment calculation needs plane information
    for (auto& face : seg_half_edge) {
      auto normal = Normal(0, 0, 0);
      const auto starting_half_edge = face->getStartingHalfEdge();
      auto current_half_edge = starting_half_edge;
      auto next_half_edge = starting_half_edge->getNextHalfEdge();
      const auto& start_location =
          starting_half_edge->getPreviousVertex()->getLocation();
      do {
        normal += crossProduct(
            current_half_edge->getVertex()->getLocation() - start_location,
            next_half_edge->getVertex()->getLocation() - start_location);
        current_half_edge = next_half_edge;
        next_half_edge = next_half_edge->getNextHalfEdge();
      } while (next_half_edge != starting_half_edge);
      normal.normalize();
      face->setPlane(Plane(normal, normal * start_location));
    }

    // Calculate volume of dodecahedron
    auto poly_vol = dodeca.calculateVolume();

    // Calculate volume of unclipped dodecahedron using AMR
    auto amr_volume_moments =
        intersectPolyhedronWithParaboloidAMR<VolumeMoments>(
            &seg_half_edge, &half_edge, aligned_paraboloid, AMR_levels);

    // Calculate volume of unclipped dodecahedron using IRL
    const auto paraboloid = Paraboloid(datum, frame, aligned_paraboloid.a(),
                                       aligned_paraboloid.b());

    auto our_volume_moments =
        getVolumeMoments<VolumeMoments>(dodeca, paraboloid);
    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;
    std::cout << "Test " << i + 1 << "/" << Ntests << std::endl;
    if (aligned_paraboloid.a() * aligned_paraboloid.b() > 0.0)
      std::cout << "ELLIPTIC" << std::endl;
    else if (aligned_paraboloid.a() * aligned_paraboloid.b() < 0.0)
      std::cout << "HYPERBOLIC" << std::endl;
    else
      std::cout << "PARABOLIC" << std::endl;
    std::cout << std::setprecision(20)
              << "M0 unclipped IRL = " << our_volume_moments.volume() / poly_vol
              << std::endl;
    std::cout << std::setprecision(20)
              << "M0 unclipped AMR = " << amr_volume_moments.volume() / poly_vol
              << std::endl;
    std::cout << "Diff AMR/IRL = "
              << std::fabs(our_volume_moments.volume() -
                           amr_volume_moments.volume()) /
                     poly_vol
              << std::endl;

    auto pt = Pt(0, 0, 0);
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      for (UnsignedIndex_t n = 0; n < 3; ++n) {
        pt[n] += frame[d][n] * amr_volume_moments.centroid().getPt()[d];
      }
    }
    pt += amr_volume_moments.volume() * datum;
    amr_volume_moments.centroid().getPt() = pt;

    std::cout << std::setprecision(20) << "M1 unclipped IRL = "
              << Pt(our_volume_moments.centroid()[0],
                    our_volume_moments.centroid()[1],
                    our_volume_moments.centroid()[2]) /
                     poly_vol
              << std::endl;
    std::cout << std::setprecision(20) << "M1 unclipped AMR = "
              << amr_volume_moments.centroid() / poly_vol << std::endl;
    std::cout << "Diff AMR/IRL = "
              << Pt(Pt(our_volume_moments.centroid()[0],
                       our_volume_moments.centroid()[1],
                       our_volume_moments.centroid()[2]) -
                    amr_volume_moments.centroid()) /
                     poly_vol
              << std::endl;
    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;

    max_error = max_error > std::fabs(our_volume_moments.volume() -
                                      amr_volume_moments.volume()) /
                                poly_vol
                    ? max_error
                    : std::fabs(our_volume_moments.volume() -
                                amr_volume_moments.volume()) /
                          poly_vol;
    rms_error +=
        std::fabs(our_volume_moments.volume() - amr_volume_moments.volume()) *
        std::fabs(our_volume_moments.volume() - amr_volume_moments.volume()) /
        poly_vol / poly_vol;
  }
  rms_error = sqrt(rms_error / static_cast<double>(Ntests));

  std::cout << "Max volume error = " << max_error << std::endl;
  std::cout << "RMS volume error = " << rms_error << std::endl;
  std::cout << "-------------------------------------------------------------"
               "---------------------------------------------------------"
            << std::endl;

  EXPECT_NEAR(max_error, 0.0, 1.0e-14);
}

TEST(ParaboloidIntersection, SISCPaperFig5) {
  using VolumeAndSuface =
      AddSurfaceOutput<Volume, ParaboloidParametrizedSurfaceOutput>;

  // Defining elliptic paraboloic
  AlignedParaboloid aligned_paraboloid({1.0, 1.0});
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Paraboloid paraboloid(datum, frame, aligned_paraboloid.a(),
                        aligned_paraboloid.b());

  // Constructing cells for each subfigure
  auto cubes = std::array<RectangularCuboid, 3>(
      {RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, -0.5),
                                          Pt(1.0, 1.0, +0.5)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, -1.5),
                                          Pt(1.0, 1.0, -0.5)),
       RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, -2.5),
                                          Pt(1.0, 1.0, -1.5))});
  std::array<std::string, 3> surface_filenames(
      {"surface_a", "surface_b", "surface_c"});
  std::array<std::string, 3> clipped_faces_filenames(
      {"_cube_a", "_cube_b", "_cube_c"});

  // Compute moments and return parametrized surface
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(cubes[i], paraboloid);
    auto temp_param_surface = temp_surface_and_moments.getSurface();
    auto temp_tri_surface = temp_param_surface.triangulate(0.1);
    temp_tri_surface.write(surface_filenames[i]);
  }

  // Generate approximate triangulation of clipped polyhedron using AMR
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    HalfEdgePolyhedronQuadratic<Pt> half_edge;
    cubes[i].setHalfEdgeVersion(&half_edge);
    auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
    auto dummy_volume = intersectPolyhedronWithParaboloidAMR<Volume>(
        &seg_half_edge, &half_edge, aligned_paraboloid, 10,
        clipped_faces_filenames[i]);
  }

  SUCCEED();
}

TEST(ParaboloidIntersection, SISCPaperFig6) {
  using VolumeMomentsAndSuface =
      AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

  AlignedParaboloid aligned_paraboloid({1.0, 1.0});  // DO NOT CHANGE
  Pt datum(0, 0, 0);
  ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  Paraboloid paraboloid(datum, frame, aligned_paraboloid.a(),
                        aligned_paraboloid.b());

  //////////////////////////////// YOU CAN CHANGE THESE PARAMETERS
  int Ntests = 3001;  // Number of tests
  double h = 1.0;     // Edge length of the cube
  ////////////////////////////////

  double max_volume_error = 0.0, rms_volume_error = 0.0;
  double max_surface_error = 0.0, rms_surface_error = 0.0;

  std::ofstream myfile;
  myfile.open("fig6.txt");
  myfile << "k m0p m0p_exact m0p_error m1xp m1xp_exact m1xp_error m1zp "
            "m1zp_exact m1zp_error m0s m0s_exact m0s_error\n";
  myfile.close();

  for (UnsignedIndex_t i = 0; i < Ntests; i++) {
    // Create and translate unit cube
    double k = (2.0 * h * h + h) * static_cast<double>(i) /
               static_cast<double>(Ntests - 1);
    RectangularCuboid cube =
        RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, -k), Pt(h, h, h - k));

    // Compute cube volume
    auto poly_vol = cube.calculateVolume();

    // Compute moments and return parametrized surface
    auto our_surface_and_moments =
        getVolumeMoments<VolumeMomentsAndSuface>(cube, paraboloid);
    auto our_moments = our_surface_and_moments.getMoments();
    auto our_surface = our_surface_and_moments.getSurface();
    auto our_surface_area = our_surface.getSurfaceArea();

    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;
    std::cout << "Test " << i + 1 << "/" << Ntests << std::endl;

    // Compute exact value for verification purposes
    auto exact_volume = ExactM0TranslatingCube(k, h);
    auto exact_m1x = ExactM1xTranslatingCube(k, h);
    auto exact_m1z = ExactM1zTranslatingCube(k, h);
    auto exact_surface_area = ExactS0TranslatingCube(k, h);
    auto exact_centroid =
        Pt(exact_m1x, exact_m1x, exact_m1z) / safelyEpsilon(exact_volume);
    auto our_centroid =
        our_moments.centroid() / safelyEpsilon(our_moments.volume());
    std::cout << std::setprecision(20)
              << "Surface EXACT  = " << exact_surface_area << std::endl;
    std::cout << std::setprecision(20)
              << "Surface IRL    = " << our_surface_area << std::endl;
    std::cout << std::setprecision(20)
              << "Vfrac unclipped EX  = " << exact_volume / poly_vol
              << std::endl;
    std::cout << std::setprecision(20)
              << "Vfrac unclipped IRL = " << our_moments.volume() / poly_vol
              << std::endl;
    std::cout << std::setprecision(20)
              << "Centroid unclipped EX  = " << exact_centroid << std::endl;
    std::cout << std::setprecision(20)
              << "Centroid unclipped IRL = " << our_centroid << std::endl;
    std::cout << "Diff Surface EX/IRL = "
              << std::fabs(our_surface_area - exact_surface_area) /
                     std::pow(poly_vol, 2.0 / 3.0)
              << std::endl;
    std::cout << "Diff Vfrac EX/IRL   = "
              << std::fabs(our_moments.volume() - exact_volume) / poly_vol
              << std::endl;
    std::cout << "Diff centroid EX/IRL   = "
              << Pt(exact_centroid - our_centroid) /
                     std::pow(static_cast<double>(poly_vol), 1.0 / 3.0)
              << std::endl;
    std::cout << "-------------------------------------------------------------"
                 "---------------------------------------------------------"
              << std::endl;

    myfile.open("fig6.txt", std::ios::app);
    // myfile << "k m0_err m1x_err m1z_err m0s_err\n";
    myfile << std::scientific << std::setprecision(20) << k << " "
           << our_moments.volume() << " " << exact_volume << " "
           << std::fabs(our_moments.volume() - exact_volume) << " "
           << our_moments.centroid()[0] << " " << exact_m1x << " "
           << std::fabs(our_moments.centroid()[0] - exact_m1x) << " "
           << our_moments.centroid()[2] << " " << exact_m1z << " "
           << std::fabs(our_moments.centroid()[2] - exact_m1z) << " "
           << our_surface_area << " " << exact_surface_area << " "
           << std::fabs(our_surface_area - exact_surface_area) << "\n";
    myfile.close();

    max_volume_error =
        max_volume_error >
                std::fabs(our_moments.volume() - exact_volume) / poly_vol
            ? max_volume_error
            : std::fabs(our_moments.volume() - exact_volume) / poly_vol;
    max_surface_error =
        max_surface_error > std::fabs(our_surface_area - exact_surface_area) /
                                std::pow(poly_vol, 2.0 / 3.0)
            ? max_surface_error
            : std::fabs(our_surface_area - exact_surface_area) /
                  std::pow(poly_vol, 2.0 / 3.0);
    rms_volume_error += std::fabs(our_moments.volume() - exact_volume) *
                        std::fabs(our_moments.volume() - exact_volume) /
                        poly_vol / poly_vol;
    rms_surface_error += std::fabs(our_surface_area - exact_surface_area) *
                         std::fabs(our_surface_area - exact_surface_area) /
                         std::pow(poly_vol, 4.0 / 3.0);
  }
  rms_volume_error = sqrt(rms_volume_error / static_cast<double>(Ntests));
  rms_surface_error = sqrt(rms_surface_error / static_cast<double>(Ntests));

  std::cout << "Max surface error = " << max_surface_error << std::endl;
  std::cout << "RMS surface error = " << rms_surface_error << std::endl;
  std::cout << "Max volume error  = " << max_volume_error << std::endl;
  std::cout << "RMS volume error  = " << rms_volume_error << std::endl;
  std::cout << "-------------------------------------------------------------"
               "---------------------------------------------------------"
            << std::endl;

  EXPECT_NEAR(max_surface_error, 0.0, 1.0e-14);
  EXPECT_NEAR(max_volume_error, 0.0, 1.0e-14);
}

TEST(ParaboloidIntersection, ProgressiveDistanceSolver) {
  AlignedParaboloid aligned_paraboloid;
  aligned_paraboloid.a() = 10.0;
  aligned_paraboloid.b() = 20.0;
  std::array<double, 3> angles{{M_PI / 10.0, M_PI / 5.0}};
  // std::array<double, 3> angles{{0.0, 0.0}};
  std::array<double, 3> translations{{0.0, 0.0, 0.0}};
  ReferenceFrame frame(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                       Normal(0.0, 0.0, 1.0));
  UnitQuaternion x_rotation(angles[0], frame[0]);
  UnitQuaternion y_rotation(angles[1], frame[1]);
  UnitQuaternion z_rotation(angles[2], frame[2]);
  auto total_rotation = x_rotation * y_rotation * z_rotation;
  total_rotation.normalize();
  frame = total_rotation * frame;
  auto datum = -Pt::fromArray(translations);
  Paraboloid paraboloid(datum, frame, aligned_paraboloid.a(),
                        aligned_paraboloid.b());

  auto cube = RectangularCuboid::fromBoundingPts(Pt(-1.0, -1.0, -1.0),
                                                 Pt(1.0, 1.0, 1.0));
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  std::uniform_real_distribution<double> random_vfrac(0.0, 1.0);

  double max_error = 0.0;
  double tolerance = 1.0e-14;
  int Ntest = 1e4;
  for (int i = 0; i < Ntest; ++i) {
    double vfrac_required = random_vfrac(eng);
    ProgressiveDistanceSolverParaboloid<RectangularCuboid> solver(
        cube, vfrac_required, tolerance, paraboloid);
    double distance = solver.getDistance();
    Paraboloid new_paraboloid(datum + distance * frame[2], frame,
                              aligned_paraboloid.a(), aligned_paraboloid.b());
    double error =
        std::fabs(getVolumeFraction(cube, new_paraboloid) - vfrac_required);
    max_error = std::max(max_error, error);
  }

  std::cout << std::setprecision(20) << Ntest
            << " Tests volume fraction matching -- max error = " << max_error
            << std::endl;

  EXPECT_NEAR(max_error, 0.0, 10.0 * tolerance);
}

TEST(ParaboloidIntersection, PtQuad) {
  auto pt_double = PtBase<double>(1.0, 1.0, 1.0);
  pt_double += PtBase<double>(1.0e-18, 1.0e-18, 1.0e-18);
  auto pt_quad = PtBase<Quad_t>(1.0q, 1.0q, 1.0q);
  pt_quad += PtBase<Quad_t>(1.0e-18q, 1.0e-18q, 1.0e-18q);
  auto pt_quad_conv = PtBase<Quad_t>(static_cast<Quad_t>(pt_double[0]),
                                     static_cast<Quad_t>(pt_double[1]),
                                     static_cast<Quad_t>(pt_double[2]));

  auto testsum0 = pt_double + Pt(1.0, 1.0, 1.0);

  std::cout << "   DBL_EPSILON = " << DBL_EPSILON << std::endl;
  std::cout << "FLT128_EPSILON = " << static_cast<double>(FLT128_EPSILON)
            << std::endl;
  std::cout << "        Pt double = " << pt_double << std::endl;
  std::cout << "        Pt double = " << pt_double << std::endl;
  std::cout << "          Pt quad = " << pt_quad << std::endl;
  std::cout << "Pt quad converted = " << pt_quad_conv << std::endl;

  auto test1 = pow(10.0q, 1.0q / 3.0q);
  auto test2 = pow(10.0q, 1.0 / 3.0);
  auto test3 = pow(10.0, 1.0 / 3.0);

  std::cout << std::setprecision(15);
  std::cout << "      10^1/3 quad/quad = " << static_cast<double>(test1)
            << std::endl;
  std::cout << "    10^1/3 quad/double = " << static_cast<double>(test2)
            << std::endl;
  std::cout << "  10^1/3 double/double = " << test3 << std::endl;

  const Quad_t alpha = 0.5q, beta = 0.5q;
  const Quad_t exact_volume = 1.0q / 2.0q - (alpha + beta) / 12.0q;
  const auto exact_m1 = PtBase<Quad_t>(
      0.0q, 0.0q,
      (9.0q * alpha * alpha + 9.0q * beta * beta + 10.0q * alpha * beta) /
              1440.0q -
          1.0q / 8.0q);
  std::cout << "   Volume exact = " << static_cast<double>(exact_volume)
            << std::endl;
  std::cout << " Centroid exact = " << exact_m1 << std::endl;

  Quad_t max_dp_error = 0.0q, max_qp_error = 0.0q;
  {
    // Create unit cube
    const auto bottom_corner = Pt(-0.5, -0.5, -0.5);
    const auto top_corner = Pt(0.5, 0.5, 0.5);
    const auto cell =
        StoredRectangularCuboid<Pt>::fromBoundingPts(bottom_corner, top_corner);
    // Reference frame and coefficients of paraboloid
    const auto origin = Pt(0.0, 0.0, 0.0);
    const auto frame = ReferenceFrame(
        Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0), Normal(0.0, 0.0, 1.0));
    // Create paraboloid object
    const auto paraboloid = Paraboloid(
        origin, frame, static_cast<double>(alpha), static_cast<double>(beta));
    // Compute moments of cell clipped by paraboloid
    const auto first_moments =
        getVolumeMoments<VolumeMoments>(cell, paraboloid);

    std::cout << std::setprecision(20) << std::scientific
              << "   Volume double = " << first_moments.volume() << std::endl;
    std::cout << " Centroid double = " << first_moments.centroid() << std::endl;
    const Quad_t volume_error =
        fabsq(static_cast<Quad_t>(first_moments.volume()) - exact_volume);
    const auto m1_error =
        PtBase<Quad_t>(PtBase<Quad_t>(first_moments.centroid()) - exact_m1);
    max_dp_error = maximum(max_dp_error, volume_error);
    max_dp_error = maximum(max_dp_error, m1_error[0]);
    max_dp_error = maximum(max_dp_error, m1_error[1]);
    max_dp_error = maximum(max_dp_error, m1_error[2]);
    std::cout << "           Error = " << static_cast<double>(volume_error)
              << "   " << m1_error << std::endl;
  }
  {  // Create unit cube
    const auto bottom_corner = PtBase<Quad_t>(-0.5q, -0.5q, -0.5q);
    const auto top_corner = PtBase<Quad_t>(0.5q, 0.5q, 0.5q);
    const auto cell = StoredRectangularCuboid<PtBase<Quad_t>>::fromBoundingPts(
        bottom_corner, top_corner);
    // Reference frame and coefficients of paraboloid
    const auto origin = PtBase<Quad_t>(0, 0, 0);
    const auto frame = ReferenceFrameBase<Quad_t>(NormalBase<Quad_t>(1, 0, 0),
                                                  NormalBase<Quad_t>(0, 1, 0),
                                                  NormalBase<Quad_t>(0, 0, 1));
    // Create paraboloid object
    const auto paraboloid = ParaboloidBase<Quad_t>(origin, frame, alpha, beta);
    // Compute moments of cell clipped by paraboloid
    const auto first_moments =
        getVolumeMoments<VolumeMomentsBase<Quad_t>>(cell, paraboloid);
    std::cout << "     Volume quad = " << first_moments.volume() << std::endl;
    std::cout << "   Centroid quad = " << first_moments.centroid() << std::endl;
    const Quad_t volume_error =
        fabsq(static_cast<Quad_t>(first_moments.volume()) - exact_volume);
    const auto m1_error =
        PtBase<Quad_t>(PtBase<Quad_t>(first_moments.centroid()) - exact_m1);
    max_qp_error = maximum(max_qp_error, volume_error);
    max_qp_error = maximum(max_qp_error, m1_error[0]);
    max_qp_error = maximum(max_qp_error, m1_error[1]);
    max_qp_error = maximum(max_qp_error, m1_error[2]);
    std::cout << "           Error = " << static_cast<double>(volume_error)
              << "   " << m1_error << std::endl;
  }

  double max_error = static_cast<double>(maximum(max_dp_error, max_qp_error));
  EXPECT_NEAR(max_error, 0.0, 1.0e-14);
}

}  // namespace