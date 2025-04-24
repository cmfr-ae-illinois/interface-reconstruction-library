// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "tools/src/paraboloid_testing/amr_run.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <random>
#include <utility>
#include <vector>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection_amr.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"

using BasePolyhedron = IRL::StoredGeneralPolyhedron<IRL::Pt, 20>;

namespace details {
std::pair<BasePolyhedron, IRL::PolyhedronConnectivity*> getGeometry(
    const std::string& a_geometry);

}

void runTranslatingCube(const double a_side_length,
                        const std::vector<std::vector<double>>& a_range_vector,
                        const std::vector<std::size_t>& a_step_vector,
                        std::fstream& a_out_file) {
  IRL::AlignedParaboloid aligned_paraboloid;
  aligned_paraboloid.a() = 1.0;  // DO NOT CHANGE
  aligned_paraboloid.b() = 1.0;  // DO NOT CHANGE
  const double h = a_side_length;

  std::size_t number_of_tests = 0;
  for (const auto& steps : a_step_vector) {
    number_of_tests += steps;
  }

  std::cout << "Total cases to run: " << number_of_tests << std::endl;

  IRL::HalfEdgePolyhedronQuadratic<IRL::Pt> half_edge;
  for (std::size_t r = 0; r < a_range_vector.size(); ++r) {
    const double step = (a_range_vector[r][1] - a_range_vector[r][0]) /
                        static_cast<double>(a_step_vector[r] - 1);
    for (std::size_t i = 0; i < a_step_vector[r]; ++i) {
      // double k = (2.0 * h * h + h) *
      //            (a_range_vector[r][0] + static_cast<double>(i) * step);
      double k = static_cast<double>(i) * step;
      IRL::RectangularCuboid cube = IRL::RectangularCuboid::fromBoundingPts(
          IRL::Pt(0.0, 0.0, -k), IRL::Pt(h, h, h - k));

      cube.setHalfEdgeVersion(&half_edge);
      auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

      for (auto& face : seg_half_edge) {
        auto normal = IRL::Normal(0.0, 0.0, 0.0);
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
        face->setPlane(IRL::Plane(normal, normal * start_location));
      }

      double exact_volume = (std::pow(k, 2.) * M_PI) / 8.;
      double exact_m1x = (2. * std::pow(k, 2.5)) / 15.;
      double exact_m1z = -0.08333333333333333 * (std::pow(k, 3.) * M_PI);
      double exact_surface_area =
          ((-1. + std::sqrt(1. + 4. * k) + 4. * k * std::sqrt(1. + 4. * k)) *
           M_PI) /
          24.;
      if (k > h) {
        exact_volume -= (std::pow(h - k, 2.) * M_PI) / 8.;
        exact_surface_area -= ((-1. + std::sqrt(1. - 4. * h + 4. * k) -
                                4. * h * std::sqrt(1. - 4. * h + 4. * k) +
                                4. * k * std::sqrt(1. - 4. * h + 4. * k)) *
                               M_PI) /
                              24.;
        exact_m1x -= (2. * std::pow(-h + k, 2.5)) / 15.;
        exact_m1z -= (std::pow(h - k, 3.) * M_PI) / 12.;
      }
      if (k > h * h) {
        exact_volume -=
            (8. * std::pow(h, 3.) * std::sqrt(-std::pow(h, 2.) + k) -
             20. * h * k * std::sqrt(-std::pow(h, 2.) + k) +
             3. * std::pow(k, 2.) * M_PI -
             6. * std::pow(k, 2.) *
                 std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) +
             6. * std::pow(k, 2.) *
                 std::atan(std::sqrt(-1. + k / std::pow(h, 2.)))) /
            24.;
        exact_surface_area -=
            (-8. * h * std::sqrt(-((std::pow(h, 2.) - k) * (1. + 4. * k))) -
             M_PI + std::sqrt(1. + 4. * k) * M_PI +
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
             2. * std::atan(std::sqrt((-std::pow(h, 2.) + k) * (1. + 4. * k)) /
                            h) -
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
        exact_m1x -= (-3. * std::pow(h, 5.) + 10. * std::pow(h, 3.) * k -
                      15. * h * std::pow(k, 2.) + 8. * std::pow(k, 2.5) +
                      8. * std::pow(-std::pow(h, 2.) + k, 2.5)) /
                     60.;
        exact_m1z -=
            (-16. * std::pow(h, 5.) * std::sqrt(-std::pow(h, 2.) + k) -
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
        exact_volume +=
            (2. * h *
                 (-4. * std::pow(h, 3.) + 6. * h * k +
                  2. * std::pow(h, 2.) * std::sqrt(-std::pow(h, 2.) + k) -
                  5. * k * std::sqrt(-std::pow(h, 2.) + k)) -
             3. * std::pow(k, 2.) *
                 std::atan(h / std::sqrt(-std::pow(h, 2.) + k)) +
             3. * std::pow(k, 2.) *
                 std::atan(std::sqrt(-1. + k / std::pow(h, 2.)))) /
            12.;
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
             std::atan(4. * h *
                       std::sqrt((-std::pow(h, 2.) + k) / (1. + 4. * k))) +
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
        exact_m1x += (-28. * std::pow(h, 5.) + 40. * std::pow(h, 3.) * k -
                      15. * h * std::pow(k, 2.) +
                      8. * std::pow(-std::pow(h, 2.) + k, 2.5)) /
                     60.;
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
        exact_volume +=
            (20. * std::pow(h, 2.) * std::sqrt(-h - std::pow(h, 2.) + k) +
             8. * std::pow(h, 3.) * std::sqrt(-h - std::pow(h, 2.) + k) -
             20. * h * k * std::sqrt(-h - std::pow(h, 2.) + k) +
             3. * std::pow(h, 2.) * M_PI - 6. * h * k * M_PI +
             3. * std::pow(k, 2.) * M_PI +
             6. * std::pow(h - k, 2.) *
                 std::atan(std::sqrt(
                     -((h + std::pow(h, 2.) - k) / std::pow(h, 2.)))) -
             6. * std::pow(h - k, 2.) *
                 std::atan(h / std::sqrt(-h - std::pow(h, 2.) + k))) /
            24.;
        exact_surface_area +=
            (-8. * h *
                 std::sqrt(4. * std::pow(h, 3.) +
                           std::pow(h, 2.) * (3. - 4. * k) + k * (1. + 4. * k) -
                           h * (1. + 8. * k)) -
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
        exact_m1x += (8. * std::pow(-h - std::pow(h, 2.) + k, 2.5) +
                      5. * h * (h + std::pow(h, 2.) - k) *
                          (-3. * h + std::pow(h, 2.) + 3. * k) +
                      8. * (-std::pow(h, 5.) + std::pow(-h + k, 2.5))) /
                     60.;
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
      auto exact_centroid = IRL::Pt(exact_m1x, exact_m1x, exact_m1z);

      // Write out k value and results
      // Total of 6 doubles
      a_out_file.write(reinterpret_cast<char*>(&k), sizeof(k));

      a_out_file.write(reinterpret_cast<char*>(&exact_volume),
                       sizeof(exact_volume));

      for (std::size_t d = 0; d < 3; ++d) {
        a_out_file.write(reinterpret_cast<char*>(&exact_centroid[d]),
                         sizeof(exact_centroid[d]));
      }

      a_out_file.write(reinterpret_cast<char*>(&exact_surface_area),
                       sizeof(exact_surface_area));
    }
  }
}

void runRandomSweep(const std::string& a_geometry,
                    const std::vector<std::vector<double>>& a_translation,
                    const std::vector<std::vector<double>>& a_rotation,
                    const std::vector<std::vector<double>>& a_coefficient,
                    const bool a_fix_to_paraboloid,
                    const std::size_t a_number_of_tests,
                    std::fstream& a_out_file) {
  /* MPI setup */
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::size_t count_printed = 0;

  const std::size_t number_of_tests_per_core =
      size == 1 ? a_number_of_tests : a_number_of_tests / (size - 1);
  const std::size_t number_of_tests_this_core =
      (rank == 0 && size > 1) ? a_number_of_tests % (size - 1)
                              : number_of_tests_per_core;
  const std::size_t batch_size = number_of_tests_per_core < RUN_BATCH_SIZE
                                     ? number_of_tests_per_core
                                     : RUN_BATCH_SIZE;
  if (rank == 0) {
    std::cout << "Rank 0 has " << number_of_tests_this_core
              << " evals; others have " << number_of_tests_per_core << " evals"
              << std::endl;
    std::cout << "Batch size = " << batch_size << std::endl;
  }
  std::vector<double> batch_data;
  batch_data.resize(12 * batch_size);

  /* Random sweep */
  auto geometry_and_connectivity = details::getGeometry(a_geometry);
  auto& geometry = geometry_and_connectivity.first;
  auto& connectivity = geometry_and_connectivity.second;

  std::random_device rd;
  int seed = static_cast<int>(rd());
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::mt19937_64 eng(static_cast<unsigned int>(seed) + rank);

  std::vector<std::uniform_real_distribution<double>> random_rotations(3);
  std::vector<std::uniform_real_distribution<double>> random_translations(3);
  for (std::size_t d = 0; d < 3; ++d) {
    random_rotations[d] = std::uniform_real_distribution<double>(
        a_rotation[0][d], a_rotation[1][d]);
    random_translations[d] = std::uniform_real_distribution<double>(
        a_translation[0][d], a_translation[1][d]);
  }

  std::vector<std::uniform_real_distribution<double>> random_coeffs(2);
  for (std::size_t d = 0; d < 2; ++d) {
    random_coeffs[d] = std::uniform_real_distribution<double>(
        a_coefficient[0][d], a_coefficient[1][d]);
  }

  IRL::HalfEdgePolyhedronQuadratic<IRL::Pt> half_edge;
  IRL::ReferenceFrame original_frame(IRL::Normal(1.0, 0.0, 0.0),
                                     IRL::Normal(0.0, 1.0, 0.0),
                                     IRL::Normal(0.0, 0.0, 1.0));

  std::size_t update_frequency =
      std::min(static_cast<std::size_t>(1e3), number_of_tests_per_core / 10);
  for (std::size_t i = 0;
       i < std::max(number_of_tests_per_core, number_of_tests_this_core); ++i) {
    std::size_t j = i % batch_size;
    if ((size == 1 || rank == 1) && (i % update_frequency == 0)) {
      std::cout << static_cast<int>(static_cast<double>(i) /
                                    static_cast<double>(
                                        std::max(number_of_tests_per_core,
                                                 number_of_tests_this_core)) *
                                    100.0)
                << "% done (" << i * size << " evals)\n";
    }
    std::array<double, 3> angles{{random_rotations[0](eng),
                                  random_rotations[1](eng),
                                  random_rotations[2](eng)}};
    IRL::Pt translations(random_translations[0](eng),
                         random_translations[1](eng),
                         random_translations[2](eng));
    IRL::AlignedParaboloid aligned_paraboloid(
        std::array<double, 2>({random_coeffs[0](eng), random_coeffs[1](eng)}));

    if (i < number_of_tests_per_core) {
      auto geom = geometry;

      if (a_fix_to_paraboloid) translations[2] = 0.0;

      IRL::UnitQuaternion x_rotation(angles[0], original_frame[0]);
      IRL::UnitQuaternion y_rotation(angles[1], original_frame[1]);
      IRL::UnitQuaternion z_rotation(angles[2], original_frame[2]);
      auto frame = x_rotation * y_rotation * z_rotation * original_frame;
      for (auto& vertex : geom) {
        IRL::Pt tmp_pt = vertex + translations;
        for (std::size_t d = 0; d < 3; ++d) {
          vertex[d] = frame[d] * tmp_pt;
        }
      }

      if (a_fix_to_paraboloid) {
        double local_space_translation = 0.0;
        auto& vertex = geom[0];
        local_space_translation =
            -aligned_paraboloid.a() * vertex[0] * vertex[0] -
            aligned_paraboloid.b() * vertex[1] * vertex[1] - vertex[2];
        for (auto& vertex : geom) {
          vertex[2] += local_space_translation;
        }
        translations += local_space_translation * frame[2];
      }

      geom.setHalfEdgeVersion(&half_edge);
      auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

      for (auto& face : seg_half_edge) {
        auto normal = IRL::Normal(0.0, 0.0, 0.0);
        const auto starting_half_edge = face->getStartingHalfEdge();
        auto current_half_edge = starting_half_edge;
        auto next_half_edge = starting_half_edge->getNextHalfEdge();
        const auto& start_location =
            starting_half_edge->getPreviousVertex()->getLocation();
        do {
          normal += IRL::crossProduct(
              current_half_edge->getVertex()->getLocation() - start_location,
              next_half_edge->getVertex()->getLocation() - start_location);
          current_half_edge = next_half_edge;
          next_half_edge = next_half_edge->getNextHalfEdge();
        } while (next_half_edge != starting_half_edge);
        normal.normalize();
        face->setPlane(IRL::Plane(normal, normal * start_location));
      }

      auto volume_moments =
          IRL::intersectPolyhedronWithParaboloidAMR<IRL::VolumeMoments>(
              &seg_half_edge, &half_edge, aligned_paraboloid, 17);

      // Move centroid to global reference frame
      // volume_moments.normalizeByVolume();
      auto centroid = IRL::Pt(0.0, 0.0, 0.0);
      for (std::size_t d = 0; d < 3; ++d) {
        for (std::size_t n = 0; n < 3; ++n) {
          centroid[n] += frame[d][n] * volume_moments.centroid()[d];
        }
      }
      // centroid -= translations;
      centroid -= translations * volume_moments.volume();
      volume_moments.centroid() = centroid;

      // Store data before sending to proc 0
      for (std::size_t d = 0; d < 3; ++d) {
        batch_data[12 * j + d] = translations[d];
        batch_data[12 * j + 3 + d] = angles[d];
        batch_data[12 * j + 9 + d] = centroid[d];
      }
      batch_data[12 * j + 6] = aligned_paraboloid.a();
      batch_data[12 * j + 7] = aligned_paraboloid.b();
      batch_data[12 * j + 8] = volume_moments.volume();
    }

    // If root: write data
    if (rank == 0 && i < number_of_tests_this_core &&
        (j == batch_size - 1 || i == number_of_tests_this_core - 1)) {
      std::size_t true_batch_size = j + 1;
      a_out_file.write(reinterpret_cast<char*>(&(batch_data.data()[0])),
                       12 * true_batch_size * sizeof(double));
      count_printed += true_batch_size;
    }

    // Else: Send back to proc 0
    if (i < number_of_tests_per_core &&
        (j == batch_size - 1 || i == number_of_tests_per_core - 1)) {
      std::size_t true_batch_size = j + 1;
      if (rank == 0) {
        for (std::size_t r = 1; r < size; ++r) {
          MPI_Status status;
          MPI_Recv(batch_data.data(), 12 * true_batch_size, MPI_DOUBLE, r,
                   1234 + r, MPI_COMM_WORLD, &status);
          a_out_file.write(reinterpret_cast<char*>(&(batch_data.data()[0])),
                           12 * true_batch_size * sizeof(double));
          count_printed += true_batch_size;
        }
      } else {
        MPI_Send(batch_data.data(), 12 * true_batch_size, MPI_DOUBLE, 0,
                 1234 + rank, MPI_COMM_WORLD);
      }
    }
  }
  if (rank == 0) {
    std::cout << "100% done" << std::endl;
    if (count_printed != a_number_of_tests) {
      std::cout << "ERROR: printed " << count_printed << " results instead of "
                << a_number_of_tests << "!" << std::endl;
    } else {
      std::cout << "SUCCESS: printed " << count_printed << " results !"
                << std::endl;
    }
  }
  delete connectivity;
}

void runOrganizedSweep(const std::string& a_geometry,
                       const std::vector<std::vector<double>>& a_translation,
                       const std::vector<std::size_t>& a_translation_steps,
                       const std::vector<std::vector<double>>& a_rotation,
                       const std::vector<std::size_t>& a_rotation_steps,
                       const std::vector<std::vector<double>>& a_coefficient,
                       const std::vector<std::size_t>& a_coefficient_steps,
                       std::fstream& a_out_file) {
  auto geometry_and_connectivity = details::getGeometry(a_geometry);
  auto& geometry = geometry_and_connectivity.first;
  auto& connectivity = geometry_and_connectivity.second;

  std::array<double, 3> translation_step_size;
  std::array<double, 3> rotation_step_size;
  for (std::size_t d = 0; d < 3; ++d) {
    translation_step_size[d] = (a_translation[1][d] - a_translation[0][d]) /
                               static_cast<double>(a_translation_steps[d] - 1);
    rotation_step_size[d] = (a_rotation[1][d] - a_rotation[0][d]) /
                            static_cast<double>(a_rotation_steps[d] - 1);
  }

  std::array<double, 2> coefficient_step_size;
  for (std::size_t d = 0; d < 2; ++d) {
    coefficient_step_size[d] = (a_coefficient[1][d] - a_coefficient[0][d]) /
                               static_cast<double>(a_coefficient_steps[d] - 1);
  }

  IRL::HalfEdgePolyhedronQuadratic<IRL::Pt> half_edge;
  IRL::ReferenceFrame original_frame(IRL::Normal(1.0, 0.0, 0.0),
                                     IRL::Normal(0.0, 1.0, 0.0),
                                     IRL::Normal(0.0, 0.0, 1.0));

  std::size_t cases_ran = 0;
  const std::size_t number_of_tests =
      a_translation_steps[0] * a_translation_steps[1] * a_translation_steps[2] *
      a_rotation_steps[0] * a_rotation_steps[1] * a_rotation_steps[2] *
      a_coefficient_steps[0] * a_coefficient_steps[1];

  std::cout << "Total cases to run: " << number_of_tests << std::endl;
  std::size_t update_frequency =
      std::min(static_cast<std::size_t>(1000), number_of_tests / 10);
  for (std::size_t tx = 0; tx < a_translation_steps[0]; ++tx) {
    for (std::size_t ty = 0; ty < a_translation_steps[1]; ++ty) {
      for (std::size_t tz = 0; tz < a_translation_steps[2]; ++tz) {
        IRL::Pt translations(
            a_translation[0][0] +
                static_cast<double>(tx) * translation_step_size[0],
            a_translation[0][1] +
                static_cast<double>(ty) * translation_step_size[1],
            a_translation[0][2] +
                static_cast<double>(tz) * translation_step_size[2]);
        for (std::size_t rx = 0; rx < a_rotation_steps[0]; ++rx) {
          for (std::size_t ry = 0; ry < a_rotation_steps[1]; ++ry) {
            for (std::size_t rz = 0; rz < a_rotation_steps[2]; ++rz) {
              std::array<double, 3> angles{
                  {a_rotation[0][0] +
                       static_cast<double>(rx) * rotation_step_size[0],
                   a_rotation[0][1] +
                       static_cast<double>(ry) * rotation_step_size[1],
                   a_rotation[0][2] +
                       static_cast<double>(rz) * rotation_step_size[2]}};
              for (std::size_t ca = 0; ca < a_coefficient_steps[0]; ++ca) {
                for (std::size_t cb = 0; cb < a_coefficient_steps[1]; ++cb) {
                  ++cases_ran;
                  if (cases_ran % update_frequency == 0) {
                    std::cout
                        << static_cast<int>(
                               static_cast<double>(cases_ran) /
                               static_cast<double>(number_of_tests) * 100.0)
                        << "% done\n";
                  }

                  IRL::AlignedParaboloid aligned_paraboloid(
                      std::array<double, 2>(
                          {a_coefficient[0][0] + static_cast<double>(ca) *
                                                     coefficient_step_size[0],
                           a_coefficient[0][1] +
                               static_cast<double>(cb) *
                                   coefficient_step_size[1]}));
                  auto geom = geometry;

                  IRL::UnitQuaternion x_rotation(angles[0], original_frame[0]);
                  IRL::UnitQuaternion y_rotation(angles[1], original_frame[1]);
                  IRL::UnitQuaternion z_rotation(angles[2], original_frame[2]);
                  auto frame =
                      x_rotation * y_rotation * z_rotation * original_frame;
                  for (auto& vertex : geom) {
                    IRL::Pt tmp_pt = vertex + translations;
                    for (std::size_t d = 0; d < 3; ++d) {
                      vertex[d] = frame[d] * tmp_pt;
                    }
                  }

                  geom.setHalfEdgeVersion(&half_edge);
                  auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

                  for (auto& face : seg_half_edge) {
                    auto normal = IRL::Normal(0.0, 0.0, 0.0);
                    const auto starting_half_edge = face->getStartingHalfEdge();
                    auto current_half_edge = starting_half_edge;
                    auto next_half_edge = starting_half_edge->getNextHalfEdge();
                    const auto& start_location =
                        starting_half_edge->getPreviousVertex()->getLocation();
                    do {
                      normal += crossProduct(
                          current_half_edge->getVertex()->getLocation() -
                              start_location,
                          next_half_edge->getVertex()->getLocation() -
                              start_location);
                      current_half_edge = next_half_edge;
                      next_half_edge = next_half_edge->getNextHalfEdge();
                    } while (next_half_edge != starting_half_edge);
                    normal.normalize();
                    face->setPlane(IRL::Plane(normal, normal * start_location));
                  }

                  auto volume_moments =
                      IRL::intersectPolyhedronWithParaboloidAMR<
                          IRL::VolumeMoments>(&seg_half_edge, &half_edge,
                                              aligned_paraboloid, 17);

                  // Move centroid to global reference frame
                  // volume_moments.normalizeByVolume();
                  auto centroid = IRL::Pt(0.0, 0.0, 0.0);
                  for (std::size_t d = 0; d < 3; ++d) {
                    for (std::size_t n = 0; n < 3; ++n) {
                      centroid[n] += frame[d][n] * volume_moments.centroid()[d];
                    }
                  }
                  // centroid -= translations;
                  centroid -= translations * volume_moments.volume();
                  volume_moments.centroid() = centroid;

                  // Write case and result to file
                  // Translations, rotations, then coefficients (8 doubles
                  // total)
                  for (std::size_t d = 0; d < 3; ++d) {
                    a_out_file.write(reinterpret_cast<char*>(&translations[d]),
                                     sizeof(translations[d]));
                  }
                  for (std::size_t d = 0; d < 3; ++d) {
                    a_out_file.write(reinterpret_cast<char*>(&angles[d]),
                                     sizeof(angles[d]));
                  }
                  a_out_file.write(
                      reinterpret_cast<char*>(&aligned_paraboloid.a()),
                      sizeof(aligned_paraboloid.a()));
                  a_out_file.write(
                      reinterpret_cast<char*>(&aligned_paraboloid.b()),
                      sizeof(aligned_paraboloid.b()));
                  // Volume and centroid (four doubles total)
                  a_out_file.write(
                      reinterpret_cast<char*>(&volume_moments.volume()),
                      sizeof(volume_moments.volume()));
                  for (std::size_t d = 0; d < 3; ++d) {
                    a_out_file.write(
                        reinterpret_cast<char*>(&volume_moments.centroid()[d]),
                        sizeof(volume_moments.centroid()[d]));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  std::cout << "100% done" << std::endl;
  delete connectivity;
}

namespace details {
std::pair<BasePolyhedron, IRL::PolyhedronConnectivity*> getGeometry(
    const std::string& a_geometry) {
  std::vector<IRL::Pt> vertices;
  IRL::PolyhedronConnectivity* connectivity = nullptr;
  if (a_geometry == "tet") {
    vertices = std::vector<IRL::Pt>{
        {IRL::Pt(1.0 / std::sqrt(3.0), 0.0, -1.0 / std::sqrt(6.0)),
         IRL::Pt(-std::sqrt(3.0) / 6.0, 0.5, -1.0 / std::sqrt(6.0)),
         IRL::Pt(-std::sqrt(3.0) / 6.0, -0.5, -1.0 / std::sqrt(6.0)),
         IRL::Pt(0.0, 0.0, 1.0 / std::sqrt(6.0))}};
    std::array<std::array<IRL::UnsignedIndex_t, 3>, 4> face_mapping{
        {{0, 1, 3}, {1, 2, 3}, {1, 0, 2}, {0, 3, 2}}};
    connectivity = new IRL::PolyhedronConnectivity(face_mapping);
  } else if (a_geometry == "cube") {
    vertices = std::vector<IRL::Pt>{
        {IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(-0.5, 0.5, -0.5),
         IRL::Pt(-0.5, 0.5, 0.5), IRL::Pt(-0.5, -0.5, 0.5),
         IRL::Pt(0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, -0.5),
         IRL::Pt(0.5, 0.5, 0.5), IRL::Pt(0.5, -0.5, 0.5)}};
    std::array<std::array<IRL::UnsignedIndex_t, 4>, 6> face_mapping{
        {{3, 2, 1, 0},
         {0, 1, 5, 4},
         {4, 5, 6, 7},
         {2, 3, 7, 6},
         {1, 2, 6, 5},
         {3, 0, 4, 7}}};
    connectivity = new IRL::PolyhedronConnectivity(face_mapping);
  } else if (a_geometry == "dodecahedron") {
    const double tau = (sqrt(5.0) + 1.0) / 2.0;
    std::array<IRL::Pt, 12> M{{IRL::Pt(0.0, tau, 1.0), IRL::Pt(0.0, -tau, 1.0),
                               IRL::Pt(0.0, tau, -1.0),
                               IRL::Pt(0.0, -tau, -1.0), IRL::Pt(1.0, 0.0, tau),
                               IRL::Pt(-1.0, 0.0, tau), IRL::Pt(1.0, 0.0, -tau),
                               IRL::Pt(-1.0, 0.0, -tau), IRL::Pt(tau, 1.0, 0.0),
                               IRL::Pt(-tau, 1.0, 0.0), IRL::Pt(tau, -1.0, 0.0),
                               IRL::Pt(-tau, -1.0, 0.0)}};
    vertices = std::vector<IRL::Pt>{{(1.0 / 3.0) * (M[0] + M[8] + M[2]),
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
    const double scale = 0.5;
    for (auto& pt : vertices) {
      pt *= scale;
    }
    std::array<std::array<IRL::UnsignedIndex_t, 5>, 12> face_mapping{
        {{5, 4, 3, 2, 1},
         {1, 2, 8, 7, 6},
         {2, 3, 10, 9, 8},
         {3, 4, 12, 11, 10},
         {4, 5, 14, 13, 12},
         {5, 1, 6, 15, 14},
         {16, 17, 18, 19, 20},
         {16, 20, 7, 8, 9},
         {9, 10, 11, 17, 16},
         {11, 12, 13, 18, 17},
         {6, 7, 20, 19, 15},
         {13, 14, 15, 19, 18}}};
    for (int i = 0; i < 12; ++i) {
      for (int j = 0; j < 5; ++j) {
        --face_mapping[i][j];
      }
    }

    connectivity = new IRL::PolyhedronConnectivity(face_mapping);
  } else if (a_geometry == "cube_hole") {
    vertices = std::vector<IRL::Pt>{
        {IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(-0.5, 0.5, -0.5),
         IRL::Pt(-0.5, 0.5, 0.5), IRL::Pt(-0.5, -0.5, 0.5),
         IRL::Pt(0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, -0.5),
         IRL::Pt(0.5, 0.5, 0.5), IRL::Pt(0.5, -0.5, 0.5),
         IRL::Pt(-0.25, 0.5, -0.25), IRL::Pt(-0.25, 0.5, 0.25),
         IRL::Pt(0.25, 0.5, 0.25), IRL::Pt(0.25, 0.5, -0.25),
         IRL::Pt(-0.25, -0.5, -0.25), IRL::Pt(-0.25, -0.5, 0.25),
         IRL::Pt(0.25, -0.5, 0.25), IRL::Pt(0.25, -0.5, -0.25)}};

    std::vector<std::vector<IRL::UnsignedIndex_t>> face_mapping(12);
    face_mapping[0] = std::vector<IRL::UnsignedIndex_t>{{3, 2, 1, 0}};
    face_mapping[1] = std::vector<IRL::UnsignedIndex_t>{{0, 1, 5, 4}};
    face_mapping[2] = std::vector<IRL::UnsignedIndex_t>{{4, 5, 6, 7}};
    face_mapping[3] = std::vector<IRL::UnsignedIndex_t>{{6, 2, 3, 7}};
    face_mapping[4] = std::vector<IRL::UnsignedIndex_t>{{1, 2, 6, 10, 9, 8}};
    face_mapping[5] = std::vector<IRL::UnsignedIndex_t>{{6, 5, 1, 8, 11, 10}};
    face_mapping[6] = std::vector<IRL::UnsignedIndex_t>{{0, 4, 7, 14, 15, 12}};
    face_mapping[7] = std::vector<IRL::UnsignedIndex_t>{{7, 3, 0, 12, 13, 14}};
    face_mapping[8] = std::vector<IRL::UnsignedIndex_t>{{8, 9, 13, 12}};
    face_mapping[9] = std::vector<IRL::UnsignedIndex_t>{{9, 10, 14, 13}};
    face_mapping[10] = std::vector<IRL::UnsignedIndex_t>{{10, 11, 15, 14}};
    face_mapping[11] = std::vector<IRL::UnsignedIndex_t>{{11, 8, 12, 15}};

    connectivity = new IRL::PolyhedronConnectivity(face_mapping);

  } else if (a_geometry == "cube_hole_convex") {
    vertices = std::vector<IRL::Pt>{
        {IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(-0.5, 0.5, -0.5),
         IRL::Pt(-0.5, 0.5, 0.5), IRL::Pt(-0.5, -0.5, 0.5),
         IRL::Pt(0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, -0.5),
         IRL::Pt(0.5, 0.5, 0.5), IRL::Pt(0.5, -0.5, 0.5),
         IRL::Pt(-0.25, 0.5, -0.25), IRL::Pt(-0.25, 0.5, 0.25),
         IRL::Pt(0.25, 0.5, 0.25), IRL::Pt(0.25, 0.5, -0.25),
         IRL::Pt(-0.25, -0.5, -0.25), IRL::Pt(-0.25, -0.5, 0.25),
         IRL::Pt(0.25, -0.5, 0.25), IRL::Pt(0.25, -0.5, -0.25)}};

    std::vector<std::vector<IRL::UnsignedIndex_t>> face_mapping(16);
    face_mapping[0] = std::vector<IRL::UnsignedIndex_t>{{3, 2, 1, 0}};
    face_mapping[1] = std::vector<IRL::UnsignedIndex_t>{{0, 1, 5, 4}};
    face_mapping[2] = std::vector<IRL::UnsignedIndex_t>{{4, 5, 6, 7}};
    face_mapping[3] = std::vector<IRL::UnsignedIndex_t>{{6, 2, 3, 7}};
    face_mapping[4] = std::vector<IRL::UnsignedIndex_t>{{1, 2, 9, 8}};
    face_mapping[5] = std::vector<IRL::UnsignedIndex_t>{{2, 6, 10, 9}};
    face_mapping[6] = std::vector<IRL::UnsignedIndex_t>{{6, 5, 11, 10}};
    face_mapping[7] = std::vector<IRL::UnsignedIndex_t>{{5, 1, 8, 11}};
    face_mapping[8] = std::vector<IRL::UnsignedIndex_t>{{0, 4, 15, 12}};
    face_mapping[9] = std::vector<IRL::UnsignedIndex_t>{{4, 7, 14, 15}};
    face_mapping[10] = std::vector<IRL::UnsignedIndex_t>{{7, 3, 13, 14}};
    face_mapping[11] = std::vector<IRL::UnsignedIndex_t>{{3, 0, 12, 13}};
    face_mapping[12] = std::vector<IRL::UnsignedIndex_t>{{8, 9, 13, 12}};
    face_mapping[13] = std::vector<IRL::UnsignedIndex_t>{{9, 10, 14, 13}};
    face_mapping[14] = std::vector<IRL::UnsignedIndex_t>{{10, 11, 15, 14}};
    face_mapping[15] = std::vector<IRL::UnsignedIndex_t>{{11, 8, 12, 15}};

    connectivity = new IRL::PolyhedronConnectivity(face_mapping);
  } else if (a_geometry == "bunny") {
    std::ifstream myfile("../vtk_data/bunny.vtk");
    std::string line;
    char space_char = ' ';
    std::vector<std::array<IRL::UnsignedIndex_t, 3>> face_mapping{};
    if (myfile.is_open()) {
      std::vector<std::string> words{};
      bool read_points = false;
      bool read_connectivity = false;
      while (getline(myfile, line)) {
        std::stringstream line_stream(line);
        std::string split_line;
        while (std::getline(line_stream, split_line, space_char)) {
          words.push_back(split_line);
        }

        if (read_points) {
          for (unsigned i = 0; i < words.size() / 3; i++) {
            vertices.push_back(IRL::Pt(std::stod(words[3 * i + 0]),
                                       std::stod(words[3 * i + 1]),
                                       std::stod(words[3 * i + 2])));
          }
        }
        if (read_connectivity) {
          for (unsigned i = 0; i < words.size() / 3; i++) {
            face_mapping.push_back(std::array<IRL::UnsignedIndex_t, 3>{
                static_cast<IRL::UnsignedIndex_t>(std::stoi(words[3 * i + 0])),
                static_cast<IRL::UnsignedIndex_t>(std::stoi(words[3 * i + 1])),
                static_cast<IRL::UnsignedIndex_t>(
                    std::stoi(words[3 * i + 2]))});
          }
        }
        if (words.size() == 0 || words[0] == "METADATA") {
          read_points = false;
          read_connectivity = false;
        } else if (words[0] == "POINTS") {
          read_points = true;
        } else if (words[0] == "CONNECTIVITY") {
          read_connectivity = true;
        } else if (words[0] == "CELL_DATA") {
          read_connectivity = false;
        }
        words.clear();
      }
      myfile.close();
    }
    std::cout << "The bunny has " << vertices.size() << " vertices!"
              << std::endl;
    std::cout << "The bunny has " << face_mapping.size() << " triangles!"
              << std::endl;
    connectivity = new IRL::PolyhedronConnectivity(face_mapping);
  } else {
    std::cout << "Unkown geometry type \"" + a_geometry + "\"" << std::endl;
    std::exit(-1);
  }

  auto polyhedron = BasePolyhedron(vertices, connectivity);

  // Normalize volume
  const double normalization_factor =
      1.0 / std::pow(polyhedron.calculateVolume(), 1.0 / 3.0);
  for (auto& vertex : polyhedron) {
    for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
      vertex[d] *= normalization_factor;
    }
  }
  // Place centroid at origin
  const auto centroid = polyhedron.calculateCentroid();
  for (auto& vertex : polyhedron) {
    for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
      vertex[d] -= centroid[d];
    }
  }

  return std::make_pair(polyhedron, connectivity);
}

}  // namespace details
