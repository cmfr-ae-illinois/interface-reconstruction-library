// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_IRL_RUN_TPP_
#define SRC_IRL_RUN_TPP_

#include <omp.h>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection_amr.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/dodecahedron.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/helpers/mymath.h"

#include "tools/src/paraboloid_testing/bunny.h"
#include "tools/src/paraboloid_testing/cube_hole.h"
#include "tools/src/paraboloid_testing/my_dodecahedron.h"

std::array<double, 10> confirmTranslatingCubeResult(
    const double a_side_length, const std::size_t a_run_size,
    std::fstream& a_results_file) {
  static constexpr std::size_t DOUBLES_IN_ONE_RUN = 6;
  const double h = a_side_length;

  IRL::ReferenceFrame original_frame(IRL::Normal(1.0, 0.0, 0.0),
                                     IRL::Normal(0.0, 1.0, 0.0),
                                     IRL::Normal(0.0, 0.0, 1.0));

  std::size_t update_frequency =
      std::min(static_cast<std::size_t>(1000), a_run_size / 10);

  std::size_t cases_run = 0;
  std::vector<double> read_data;
  std::array<double, 10> results{
      {0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0}};

  std::cout << "k m0 m0_err m1x m1x_err m1z m1z_err m0s m0s_err" << std::endl;
  while (cases_run < a_run_size) {
    const std::size_t batch_size = cases_run + RUN_BATCH_SIZE <= a_run_size
                                       ? RUN_BATCH_SIZE
                                       : a_run_size - cases_run;
    read_data.resize(batch_size * DOUBLES_IN_ONE_RUN);
    a_results_file.read(reinterpret_cast<char*>(read_data.data()),
                        sizeof(double) * batch_size * DOUBLES_IN_ONE_RUN);
    for (std::size_t i = 0; i < batch_size; ++i) {
      const double* start = read_data.data() + i * DOUBLES_IN_ONE_RUN;
      const double k = start[0];

      IRL::RectangularCuboid cube = IRL::RectangularCuboid::fromBoundingPts(
          IRL::Pt(0.0, 0.0, -k), IRL::Pt(h, h, h - k));
      IRL::Paraboloid paraboloid(IRL::Pt(0.0, 0.0, 0.0), original_frame, 1.0,
                                 1.0);

      const double start_time = omp_get_wtime();
      auto volume_and_surface = IRL::getVolumeMoments<IRL::AddSurfaceOutput<
          IRL::VolumeMoments, IRL::ParaboloidParametrizedSurfaceOutput>>(cube,
                                                               paraboloid);
      const double end_time = omp_get_wtime();
      results[9] += end_time - start_time;
      auto& moments = volume_and_surface.getMoments();
      auto& surface = volume_and_surface.getSurface();

      const double volume_err = std::fabs(moments.volume() - start[1]);
      results[0] += volume_err;
      results[1] += volume_err * volume_err;
      results[2] = std::max(results[2], volume_err);

      auto amr_centroid = IRL::Pt::fromRawDoublePointer(start + 2);
      const double err_x = std::fabs(moments.centroid()[0] - amr_centroid[0]);
      const double err_y = std::fabs(moments.centroid()[1] - amr_centroid[1]);
      const double err_z = std::fabs(moments.centroid()[2] - amr_centroid[2]);
      results[3] += err_x + err_y + err_z;
      results[4] += err_x * err_x + err_y * err_y + err_z * err_z;
      results[5] = std::max(results[5], err_x);
      results[5] = std::max(results[5], err_y);
      results[5] = std::max(results[5], err_z);

      const double err_surf = std::fabs(surface.getSurfaceArea() - start[5]);
      results[6] += err_surf;
      results[7] += err_surf * err_surf;
      results[8] = std::max(results[8], err_surf);
    }
    cases_run += batch_size;
  }
  const double case_count = static_cast<double>(a_run_size);
  results[0] /= case_count;
  results[1] = std::sqrt(results[1]) / case_count;
  results[3] /= case_count;
  results[4] = std::sqrt(results[4]) / case_count;
  results[6] /= case_count;
  results[7] = std::sqrt(results[7]) / case_count;
  results[9] = results[9] / case_count;

  return results;
}

struct TetCreator {
 public:
  static IRL::Tet create(void) {
    IRL::Tet tet =
        IRL::Tet({IRL::Pt(-std::sqrt(3.0) / 6.0, 0.5, -1.0 / std::sqrt(6.0)),
                  IRL::Pt(-std::sqrt(3.0) / 6.0, -0.5, -1.0 / std::sqrt(6.0)),
                  IRL::Pt(0.0, 0.0, 1.0 / std::sqrt(6.0)),
                  IRL::Pt(1.0 / std::sqrt(3.0), 0.0, -1.0 / std::sqrt(6.0))});
    const double normalization_factor =
        1.0 / std::pow(tet.calculateVolume(), 1.0 / 3.0);
    for (auto& vertex : tet) {
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] *= normalization_factor;
      }
    }
    const auto centroid = tet.calculateCentroid();
    for (auto& vertex : tet) {
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] -= centroid[d];
      }
    }
    return tet;
  }
};

struct CubeCreator {
 public:
  static IRL::RectangularCuboid create(void) {
    return IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(-0.5, -0.5, -0.5),
                                                   IRL::Pt(0.5, 0.5, 0.5));
  }
};

struct DodecahedronCreator {
 public:
  static IRL::MyDodecahedron create(void) {
    double tau = (sqrt(5.0) + 1.0) / 2.0;
    std::array<IRL::Pt, 12> M{{IRL::Pt(0.0, tau, 1.0), IRL::Pt(0.0, -tau, 1.0),
                               IRL::Pt(0.0, tau, -1.0),
                               IRL::Pt(0.0, -tau, -1.0), IRL::Pt(1.0, 0.0, tau),
                               IRL::Pt(-1.0, 0.0, tau), IRL::Pt(1.0, 0.0, -tau),
                               IRL::Pt(-1.0, 0.0, -tau), IRL::Pt(tau, 1.0, 0.0),
                               IRL::Pt(-tau, 1.0, 0.0), IRL::Pt(tau, -1.0, 0.0),
                               IRL::Pt(-tau, -1.0, 0.0)}};
    std::array<IRL::Pt, 20> vertex_list{{(1.0 / 3.0) * (M[0] + M[8] + M[2]),
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
    for (auto& pt : vertex_list) {
      pt *= scale;
    }

    IRL::MyDodecahedron dodeca = IRL::MyDodecahedron::fromRawPtPointer(
        vertex_list.size(), vertex_list.data());
    const double normalization_factor =
        1.0 / std::pow(dodeca.calculateVolume(), 1.0 / 3.0);
    for (auto& vertex : dodeca) {
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] *= normalization_factor;
      }
    }
    const auto centroid = dodeca.calculateCentroid();
    for (auto& vertex : dodeca) {
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] -= centroid[d];
      }
    }
    return dodeca;
  }
};

struct CubeHoleCreator {
 public:
  static IRL::CubeHole create(void) {
    IRL::CubeHole cubehole =
        IRL::CubeHole({IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(-0.5, 0.5, -0.5),
                       IRL::Pt(-0.5, 0.5, 0.5), IRL::Pt(-0.5, -0.5, 0.5),
                       IRL::Pt(0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, -0.5),
                       IRL::Pt(0.5, 0.5, 0.5), IRL::Pt(0.5, -0.5, 0.5),
                       IRL::Pt(-0.25, 0.5, -0.25), IRL::Pt(-0.25, 0.5, 0.25),
                       IRL::Pt(0.25, 0.5, 0.25), IRL::Pt(0.25, 0.5, -0.25),
                       IRL::Pt(-0.25, -0.5, -0.25), IRL::Pt(-0.25, -0.5, 0.25),
                       IRL::Pt(0.25, -0.5, 0.25), IRL::Pt(0.25, -0.5, -0.25)});
    const double normalization_factor =
        1.0 / std::pow(cubehole.calculateVolume(), 1.0 / 3.0);
    for (auto& vertex : cubehole) {
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] *= normalization_factor;
      }
    }
    const auto centroid = cubehole.calculateCentroid();
    for (auto& vertex : cubehole) {
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] -= centroid[d];
      }
    }
    return cubehole;
  }
};

struct BunnyCreator {
 public:
  static IRL::Bunny create(void) {
    std::ifstream myfile("../vtk_data/bunny.vtk");
    std::string line;
    char space_char = ' ';
    std::vector<IRL::Pt> vertices;
    if (myfile.is_open()) {
      std::vector<std::string> words{};
      bool read_points = false;
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
        if (words.size() == 0 || words[0] == "METADATA") {
          read_points = false;
        } else if (words[0] == "POINTS") {
          read_points = true;
        }
        words.clear();
      }
      myfile.close();
    }

    auto bunny = IRL::Bunny::fromRawPtPointer(vertices.size(), vertices.data());
    const double normalization_factor =
        1.0 / std::pow(bunny.calculateVolume(), 1.0 / 3.0);
    for (auto& vertex : bunny) {
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] *= normalization_factor;
      }
    }
    const auto centroid = bunny.calculateCentroid();
    for (auto& vertex : bunny) {
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        vertex[d] -= centroid[d];
      }
    }
    return bunny;
  }
};

inline std::array<double, 10> confirmParameterSweepResult(
    const std::string& a_geometry, const std::size_t a_run_size,
    std::fstream& a_results_file) {
  if (a_geometry == "tet") {
    return performParameterSweep<TetCreator>(a_run_size, a_results_file);
  } else if (a_geometry == "cube") {
    return performParameterSweep<CubeCreator>(a_run_size, a_results_file);
  } else if (a_geometry == "dodecahedron") {
    return performParameterSweep<DodecahedronCreator>(a_run_size,
                                                      a_results_file);
  } else if (a_geometry == "cube_hole") {
    return performParameterSweep<CubeHoleCreator>(a_run_size, a_results_file);
  } else if (a_geometry == "bunny") {
    return performParameterSweep<BunnyCreator>(a_run_size, a_results_file);
  } else {
    std::cout << "Unkown geometry type \"" + a_geometry + "\"" << std::endl;
    std::exit(-1);
  }
}

template <class GeometryType>
std::array<double, 10> performParameterSweep(const std::size_t a_run_size,
                                             std::fstream& a_results_file) {
  static constexpr std::size_t DOUBLES_IN_ONE_RUN = 12;
  std::cout << "Starting parameter sweep" << std::endl;
  auto geom = GeometryType::create();
  IRL::ReferenceFrame original_frame(IRL::Normal(1.0, 0.0, 0.0),
                                     IRL::Normal(0.0, 1.0, 0.0),
                                     IRL::Normal(0.0, 0.0, 1.0));

  std::size_t update_frequency =
      std::min(static_cast<std::size_t>(500000), a_run_size / 10);

  std::size_t cases_run = 0;
  std::vector<double> read_data;
  std::array<double, 10> results{
      {0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0}};
  double cases_elliptic = 0.0, cases_hyperbolic = 0.0, cases_parabolic = 0.0;
  double time_elliptic = 0.0, time_hyperbolic = 0.0, time_parabolic = 0.0;
  while (cases_run < a_run_size) {
    const std::size_t batch_size = cases_run + RUN_BATCH_SIZE <= a_run_size
                                       ? RUN_BATCH_SIZE
                                       : a_run_size - cases_run;
    read_data.resize(batch_size * DOUBLES_IN_ONE_RUN);
    a_results_file.read(reinterpret_cast<char*>(read_data.data()),
                        sizeof(double) * batch_size * DOUBLES_IN_ONE_RUN);

    for (std::size_t i = 0; i < batch_size; ++i) {
      if ((cases_run + i) % update_frequency == 0) {
        std::cout << static_cast<int>(static_cast<double>(cases_run + i) /
                                      static_cast<double>(a_run_size) * 100.0)
                  << "% done\n";
      }

      const double* start = read_data.data() + i * DOUBLES_IN_ONE_RUN;

      IRL::UnitQuaternion x_rotation(start[3], original_frame[0]);
      IRL::UnitQuaternion y_rotation(start[4], original_frame[1]);
      IRL::UnitQuaternion z_rotation(start[5], original_frame[2]);
      auto frame = x_rotation * y_rotation * z_rotation * original_frame;

      IRL::Paraboloid paraboloid(-IRL::Pt::fromRawDoublePointer(start), frame,
                                 start[6], start[7]);

      double start_time = omp_get_wtime();
      auto moments =
          IRL::getVolumeMoments<IRL::VolumeMoments>(geom, paraboloid);
      double end_time = omp_get_wtime();
      results[9] += end_time - start_time;

      if (start[6] * start[7] > 0.0) {
        time_elliptic += end_time - start_time;
        cases_elliptic += 1.0;
      } else if (start[6] * start[7] < 0.0) {
        time_hyperbolic += end_time - start_time;
        cases_hyperbolic += 1.0;
      } else {
        time_parabolic += end_time - start_time;
        cases_parabolic += 1.0;
      }

      const double volume_err = std::fabs(moments.volume() - start[8]);
      results[0] += volume_err;
      results[1] += volume_err * volume_err;
      results[2] = std::max(results[2], volume_err);
      const double err_x = std::fabs(moments.centroid()[0] - start[9]);
      const double err_y = std::fabs(moments.centroid()[1] - start[10]);
      const double err_z = std::fabs(moments.centroid()[2] - start[11]);
      results[3] += err_x + err_y + err_z;
      results[4] += err_x * err_x + err_y * err_y + err_z * err_z;
      results[5] = std::max(results[5], err_x);
      results[5] = std::max(results[5], err_y);
      results[5] = std::max(results[5], err_z);
    }
    cases_run += batch_size;
  }
  std::cout << "100% done\n";
  std::cout << "Elliptic cases   -> "
            << 100.0 * cases_elliptic /
                   (cases_elliptic + cases_hyperbolic + cases_parabolic)
            << " %" << std::endl;
  std::cout << "Elliptic time    -> "
            << 100.0 * time_elliptic /
                   (time_elliptic + time_hyperbolic + time_parabolic)
            << " %" << std::endl;
  std::cout << "Hyperbolic cases -> "
            << 100.0 * cases_hyperbolic /
                   (cases_elliptic + cases_hyperbolic + cases_parabolic)
            << " %" << std::endl;
  std::cout << "Hyperbolic time  -> "
            << 100.0 * time_hyperbolic /
                   (time_elliptic + time_hyperbolic + time_parabolic)
            << " %" << std::endl;
  std::cout << "Parabolic cases  -> "
            << 100.0 * cases_parabolic /
                   (cases_elliptic + cases_hyperbolic + cases_parabolic)
            << " %" << std::endl;
  std::cout << "Parabolic time   -> "
            << 100.0 * time_parabolic /
                   (time_elliptic + time_hyperbolic + time_parabolic)
            << " %" << std::endl;
  const double case_count = static_cast<double>(a_run_size);
  results[0] /= case_count;
  results[1] = std::sqrt(results[1]) / case_count;
  results[3] /= 3 * case_count;
  results[4] = std::sqrt(results[4]) / (3 * case_count);
  results[9] = results[9] / case_count;

  return results;
}

#endif  // SRC_IRL_RUN_TPP_
