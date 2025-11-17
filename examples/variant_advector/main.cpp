// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <mpi.h>
#include <chrono>
#include <iostream>
#include <string>

#include "examples/variant_advector/reconstruction_types.h"
#include "examples/variant_advector/solver.h"
#include "examples/variant_advector/vof_advection.h"

#include "examples/variant_advector/deformation_3d.h"
#include "examples/variant_advector/rotation_3d.h"
#include "examples/variant_advector/translation_3d.h"

static int startSimulation(const std::string& a_case_name,
                           const std::string& a_advection_method,
                           const std::string& a_reconstruction_method,
                           const double a_max_cfl, const double a_final_time,
                           const int a_viz_frequency, const int a_nx);

int main(int argc, char* argv[]) {
  if (argc != 8) {
    std::cout << "Incorrect amount of command line arguments supplied. \n";
    std::exit(-1);
  }

  // Initialize MPI
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string case_name = argv[1];
  std::string advection_method = argv[2];
  std::string reconstruction_method = argv[3];
  double max_cfl = std::stod(argv[4]);
  double final_time = std::stod(argv[5]);
  int viz_frequency = std::stoi(argv[6]);
  int Nx = std::stoi(argv[7]);

  auto start = std::chrono::system_clock::now();
  startSimulation(case_name, advection_method, reconstruction_method, max_cfl,
                  final_time, viz_frequency, Nx);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = end - start;
  if (rank == 0) {
    printf("Total run time: %20f s\n\n", runtime.count());
  }

  MPI_Finalize();

  return 0;
}

static int startSimulation(const std::string& a_case_name,
                           const std::string& a_advection_method,
                           const std::string& a_reconstruction_method,
                           const double a_max_cfl, const double a_final_time,
                           const int a_viz_frequency, const int a_nx) {
  if (a_case_name == "Translation3D") {
    return runSimulation<Translation3D>(a_case_name, a_advection_method,
                                        a_reconstruction_method, a_max_cfl,
                                        a_final_time, a_viz_frequency, a_nx);
  } else if (a_case_name == "Deformation3D") {
    return runSimulation<Deformation3D>(a_case_name, a_advection_method,
                                        a_reconstruction_method, a_max_cfl,
                                        a_final_time, a_viz_frequency, a_nx);
  } else if (a_case_name == "Rotation3D") {
    return runSimulation<Rotation3D>(a_case_name, a_advection_method,
                                     a_reconstruction_method, a_max_cfl,
                                     a_final_time, a_viz_frequency, a_nx);
  } else {
    std::cout << "Unknown case : " << a_case_name << '\n';
    std::cout
        << "Valid entries are: Translation3D, Deformation3D, Rotation3D. \n";
    std::exit(-1);
  }
  return -1;
}
