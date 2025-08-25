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

#include "examples/sphere_initialization/reconstruction_types.h"
#include "examples/sphere_initialization/solver.h"

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout
        << "Incorrect amount of command line arguments supplied. Arguments are "
           "reconstruction method and number of cells. \n";
    std::exit(-1);
  }

  // Initialize MPI
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string reconstruction_method = argv[1];
  int Nx = std::stoi(argv[2]);

  auto start = std::chrono::system_clock::now();
  reconstructSurface(reconstruction_method, Nx);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = end - start;
  if (rank == 0) {
    printf("Total run time: %20f s\n\n", runtime.count());
  }

  MPI_Finalize();

  return 0;
}
