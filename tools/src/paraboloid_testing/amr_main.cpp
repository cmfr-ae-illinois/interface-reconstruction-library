// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <omp.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include <nlohmann/json.hpp>

#include "tools/src/paraboloid_testing/amr_run.h"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "Must supply name of input file on command line" << std::endl;
    std::exit(-1);
  }

  std::ifstream input_file(argv[1]);
  if (!input_file.good()) {
    std::cout << "Problem reading input file" << std::endl;
    std::exit(-1);
  }

  auto input = nlohmann::json::parse(input_file, nullptr, true, true);

  const auto test_name = input["test_name"].get<std::string>();
  const auto file_name = input["out_file_name"].get<std::string>();
  std::fstream file;

  /* MPI setup */
  int rank, size;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double start_time = omp_get_wtime();

  if (rank == 0) {
    file.open(file_name, file.binary | file.trunc | file.out);
    std::string json_string = input.dump();
    std::size_t string_length = json_string.size();
    file.write(reinterpret_cast<char*>(&string_length), sizeof(std::size_t));
    file.write(json_string.c_str(), string_length);
  }

  if (test_name == "face_only") {
  } else if (test_name == "basic_cube") {
  } else if (test_name == "translating_cube") {
    const auto side_length =
        input["TranslatingCube"]["side_length"].get<double>();
    const auto range_vector = input["TranslatingCube"]["ranges"]
                                  .get<std::vector<std::vector<double>>>();
    const auto range_steps =
        input["TranslatingCube"]["steps"].get<std::vector<std::size_t>>();
    if (range_vector.size() != range_steps.size()) {
      std::cout << "A number of steps must be given for each range in the "
                   "supplied \"ranges\" vector.\n";
      std::exit(-1);
    }
    for (auto& range : range_vector) {
      if (range.size() != 2) {
        std::cout << "Each provided range in the \"ranges\" vector must be "
                     "exactly two entries.\n";
        std::exit(-1);
      }
    }

    runTranslatingCube(side_length, range_vector, range_steps, file);

  } else if (test_name == "parameter_sweep") {
    const auto geometry = input["Sweep"]["geometry"].get<std::string>();

    const auto translation =
        input["Sweep"]["translation"].get<std::vector<std::vector<double>>>();
    if (translation.size() != 2 || translation[0].size() != 3 ||
        translation[1].size() != 3) {
      std::cout << "[\"Sweep\"][\"translation\"] must specify a starting "
                   "and ending point as two vectors, such as "
                   "[[0.0,0.0,0.0],[1.0,0.0,1.0]] \n";
      std::exit(-1);
    }

    const auto rotation =
        input["Sweep"]["rotation"].get<std::vector<std::vector<double>>>();
    if (rotation.size() != 2 || rotation[0].size() != 3 ||
        rotation[1].size() != 3) {
      std::cout << "[\"Sweep\"][\"rotation\"] must specify a starting "
                   "and ending rotation (in radians) as two vectors, such as "
                   "[[-3.14,-3.14,-3.14],[3.14,3.14,3.14]] \n";
      std::exit(-1);
    }

    const auto coefficient =
        input["Sweep"]["coefficient"].get<std::vector<std::vector<double>>>();
    if (coefficient.size() != 2 || coefficient[0].size() != 2 ||
        coefficient[1].size() != 2) {
      std::cout << "[\"Sweep\"][\"coefficient\"] must specify a starting "
                   "and ending rotation for the a and b value of the "
                   "paraboloid, such as "
                   "[[-5.0,-5.0],[5.0,5.0]] \n";
      std::exit(-1);
    }

    const auto random = input["Sweep"]["random"].get<bool>();
    if (random) {
      const auto fix_to_paraboloid =
          input["Sweep"]["fix_to_paraboloid"].get<bool>();
      const auto number_of_tests =
          input["Sweep"]["number_of_tests"].get<std::size_t>();
      runRandomSweep(geometry, translation, rotation, coefficient,
                     fix_to_paraboloid, number_of_tests, file);

    } else {
      const auto translation_steps =
          input["Sweep"]["translation_steps"].get<std::vector<std::size_t>>();
      const auto rotation_steps =
          input["Sweep"]["rotation_steps"].get<std::vector<std::size_t>>();
      const auto coefficient_steps =
          input["Sweep"]["coefficient_steps"].get<std::vector<std::size_t>>();

      if (translation_steps.size() != 3 || rotation_steps.size() != 3 ||
          coefficient_steps.size() != 2) {
        std::cout << "The step counts for translation, rotation, and "
                     "coefficient must be 3, 3, and 2, respectively.\n";
        std::exit(-1);
      }
      runOrganizedSweep(geometry, translation, translation_steps, rotation,
                        rotation_steps, coefficient, coefficient_steps, file);
    }
  } else {
    std::cout << "Unkown test type \"" + test_name + "\"" << std::endl;
    std::exit(-1);
  }

  double end_time = omp_get_wtime();

  if (rank == 0) {
    file.close();
    std::cout << "Results generated in " << end_time - start_time << " seconds"
              << std::endl;
  }

  MPI_Finalize();
}
