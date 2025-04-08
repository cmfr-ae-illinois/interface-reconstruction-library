// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <nlohmann/json.hpp>

#include "tools/src/paraboloid_testing/irl_run.h"

static constexpr std::size_t MAX_INPUT_STRING_SIZE = 2000;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout
        << "Must supply name for file to write results to as first argument.\n";
    std::exit(-1);
  }
  if (argc < 3) {
    std::cout
        << "Must supply list of results files generated from AMR executable"
        << std::endl;
    std::exit(-1);
  }

  std::vector<std::string> file_names(argc - 2);
  for (int i = 2; i < argc; ++i) {
    file_names[i - 2] = std::string(argv[i]);
  }

  auto out_file = fopen(argv[1], "w");
  for (const auto& file_name : file_names) {
    std::cout << "Confirming results from " << file_name << std::endl;
    std::fstream file(file_name, file.binary | file.in);

    std::size_t input_length;
    file.read(reinterpret_cast<char*>(&input_length), sizeof(std::size_t));
    std::string input_file;
    input_file.resize(input_length);
    file.read(input_file.data(), input_length);
    auto input = nlohmann::json::parse(input_file, nullptr, true, true);
    std::cout << "Input file: " << input_file << std::endl;

    const auto test_name = input["test_name"].get<std::string>();
    // Functions should return L1, L2, Linf, and total time
    std::size_t total_cases = 0;
    std::array<double, 10> results;
    if (test_name == "face_only") {
    } else if (test_name == "basic_cube") {
    } else if (test_name == "translating_cube") {
      const auto side_length =
          input["TranslatingCube"]["side_length"].get<double>();
      const auto range_steps =
          input["TranslatingCube"]["steps"].get<std::vector<std::size_t>>();
      total_cases = 0;
      for (auto& steps : range_steps) {
        total_cases += steps;
      }
      results = confirmTranslatingCubeResult(side_length, total_cases, file);

    } else if (test_name == "parameter_sweep") {
      auto geometry = input["Sweep"]["geometry"].get<std::string>();

      if (input["Sweep"]["random"].get<bool>()) {
        total_cases = input["Sweep"]["number_of_tests"].get<std::size_t>();
      } else {
        const auto translation_steps =
            input["Sweep"]["translation_steps"].get<std::vector<std::size_t>>();
        const auto rotation_steps =
            input["Sweep"]["rotation_steps"].get<std::vector<std::size_t>>();
        const auto coefficient_steps =
            input["Sweep"]["coefficient_steps"].get<std::vector<std::size_t>>();
        total_cases = translation_steps[0] * translation_steps[1] *
                      translation_steps[2] * rotation_steps[0] *
                      rotation_steps[1] * rotation_steps[2] *
                      coefficient_steps[0] * coefficient_steps[1];
      }
      std::cout << "Checking parameter sweep for geometry " << geometry
                << " with " << total_cases << " cases" << std::endl;
      results = confirmParameterSweepResult(geometry, total_cases, file);
    } else {
      std::cout << "Unknown test type \"" + test_name + "\"" << std::endl;
      std::exit(-1);
    }
    fprintf(out_file, "%s %10i", file_name.c_str(),
            static_cast<int>(total_cases));
    for (const auto& entry : results) {
      fprintf(out_file, " %15.8E", entry);
    }
    fprintf(out_file, "\n");
    std::cout << '\n';
  }
  fclose(out_file);
}
