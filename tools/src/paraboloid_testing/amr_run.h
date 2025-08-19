// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_AMR_RUN_H_
#define SRC_AMR_RUN_H_

#include <mpi.h>
#include <array>
#include <fstream>
#include <string>
#include <vector>

static constexpr std::size_t RUN_BATCH_SIZE = 1000;

// Writes out a binary file with 6 doubles per case, organized as
// 1. k value (determines paraboloid positioning for case)
// 2. volume
// 3. centroid x
// 4. centroid y
// 5. centroid z
// 6. surface area
void runTranslatingCube(const double a_side_length,
                        const std::vector<std::vector<double>>& a_range_vector,
                        const std::vector<std::size_t>& a_step_vector,
                        std::fstream& a_out_file);

// Writes out a binary file with 12 doubles per case, organized as
// 1. translation x
// 2. translation y
// 3. translation z
// 4. Rotation amount (in radians) about x axis
// 5. Rotation amount (in radians) about y axis
// 6. Rotation amount (in radians) about z axis
// 7. a coefficient for paraboloid
// 8. b coefficient for paraboloid
// 9. volume
// 10. centroid x
// 11. centroid y
// 12. centroid z
void runRandomSweep(const std::string& a_geometry,
                    const std::vector<std::vector<double>>& a_translation,
                    const std::vector<std::vector<double>>& a_rotation,
                    const std::vector<std::vector<double>>& a_coefficient,
                    const bool a_fix_to_paraboloid,
                    const std::size_t a_number_of_tests,
                    std::fstream& a_out_file);

// Writes out a binary file with 12 doubles per case, organized as
// 1. translation x
// 2. translation y
// 3. translation z
// 4. Rotation amount (in radians) about x axis
// 5. Rotation amount (in radians) about y axis
// 6. Rotation amount (in radians) about z axis
// 7. a coefficient for paraboloid
// 8. b coefficient for paraboloid
// 9. volume
// 10. centroid x
// 11. centroid y
// 12. centroid z
void runOrganizedSweep(const std::string& a_geometry,
                       const std::vector<std::vector<double>>& a_translation,
                       const std::vector<std::size_t>& a_translation_steps,
                       const std::vector<std::vector<double>>& a_rotation,
                       const std::vector<std::size_t>& a_rotation_steps,
                       const std::vector<std::vector<double>>& a_coefficient,
                       const std::vector<std::size_t>& a_coefficient_steps,
                       std::fstream& a_out_file);

#endif  // SRC_AMR_RUN_H_
