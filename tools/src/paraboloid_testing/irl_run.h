// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_IRL_RUN_H_
#define SRC_IRL_RUN_H_

#include <array>
#include <fstream>
#include <string>
#include <vector>

static constexpr std::size_t RUN_BATCH_SIZE = 1000;

std::array<double, 10> confirmTranslatingCubeResult(
    const double a_side_length, const std::size_t a_run_size,
    std::fstream& a_results_file);

std::array<double, 10> confirmParameterSweepResult(
    const std::string& a_geometry, const std::size_t a_run_size,
    std::fstream& a_results_file);

template <class GeometryType>
std::array<double, 10> performParameterSweep(const std::size_t a_run_size,
                                             std::fstream& a_results_file);

#include "tools/src/paraboloid_testing/irl_run.tpp"

#endif  // SRC_IRL_RUN_H_
