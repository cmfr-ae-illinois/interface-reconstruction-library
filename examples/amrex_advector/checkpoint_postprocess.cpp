// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <string>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include "examples/amrex_advector/amrcore_advector.h"

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv);

  {
    BL_PROFILE("checkpoint_postprocess()");

    std::string initial_checkpoint;
    std::string final_checkpoint;
    {
      amrex::ParmParse pp("postprocess");
      pp.get("initial_checkpoint", initial_checkpoint);
      pp.get("final_checkpoint", final_checkpoint);
    }

    const auto start_total = amrex::second();

    AmrCoreAdv amr_core_adv;
    amr_core_adv.PostProcessCheckpointPair(initial_checkpoint,
                                           final_checkpoint);

    auto total_time = amrex::second() - start_total;
    amrex::ParallelDescriptor::ReduceRealMax(
        total_time, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "\nCheckpoint postprocess time: " << total_time << "\n";
  }

  amrex::Finalize();
}
