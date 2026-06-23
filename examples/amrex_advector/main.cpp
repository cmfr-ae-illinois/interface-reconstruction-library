// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <mpi.h>
#include <chrono>
#include <iostream>
#include <string>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include "examples/amrex_advector/amrcore_advector.h"

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv);

  {
    // timer for profiling
    BL_PROFILE("main()");

    // wallclock time
    const auto strt_total = amrex::second();

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    AmrCoreAdv amr_core_adv;

    // initialize AMR data
    amr_core_adv.InitData();

    // advance solution to final time
    amr_core_adv.Evolve();

    // wallclock time
    auto end_total = amrex::second() - strt_total;

    if (amr_core_adv.Verbose()) {
      // print wallclock time
      amrex::ParallelDescriptor::ReduceRealMax(
          end_total, amrex::ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "\nReconstruction Time: " << amr_core_adv.RecTime()
                     << '\n';
      amrex::Print() << "     Advection Time: " << amr_core_adv.AdvTime()
                     << '\n';
      amrex::Print() << "         Total Time: " << end_total << '\n';
    }
  }

  amrex::Finalize();
}