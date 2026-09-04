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
      amrex::Print() << "\n               Reconstruction Time: "
                     << std::scientific << std::setprecision(2)
                     << amr_core_adv.RecTime() << '\n';
      double max_rec_time_loop = amr_core_adv.RecLoopTime();
      amrex::ParallelDescriptor::ReduceRealMax(
          max_rec_time_loop, amrex::ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "      Max reconstruction Loop Time: "
                     << std::scientific << std::setprecision(2)
                     << max_rec_time_loop << '\n';
      const double rec_time_per_mixed_cell = amr_core_adv.RecTimePerMixedCell();
      amrex::Print() << "Reconstruction Time per Mixed cell: "
                     << std::scientific << std::setprecision(2)
                     << rec_time_per_mixed_cell << '\n';
      const double niter_mof2_per_mixed_cell =
          amr_core_adv.MOF2IterPerMixedCell();
      if (niter_mof2_per_mixed_cell > 0.0) {
        amrex::Print() << " MOF2 Num. of Iter. per Mixed cell: " << std::fixed
                       << std::setprecision(2) << niter_mof2_per_mixed_cell
                       << '\n';
      }
      amrex::Print() << "                    Advection Time: "
                     << std::scientific << std::setprecision(2)
                     << amr_core_adv.AdvTime() << '\n';
      amrex::Print() << "                   Total Wall Time: "
                     << std::scientific << std::setprecision(2) << end_total
                     << "\n\n";
    }
  }

  amrex::Finalize();
}
