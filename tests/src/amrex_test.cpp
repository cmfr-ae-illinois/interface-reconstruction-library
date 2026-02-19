// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/amrex/sepunion_multifab.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(AMReX, PeriodicGhostUpdate) {
  int argc = 0;
  char** argv = nullptr;
  amrex::Initialize(argc, argv);

  // Create box array and distribution mapping
  const int nx = 32, ny = 48, nz = 64;
  amrex::IntVect dom_lo(0, 0, 0), dom_hi(nx - 1, ny - 1, nz - 1);
  amrex::Box domain(dom_lo, dom_hi);
  amrex::Array<int, 3> is_periodic{1, 1, 1};
  amrex::RealBox real_box({0., 0., 0.}, {1., 1., 1.});
  amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian,
                       is_periodic);
  amrex::BoxArray ba(domain);
  ba.maxSize(16);
  amrex::DistributionMapping dm(ba);

  int ncomp = 1;  // one SeparatorUnion per cell
  int ngrow = 2;  // 2 layer of ghost cells

  // Create MultiFab of SeparatorUnions
  amrex::SepUnionMultiFab sepu_fab(ba, dm, ncomp, ngrow);

  // Fill in MultiFab with Paraboloids, with datum = (i,j,k)
  for (amrex::MFIter mfi(sepu_fab); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    auto& fab = sepu_fab[mfi];

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          SeparatorUnion& cell = fab(amrex::IntVect(i, j, k), 0);
          cell = Paraboloid();
          cell.getParaboloid().setDatum(Pt(i, j, k));
        }
      }
    }
  }

  // Update ghost layers (default is COPY)
  sepu_fab.FillBoundary(geom.periodicity());

  // Verify that ghost layers have been updating correctly
  for (amrex::MFIter mfi(sepu_fab); mfi.isValid(); ++mfi) {
    auto& fab = sepu_fab[mfi];
    const amrex::Box& bx = fab.box();

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          // Compute indices of parent cell
          const int ii = ((i % nx) + nx) % nx;
          const int jj = ((j % ny) + ny) % ny;
          const int kk = ((k % nz) + nz) % nz;
          SeparatorUnion& cell = fab(amrex::IntVect(i, j, k), 0);
          const auto datum = cell.getParaboloid().getDatum();
          EXPECT_EQ(cell.type(), SeparatorUnion::SeparatorType::Paraboloid);
          EXPECT_NEAR(datum[0], static_cast<double>(ii), 1.0e-15);
          EXPECT_NEAR(datum[1], static_cast<double>(jj), 1.0e-15);
          EXPECT_NEAR(datum[2], static_cast<double>(kk), 1.0e-15);
        }
      }
    }
  }

  amrex::Finalize();
}
}  // namespace
