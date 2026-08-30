// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <array>
#include <sstream>
#include <string>

#include "examples/implicit_surface_reconstruction/binary.h"
#include "examples/implicit_surface_reconstruction/surface_select.h"

namespace {

using ExactMomentData =
    Data<std::pair<IRL::GeneralMoments3D<2>, IRL::GeneralSurfaceMoments3D<2>>>;

int numberOfMomentComponents(const int moment_order) {
  if (moment_order < 0 || moment_order > 2) {
    std::ostringstream oss;
    oss << "moment order must be 0, 1, or 2; got " << moment_order;
    amrex::Abort(oss.str());
  }
  return (moment_order + 1) * (moment_order + 2) * (moment_order + 3) / 6;
}

amrex::Geometry makeGeometry(const BasicMesh& mesh) {
  const int ncell = mesh.getNx();
  const amrex::IntVect domain_lo(0, 0, 0);
  const amrex::IntVect domain_hi(ncell - 1, ncell - 1, ncell - 1);
  const amrex::Box domain(domain_lo, domain_hi);
  const amrex::RealBox real_box(
      {AMREX_D_DECL(mesh.x(mesh.imin()), mesh.y(mesh.jmin()),
                    mesh.z(mesh.kmin()))},
      {AMREX_D_DECL(mesh.x(mesh.imax() + 1), mesh.y(mesh.jmax() + 1),
                    mesh.z(mesh.kmax() + 1))});
  const int coord = 0;
  const std::array<int, AMREX_SPACEDIM> periodic{AMREX_D_DECL(1, 1, 1)};
  return amrex::Geometry(domain, &real_box, coord, periodic.data());
}

BasicMesh makeMesh(const int ncell, const std::string& shape) {
  BasicMesh mesh(ncell, ncell, ncell, 0);
  auto surface = makeSurface(shape, mesh);
  (void)surface;
  return mesh;
}

void fillExactMomentMultiFabs(const ExactMomentData& exact_moments,
                              amrex::MultiFab& exact_volume_moments,
                              amrex::MultiFab& exact_surface_moments) {
  const int volume_ncomp = exact_volume_moments.nComp();
  const int surface_ncomp = exact_surface_moments.nComp();

  for (amrex::MFIter mfi(exact_volume_moments); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    auto volume_arr = exact_volume_moments.array(mfi);
    auto surface_arr = exact_surface_moments.array(mfi);
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    amrex::ParallelFor(
        bx, [=, &exact_moments] AMREX_GPU_DEVICE(int i, int j, int k) {
          const auto& volume_moments = exact_moments(i, j, k).first;
          const auto& surface_moments = exact_moments(i, j, k).second;

          for (int n = 0; n < volume_ncomp; ++n) {
            volume_arr(i, j, k, n) = volume_moments[n];
          }
          for (int n = 0; n < surface_ncomp; ++n) {
            surface_arr(i, j, k, n) = surface_moments[n];
          }
        });
  }
}

void recenterMomentMultiFab(amrex::MultiFab& moments,
                            const amrex::Geometry& geom) {
  const int ncomp = moments.nComp();
  if (ncomp != 1 && ncomp != 4 && ncomp != 10) {
    std::ostringstream oss;
    oss << "Moment MultiFab must have 1, 4, or 10 components; got " << ncomp;
    amrex::Abort(oss.str());
  }

  if (ncomp == 1) return;

  const auto dx = geom.CellSizeArray();
  const auto problo = geom.ProbLoArray();

  for (amrex::MFIter mfi(moments); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    auto arr = moments.array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const double xc = problo[0] + (static_cast<double>(i) + 0.5) * dx[0];
      const double yc = problo[1] + (static_cast<double>(j) + 0.5) * dx[1];
      const double zc = problo[2] + (static_cast<double>(k) + 0.5) * dx[2];

      const double m0 = arr(i, j, k, 0);
      const double mx = arr(i, j, k, 1);
      const double my = arr(i, j, k, 2);
      const double mz = arr(i, j, k, 3);

      arr(i, j, k, 1) = mx - m0 * xc;
      arr(i, j, k, 2) = my - m0 * yc;
      arr(i, j, k, 3) = mz - m0 * zc;

      if (ncomp == 4) return;

      const double mxx = arr(i, j, k, 4);
      const double mxy = arr(i, j, k, 5);
      const double mxz = arr(i, j, k, 6);
      const double myy = arr(i, j, k, 7);
      const double myz = arr(i, j, k, 8);
      const double mzz = arr(i, j, k, 9);

      arr(i, j, k, 4) = mxx - 2.0 * xc * mx + m0 * xc * xc;
      arr(i, j, k, 5) = mxy - xc * my - yc * mx + m0 * xc * yc;
      arr(i, j, k, 6) = mxz - xc * mz - zc * mx + m0 * xc * zc;
      arr(i, j, k, 7) = myy - 2.0 * yc * my + m0 * yc * yc;
      arr(i, j, k, 8) = myz - yc * mz - zc * my + m0 * yc * zc;
      arr(i, j, k, 9) = mzz - 2.0 * zc * mz + m0 * zc * zc;
    });
  }
}

void printMomentSums(const amrex::MultiFab& moments, const std::string& name) {
  const int ncomp = moments.nComp();
  std::vector<double> sums(ncomp, 0.0);

  for (amrex::MFIter mfi(moments); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const auto arr = moments.const_array(mfi);
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          for (int n = 0; n < ncomp; ++n) {
            sums[n] += arr(i, j, k, n);
          }
        }
      }
    }
  }

  for (int n = 0; n < ncomp; ++n) {
    amrex::ParallelDescriptor::ReduceRealSum(
        sums[n], amrex::ParallelDescriptor::IOProcessorNumber());
  }

  amrex::Print() << "\n" << name << " moment sums before recentering:\n";
  for (int n = 0; n < ncomp; ++n) {
    amrex::Print() << "  component " << n << " = " << std::scientific
                   << std::setprecision(16) << sums[n] << "\n";
  }
}

void printUsage() {
  amrex::Print() << "amrex_reconstruction_convergence inputs:\n"
                 << "  binary_file = path/to/exact_moments.bin\n"
                 << "  shape = sphere|ellipsoid|genus|orthocircle\n"
                 << "  nx_fine = fine binary resolution\n"
                 << "Optional:\n"
                 << "  factor = 1\n"
                 << "  max_grid_size = 32\n"
                 << "  moment_order = 2\n"
                 << "  volume_moment_order = moment_order\n"
                 << "  surface_moment_order = moment_order\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv);

  {
    std::string binary_file;
    std::string shape;
    int nx_fine = -1;
    int factor = 1;
    int max_grid_size = 32;
    int moment_order = 2;
    int volume_moment_order = -1;
    int surface_moment_order = -1;

    amrex::ParmParse pp;
    pp.query("binary_file", binary_file);
    pp.query("shape", shape);
    pp.query("nx_fine", nx_fine);
    pp.query("factor", factor);
    pp.query("max_grid_size", max_grid_size);
    pp.query("moment_order", moment_order);
    pp.query("volume_moment_order", volume_moment_order);
    pp.query("surface_moment_order", surface_moment_order);

    if (volume_moment_order < 0) volume_moment_order = moment_order;
    if (surface_moment_order < 0) surface_moment_order = moment_order;

    if (binary_file.empty() || shape.empty() || nx_fine <= 0 || factor <= 0 ||
        nx_fine % factor != 0) {
      printUsage();
      amrex::Abort(
          "Provide binary_file, positive nx_fine, and factor dividing "
          "nx_fine.");
    }

    // extracting exact moment data for choice of resolution
    const int ncell = nx_fine / factor;
    BasicMesh mesh = makeMesh(ncell, shape);
    ExactMomentData exact_moments(&mesh);
    coarsenMomentsFromBinary<2, 2>(binary_file, factor, &exact_moments);

    // creating multifab for exact volumetric and surface moments
    const amrex::Geometry geom = makeGeometry(mesh);
    amrex::BoxArray ba(geom.Domain());
    ba.maxSize(max_grid_size);
    const amrex::DistributionMapping dm(ba);
    const int volume_ncomp = numberOfMomentComponents(volume_moment_order);
    const int surface_ncomp = numberOfMomentComponents(surface_moment_order);
    amrex::MultiFab exact_volume_moments(ba, dm, volume_ncomp, 0);
    amrex::MultiFab exact_surface_moments(ba, dm, surface_ncomp, 0);
    exact_volume_moments.setVal(0.0);
    exact_surface_moments.setVal(0.0);

    // filling multifab with exact moments
    fillExactMomentMultiFabs(exact_moments, exact_volume_moments,
                             exact_surface_moments);

    printMomentSums(exact_volume_moments, "Exact volume");
    printMomentSums(exact_surface_moments, "Exact surface");

    // recentering moments about each cell center
    recenterMomentMultiFab(exact_volume_moments, geom);
    recenterMomentMultiFab(exact_surface_moments, geom);

    amrex::Print() << "Filled and recentered exact moment MultiFabs on "
                   << ncell << "^3 mesh\n"
                   << "  volume components = " << volume_ncomp << "\n"
                   << "  surface components = " << surface_ncomp << "\n"
                   << "  domain = [" << mesh.x(mesh.imin()) << ", "
                   << mesh.x(mesh.imax() + 1) << "] x [" << mesh.y(mesh.jmin())
                   << ", " << mesh.y(mesh.jmax() + 1) << "] x ["
                   << mesh.z(mesh.kmin()) << ", " << mesh.z(mesh.kmax() + 1)
                   << "]\n";
  }

  amrex::Finalize();
  return 0;
}
