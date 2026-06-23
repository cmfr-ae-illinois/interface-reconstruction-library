// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <mpi.h>
#include <iostream>

#include "examples/variant_advector/vof_advection.h"

void resetMoments(
    const Data<IRL::LocalizedSeparatorVariantLink>& a_link_localized_interface,
    Data<IRL::VolumeMoments>* a_liq_moments,
    Data<IRL::VolumeMoments>* a_gas_moments) {
  using MomentType = IRL::SeparatedMoments<IRL::VolumeMoments>;
  const BasicMesh& mesh = a_link_localized_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto interface =
            a_link_localized_interface(i, j, k).getNextReconstruction();
        const auto moments = IRL::getVolumeMoments<MomentType>(cell, interface);
        (*a_liq_moments)(i, j, k).centroid() = moments[0].centroid();
        (*a_gas_moments)(i, j, k).centroid() = moments[1].centroid();
      }
    }
  }
}

std::array<int, 3> getIndexFromTag(const BasicMesh& a_mesh,
                                   const IRL::UnsignedIndex_t a_tag) {
  std::array<int, 3> indices;
  auto int_tag = static_cast<int64_t>(a_tag);
  indices[0] = static_cast<int>(int_tag / (a_mesh.getNzo() * a_mesh.getNyo()));
  indices[1] = static_cast<int>(
      (int_tag - indices[0] * (a_mesh.getNzo() * a_mesh.getNyo())) /
      a_mesh.getNzo());
  indices[2] =
      static_cast<int>(int_tag - a_mesh.getNzo() * indices[1] -
                       (indices[0] * (a_mesh.getNzo() * a_mesh.getNyo())));

  return {indices[0] + a_mesh.imino(), indices[1] + a_mesh.jmino(),
          indices[2] + a_mesh.kmino()};
}

void advectVOF(
    const std::string& a_advection_method,
    const std::string& a_reconstruction_method, const double a_dt,
    const double a_time, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface,
    Data<IRL::VolumeMoments>* a_liq_moments,
    Data<IRL::VolumeMoments>* a_gas_moments) {
  if (a_advection_method == "FullLagrangianCorrected" ||
      a_advection_method == "Default") {
    FullLagrangianCorrected::advectVOF(
        a_reconstruction_method, a_dt, a_time, a_U, a_V, a_W,
        a_link_localized_interface, a_liq_moments, a_gas_moments);
  } else {
    std::cout << "Unknown advection method of : " << a_advection_method << '\n';
    std::cout
        << "Value entries are:\n    FullLagrangianCorrected (Default). \n";
    std::exit(-1);
  }
}

const IRL::Polyhedron24 constructCorrectedPreimage(
    const int i, const int j, const int k, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W) {
  const BasicMesh& mesh = a_U.getMesh();
  std::array<IRL::Pt, 14> cell, preimage;
  IRL::CappedDodecahedron flux;
  // Initialize cell corners
  cell[0] = IRL::Pt(mesh.x(i + 1), mesh.y(j), mesh.z(k + 1));
  cell[1] = IRL::Pt(mesh.x(i + 1), mesh.y(j), mesh.z(k));
  cell[2] = IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k));
  cell[3] = IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
  cell[4] = IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k + 1));
  cell[5] = IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k));
  cell[6] = IRL::Pt(mesh.x(i), mesh.y(j + 1), mesh.z(k));
  cell[7] = IRL::Pt(mesh.x(i), mesh.y(j + 1), mesh.z(k + 1));
  // Compute the back projected cell corners
  for (int n = 0; n < 8; ++n) {
    preimage[n] = project_vertex(cell[n], -a_dt, a_U, a_V, a_W);
  }
  // Compute soleinoidal flux volumes
  std::array<double, 6> flux_volumes;
  flux_volumes[0] =
      0.5 * a_dt * (a_U(i - 1, j, k) + a_U(i, j, k)) * mesh.dy() * mesh.dz();
  flux_volumes[1] =
      0.5 * a_dt * (a_U(i, j, k) + a_U(i + 1, j, k)) * mesh.dy() * mesh.dz();
  flux_volumes[2] =
      0.5 * a_dt * (a_V(i, j - 1, k) + a_V(i, j, k)) * mesh.dx() * mesh.dz();
  flux_volumes[3] =
      0.5 * a_dt * (a_V(i, j, k) + a_V(i, j + 1, k)) * mesh.dx() * mesh.dz();
  flux_volumes[4] =
      0.5 * a_dt * (a_W(i, j, k - 1) + a_W(i, j, k)) * mesh.dx() * mesh.dy();
  flux_volumes[5] =
      0.5 * a_dt * (a_W(i, j, k) + a_W(i, j, k + 1)) * mesh.dx() * mesh.dy();
  // Create face flux hexahedra to compute correction
  for (int f = 0; f < 6; f++) {
    for (int i = 0; i < 4; i++) {
      flux[i] = cell[flux_id_table[f][i]];
      flux[i + 4] = preimage[flux_id_table[f][i]];
    }
    flux[8] = project_vertex(0.25 * (flux[0] + flux[1] + flux[2] + flux[3]),
                             -a_dt, a_U, a_V, a_W);
    flux.adjustCapToMatchVolume(flux_volumes[f]);
    preimage[face_center_id_table[f]] = flux[8];
  }
  return IRL::Polyhedron24::fromRawPtPointer(14, preimage.data());
};

void FullLagrangianCorrected::advectVOF(
    const std::string& a_reconstruction_method, const double a_dt,
    const double a_time, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface,
    Data<IRL::VolumeMoments>* a_liq_moments,
    Data<IRL::VolumeMoments>* a_gas_moments) {
  const BasicMesh& mesh = a_liq_moments->getMesh();

  // Reset centroid locationt to enforce realizability and compatibility with
  // the reconstructed interface
  resetMoments(*a_link_localized_interface, a_liq_moments, a_gas_moments);

  // Calculate maximum CFL number (will be needed to compute band)
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double CFLx = std::abs(a_U(i, j, k)) * a_dt / mesh.dx();
        const double CFLy = std::abs(a_V(i, j, k)) * a_dt / mesh.dy();
        const double CFLz = std::abs(a_W(i, j, k)) * a_dt / mesh.dz();
        CFL = std::fmax(CFL, std::fmax(CFLx, std::fmax(CFLy, CFLz)));
      }
    }
  }

  // Initialize band with mixed cells
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        band(i, j, k) = 0;
        const double liquid_volume_fraction =
            (*a_liq_moments)(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          band(i, j, k) = 1;
        }
      }
    }
  }
  band.updateBorder();

  // Extend band based on CFL number
  const int nlayers = static_cast<int>(std::ceil(CFL));
  for (int n = 0; n < nlayers; ++n) {
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
          if (band(i, j, k) == 0) {
            for (int ii = -1; ii <= 1; ++ii) {
              for (int jj = -1; jj <= 1; ++jj) {
                for (int kk = -1; kk <= 1; ++kk) {
                  if (band(i + ii, j + jj, k + kk) == n + 1) {
                    band(i, j, k) = n + 2;
                  }
                }
              }
            }
          }
        }
      }
    }
    band.updateBorder();
  }

  int nmixed_global = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (band(i, j, k) > 0) {
          nmixed_global++;
        }
      }
    }
  }

  // Share mixed cells between processors
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int nmixed_local = std::max(nmixed_global / size, 1);
  std::vector<int> proc_offset(size + 1);
  proc_offset[0] = 0;
  for (int r = 0; r < size; r++)
    proc_offset[r + 1] = proc_offset[r] + nmixed_local;
  proc_offset[size] = nmixed_global;
  for (int r = 1; r < size + 1; r++)
    proc_offset[r] = std::min(proc_offset[r], nmixed_global);
  nmixed_local = proc_offset[rank + 1] - proc_offset[rank];

  const int size_moments = 2 * 4;  // M0 and M1 for gas + liquid
  std::vector<double> moments_local(size_moments * nmixed_local);
  std::vector<double> moments_global(size_moments * nmixed_global);

  int count = 0, count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (band(i, j, k) > 0) {
          if (count >= proc_offset[rank] && count < proc_offset[rank + 1]) {
            // Construct corrected preimage
            const auto preimage =
                constructCorrectedPreimage(i, j, k, a_dt, a_U, a_V, a_W);
            // Now perform the actual cutting.
            const auto moments = IRL::getVolumeMoments<
                IRL::SeparatedMoments<IRL::VolumeMoments>>(
                preimage, (*a_link_localized_interface)(i, j, k));
            IRL::Pt l_centroid =
                moments[0].centroid() *
                (1.0 / IRL::safelyEpsilon(moments[0].volume()));
            IRL::Pt g_centroid =
                moments[1].centroid() *
                (1.0 / IRL::safelyEpsilon(moments[1].volume()));
            l_centroid = project_vertex(l_centroid, a_dt, a_U, a_V, a_W);
            g_centroid = project_vertex(g_centroid, a_dt, a_U, a_V, a_W);
            moments_local[count_local++] = moments[0].volume();
            moments_local[count_local++] = moments[1].volume();
            for (int d = 0; d < 3; d++) {
              moments_local[count_local++] =
                  l_centroid[d] * moments[0].volume();
              moments_local[count_local++] =
                  g_centroid[d] * moments[1].volume();
            }
          }
          count++;
        }
      }
    }
  }

  std::vector<int> proc_count(size);
  for (int r = 0; r < size; r++) {
    proc_count[r] = size_moments * (proc_offset[r + 1] - proc_offset[r]);
    proc_offset[r] = size_moments * proc_offset[r];
  }
  MPI_Allgatherv(moments_local.data(), moments_local.size(), MPI_DOUBLE,
                 moments_global.data(), proc_count.data(), proc_offset.data(),
                 MPI_DOUBLE, MPI_COMM_WORLD);

  count = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (band(i, j, k) > 0) {
          (*a_liq_moments)(i, j, k).volume() = moments_global[count++];
          (*a_gas_moments)(i, j, k).volume() = moments_global[count++];
          for (int d = 0; d < 3; d++) {
            (*a_liq_moments)(i, j, k).centroid()[d] = moments_global[count++];
            (*a_gas_moments)(i, j, k).centroid()[d] = moments_global[count++];
          }
        }
      }
    }
  }

  a_liq_moments->updateBorder();
  a_gas_moments->updateBorder();
  correctMomentLocations(a_liq_moments, a_gas_moments);
}

void correctMomentLocations(Data<IRL::VolumeMoments>* a_liq_moments,
                            Data<IRL::VolumeMoments>* a_gas_moments) {
  const BasicMesh& mesh = (*a_liq_moments).getMesh();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liq_moments)(i, j, k).centroid()[0] -=
            (*a_liq_moments)(i, j, k).volume() * mesh.lx();
        (*a_gas_moments)(i, j, k).centroid()[0] -=
            (*a_gas_moments)(i, j, k).volume() * mesh.lx();
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liq_moments)(i, j, k).centroid()[0] +=
            (*a_liq_moments)(i, j, k).volume() * mesh.lx();
        (*a_gas_moments)(i, j, k).centroid()[0] +=
            (*a_gas_moments)(i, j, k).volume() * mesh.lx();
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liq_moments)(i, j, k).centroid()[1] -=
            (*a_liq_moments)(i, j, k).volume() * mesh.ly();
        (*a_gas_moments)(i, j, k).centroid()[1] -=
            (*a_gas_moments)(i, j, k).volume() * mesh.ly();
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liq_moments)(i, j, k).centroid()[1] +=
            (*a_liq_moments)(i, j, k).volume() * mesh.ly();
        (*a_gas_moments)(i, j, k).centroid()[1] +=
            (*a_gas_moments)(i, j, k).volume() * mesh.ly();
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        (*a_liq_moments)(i, j, k).centroid()[2] -=
            (*a_liq_moments)(i, j, k).volume() * mesh.lz();
        (*a_gas_moments)(i, j, k).centroid()[2] -=
            (*a_gas_moments)(i, j, k).volume() * mesh.lz();
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        (*a_liq_moments)(i, j, k).centroid()[2] +=
            (*a_liq_moments)(i, j, k).volume() * mesh.lz();
        (*a_gas_moments)(i, j, k).centroid()[2] +=
            (*a_gas_moments)(i, j, k).volume() * mesh.lz();
      }
    }
  }
}