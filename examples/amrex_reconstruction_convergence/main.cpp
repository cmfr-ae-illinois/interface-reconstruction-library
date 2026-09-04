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
#include <AMReX_Utility.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "irl/amrex/sepunion_multifab.h"
#include "irl/generic_cutting/generic_cutting_definitions.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"

static constexpr int comp_vf = 0;
static constexpr int comp_m0 = 1;
static constexpr int comp_m1_l = 2;
static constexpr int comp_m1_g = 5;
static constexpr int comp_m2_l = 8;
static constexpr int comp_m2_g = 14;

struct InterfaceScalarField {
  std::string name;
  amrex::MultiFab polygon_scalar_data;
  amrex::MultiFab paraboloid_scalar_data;
  std::vector<double> flattened_polygon_scalar_data;
  std::vector<double> flattened_paraboloid_scalar_data;

  InterfaceScalarField() = default;

  InterfaceScalarField(const std::string& a_name, const amrex::BoxArray& a_ba,
                       const amrex::DistributionMapping& a_dm,
                       const int a_ngrow = 0)
      : name(a_name),
        polygon_scalar_data(a_ba, a_dm, 1, a_ngrow),
        paraboloid_scalar_data(a_ba, a_dm, 1, a_ngrow) {
    polygon_scalar_data.setVal(0.0);
    paraboloid_scalar_data.setVal(0.0);
  }

  void clearFlattenedData() {
    flattened_polygon_scalar_data.clear();
    flattened_paraboloid_scalar_data.clear();
  }
};

#include "examples/amrex_advector/reconstruction_cf.h"
#include "examples/amrex_advector/reconstruction_elvira.h"
#include "examples/amrex_advector/reconstruction_hybrid.h"
#include "examples/amrex_advector/reconstruction_hybrid2.h"
#include "examples/amrex_advector/reconstruction_ivf.h"
#include "examples/amrex_advector/reconstruction_lvira.h"
#include "examples/amrex_advector/reconstruction_mof1.h"
#include "examples/amrex_advector/reconstruction_mof2.h"
#include "examples/amrex_advector/reconstruction_plicnet.h"
#include "examples/amrex_advector/reconstruction_pu.h"
#include "examples/amrex_advector/reconstruction_vf.h"
#include "examples/amrex_advector/reconstruction_vf2.h"
#include "examples/implicit_surface_reconstruction/binary.h"
#include "examples/implicit_surface_reconstruction/surface_select.h"

namespace {

constexpr int NUM_MOMENTS = 10;
using MomentArray = std::array<double, NUM_MOMENTS>;

enum class MomentType { Volume, Surface };

struct MomentErrorNorms {
  double m0_linf = 0.0;
  double m0_l1 = 0.0;
  double m0_l2 = 0.0;
  double m1_linf = 0.0;
  double m1_l1 = 0.0;
  double m1_l2 = 0.0;
  double m2_linf = 0.0;
  double m2_l1 = 0.0;
  double m2_l2 = 0.0;
};

int numberOfMomentComponents(const int moment_order) {
  if (moment_order < 0 || moment_order > 2) {
    std::ostringstream oss;
    oss << "moment order must be 0, 1, or 2; got " << moment_order;
    amrex::Abort(oss.str());
  }
  return (moment_order + 1) * (moment_order + 2) * (moment_order + 3) / 6;
}

int reconstructionRequiredVolumeOrder(const std::string& method) {
  if (method == "mof2" || method == "supermof2") return 2;
  if (method == "mof" || method == "mof1" || method == "plicnet") return 1;
  return 0;
}

int advectorMomentComponents(const int volume_order) {
  if (volume_order >= 2) return 20;
  if (volume_order >= 1) return 8;
  return 2;
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

IRL::RectangularCuboid makeCell(const amrex::Geometry& geom, const int i,
                                const int j, const int k) {
  const auto dx = geom.CellSizeArray();
  const auto problo = geom.ProbLoArray();
  const double x = problo[0] + static_cast<double>(i) * dx[0];
  const double y = problo[1] + static_cast<double>(j) * dx[1];
  const double z = problo[2] + static_cast<double>(k) * dx[2];
  return IRL::RectangularCuboid::fromBoundingPts(
      IRL::Pt(x, y, z), IRL::Pt(x + dx[0], y + dx[1], z + dx[2]));
}

void copyVolumeMomentsToArray(const IRL::VolumeMoments& moments,
                              MomentArray& result) {
  result[0] = moments.volume();
  result[1] = moments.centroid()[0];
  result[2] = moments.centroid()[1];
  result[3] = moments.centroid()[2];
}

template <class MomentsType>
void copyGeneralMomentsToArray(const MomentsType& moments,
                               MomentArray& result) {
  for (IRL::UnsignedIndex_t m = 0; m < moments.size(); ++m) {
    result[m] = moments[m];
  }
}

MomentArray getFullCellVolumeMoments(const IRL::RectangularCuboid& cell,
                                     const int max_order) {
  MomentArray result{};
  if (max_order == 0) {
    result[0] = IRL::getVolumeMoments<IRL::Volume>(cell);
  } else if (max_order == 1) {
    copyVolumeMomentsToArray(IRL::getVolumeMoments<IRL::VolumeMoments>(cell),
                             result);
  } else {
    copyGeneralMomentsToArray(
        IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell), result);
  }
  return result;
}

void coarsenMomentTypeFromBinaryToMultiFab(const std::string& binary_file,
                                           const int factor,
                                           const BasicMesh& mesh,
                                           const MomentType moment_type,
                                           amrex::MultiFab& exact_moments) {
  const int ncomp = exact_moments.nComp();
  if (ncomp > NUM_MOMENTS) {
    amrex::Abort("Only moment orders up to 2 are supported.");
  }

  std::ifstream stream(binary_file, std::ios::binary);
  if (!stream) {
    amrex::Abort("Cannot open binary moment file: " + binary_file);
  }

  const auto header = sparse_moment_io::readHeader<2, 2>(&stream);
  const InsideCellMask mask = sparse_moment_io::readMask(&stream, header);
  const auto records = sparse_moment_io::readMixedCells<2, 2>(&stream, header);

  if (header.nx != static_cast<std::uint32_t>(factor * mesh.getNx()) ||
      header.ny != static_cast<std::uint32_t>(factor * mesh.getNy()) ||
      header.nz != static_cast<std::uint32_t>(factor * mesh.getNz())) {
    amrex::Abort(
        "Fine binary dimensions must equal factor times coarse dimensions.");
  }

  exact_moments.setVal(0.0);

  const double fine_dx = mesh.lx() / static_cast<double>(header.nx);
  const double fine_dy = mesh.ly() / static_cast<double>(header.ny);
  const double fine_dz = mesh.lz() / static_cast<double>(header.nz);
  const double x_lower = mesh.x(mesh.imin());
  const double y_lower = mesh.y(mesh.jmin());
  const double z_lower = mesh.z(mesh.kmin());

  for (amrex::MFIter mfi(exact_moments); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    auto moment_arr = exact_moments.array(mfi);
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    if (moment_type == MomentType::Volume) {
      // Full liquid fine cells are stored only in the mask, so rebuild their
      // moments geometrically and accumulate them into the owning coarse cell.
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
          for (int i = lo.x; i <= hi.x; ++i) {
            for (int fi = i * factor; fi < (i + 1) * factor; ++fi) {
              for (int fj = j * factor; fj < (j + 1) * factor; ++fj) {
                for (int fk = k * factor; fk < (k + 1) * factor; ++fk) {
                  const std::size_t fine_index = getLinearCellIndex(
                      fi, fj, fk, static_cast<int>(header.ny),
                      static_cast<int>(header.nz));
                  if (!mask.get(fine_index)) continue;

                  const IRL::Pt lower(x_lower + fi * fine_dx,
                                      y_lower + fj * fine_dy,
                                      z_lower + fk * fine_dz);
                  const IRL::Pt upper(lower.x() + fine_dx, lower.y() + fine_dy,
                                      lower.z() + fine_dz);
                  const auto fine_cell =
                      IRL::RectangularCuboid::fromBoundingPts(lower, upper);
                  const auto moments =
                      IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
                          fine_cell);
                  for (int n = 0; n < ncomp; ++n) {
                    moment_arr(i, j, k, n) += moments[n];
                  }
                }
              }
            }
          }
        }
      }
    }

    // Mixed cells carry explicit exact moments in the binary records.
    for (const auto& record : records) {
      int fine_i, fine_j, fine_k;
      getCellIndicesFromLinearIndex(
          record.linear_index, static_cast<int>(header.ny),
          static_cast<int>(header.nz), &fine_i, &fine_j, &fine_k);
      const int i = fine_i / factor;
      const int j = fine_j / factor;
      const int k = fine_k / factor;
      if (!bx.contains(amrex::IntVect(AMREX_D_DECL(i, j, k)))) continue;

      const double* source =
          moment_type == MomentType::Volume ? record.volume : record.surface;
      for (int n = 0; n < ncomp; ++n) {
        moment_arr(i, j, k, n) += source[n];
      }
    }
  }
}

void initializeSepUnionMultiFab(amrex::SepUnionMultiFab& mf) {
  for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    auto arr = mf.array(mfi);
    const amrex::Box& bx = mfi.growntilebox();
    const int ncomp = mf.nComp();
    amrex::ParallelFor(bx, ncomp,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                         arr(i, j, k, n) = IRL::SeparatorUnion();
                       });
  }
}

void fillAdvectorMomentsFromExactVolume(const amrex::MultiFab& exact_volume,
                                        const amrex::Geometry& geom,
                                        amrex::MultiFab& adv_moments) {
  const int exact_ncomp = exact_volume.nComp();
  const int adv_ncomp = adv_moments.nComp();
  const auto dx = geom.CellSizeArray();
  const double cell_vol = dx[0] * dx[1] * dx[2];

  adv_moments.setVal(0.0);

  for (amrex::MFIter mfi(exact_volume); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const auto exact = exact_volume.const_array(mfi);
    auto adv = adv_moments.array(mfi);
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          const double m0 = exact(i, j, k, 0);
          adv(i, j, k, comp_vf) = m0 / cell_vol;
          adv(i, j, k, comp_m0) = m0;

          MomentArray full = getFullCellVolumeMoments(makeCell(geom, i, j, k),
                                                      adv_ncomp >= 20  ? 2
                                                      : adv_ncomp >= 8 ? 1
                                                                       : 0);
          if (adv_ncomp >= 8) {
            if (exact_ncomp < 4) {
              amrex::Abort(
                  "Selected reconstruction requires first volume moments.");
            }
            for (int n = 0; n < 3; ++n) {
              const double liquid = exact(i, j, k, 1 + n);
              adv(i, j, k, comp_m1_l + n) = liquid;
              adv(i, j, k, comp_m1_g + n) = full[1 + n] - liquid;
            }
          }
          if (adv_ncomp >= 20) {
            if (exact_ncomp < 10) {
              amrex::Abort(
                  "Selected reconstruction requires second volume moments.");
            }
            for (int n = 0; n < 6; ++n) {
              const double liquid = exact(i, j, k, 4 + n);
              adv(i, j, k, comp_m2_l + n) = liquid;
              adv(i, j, k, comp_m2_g + n) = full[4 + n] - liquid;
            }
          }
        }
      }
    }
  }

  adv_moments.FillBoundary(geom.periodicity());
}

void reconstructInterface(const std::string& method,
                          amrex::SepUnionMultiFab& interface,
                          amrex::SepUnionMultiFab& interface_with_ghost,
                          const amrex::MultiFab& moments,
                          const amrex::Geometry& geom,
                          amrex::Real* reconstruction_loop_time) {
  if (method == "elvira" || method == "default") {
    ELVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom,
                              nullptr, reconstruction_loop_time);
  } else if (method == "lvira") {
    LVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom,
                             nullptr, reconstruction_loop_time);
  } else if (method == "plicnet") {
    PLICNet::GetReconstruction(interface, interface_with_ghost, moments, geom,
                               nullptr, reconstruction_loop_time);
  } else if (method == "mof" || method == "mof1") {
    MOF1::GetReconstruction(interface, interface_with_ghost, moments, geom,
                            nullptr, reconstruction_loop_time);
  } else if (method == "vf") {
    VF::GetReconstruction(interface, interface_with_ghost, moments, geom,
                          nullptr, reconstruction_loop_time);
  } else if (method == "vf2") {
    VF2::GetReconstruction(interface, interface_with_ghost, moments, geom,
                           nullptr, reconstruction_loop_time);
  } else if (method == "ivf") {
    iVF::GetReconstruction(interface, interface_with_ghost, moments, geom,
                           nullptr, reconstruction_loop_time);
  } else if (method == "pu") {
    PU::GetReconstruction(interface, interface_with_ghost, moments, geom,
                          nullptr, reconstruction_loop_time);
  } else if (method == "cf") {
    CF::GetReconstruction(interface, interface_with_ghost, moments, geom,
                          nullptr, reconstruction_loop_time);
  } else if (method == "mof2") {
    MOF2::GetReconstruction(interface, interface_with_ghost, moments, geom,
                            nullptr, reconstruction_loop_time);
  } else if (method == "supermof2") {
    SuperMOF2::GetReconstruction(interface, interface_with_ghost, moments, geom,
                                 nullptr, reconstruction_loop_time);
  } else if (method == "hybrid") {
    HYBRID::GetReconstruction(interface, interface_with_ghost, moments, geom,
                              nullptr, reconstruction_loop_time);
  } else if (method == "hybrid2") {
    HYBRID2::GetReconstruction(interface, interface_with_ghost, moments, geom,
                               nullptr);
  } else {
    amrex::Abort("Unknown reconstruction method: " + method);
  }
}

amrex::Long countLocalMixedCells(const amrex::MultiFab& moments) {
  amrex::Long local_mixed_cells = 0;

  for (amrex::MFIter mfi(moments); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto moments_array = moments.const_array(mfi);
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          const double vf = moments_array(i, j, k, comp_vf);
          if (vf > IRL::global_constants::VF_LOW &&
              vf < IRL::global_constants::VF_HIGH) {
            ++local_mixed_cells;
          }
        }
      }
    }
  }

  return local_mixed_cells;
}

MomentArray recenteredMoments(const MomentArray& moments, const IRL::Pt& xc) {
  MomentArray centered{};

  const double x = xc[0];
  const double y = xc[1];
  const double z = xc[2];
  const double M0 = moments[0];
  const double Mx = moments[1];
  const double My = moments[2];
  const double Mz = moments[3];

  centered[0] = M0;
  centered[1] = Mx - M0 * x;
  centered[2] = My - M0 * y;
  centered[3] = Mz - M0 * z;
  centered[4] = moments[4] - 2.0 * x * Mx + M0 * x * x;
  centered[5] = moments[5] - x * My - y * Mx + M0 * x * y;
  centered[6] = moments[6] - x * Mz - z * Mx + M0 * x * z;
  centered[7] = moments[7] - 2.0 * y * My + M0 * y * y;
  centered[8] = moments[8] - y * Mz - z * My + M0 * y * z;
  centered[9] = moments[9] - 2.0 * z * Mz + M0 * z * z;

  return centered;
}

MomentArray calculatePolygonSurfaceMoments(const IRL::Polygon& polygon,
                                           const int max_order) {
  MomentArray moments{};
  const int n = polygon.getNumberOfVertices();
  if (n < 3) return moments;

  const IRL::Pt& v0 = polygon[0];
  for (int iv = 1; iv < n - 1; ++iv) {
    const IRL::Pt& v1 = polygon[iv];
    const IRL::Pt& v2 = polygon[iv + 1];

    IRL::Polygon triangle;
    triangle.addVertex(v0);
    triangle.addVertex(v1);
    triangle.addVertex(v2);
    triangle.calculateAndSetPlaneOfExistence();

    const double area = triangle.calculateVolume();
    moments[0] += area;
    if (max_order == 0) continue;

    const IRL::Pt centroid = (v0 + v1 + v2) * (1.0 / 3.0);
    moments[1] += centroid[0] * area;
    moments[2] += centroid[1] * area;
    moments[3] += centroid[2] * area;
    if (max_order == 1) continue;

    const auto a = v1 - v0;
    const auto b = v2 - v0;
    const double factor = 2.0 * area;

    auto accumulate = [&](const IRL::Pt& u, const IRL::Pt& v,
                          const double scale) {
      moments[4] += scale * u[0] * v[0];
      moments[5] += scale * u[0] * v[1];
      moments[6] += scale * u[0] * v[2];
      moments[7] += scale * u[1] * v[1];
      moments[8] += scale * u[1] * v[2];
      moments[9] += scale * u[2] * v[2];
    };

    accumulate(v0, v0, factor * 0.5);
    accumulate(v0, a, factor / 6.0);
    accumulate(a, v0, factor / 6.0);
    accumulate(v0, b, factor / 6.0);
    accumulate(b, v0, factor / 6.0);
    accumulate(a, a, factor / 12.0);
    accumulate(b, b, factor / 12.0);
    accumulate(a, b, factor / 24.0);
    accumulate(b, a, factor / 24.0);
  }

  return moments;
}

MomentArray getReconstructedVolumeMoments(const IRL::RectangularCuboid& cell,
                                          const IRL::SeparatorUnion& interface,
                                          const double vf,
                                          const int max_order) {
  MomentArray result{};
  if (vf <= IRL::global_constants::VF_LOW) return result;
  if (vf >= IRL::global_constants::VF_HIGH) {
    return getFullCellVolumeMoments(cell, max_order);
  }

  if (interface.type() == IRL::SeparatorUnion::SeparatorType::OnePlane) {
    const auto separator =
        IRL::PlanarSeparator::fromOnePlane(interface.getPlane());
    if (max_order == 0) {
      result[0] = IRL::getVolumeMoments<IRL::Volume>(cell, separator);
    } else if (max_order == 1) {
      copyVolumeMomentsToArray(
          IRL::getVolumeMoments<IRL::VolumeMoments>(cell, separator), result);
    } else {
      copyGeneralMomentsToArray(
          IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, separator),
          result);
    }
    return result;
  }

  if (interface.type() == IRL::SeparatorUnion::SeparatorType::Paraboloid) {
    const auto paraboloid = interface.getParaboloid();
    if (max_order == 0) {
      result[0] = IRL::getVolumeMoments<IRL::Volume>(cell, paraboloid);
    } else if (max_order == 1) {
      copyVolumeMomentsToArray(
          IRL::getVolumeMoments<IRL::VolumeMoments>(cell, paraboloid), result);
    } else {
      copyGeneralMomentsToArray(
          IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid),
          result);
    }
    return result;
  }

  amrex::Abort("Mixed cell has unsupported interface type.");
  return result;
}

MomentArray getReconstructedSurfaceMoments(const IRL::RectangularCuboid& cell,
                                           const IRL::SeparatorUnion& interface,
                                           const double vf,
                                           const int max_order) {
  MomentArray result{};
  if (vf <= IRL::global_constants::VF_LOW ||
      vf >= IRL::global_constants::VF_HIGH) {
    return result;
  }

  if (interface.type() == IRL::SeparatorUnion::SeparatorType::OnePlane) {
    const auto separator =
        IRL::PlanarSeparator::fromOnePlane(interface.getPlane());
    IRL::Polygon polygon = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
        cell, separator, separator[0]);
    if (polygon.getNumberOfVertices() > 2) {
      polygon.calculateAndSetPlaneOfExistence();
      result = calculatePolygonSurfaceMoments(polygon, max_order);
    }
    return result;
  }

  if (interface.type() == IRL::SeparatorUnion::SeparatorType::Paraboloid) {
    using VolumeMomentsAndSurface =
        IRL::AddSurfaceOutput<IRL::Volume,
                              IRL::ParaboloidParametrizedSurfaceOutput>;
    const auto paraboloid = interface.getParaboloid();
    auto volume_and_surface =
        IRL::getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid);
    auto surface = volume_and_surface.getSurface();
    if (max_order == 0) {
      copyGeneralMomentsToArray(surface.getSurfaceMoments<0>(), result);
    } else if (max_order == 1) {
      copyGeneralMomentsToArray(surface.getSurfaceMoments<1>(), result);
    } else {
      copyGeneralMomentsToArray(surface.getSurfaceMoments<2>(), result);
    }
    return result;
  }

  amrex::Abort("Mixed cell has unsupported surface interface type.");
  return result;
}

IRL::MixedPolygonBezierSurface buildInterfaceSurface(
    const amrex::SepUnionMultiFab& interface, const amrex::MultiFab& moments,
    const amrex::Geometry& geom) {
  IRL::MixedPolygonBezierSurface surface;
  const auto problo = geom.ProbLoArray();
  const auto dx = geom.CellSizeArray();

  for (amrex::MFIter mfi(moments, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const auto moments_fab = moments.const_array(mfi);
    const auto interface_fab = interface.const_array(mfi);
    const amrex::Box& box = mfi.tilebox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int k = lo.z; k <= hi.z; ++k) {
      const double z = problo[2] + k * dx[2];
      for (int j = lo.y; j <= hi.y; ++j) {
        const double y = problo[1] + j * dx[1];
        for (int i = lo.x; i <= hi.x; ++i) {
          const double vf = moments_fab(i, j, k, comp_vf);
          if (vf <= IRL::global_constants::VF_LOW ||
              vf >= IRL::global_constants::VF_HIGH) {
            continue;
          }

          const double x = problo[0] + i * dx[0];
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(x, y, z), IRL::Pt(x + dx[0], y + dx[1], z + dx[2]));

          if (interface_fab(i, j, k).type() ==
              IRL::SeparatorUnion::SeparatorType::OnePlane) {
            const auto planar_sep = IRL::PlanarSeparator::fromOnePlane(
                interface_fab(i, j, k).getPlane());
            const auto polygon =
                IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                    cell, planar_sep, planar_sep[0]);
            if (polygon.getNumberOfVertices() > 2) {
              surface.addPolygon(polygon);
            }
          } else if (interface_fab(i, j, k).type() ==
                     IRL::SeparatorUnion::SeparatorType::Paraboloid) {
            using VolumeAndSurface =
                IRL::AddSurfaceOutput<IRL::Volume,
                                      IRL::ParaboloidParametrizedSurfaceOutput>;
            const auto paraboloid = interface_fab(i, j, k).getParaboloid();
            auto volume_and_surface =
                IRL::getVolumeMoments<VolumeAndSurface>(cell, paraboloid);
            const double area =
                volume_and_surface.getSurface().getSurfaceArea();
            const double cell_dx = std::cbrt(dx[0] * dx[1] * dx[2]);
            if (area > 1.0e-4 * cell_dx * cell_dx &&
                area < 10.0 * cell_dx * cell_dx) {
              surface.addSurface(volume_and_surface.getSurface()
                                     .getQuadraticBezierTriangleApprox());
            }
          }
        }
      }
    }
  }

  return surface;
}

std::string joinPath(const std::string& directory, const std::string& name) {
  if (directory.empty() || directory == "." || name.empty() ||
      name.front() == '/') {
    return name;
  }
  if (directory.back() == '/') {
    return directory + name;
  }
  return directory + "/" + name;
}

void outputReconstructedInterface(const amrex::SepUnionMultiFab& interface,
                                  const amrex::MultiFab& moments,
                                  const amrex::Geometry& geom,
                                  const std::string& output_dir,
                                  const std::string& output_name) {
  IRL::MixedPolygonBezierSurface surface =
      buildInterfaceSurface(interface, moments, geom);

  std::string name = output_name;
  const std::string suffix = ".vtu";
  if (name.size() >= suffix.size() &&
      name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
    name.resize(name.size() - suffix.size());
  }

  const std::string file_base = joinPath(output_dir, name);
  const std::string filename = file_base + ".vtu";
  const int rank = amrex::ParallelDescriptor::MyProc();
  const int size = amrex::ParallelDescriptor::NProcs();

  const auto& points = surface.getPointList();
  const auto& polygons = surface.getPolygonList();
  const auto& bezier_triangles = surface.getBezierTriangleList();
  const int number_of_points = static_cast<int>(points.size());
  const int number_of_polygons = static_cast<int>(polygons.first.size());
  const int number_of_triangles = static_cast<int>(bezier_triangles.size());

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::UtilCreateDirectory(output_dir, 0755);
    FILE* file = std::fopen(filename.c_str(), "w");
    if (file == nullptr) {
      amrex::FileOpenFailed(filename);
    }
    std::fprintf(file, "<?xml version=\"1.0\"?>\n");
    std::fprintf(file,
                 "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
                 "byte_order=\"LittleEndian\">\n");
    std::fprintf(file, "  <UnstructuredGrid>\n");
    std::fclose(file);
  }
  amrex::ParallelDescriptor::Barrier();

  for (int r = 0; r < size; ++r) {
    if (rank == r) {
      FILE* file = std::fopen(filename.c_str(), "a");
      if (file == nullptr) {
        amrex::FileOpenFailed(filename);
      }
      std::fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
                   number_of_points, number_of_triangles + number_of_polygons);
      std::fprintf(file,
                   "<Points>\n<DataArray type=\"Float64\" "
                   "NumberOfComponents=\"3\">\n");
      for (IRL::UnsignedIndex_t i = 0; i < points.size(); ++i) {
        std::fprintf(file, "%15.8E %15.8E %15.8E ", std::get<0>(points[i])[0],
                     std::get<0>(points[i])[1], std::get<0>(points[i])[2]);
      }
      std::fprintf(file, "\n</DataArray>\n</Points>\n");
      std::fprintf(file,
                   "<PointData RationalWeights=\"RationalWeights\">\n"
                   "<DataArray type=\"Float64\" Name=\"RationalWeights\" "
                   "format=\"ascii\">\n");
      for (IRL::UnsignedIndex_t i = 0; i < points.size(); ++i) {
        std::fprintf(file, "%15.8E ", std::get<1>(points[i]));
      }
      std::fprintf(file, "\n</DataArray>\n</PointData>\n");

      std::fprintf(file, "<Cells>\n");
      std::fprintf(file,
                   "<DataArray type=\"Int64\" Name=\"connectivity\" "
                   "format=\"ascii\">\n");
      for (IRL::UnsignedIndex_t i = 0; i < bezier_triangles.size(); ++i) {
        for (IRL::UnsignedIndex_t j = 0; j < bezier_triangles[i].size(); ++j) {
          std::fprintf(file, "%d ", bezier_triangles[i][j]);
        }
      }
      for (IRL::UnsignedIndex_t i = 0; i < polygons.second.size(); ++i) {
        std::fprintf(file, "%d ", polygons.second[i]);
      }
      std::fprintf(file, "\n</DataArray>\n");

      std::fprintf(file,
                   "<DataArray type=\"Int64\" Name=\"offsets\" "
                   "format=\"ascii\">\n");
      IRL::UnsignedIndex_t count = 0;
      for (IRL::UnsignedIndex_t i = 0; i < bezier_triangles.size(); ++i) {
        count += bezier_triangles[i].size();
        std::fprintf(file, "%d ", count);
      }
      for (IRL::UnsignedIndex_t i = 0; i < polygons.first.size(); ++i) {
        count += polygons.first[i];
        std::fprintf(file, "%d ", count);
      }
      std::fprintf(file, "\n</DataArray>\n");
      std::fprintf(file,
                   "<DataArray type=\"UInt8\" Name=\"types\" "
                   "format=\"ascii\">\n");
      for (IRL::UnsignedIndex_t i = 0; i < bezier_triangles.size(); ++i) {
        std::fprintf(file, "76 ");
      }
      for (IRL::UnsignedIndex_t i = 0; i < polygons.first.size(); ++i) {
        std::fprintf(file, "7 ");
      }
      std::fprintf(file, "\n</DataArray>\n");
      std::fprintf(file, "</Cells>\n");
      std::fprintf(file, "</Piece>\n");
      std::fclose(file);
    }
    amrex::ParallelDescriptor::Barrier();
  }

  if (amrex::ParallelDescriptor::IOProcessor()) {
    FILE* file = std::fopen(filename.c_str(), "a");
    if (file == nullptr) {
      amrex::FileOpenFailed(filename);
    }
    std::fprintf(file, "  </UnstructuredGrid>\n</VTKFile>\n");
    std::fclose(file);
    amrex::Print() << "\nWrote reconstructed interface: " << filename << "\n";
  }
}

MomentErrorNorms computeAndPrintMomentErrorNorms(
    const amrex::MultiFab& exact_moments,
    const amrex::SepUnionMultiFab& interface_with_ghost,
    const amrex::MultiFab& adv_moments, const amrex::Geometry& geom,
    const MomentType moment_type, const int moment_order,
    const std::string& label) {
  const int ncomp = numberOfMomentComponents(moment_order);
  const auto dx = geom.CellSizeArray();
  const auto problo = geom.ProbLoArray();

  double local_linf_m0 = 0.0;
  double local_linf_m1 = 0.0;
  double local_linf_m2 = 0.0;
  double local_l1_m0 = 0.0;
  double local_l1_m1 = 0.0;
  double local_l1_m2 = 0.0;
  double local_l2_m0 = 0.0;
  double local_l2_m1 = 0.0;
  double local_l2_m2 = 0.0;

  for (amrex::MFIter mfi(exact_moments); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const auto exact = exact_moments.const_array(mfi);
    const auto moments = adv_moments.const_array(mfi);
    const auto iface = interface_with_ghost.const_array(mfi);
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          const auto cell = makeCell(geom, i, j, k);
          const double vf = moments(i, j, k, comp_vf);
          MomentArray reconstructed =
              moment_type == MomentType::Volume
                  ? getReconstructedVolumeMoments(cell, iface(i, j, k), vf,
                                                  moment_order)
                  : getReconstructedSurfaceMoments(cell, iface(i, j, k), vf,
                                                   moment_order);

          MomentArray exact_array{};
          for (int n = 0; n < ncomp; ++n) {
            exact_array[n] = exact(i, j, k, n);
          }

          const IRL::Pt center(
              problo[0] + (static_cast<double>(i) + 0.5) * dx[0],
              problo[1] + (static_cast<double>(j) + 0.5) * dx[1],
              problo[2] + (static_cast<double>(k) + 0.5) * dx[2]);
          exact_array = recenteredMoments(exact_array, center);
          reconstructed = recenteredMoments(reconstructed, center);

          const double d0 = std::abs(reconstructed[0] - exact_array[0]);
          local_linf_m0 = std::max(local_linf_m0, d0);
          local_l1_m0 += d0;
          local_l2_m0 += d0 * d0;

          if (moment_order >= 1) {
            const double dx0 = reconstructed[1] - exact_array[1];
            const double dy0 = reconstructed[2] - exact_array[2];
            const double dz0 = reconstructed[3] - exact_array[3];
            const double e1 = std::sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
            local_linf_m1 = std::max(local_linf_m1, e1);
            local_l1_m1 += e1;
            local_l2_m1 += e1 * e1;
          }

          if (moment_order >= 2) {
            const double dxx = reconstructed[4] - exact_array[4];
            const double dxy = reconstructed[5] - exact_array[5];
            const double dxz = reconstructed[6] - exact_array[6];
            const double dyy = reconstructed[7] - exact_array[7];
            const double dyz = reconstructed[8] - exact_array[8];
            const double dzz = reconstructed[9] - exact_array[9];
            const double e2_sq = dxx * dxx + dyy * dyy + dzz * dzz +
                                 2.0 * (dxy * dxy + dxz * dxz + dyz * dyz);
            const double e2 = std::sqrt(e2_sq);
            local_linf_m2 = std::max(local_linf_m2, e2);
            local_l1_m2 += e2;
            local_l2_m2 += e2_sq;
          }
        }
      }
    }
  }

  amrex::ParallelDescriptor::ReduceRealMax(
      local_linf_m0, amrex::ParallelDescriptor::IOProcessorNumber());
  amrex::ParallelDescriptor::ReduceRealSum(
      local_l1_m0, amrex::ParallelDescriptor::IOProcessorNumber());
  amrex::ParallelDescriptor::ReduceRealSum(
      local_l2_m0, amrex::ParallelDescriptor::IOProcessorNumber());
  if (moment_order >= 1) {
    amrex::ParallelDescriptor::ReduceRealMax(
        local_linf_m1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
        local_l1_m1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
        local_l2_m1, amrex::ParallelDescriptor::IOProcessorNumber());
  }
  if (moment_order >= 2) {
    amrex::ParallelDescriptor::ReduceRealMax(
        local_linf_m2, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
        local_l1_m2, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
        local_l2_m2, amrex::ParallelDescriptor::IOProcessorNumber());
  }

  const double inv_ncells = 1.0 / static_cast<double>(geom.Domain().numPts());
  const double h = dx[0];
  const double m0_scale =
      moment_type == MomentType::Volume ? 1.0 / (h * h * h) : 1.0 / (h * h);
  const double m1_scale = m0_scale / h;
  const double m2_scale = m1_scale / h;

  MomentErrorNorms norms;
  norms.m0_linf = local_linf_m0 * m0_scale;
  norms.m0_l1 = local_l1_m0 * inv_ncells * m0_scale;
  norms.m0_l2 = std::sqrt(local_l2_m0 * inv_ncells) * m0_scale;
  if (moment_order >= 1) {
    norms.m1_linf = local_linf_m1 * m1_scale;
    norms.m1_l1 = local_l1_m1 * inv_ncells * m1_scale;
    norms.m1_l2 = std::sqrt(local_l2_m1 * inv_ncells) * m1_scale;
  }
  if (moment_order >= 2) {
    norms.m2_linf = local_linf_m2 * m2_scale;
    norms.m2_l1 = local_l1_m2 * inv_ncells * m2_scale;
    norms.m2_l2 = std::sqrt(local_l2_m2 * inv_ncells) * m2_scale;
  }

  amrex::Print() << "\n" << label << " reconstructed moment error norms\n";
  amrex::Print() << "  M0 Linf = " << std::scientific << std::setprecision(16)
                 << norms.m0_linf << "\n";
  amrex::Print() << "  M0 L1   = " << std::scientific << std::setprecision(16)
                 << norms.m0_l1 << "\n";
  amrex::Print() << "  M0 L2   = " << std::scientific << std::setprecision(16)
                 << norms.m0_l2 << "\n";
  if (moment_order >= 1) {
    amrex::Print() << "  M1 Linf = " << std::scientific << std::setprecision(16)
                   << norms.m1_linf << "\n";
    amrex::Print() << "  M1 L1   = " << std::scientific << std::setprecision(16)
                   << norms.m1_l1 << "\n";
    amrex::Print() << "  M1 L2   = " << std::scientific << std::setprecision(16)
                   << norms.m1_l2 << "\n";
  }
  if (moment_order >= 2) {
    amrex::Print() << "  M2 Linf = " << std::scientific << std::setprecision(16)
                   << norms.m2_linf << "\n";
    amrex::Print() << "  M2 L1   = " << std::scientific << std::setprecision(16)
                   << norms.m2_l1 << "\n";
    amrex::Print() << "  M2 L2   = " << std::scientific << std::setprecision(16)
                   << norms.m2_l2 << "\n";
  }

  return norms;
}

double relativeMomentError(const MomentArray& reconstructed,
                           const MomentArray& exact, const int order) {
  double error_squared = 0.0;
  double exact_squared = 0.0;
  int begin = 0;
  int end = 1;
  if (order == 1) {
    begin = 1;
    end = 4;
  } else if (order == 2) {
    begin = 4;
    end = 10;
  }

  for (int n = begin; n < end; ++n) {
    const double weight =
        order == 2 && (n == 5 || n == 6 || n == 8) ? 2.0 : 1.0;
    const double difference = reconstructed[n] - exact[n];
    error_squared += weight * difference * difference;
    exact_squared += weight * exact[n] * exact[n];
  }

  const double error = std::sqrt(error_squared);
  const double denominator = std::sqrt(exact_squared);
  if (denominator > std::numeric_limits<double>::epsilon()) {
    return error / denominator;
  }
  return error <= std::numeric_limits<double>::epsilon()
             ? 0.0
             : std::numeric_limits<double>::infinity();
}

void outputMixedCellMomentErrors(
    const std::string& filename, const amrex::MultiFab& exact_volume_moments,
    const amrex::MultiFab* exact_surface_moments,
    const amrex::SepUnionMultiFab& interface_with_ghost,
    const amrex::MultiFab& adv_moments, const amrex::Geometry& geom,
    const int volume_moment_order, const int surface_moment_order,
    const bool output_volume, const bool output_surface) {
  const int number_of_columns = 1 +
                                (output_volume ? volume_moment_order + 1 : 0) +
                                (output_surface ? surface_moment_order + 1 : 0);
  const auto dx = geom.CellSizeArray();
  const auto problo = geom.ProbLoArray();
  std::vector<double> local_rows;

  for (amrex::MFIter mfi(adv_moments); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto adv = adv_moments.const_array(mfi);
    const auto exact_volume = exact_volume_moments.const_array(mfi);
    const auto iface = interface_with_ghost.const_array(mfi);
    const auto exact_surface = output_surface
                                   ? exact_surface_moments->const_array(mfi)
                                   : amrex::Array4<const amrex::Real>{};
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          const double vf = adv(i, j, k, comp_vf);
          if (vf <= IRL::global_constants::VF_LOW ||
              vf >= IRL::global_constants::VF_HIGH) {
            continue;
          }

          const auto cell = makeCell(geom, i, j, k);
          const IRL::Pt center(
              problo[0] + (static_cast<double>(i) + 0.5) * dx[0],
              problo[1] + (static_cast<double>(j) + 0.5) * dx[1],
              problo[2] + (static_cast<double>(k) + 0.5) * dx[2]);
          local_rows.push_back(vf);

          if (output_volume) {
            MomentArray exact{};
            for (int n = 0; n < numberOfMomentComponents(volume_moment_order);
                 ++n) {
              exact[n] = exact_volume(i, j, k, n);
            }
            MomentArray reconstructed = getReconstructedVolumeMoments(
                cell, iface(i, j, k), vf, volume_moment_order);
            exact = recenteredMoments(exact, center);
            reconstructed = recenteredMoments(reconstructed, center);
            for (int order = 0; order <= volume_moment_order; ++order) {
              local_rows.push_back(
                  relativeMomentError(reconstructed, exact, order));
            }
          }

          if (output_surface) {
            MomentArray exact{};
            for (int n = 0; n < numberOfMomentComponents(surface_moment_order);
                 ++n) {
              exact[n] = exact_surface(i, j, k, n);
            }
            MomentArray reconstructed = getReconstructedSurfaceMoments(
                cell, iface(i, j, k), vf, surface_moment_order);
            exact = recenteredMoments(exact, center);
            reconstructed = recenteredMoments(reconstructed, center);
            for (int order = 0; order <= surface_moment_order; ++order) {
              local_rows.push_back(
                  relativeMomentError(reconstructed, exact, order));
            }
          }
        }
      }
    }
  }

  if (local_rows.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    amrex::Abort("Too many local mixed-cell error values for MPI_Gatherv.");
  }
  const int local_count = static_cast<int>(local_rows.size());
  const int root = amrex::ParallelDescriptor::IOProcessorNumber();
  const std::vector<int> counts =
      amrex::ParallelDescriptor::Gather(local_count, root);
  std::vector<int> offsets;
  std::vector<double> global_rows;
  if (amrex::ParallelDescriptor::IOProcessor()) {
    offsets.resize(counts.size(), 0);
    for (std::size_t rank = 1; rank < counts.size(); ++rank) {
      offsets[rank] = offsets[rank - 1] + counts[rank - 1];
    }
    const int global_count = offsets.back() + counts.back();
    global_rows.resize(static_cast<std::size_t>(global_count));
  }
  amrex::ParallelDescriptor::Gatherv(local_rows.data(), local_count,
                                     global_rows.data(), counts, offsets, root);

  if (!amrex::ParallelDescriptor::IOProcessor()) return;
  std::ofstream csv(filename, std::ios::trunc);
  if (!csv) amrex::FileOpenFailed(filename);
  csv << "volume_fraction";
  if (output_volume) {
    for (int order = 0; order <= volume_moment_order; ++order) {
      csv << ",volume_M" << order << "_relative_error";
    }
  }
  if (output_surface) {
    for (int order = 0; order <= surface_moment_order; ++order) {
      csv << ",surface_M" << order << "_relative_error";
    }
  }
  csv << '\n' << std::scientific << std::setprecision(16);
  for (std::size_t offset = 0; offset < global_rows.size();
       offset += static_cast<std::size_t>(number_of_columns)) {
    csv << global_rows[offset];
    for (int column = 1; column < number_of_columns; ++column) {
      csv << ',' << global_rows[offset + static_cast<std::size_t>(column)];
    }
    csv << '\n';
  }
  amrex::Print() << "  wrote mixed-cell moment errors: " << filename << "\n";
}

void appendMomentMetricsCsv(const std::string& csv_file,
                            const std::string& shape, const std::string& method,
                            const int nx_fine, const int factor, const int nx,
                            const std::string& moment_type,
                            const int moment_order,
                            const MomentErrorNorms& norms) {
  if (!amrex::ParallelDescriptor::IOProcessor() || csv_file.empty()) return;

  const bool write_header = !std::ifstream(csv_file).good();
  std::ofstream csv(csv_file, std::ios::app);
  if (!csv) {
    amrex::FileOpenFailed(csv_file);
  }

  if (write_header) {
    csv << "shape,method,nx_fine,factor,nx,moment_type,moment_order,"
        << "M0_Linf,M0_L1,M0_L2,M1_Linf,M1_L1,M1_L2,"
        << "M2_Linf,M2_L1,M2_L2\n";
  }

  csv << shape << ',' << method << ',' << nx_fine << ',' << factor << ',' << nx
      << ',' << moment_type << ',' << moment_order << ',' << std::scientific
      << std::setprecision(16) << norms.m0_linf << ',' << norms.m0_l1 << ','
      << norms.m0_l2 << ',' << norms.m1_linf << ',' << norms.m1_l1 << ','
      << norms.m1_l2 << ',' << norms.m2_linf << ',' << norms.m2_l1 << ','
      << norms.m2_l2 << "\n";
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
  amrex::Print()
      << "amrex_reconstruction_convergence inputs:\n"
      << "  binary_file = path/to/exact_moments.bin\n"
      << "  shape = sphere|ellipsoid|genus|orthocircle\n"
      << "  nx_fine = fine binary resolution\n"
      << "  reconstruction_method = vf\n"
      << "Optional:\n"
      << "  factor = 1\n"
      << "  factors = 1 2 4 8 16\n"
      << "  max_grid_size = 32\n"
      << "  moment_order = 2\n"
      << "  volume_moment_order = moment_order\n"
      << "  surface_moment_order = moment_order\n"
      << "  do_volume = 1\n"
      << "  do_surface = 1\n"
      << "  output_interface = 0\n"
      << "  interface_output_dir = .\n"
      << "  interface_output_name = interface_<shape>_<method>_f<factor>\n"
      << "  output_csv = 0\n"
      << "  csv_file = metrics.csv\n"
      << "  output_cell_errors = 0\n"
      << "  cell_error_output_dir = .\n"
      << "  cell_error_output_name = <shape>_<method>_cell_errors\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv);

  {
    std::string binary_file;
    std::string shape;
    std::string reconstruction_method = "vf";
    int nx_fine = -1;
    int factor = 1;
    std::vector<int> factors;
    int max_grid_size = 32;
    int moment_order = 2;
    int volume_moment_order = -1;
    int surface_moment_order = -1;
    int do_volume = 1;
    int do_surface = 1;
    int output_interface = 0;
    int output_csv = 0;
    int output_cell_errors = 0;
    std::string csv_file = "metrics.csv";
    std::string interface_output_dir = ".";
    std::string interface_output_name;
    std::string cell_error_output_dir = ".";
    std::string cell_error_output_name;

    amrex::ParmParse pp;
    pp.query("binary_file", binary_file);
    pp.query("shape", shape);
    pp.query("nx_fine", nx_fine);
    pp.query("factor", factor);
    pp.queryarr("factors", factors);
    pp.query("max_grid_size", max_grid_size);
    pp.query("moment_order", moment_order);
    pp.query("volume_moment_order", volume_moment_order);
    pp.query("surface_moment_order", surface_moment_order);
    pp.query("do_volume", do_volume);
    pp.query("do_surface", do_surface);
    pp.query("output_interface", output_interface);
    pp.query("output_csv", output_csv);
    pp.query("output_cell_errors", output_cell_errors);
    pp.query("csv_file", csv_file);
    pp.query("interface_output_dir", interface_output_dir);
    pp.query("interface_output_name", interface_output_name);
    pp.query("cell_error_output_dir", cell_error_output_dir);
    pp.query("cell_error_output_name", cell_error_output_name);
    pp.query("reconstruction_method", reconstruction_method);

    amrex::ParmParse pprec("reconstruction");
    pprec.query("name", reconstruction_method);

    if (factors.empty()) {
      factors.push_back(factor);
    }
    if (volume_moment_order < 0) volume_moment_order = moment_order;
    if (surface_moment_order < 0) surface_moment_order = moment_order;
    if (!do_volume && !do_surface && !output_interface && !output_cell_errors) {
      amrex::Abort("Enable at least one output or error calculation.");
    }
    if (output_cell_errors && !do_volume && !do_surface) {
      amrex::Abort(
          "output_cell_errors requires do_volume or do_surface to be enabled.");
    }

    if (binary_file.empty() || shape.empty() || nx_fine <= 0) {
      printUsage();
      amrex::Abort("Provide binary_file, shape, and positive nx_fine.");
    }

    if (output_csv && amrex::ParallelDescriptor::IOProcessor()) {
      std::ofstream csv(csv_file, std::ios::trunc);
      if (!csv) {
        amrex::FileOpenFailed(csv_file);
      }
      csv << "shape,method,nx_fine,factor,nx,moment_type,moment_order,"
          << "M0_Linf,M0_L1,M0_L2,M1_Linf,M1_L1,M1_L2,"
          << "M2_Linf,M2_L1,M2_L2\n";
    }

    for (const int current_factor : factors) {
      if (current_factor <= 0 || nx_fine % current_factor != 0) {
        std::ostringstream oss;
        oss << "Factor must be positive and divide nx_fine; got factor = "
            << current_factor << ", nx_fine = " << nx_fine;
        amrex::Abort(oss.str());
      }

      const int ncell = nx_fine / current_factor;
      BasicMesh mesh = makeMesh(ncell, shape);
      const amrex::Geometry geom = makeGeometry(mesh);
      amrex::BoxArray ba(geom.Domain());
      ba.maxSize(max_grid_size);
      const amrex::DistributionMapping dm(ba);

      const int required_volume_order =
          std::max(volume_moment_order,
                   reconstructionRequiredVolumeOrder(reconstruction_method));
      const int adv_ncomp = advectorMomentComponents(required_volume_order);
      const int num_grow = required_volume_order >= 2 ? 2 : 1;

      amrex::Print() << "\nAMReX reconstruction convergence setup\n"
                     << "  method = " << reconstruction_method << "\n"
                     << "  mesh = " << ncell << "^3\n"
                     << "  factor = " << current_factor << "\n"
                     << "  advector moment components = " << adv_ncomp << "\n";

      amrex::SepUnionMultiFab interface(ba, dm, 1, 0);
      amrex::SepUnionMultiFab interface_with_ghost(ba, dm, 1, num_grow);
      initializeSepUnionMultiFab(interface);
      initializeSepUnionMultiFab(interface_with_ghost);

      amrex::MultiFab adv_moments(ba, dm, adv_ncomp, num_grow);

      amrex::MultiFab exact_volume_for_reconstruction(
          ba, dm, numberOfMomentComponents(required_volume_order), 0);
      coarsenMomentTypeFromBinaryToMultiFab(binary_file, current_factor, mesh,
                                            MomentType::Volume,
                                            exact_volume_for_reconstruction);
      fillAdvectorMomentsFromExactVolume(exact_volume_for_reconstruction, geom,
                                         adv_moments);
      amrex::Real reconstruction_loop_time = 0.0;
      reconstructInterface(reconstruction_method, interface,
                           interface_with_ghost, adv_moments, geom,
                           &reconstruction_loop_time);

      const amrex::Long local_mixed_cells = countLocalMixedCells(adv_moments);
      if (local_mixed_cells == 0) reconstruction_loop_time = 0.0;

      amrex::Long global_mixed_cells = local_mixed_cells;
      amrex::Real global_reconstruction_loop_time = reconstruction_loop_time;
      amrex::ParallelDescriptor::ReduceLongSum(
          global_mixed_cells, amrex::ParallelDescriptor::IOProcessorNumber());
      amrex::ParallelDescriptor::ReduceRealSum(
          global_reconstruction_loop_time,
          amrex::ParallelDescriptor::IOProcessorNumber());

      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "  global mixed cells = " << global_mixed_cells
                       << "\n"
                       << "  summed reconstruction loop time = "
                       << global_reconstruction_loop_time << " s\n";
        if (global_mixed_cells > 0) {
          amrex::Print() << "  average reconstruction time per mixed cell = "
                         << global_reconstruction_loop_time /
                                static_cast<double>(global_mixed_cells)
                         << " s\n";
        } else {
          amrex::Print() << "  average reconstruction time per mixed cell = N/A"
                         << " (no mixed cells)\n";
        }
      }

      if (output_interface) {
        std::string current_interface_name = interface_output_name;
        if (current_interface_name.empty()) {
          current_interface_name = shape + "_" + reconstruction_method + "_Nx" +
                                   std::to_string(ncell) + "_interface";
        }
        outputReconstructedInterface(interface, adv_moments, geom,
                                     interface_output_dir,
                                     current_interface_name);
      }

      if (do_volume) {
        const MomentErrorNorms volume_norms = computeAndPrintMomentErrorNorms(
            exact_volume_for_reconstruction, interface_with_ghost, adv_moments,
            geom, MomentType::Volume, volume_moment_order, "Volume");
        if (output_csv) {
          appendMomentMetricsCsv(csv_file, shape, reconstruction_method,
                                 nx_fine, current_factor, ncell, "volume",
                                 volume_moment_order, volume_norms);
        }
      }

      amrex::MultiFab exact_surface_moments;
      if (do_surface) {
        exact_surface_moments.define(
            ba, dm, numberOfMomentComponents(surface_moment_order), 0);
        coarsenMomentTypeFromBinaryToMultiFab(binary_file, current_factor, mesh,
                                              MomentType::Surface,
                                              exact_surface_moments);
        const MomentErrorNorms surface_norms = computeAndPrintMomentErrorNorms(
            exact_surface_moments, interface_with_ghost, adv_moments, geom,
            MomentType::Surface, surface_moment_order, "Surface");
        if (output_csv) {
          appendMomentMetricsCsv(csv_file, shape, reconstruction_method,
                                 nx_fine, current_factor, ncell, "surface",
                                 surface_moment_order, surface_norms);
        }
      }

      if (output_cell_errors) {
        std::string output_name = cell_error_output_name;
        if (output_name.empty()) {
          output_name = shape + "_" + reconstruction_method + "_cell_errors";
        }
        std::string filename =
            output_name + "_f" + std::to_string(current_factor) + ".csv";
        if (!cell_error_output_dir.empty() && cell_error_output_dir != ".") {
          filename = cell_error_output_dir + "/" + filename;
        }
        outputMixedCellMomentErrors(
            filename, exact_volume_for_reconstruction,
            do_surface ? &exact_surface_moments : nullptr, interface_with_ghost,
            adv_moments, geom, volume_moment_order, surface_moment_order,
            do_volume != 0, do_surface != 0);
      }
    }
  }

  amrex::Finalize();
  return 0;
}
