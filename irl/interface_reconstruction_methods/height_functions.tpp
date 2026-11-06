// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_HEIGHT_FUNCTIONS_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_HEIGHT_FUNCTIONS_TPP_

#include "irl/interface_reconstruction_methods/elvira.h"
#include "irl/interface_reconstruction_methods/elvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/lvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/lvira_optimization.h"

namespace IRL {

template <int Nx, int Ny, int Nz, class SeparatorVariantType>
inline Normal
HeightFunctions_3D<Nx, Ny, Nz, SeparatorVariantType>::lvira_normal(
    const std::array<std::array<std::array<double, Nz>, Ny>, Nx>& a_vfrac) {
  static_assert(Nx >= 3 and Ny >= 3 and Nz >= 3,
                "LVIRA stencil needs a minimum of 3 cells in each direction");
  ELVIRANeighborhood elvira_neighborhood;
  LVIRANeighborhood<RectangularCuboid> lvira_neighborhood;
  elvira_neighborhood.resize(27);
  lvira_neighborhood.resize(27);
  std::array<RectangularCuboid, 27> cells;
  std::array<double, 27> cells_vfrac;
  lvira_neighborhood.setCenterOfStencil(1 * 9 + 1 * 3 + 1);
  const int istart = (Nx % 2 == 0) ? Nx / 2 - 2 : Nx / 2 - 1;
  const int jstart = (Ny % 2 == 0) ? Ny / 2 - 2 : Ny / 2 - 1;
  const int kstart = (Nz % 2 == 0) ? Nz / 2 - 2 : Nz / 2 - 1;
  for (int kk = 0; kk < 3; ++kk) {
    for (int jj = 0; jj < 3; ++jj) {
      for (int ii = 0; ii < 3; ++ii) {
        const int neigh_id = kk * 9 + jj * 3 + ii;
        int ilo = istart + ii;
        int ihi = ilo + 1;
        if (Nx % 2 == 0) {
          if (ii > 0) ihi += 1;
          if (ii > 1) ilo += 1;
        }
        int jlo = jstart + jj;
        int jhi = jlo + 1;
        if (Ny % 2 == 0) {
          if (jj > 0) jhi += 1;
          if (jj > 1) jlo += 1;
        }
        int klo = kstart + kk;
        int khi = klo + 1;
        if (Nz % 2 == 0) {
          if (kk > 0) khi += 1;
          if (kk > 1) klo += 1;
        }
        cells[neigh_id] = RectangularCuboid::fromBoundingPts(Pt(ilo, jlo, klo),
                                                             Pt(ihi, jhi, khi));
        cells_vfrac[neigh_id] = 0.0;
        for (int k = klo; k < khi; ++k) {
          for (int j = jlo; j < jhi; ++j) {
            for (int i = ilo; i < ihi; ++i) {
              cells_vfrac[neigh_id] += a_vfrac[i][j][k];
            }
          }
        }
        cells_vfrac[neigh_id] /=
            static_cast<double>((ihi - ilo) * (jhi - jlo) * (khi - klo));
        elvira_neighborhood.setMember(&cells[neigh_id], &cells_vfrac[neigh_id],
                                      ii - 1, jj - 1, kk - 1);
        lvira_neighborhood.setMember(static_cast<UnsignedIndex_t>(neigh_id),
                                     &cells[neigh_id], &cells_vfrac[neigh_id]);
      }
    }
  }
  ELVIRA_3D elvira_system;
  LVIRA_3D<RectangularCuboid> lvira_system;
  // Solve LVIRA problem with ELVIRA as initial guess
  const PlanarSeparator separator = lvira_system.solve(
      lvira_neighborhood, elvira_system.solve(&elvira_neighborhood));
  return separator[0].normal();
}

template <int Nx, int Ny, int Nz, class SeparatorVariantType>
SeparatorVariantType
HeightFunctions_3D<Nx, Ny, Nz, SeparatorVariantType>::solve(
    const std::array<std::array<std::array<double, Nz>, Ny>, Nx>& a_vfrac,
    const Pt a_stencil_center, const double a_dx, const Normal& a_normal) {
  static_assert(Nx >= 3 and Ny >= 3 and Nz >= 3,
                "Height functions need a minimum of 3 cells in each direction");
  Normal normal = a_normal;
  if (squaredMagnitude(normal) < 0.5) {
    // Compute normal using LVIRA
    normal = this->lvira_normal(a_vfrac);
  }

  // Find dominant normal direction
  int dir = 0;
  if (std::fabs(normal[dir]) < std::fabs(normal[1])) dir = 1;
  if (std::fabs(normal[dir]) < std::fabs(normal[2])) dir = 2;

  // Compute heights
  int sign = normal[dir] >= 0.0 ? 1 : -1;
  int nvalid = this->compute_heights(a_vfrac, dir, sign);
  if (nvalid < 6) {
    const int dir1 = (dir + 1) % 3;
    const int dir2 = (dir + 2) % 3;
    if (std::fabs(normal[dir1]) > std::fabs(normal[dir2])) {
      dir = dir1;
      sign = normal[dir] >= 0.0 ? 1 : -1;
      nvalid = this->compute_heights(a_vfrac, dir, sign);
      if (nvalid < 6) {
        dir = dir2;
        sign = normal[dir] >= 0.0 ? 1 : -1;
        nvalid = this->compute_heights(a_vfrac, dir, sign);
      }
    } else {
      dir = dir2;
      sign = normal[dir] >= 0.0 ? 1 : -1;
      nvalid = this->compute_heights(a_vfrac, dir, sign);
      if (nvalid < 6) {
        dir = dir1;
        sign = normal[dir] >= 0.0 ? 1 : -1;
        nvalid = this->compute_heights(a_vfrac, dir, sign);
      }
    }
  }

  // Reconstruct paraboloid from heights
  if (nvalid >= 6) {
    const int rank = compute_coeffs(a_dx);
    if (rank >= 6) {
      const double a = coefficients_m(0), b = coefficients_m(1),
                   c = coefficients_m(2), d = coefficients_m(3),
                   e = coefficients_m(4), f = coefficients_m(5);
      const double theta = 0.5 * std::atan2(e, (safelyTiny(d - f)));
      const double cos_t = std::cos(theta), sin_t = std::sin(theta);
      const double A =
          -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
      const double B =
          -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);
      const double denom_inv = 1.0 / safelyTiny(4.0 * d * f - e * e);
      const double u = (2.0 * b * f - c * e) * denom_inv;
      const double v = -(b * e - 2.0 * d * c) * denom_inv;
      const double w = -(a + (-b * b * f + b * c * e - c * c * d) * denom_inv);

      const std::array<Normal, 3> evecs(
          {Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1)});
      const ReferenceFrame frame =
          ReferenceFrame(evecs[(dir + 1) % 3],
                         static_cast<double>(sign) * evecs[(dir + 2) % 3],
                         static_cast<double>(sign) * evecs[dir]);
      const Pt paraboloid_datum =
          a_stencil_center - u * frame[0] - v * frame[1] - w * frame[2];
      const UnitQuaternion rotation(theta, frame[2]);
      const ReferenceFrame paraboloid_frame = rotation * frame;
      if (nvalid == Nhx_m * Nhy_m) {
        status_m = HeightFunctionStatus::RETURNED_PARABOLOID_FROM_FULL_STENCIL;
      } else {
        status_m =
            HeightFunctionStatus::RETURNED_PARABOLOID_FROM_PARTIAL_STENCIL;
      }
      return Paraboloid(paraboloid_datum, paraboloid_frame, A, B);
    }
  }
  status_m = HeightFunctionStatus::RETURNED_PLANE;
  return Plane(normal, normal * a_stencil_center);
}

template <int Nx, int Ny, int Nz, class SeparatorVariantType>
inline int HeightFunctions_3D<Nx, Ny, Nz, SeparatorVariantType>::compute_coeffs(
    const double dx) {
  int nvalid = 0;
  for (int i = 0; i < Nhx_m; ++i) {
    for (int j = 0; j < Nhy_m; ++j) {
      if (valid_m(i, j)) {
        nvalid++;
      }
    }
  }
  Eigen::MatrixXd M(nvalid, 6);
  Eigen::VectorXd H(nvalid);
  int count = 0;
  for (int i = 0; i < Nhx_m; ++i) {
    for (int j = 0; j < Nhy_m; ++j) {
      if (valid_m(i, j)) {
        const double x0 =
            (static_cast<double>(i) - static_cast<double>(Nhx_m) / 2.0) * dx;
        const double y0 =
            (static_cast<double>(j) - static_cast<double>(Nhy_m) / 2.0) * dx;
        const double x1 = x0 + dx;
        const double y1 = y0 + dx;
        const double w = 1.0;  // / (1.0 + std::fabs(1 - i) + std::fabs(1 - j));
        H(count) = w * dx * heights_m(i, j);
        M(count, 0) = w * 1.0;                                         // 1
        M(count, 1) = w * (x1 * x1 - x0 * x0) / (2.0 * dx);            // x
        M(count, 2) = w * (y1 * y1 - y0 * y0) / (2.0 * dx);            // y
        M(count, 3) = w * (x1 * x1 * x1 - x0 * x0 * x0) / (3.0 * dx);  // x^2
        M(count, 4) = w * (x1 * x1 - x0 * x0) * (y1 * y1 - y0 * y0) /
                      (4.0 * dx * dx);                                 // xy
        M(count, 5) = w * (y1 * y1 * y1 - y0 * y0 * y0) / (3.0 * dx);  // y^2
        count++;
      }
    }
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  coefficients_m = svd.solve(H);
  return svd.rank();
}

template <int Nx, int Ny, int Nz, class SeparatorVariantType>
int HeightFunctions_3D<Nx, Ny, Nz, SeparatorVariantType>::compute_heights(
    const std::array<std::array<std::array<double, Nz>, Ny>, Nx>& a_vfrac,
    const int a_dir, const int a_sign) {
  if (a_dir == 0) {
    Nhx_m = (Ny % 2 == 1) ? 3 : 4;
    Nhy_m = (Nz % 2 == 1) ? 3 : 4;
  } else if (a_dir == 1) {
    Nhx_m = (Nz % 2 == 1) ? 3 : 4;
    Nhy_m = (Nx % 2 == 1) ? 3 : 4;
  } else {
    Nhx_m = (Nx % 2 == 1) ? 3 : 4;
    Nhy_m = (Ny % 2 == 1) ? 3 : 4;
  }
  heights_m.resize(Nhx_m, Nhy_m);
  valid_m.resize(Nhx_m, Nhy_m);

  const std::array<int, 3> half_columns({Nx / 2, Ny / 2, Nz / 2});
  const int t1 = (a_dir + 1) % 3;
  const int t2 = (a_dir + 2) % 3;
  const int s = a_sign >= 0 ? 1 : -1;

  int nvalid = 0;
  for (int a = 0; a < Nhx_m; ++a) {
    for (int b = 0; b < Nhy_m; ++b) {
      heights_m(a, b) = 0.0;
      valid_m(a, b) = false;

      auto column = [&](const int m) {
        int idx[3];
        idx[t1] = half_columns[t1] + a - Nhx_m / 2;
        idx[t2] = half_columns[t2] + (s > 0 ? b : Nhy_m - 1 - b) - Nhy_m / 2;
        idx[a_dir] = half_columns[a_dir] + s * m;
        return a_vfrac[idx[0]][idx[1]][idx[2]];
      };
      auto is_full = [&](const double c) { return c > 1.0 - tol_m; };
      auto is_empty = [&](const double c) { return c < tol_m; };

      int top = 0;
      bool found_top = false, consistent = true;
      bool seen_interfacial = !is_full(column(0));
      for (int m = 0; m <= half_columns[a_dir]; ++m) {
        const double c = column(m);
        if (is_empty(c)) {
          top = m;
          found_top = true;
          break;
        } else if (is_full(c)) {
          if (seen_interfacial) {
            consistent = false;
            break;
          }
        } else {
          seen_interfacial = true;
        }
      }

      int bottom = 0;
      bool found_bottom = false;
      seen_interfacial = !is_empty(column(0));
      for (int m = 0; consistent && m >= -half_columns[a_dir]; --m) {
        const double c = column(m);
        if (is_full(c)) {
          bottom = m;
          found_bottom = true;
          break;
        } else if (is_empty(c)) {
          if (seen_interfacial) {
            consistent = false;
            break;
          }
        } else {
          seen_interfacial = true;
        }
      }

      if (!consistent || !found_top || !found_bottom) {
        continue;
      }

      double sum = 0.0;
      for (int m = bottom; m <= top; ++m) {
        sum += column(m);
      }
      heights_m(a, b) = static_cast<double>(bottom) - 0.5 + sum;
      valid_m(a, b) = true;
      ++nvalid;
    }
  }
  return nvalid;
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_HEIGHT_FUNCTIONS_TPP_
