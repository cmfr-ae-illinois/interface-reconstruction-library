// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include <mpi.h>

#include "examples/paraboloid_advector/reconstruction_types.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/interface_reconstruction_methods/constrained_optimization_behavior.h"
#include "irl/interface_reconstruction_methods/elvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/plvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"
#include "irl/optimization/constrained_levenberg_marquardt.h"
#include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
#include "irl/paraboloid_reconstruction/hessian_paraboloid.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/planar_separator.h"

#include <Eigen/Dense>
#include <Eigen/QR>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/vof_advection.h"

void getReconstruction(
    const std::string& a_reconstruction_method,
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::GeneralMoments3D<2>>& a_gas_moments,
    Data<IRL::LocalizedParaboloidLink<double>>* a_localized_paraboloid_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::Paraboloid>* a_interface) {
  if (a_reconstruction_method == "Jibben") {
    Jibben::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_W, a_interface, a_localized_paraboloid_link);
  } else if (a_reconstruction_method == "CentroidFit") {
    Centroid::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                                a_W, a_interface, a_localized_paraboloid_link);
  } else if (a_reconstruction_method == "PLIC") {
    PLIC::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_W, a_interface, a_localized_paraboloid_link);
  } else if (a_reconstruction_method == "PLICSuperMOF1") {
    PLICSuperMOF1::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U,
                                     a_V, a_W, a_interface,
                                     a_localized_paraboloid_link);
  } else if (a_reconstruction_method == "PLICMOF1") {
    PLICMOF1::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                                a_W, a_interface, a_localized_paraboloid_link);
  } else if (a_reconstruction_method == "PPICMOF2") {
    PPICMOF2::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                                a_W, a_interface, a_localized_paraboloid_link);
  } else if (a_reconstruction_method == "PPICSuperMOF2") {
    PPICSuperMOF2::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U,
                                     a_V, a_W, a_interface,
                                     a_localized_paraboloid_link);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: PLIC, CentroidFit, Jibben. \n";
    std::exit(-1);
  }
}

// Wendland radial basis function
// Wendland, H. (1995). Piecewise polynomial, positive definite and
// compactly supported radial functions of minimal degree. Advances in
// Computational Mathematics, 4(1), 389â€“396.
double wgauss(const double d, const double h) {
  if (d >= h) {
    return 0.0;
  } else {
    return (1.0 + 4.0 * d / h) * std::pow(1.0 - d / h, 4.0);
  }
}

void updateReconstructionELVIRA(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    Data<IRL::PlanarSeparator>* a_liquid_gas_interface) {
  IRL::ELVIRANeighborhood neighborhood;
  const BasicMesh& mesh = a_liquid_moments.getMesh();
  neighborhood.resize(27);
  IRL::RectangularCuboid cells[27];
  std::array<double, 27> cells_vfrac;
  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
            liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, liquid_volume_fraction - 0.5);
          (*a_liquid_gas_interface)(i, j, k) =
              IRL::PlanarSeparator::fromOnePlane(
                  IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }
        // Build surrounding stencil information for ELVIRA.
        for (int kk = k - 1; kk < k + 2; ++kk) {
          for (int jj = j - 1; jj < j + 2; ++jj) {
            for (int ii = i - 1; ii < i + 2; ++ii) {
              // Reversed order, bad for cache locality but thats okay..
              cells[(kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1)] =
                  IRL::RectangularCuboid::fromBoundingPts(
                      IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                      IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
              cells_vfrac[(kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1)] =
                  a_liquid_moments(ii, jj, kk)[0] / mesh.cell_volume();
              neighborhood.setMember(
                  &cells[(kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1)],
                  &cells_vfrac[(kk - k + 1) * 9 + (jj - j + 1) * 3 +
                               (ii - i + 1)],
                  ii - i, jj - j, kk - k);
            }
          }
        }
        // Now perform actual ELVIRA and obtain interface PlanarSeparator
        (*a_liquid_gas_interface)(i, j, k) =
            reconstructionWithELVIRA3D(neighborhood);
      }
    }
  }
  // Update border with simple ghost-cell fill and correct distances for
  // assumed periodic boundary
  a_liquid_gas_interface->updateBorder();
  // correctInterfacePlaneBorders(a_liquid_gas_interface);
}

// Reconstruction with LVIRA - use input PlanarSeparator as initial guess
void updateReconstructionLVIRA(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments, const int a_nneigh,
    Data<IRL::PlanarSeparator>* a_liquid_gas_interface) {
  const BasicMesh& mesh = a_liquid_moments.getMesh();

  IRL::LVIRANeighborhood<IRL::RectangularCuboid> neighborhood;
  std::vector<IRL::RectangularCuboid> cells;
  std::vector<double> cells_vfrac;
  // std::vector<double> weights; // maybe later

  const int grid_size =
      (a_nneigh * 2 + 1) * (a_nneigh * 2 + 1) * (a_nneigh * 2 + 1);
  neighborhood.resize(grid_size);
  cells.resize(grid_size);
  cells_vfrac.resize(grid_size);

  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
            liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, liquid_volume_fraction - 0.5);
          (*a_liquid_gas_interface)(i, j, k) =
              IRL::PlanarSeparator::fromOnePlane(
                  IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }

        // Build surrounding stencil information for LVIRA.
        IRL::UnsignedIndex_t ndata = 0;
        for (int kk = k - a_nneigh; kk <= k + a_nneigh; ++kk) {
          for (int jj = j - a_nneigh; jj <= j + a_nneigh; ++jj) {
            for (int ii = i - a_nneigh; ii <= i + a_nneigh; ++ii) {
              // Trap center cell
              if (ii == i && jj == j && kk == k) {
                neighborhood.setCenterOfStencil(ndata);
              }
              cells[ndata] = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                  IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
              cells_vfrac[ndata] =
                  a_liquid_moments(ii, jj, kk)[0] / mesh.cell_volume();

              neighborhood.setMember(ndata, &cells[ndata], &cells_vfrac[ndata]);
              // Increment counter
              ++ndata;
            }
          }
        }
        auto found_planar_separator = (*a_liquid_gas_interface)(i, j, k);
        // Now perform actual LVIRA and obtain interface PlanarSeparator
        (*a_liquid_gas_interface)(i, j, k) =
            reconstructionWithLVIRA3D(neighborhood, found_planar_separator);
      }
    }
  }
}

void updatePolygon(const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
                   const Data<IRL::PlanarSeparator>& a_liquid_gas_interface,
                   Data<IRL::Polygon>* a_interface_polygon) {
  const BasicMesh& mesh = a_liquid_moments.getMesh();
  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
            liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        (*a_interface_polygon)(i, j, k) =
            IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                cell, a_liquid_gas_interface(i, j, k),
                a_liquid_gas_interface(i, j, k)[0]);
      }
    }
  }
}

std::array<double, 6> fitParaboloidToPLICHeights(
    const Data<IRL::Polygon>& a_polygon,
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const IRL::Pt& a_reference_point, const IRL::ReferenceFrame& a_frame,
    const int a_i, const int a_j, const int a_k, const int a_nneigh,
    const double a_width) {
  const BasicMesh& mesh = a_polygon.getMesh();
  const double meshsize = 1.0;  // mesh.dx();
  const int ic(a_i), jc(a_j), kc(a_k);
  IRL::UnsignedIndex_t m = 0;
  for (int k = kc - a_nneigh; k <= kc + a_nneigh; ++k) {
    for (int j = jc - a_nneigh; j <= jc + a_nneigh; ++j) {
      for (int i = ic - a_nneigh; i <= ic + a_nneigh; ++i) {
        if (a_polygon(i, j, k).getNumberOfVertices() > 0) m++;
      }
    }
  }

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, 6);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(m);
  const IRL::Pt pref = a_reference_point;
  const auto frame = a_frame;

  m = 0;
  for (int k = kc - a_nneigh; k <= kc + a_nneigh; ++k) {
    for (int j = jc - a_nneigh; j <= jc + a_nneigh; ++j) {
      for (int i = ic - a_nneigh; i <= ic + a_nneigh; ++i) {
        const IRL::UnsignedIndex_t shape =
            a_polygon(i, j, k).getNumberOfVertices();
        if (shape == 0) {
          continue;
        }
        // Local polygon normal and centroid
        IRL::Pt ploc = a_polygon(i, j, k).calculateCentroid();
        IRL::Normal nloc = a_polygon(i, j, k).calculateNormal();
        if (frame[2] * nloc <= 0.0) {
          continue;
        }
        ploc -= pref;
        ploc /= meshsize;
        const IRL::Pt tmp_pt = ploc;
        const IRL::Normal tmp_n = nloc;
        for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
          ploc[d] = frame[d] * tmp_pt;
          nloc[d] = frame[d] * tmp_n;
        }
        // Plane coefficients
        Eigen::VectorXd reconstruction_plane_coeffs(3);
        reconstruction_plane_coeffs << -(ploc * nloc) / meshsize, nloc[0],
            nloc[1];
        reconstruction_plane_coeffs /= -nloc[2];
        // Get weighting
        const double gaussianweight =  // 1.0;
            a_width <= 0.0
                ? 1.0
                : wgauss(std::sqrt(static_cast<IRL::Vec3<double>>(ploc) *
                                   static_cast<IRL::Vec3<double>>(ploc)),
                         a_width);
        const double vfrac = a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        double vfrac_weight = 1.0;
        const double limit_vfrac = 0.1;
        if (vfrac < limit_vfrac) {
          vfrac_weight = 0.5 - 0.5 * std::cos(M_PI * vfrac / limit_vfrac);
        } else if (vfrac > 1.0 - limit_vfrac) {
          vfrac_weight =
              0.5 - 0.5 * std::cos(M_PI * (1.0 - vfrac) / limit_vfrac);
        }
        double ww = 1.0;
        ww *= gaussianweight;
        ww *= vfrac_weight;
        // Integrals
        Eigen::VectorXd integrals = Eigen::VectorXd::Zero(6);
        double b_dot_sum = 0.0;
        for (IRL::UnsignedIndex_t v = 0; v < shape; ++v) {
          IRL::UnsignedIndex_t vn = (v + 1) % shape;
          IRL::Pt vert1 = a_polygon(i, j, k)[v];
          IRL::Pt vert2 = a_polygon(i, j, k)[vn];
          vert1 -= pref;
          vert2 -= pref;
          vert1 /= meshsize;
          vert2 /= meshsize;
          IRL::Pt tmp_pt1 = vert1;
          IRL::Pt tmp_pt2 = vert2;
          for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
            vert1[d] = frame[d] * tmp_pt1;
            vert2[d] = frame[d] * tmp_pt2;
          }

          const double xv = vert1[0];
          const double yv = vert1[1];
          const double xvn = vert2[0];
          const double yvn = vert2[1];

          Eigen::VectorXd integral_to_add(6);
          integral_to_add << (xv * yvn - xvn * yv) / 2.0,
              (xv + xvn) * (xv * yvn - xvn * yv) / 6.0,
              (yv + yvn) * (xv * yvn - xvn * yv) / 6.0,
              (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0,
              (yvn - yv) *
                  (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
                   2.0 * xv * xvn * yvn + xvn * xvn * yv +
                   3.0 * xvn * xvn * yvn) /
                  24.0,
              (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;
          integrals += integral_to_add;
        }
        if (ww > 0.0) {
          A.row(m) << ww * integrals(0), ww * integrals(1), ww * integrals(2),
              ww * integrals(3), ww * integrals(4), ww * integrals(5);
          b(m) = ww * integrals.head(3).dot(reconstruction_plane_coeffs);
        }
        m++;
      }
    }
  }

  // Unconstrained LS solution
  Eigen::VectorXd sol = A.completeOrthogonalDecomposition().pseudoInverse() * b;
  return std::array<double, 6>{
      {sol(0), sol(1), sol(2), sol(3), sol(4), sol(5)}};
}

std::array<double, 6> fitParaboloidToPLICHeightsOld(
    const Data<IRL::Polygon>& a_polygon,
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const IRL::Pt& a_reference_point, const IRL::ReferenceFrame& a_frame,
    const int a_i, const int a_j, const int a_k, const int a_nneigh,
    const double a_width) {
  const BasicMesh& mesh = a_polygon.getMesh();
  const double meshsize = 1.0;  // mesh.dx();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6, 6);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(6);
  const int ic(a_i), jc(a_j), kc(a_k);
  const IRL::Pt pref = a_reference_point;
  const auto frame = a_frame;

  for (int k = kc - a_nneigh; k <= kc + a_nneigh; ++k) {
    for (int j = jc - a_nneigh; j <= jc + a_nneigh; ++j) {
      for (int i = ic - a_nneigh; i <= ic + a_nneigh; ++i) {
        const IRL::UnsignedIndex_t shape =
            a_polygon(i, j, k).getNumberOfVertices();
        if (shape == 0) {
          continue;
        }
        // Local polygon normal and centroid
        IRL::Pt ploc = a_polygon(i, j, k).calculateCentroid();
        IRL::Normal nloc = a_polygon(i, j, k).calculateNormal();
        // if (frame[2] * nloc <= 0.0) {
        //   continue;
        // }
        ploc -= pref;
        ploc /= meshsize;
        const IRL::Pt tmp_pt = ploc;
        const IRL::Normal tmp_n = nloc;
        for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
          ploc[d] = frame[d] * tmp_pt;
          nloc[d] = frame[d] * tmp_n;
        }
        // Plane coefficients
        Eigen::VectorXd reconstruction_plane_coeffs(3);
        reconstruction_plane_coeffs << -(ploc * nloc) / meshsize, nloc[0],
            nloc[1];
        reconstruction_plane_coeffs /= -nloc[2];
        // Integrals
        Eigen::VectorXd integrals = Eigen::VectorXd::Zero(6);
        double b_dot_sum = 0.0;
        for (IRL::UnsignedIndex_t v = 0; v < shape; ++v) {
          IRL::UnsignedIndex_t vn = (v + 1) % shape;
          IRL::Pt vert1 = a_polygon(i, j, k)[v];
          IRL::Pt vert2 = a_polygon(i, j, k)[vn];
          vert1 -= pref;
          vert2 -= pref;
          vert1 /= meshsize;
          vert2 /= meshsize;
          IRL::Pt tmp_pt1 = vert1;
          IRL::Pt tmp_pt2 = vert2;
          for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
            vert1[d] = frame[d] * tmp_pt1;
            vert2[d] = frame[d] * tmp_pt2;
          }

          const double xv = vert1[0];
          const double yv = vert1[1];
          const double xvn = vert2[0];
          const double yvn = vert2[1];

          Eigen::VectorXd integral_to_add(6);
          integral_to_add << (xv * yvn - xvn * yv) / 2.0,
              (xv + xvn) * (xv * yvn - xvn * yv) / 6.0,
              (yv + yvn) * (xv * yvn - xvn * yv) / 6.0,
              (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0,
              (yvn - yv) *
                  (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
                   2.0 * xv * xvn * yvn + xvn * xvn * yv +
                   3.0 * xvn * xvn * yvn) /
                  24.0,
              (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;
          integrals += integral_to_add;
        }
        b_dot_sum += integrals.head(3).dot(reconstruction_plane_coeffs);

        // Get weighting
        const double gaussianweight =  // 1.0;
            a_width <= 0.0
                ? 1.0
                : wgauss(std::sqrt(static_cast<IRL::Vec3<double>>(ploc) *
                                   static_cast<IRL::Vec3<double>>(ploc)),
                         a_width);
        const double vfrac = a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        double vfrac_weight = 1.0;
        const double limit_vfrac = 0.1;
        if (vfrac < limit_vfrac) {
          vfrac_weight = 0.5 - 0.5 * std::cos(M_PI * vfrac / limit_vfrac);
        } else if (vfrac > 1.0 - limit_vfrac) {
          vfrac_weight =
              0.5 - 0.5 * std::cos(M_PI * (1.0 - vfrac) / limit_vfrac);
        }
        double ww = 1.0;
        ww *= gaussianweight;
        ww *= vfrac_weight;

        if (ww > 0.0) {
          A += ww * integrals * integrals.transpose();
          b += ww * integrals * b_dot_sum;
        }
      }
    }
  }
  Eigen::VectorXd sol = A.colPivHouseholderQr().solve(b);
  return std::array<double, 6>{
      {sol(0), sol(1), sol(2), sol(3), sol(4), sol(5)}};
}

std::array<double, 6> fitParaboloidToCentroids(
    const Data<IRL::Polygon>& a_polygon,
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const IRL::Pt& a_reference_point, const IRL::ReferenceFrame& a_frame,
    const int a_i, const int a_j, const int a_k, const int a_nneigh,
    const double a_width) {
  const BasicMesh& mesh = a_polygon.getMesh();
  const double meshsize = 1.0;  // mesh.dx();
  // const int ncells = std::pow(2 * a_nneigh + 1, 3);
  const int ncells =
      (2 * a_nneigh + 1) * (2 * a_nneigh + 1) * (2 * a_nneigh + 1);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(ncells, 6);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(ncells);
  const int ic(a_i), jc(a_j), kc(a_k);
  const IRL::Pt pref = a_reference_point;
  const auto frame = a_frame;

  IRL::UnsignedIndex_t ndata = 0;
  for (int k = kc - a_nneigh; k <= kc + a_nneigh; ++k) {
    for (int j = jc - a_nneigh; j <= jc + a_nneigh; ++j) {
      for (int i = ic - a_nneigh; i <= ic + a_nneigh; ++i) {
        if (a_polygon(i, j, k).getNumberOfVertices() == 0) {
          continue;
        }
        IRL::Pt ploc = a_polygon(i, j, k).calculateCentroid();
        const IRL::Normal nloc = a_polygon(i, j, k).calculateNormal();
        const double surf = a_polygon(i, j, k).calculateAbsoluteVolume() /
                            (meshsize * meshsize);
        const double normalproj = std::max(frame[2] * nloc, 0.0);
        // if (normalproj <= 0.0) continue;
        ploc -= pref;
        ploc /= meshsize;
        const IRL::Pt tmp_pt = ploc;
        for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
          ploc[d] = frame[d] * tmp_pt;
        }
        const double gaussianweight =
            // 1.0;
            a_width <= 0.0
                ? 1.0
                : wgauss(std::sqrt(static_cast<IRL::Vec3<double>>(ploc) *
                                   static_cast<IRL::Vec3<double>>(ploc)),
                         a_width);

        const double vfrac = a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        double vfrac_weight = 1.0;
        const double limit_vfrac = 0.1;
        if (vfrac < limit_vfrac) {
          vfrac_weight = 0.5 - 0.5 * std::cos(M_PI * vfrac / limit_vfrac);
        } else if (vfrac > 1.0 - limit_vfrac) {
          vfrac_weight =
              0.5 - 0.5 * std::cos(M_PI * (1.0 - vfrac) / limit_vfrac);
        }
        double ww = 1.0;
        ww *= normalproj;
        ww *= surf;
        ww *= gaussianweight;
        ww *= vfrac_weight;

        if (ww > 0.0) {
          // Store least squares matrix and RHS
          A(ndata, 0) = std::sqrt(ww);
          A(ndata, 1) = 0.0;
          A(ndata, 2) = 0.0;
          A(ndata, 3) = std::sqrt(ww) * ploc[0] * ploc[0];
          A(ndata, 4) = std::sqrt(ww) * ploc[0] * ploc[1];
          A(ndata, 5) = std::sqrt(ww) * ploc[1] * ploc[1];
          b(ndata) = std::sqrt(ww) * ploc[2];
          // Increment counter
          ++ndata;
        }
      }
    }
  }
  A.conservativeResize(ndata, Eigen::NoChange);
  b.conservativeResize(ndata, Eigen::NoChange);
  Eigen::VectorXd sol = A.colPivHouseholderQr().solve(b);
  return std::array<double, 6>{
      // {sol(0), sol(1), sol(2), sol(3), sol(4), sol(5)}};
      {sol(0), 0.0, 0.0, sol(3), sol(4), sol(5)}};
}

void PLICMOF1::getReconstruction(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::GeneralMoments3D<2>>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::Paraboloid>* a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  const BasicMesh& mesh = a_liquid_moments.getMesh();

  const double cell_volume = mesh.cell_volume();
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / cell_volume;
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
        } else {
          auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
              IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          const IRL::Pt cell_centroid =
              0.5 * (IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)) +
                     IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          IRL::Pt liquid_centroid = IRL::Pt(a_liquid_moments(i, j, k)[1],
                                            a_liquid_moments(i, j, k)[2],
                                            a_liquid_moments(i, j, k)[3]);
          IRL::Pt gas_centroid = cell_volume * cell_centroid - liquid_centroid;
          liquid_centroid *= 1.0 / a_liquid_moments(i, j, k)[0];
          gas_centroid *= 1.0 / (cell_volume - a_liquid_moments(i, j, k)[0]);
          IRL::SeparatedMoments<IRL::VolumeMoments> svm(
              IRL::VolumeMoments(a_liquid_moments(i, j, k)[0], liquid_centroid),
              IRL::VolumeMoments(cell_volume - a_liquid_moments(i, j, k)[0],
                                 gas_centroid));
          auto planar_separator =
              IRL::reconstructionWithMOF3D(cell, svm, 0.5, 0.5);

          const IRL::Normal normal = planar_separator[0].normal();
          IRL::ReferenceFrame fit_frame;
          int largest_dir = 0;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
            largest_dir = 1;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
            largest_dir = 2;
          if (largest_dir == 0)
            fit_frame[0] = crossProduct(normal, IRL::Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 1)
            fit_frame[0] = crossProduct(normal, IRL::Normal(0.0, 0.0, 1.0));
          else
            fit_frame[0] = crossProduct(normal, IRL::Normal(1.0, 0.0, 0.0));
          fit_frame[0].normalize();
          fit_frame[1] = crossProduct(normal, fit_frame[0]);
          fit_frame[2] = normal;
          IRL::Paraboloid paraboloid =
              IRL::Paraboloid(cell_centroid, fit_frame, 1.0e-3, -1.0e-3);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                              paraboloid);

          if (solver_distance.getDistance() == -DBL_MAX) {
            (*a_interface)(i, j, k) =
                IRL::Paraboloid(cell_centroid, fit_frame, 1.0e-3, -1.0e-3);
          } else {
            auto new_datum =
                IRL::Pt(paraboloid.getDatum() +
                        solver_distance.getDistance() * fit_frame[2]);
            paraboloid.setDatum(new_datum);
            (*a_interface)(i, j, k) = paraboloid;
          }
        }
      }
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void RecenterMoments(IRL::GeneralMoments3D<2>* moments, const IRL::Pt& center) {
  const double M0 = (*moments)[0];
  const Eigen::Matrix<double, 3, 1> M1{(*moments)[1], (*moments)[2],
                                       (*moments)[3]};
  const Eigen::Matrix<double, 3, 3> M2{
      {(*moments)[4], (*moments)[5], (*moments)[6]},
      {(*moments)[5], (*moments)[7], (*moments)[8]},
      {(*moments)[6], (*moments)[8], (*moments)[9]}};
  const Eigen::Matrix<double, 3, 1> C{center[0], center[1], center[2]};
  auto M1_final = M1 - M0 * C;
  auto M2_final =
      M2 - M1 * C.transpose() - C * M1.transpose() + M0 * C * C.transpose();
  (*moments)[1] = M1_final(0);
  (*moments)[2] = M1_final(1);
  (*moments)[3] = M1_final(2);
  (*moments)[4] = M2_final(0, 0);
  (*moments)[5] = M2_final(0, 1);
  (*moments)[6] = M2_final(0, 2);
  (*moments)[7] = M2_final(1, 1);
  (*moments)[8] = M2_final(1, 2);
  (*moments)[9] = M2_final(2, 2);

  // if (M2_final(0, 0) < 0.0) {
  //   std::cout << "WARNING: c M2xx = " << M2_final(0, 0) << "; M0 = " << M0
  //             << std::endl;
  // }
  // if (M2_final(1, 1) < 0.0) {
  //   std::cout << "WARNING: c M2yy = " << M2_final(1, 1) << "; M0 = " << M0
  //             << std::endl;
  // }
  // if (M2_final(2, 2) < 0.0) {
  //   std::cout << "WARNING: c M2zz = " << M2_final(2, 2) << "; M0 = " << M0
  //             << std::endl;
  // }
  // Recenter second moments
  // (*moments)[4] +=
  //     (*moments)[0] * center[0] * center[0] - 2.0 * (*moments)[1] *
  //     center[0];
  // (*moments)[5] += (*moments)[0] * center[0] * center[1] -
  //                  (*moments)[1] * center[1] - (*moments)[2] * center[0];
  // (*moments)[6] += (*moments)[0] * center[0] * center[2] -
  //                  (*moments)[1] * center[2] - (*moments)[3] * center[0];
  // (*moments)[7] +=
  //     (*moments)[0] * center[1] * center[1] - 2.0 * (*moments)[2] *
  //     center[1];
  // (*moments)[8] += (*moments)[0] * center[1] * center[2] -
  //                  (*moments)[2] * center[2] - (*moments)[3] * center[1];
  // (*moments)[9] +=
  //     (*moments)[0] * center[2] * center[2] - 2.0 * (*moments)[3] *
  //     center[2];

  // // Recenter first moments
  // (*moments)[1] -= (*moments)[0] * center[0];
  // (*moments)[2] -= (*moments)[0] * center[1];
  // (*moments)[3] -= (*moments)[0] * center[2];
}

struct MOF1Functor {
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  typedef double Scalar;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // Variables
  int m_inputs, m_values;
  double m_cell_volume;
  IRL::RectangularCuboid m_cell;
  IRL::VolumeMoments m_moments;
  IRL::ReferenceFrame m_frame;

  // Constructor
  MOF1Functor(int inputs, int values, const IRL::RectangularCuboid& cell,
              const IRL::VolumeMoments& moments)
      : m_inputs(inputs),
        m_values(values),
        m_cell(cell),
        m_cell_volume(cell.calculateVolume()),
        m_moments(moments) {}

  void setframe(const IRL::Normal& guess_direction) {
    IRL::Normal normal = guess_direction;
    normal.normalize();
    int largest_dir = 0;
    if (std::fabs(normal[largest_dir]) < std::fabs(normal[1])) largest_dir = 1;
    if (std::fabs(normal[largest_dir]) < std::fabs(normal[2])) largest_dir = 2;
    if (largest_dir == 0)
      m_frame[0] = crossProduct(normal, IRL::Normal(0.0, 1.0, 0.0));
    else if (largest_dir == 1)
      m_frame[0] = crossProduct(normal, IRL::Normal(0.0, 0.0, 1.0));
    else
      m_frame[0] = crossProduct(normal, IRL::Normal(1.0, 0.0, 0.0));
    m_frame[0].normalize();
    m_frame[1] = crossProduct(normal, m_frame[0]);
    m_frame[2] = normal;
  }

  const IRL::PlanarSeparator getseparator(const Eigen::VectorXd& x) const {
    const auto rotation = IRL::UnitQuaternion(x(1), m_frame[1]) *
                          IRL::UnitQuaternion(x(0), m_frame[0]);
    const auto newframe = rotation * m_frame;
    const double guess_distance = m_frame[2] * m_moments.centroid() /
                                  IRL::safelyEpsilon(m_moments.volume());
    IRL::PlanarSeparator separator = IRL::PlanarSeparator::fromOnePlane(
        IRL::Plane(newframe[2], guess_distance));
    IRL::setDistanceToMatchVolumeFractionPartialFill(
        m_cell, m_moments.volume() / IRL::safelyEpsilon(m_cell_volume),
        &separator);
    return separator;
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto moments = IRL::getVolumeMoments<IRL::VolumeMoments>(
        m_cell, this->getseparator(x));
    fvec(0) = moments.volume() - m_moments.volume();
    fvec(1) = moments.centroid()[0] - m_moments.centroid()[0];
    fvec(2) = moments.centroid()[1] - m_moments.centroid()[1];
    fvec(3) = moments.centroid()[2] - m_moments.centroid()[2];
    return 0;
  }

  // Get methods for users to determine function input and output dimensions
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

void PLICSuperMOF1::getReconstruction(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::GeneralMoments3D<2>>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::Paraboloid>* a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  const BasicMesh& mesh = a_liquid_moments.getMesh();

  const double cell_volume = mesh.cell_volume();
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / cell_volume;
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
        } else {
          auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
              IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          auto super_cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i - 1), mesh.y(j - 1), mesh.z(k - 1)),
              IRL::Pt(mesh.x(i + 2), mesh.y(j + 2), mesh.z(k + 2)));
          const IRL::Pt cell_centroid =
              0.5 * (IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)) +
                     IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          double liquid_volume = 0.0;
          IRL::Pt liquid_centroid = IRL::Pt(0, 0, 0);
          for (int ii = -1; ii <= 1; ++ii) {
            for (int jj = -1; jj <= 1; ++jj) {
              for (int kk = -1; kk <= 1; ++kk) {
                liquid_volume += a_liquid_moments(i + ii, j + jj, k + kk)[0];
                liquid_centroid +=
                    IRL::Pt(a_liquid_moments(i + ii, j + jj, k + kk)[1],
                            a_liquid_moments(i + ii, j + jj, k + kk)[2],
                            a_liquid_moments(i + ii, j + jj, k + kk)[3]);
              }
            }
          }
          IRL::Pt gas_centroid =
              27.0 * cell_volume * cell_centroid - liquid_centroid;
          liquid_centroid *= 1.0 / liquid_volume;
          gas_centroid *= 1.0 / (27.0 * cell_volume - liquid_volume);
          IRL::SeparatedMoments<IRL::VolumeMoments> svm(
              IRL::VolumeMoments(liquid_volume, liquid_centroid),
              IRL::VolumeMoments(27.0 * cell_volume - liquid_volume,
                                 gas_centroid));
          // auto planar_separator =
          //     IRL::reconstructionWithMOF3D(super_cell, svm, 0.5, 0.5);

          ///////////////////////// EIGEN TEST
          auto centroid_line = IRL::Normal(gas_centroid - liquid_centroid);
          const auto moments_for_fit = IRL::VolumeMoments(
              liquid_volume, liquid_volume * liquid_centroid);
          MOF1Functor myMOF1Functor(2, 4, super_cell, moments_for_fit);
          myMOF1Functor.setframe(centroid_line);
          Eigen::NumericalDiff<MOF1Functor> NDMOF1Functor(myMOF1Functor);
          Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF1Functor>, double>
              MOF1LevenbergMarquardt(NDMOF1Functor);
          MOF1LevenbergMarquardt.parameters.ftol = 1e-15;
          MOF1LevenbergMarquardt.parameters.xtol = 1e-15;
          MOF1LevenbergMarquardt.parameters.maxfev = 100;  // Max iterations
          Eigen::VectorXd x(2);
          x(0) = 0.0;
          x(1) = 0.0;
          MOF1LevenbergMarquardt.minimize(x);
          const auto planar_separator = myMOF1Functor.getseparator(x);
          // std::cout << "              Eigen separator: "
          //           << planar_separator[0].normal() << " "
          //           << planar_separator[0].distance() << std::endl;
          // std::cout << "                IRL separator: "
          //           << planar_separator[0].normal() << " "
          //           << planar_separator[0].distance() << std::endl;
          ///////////////////////// EIGEN TEST

          const IRL::Normal normal = planar_separator[0].normal();
          IRL::ReferenceFrame fit_frame;
          int largest_dir = 0;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
            largest_dir = 1;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
            largest_dir = 2;
          if (largest_dir == 0)
            fit_frame[0] = crossProduct(normal, IRL::Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 1)
            fit_frame[0] = crossProduct(normal, IRL::Normal(0.0, 0.0, 1.0));
          else
            fit_frame[0] = crossProduct(normal, IRL::Normal(1.0, 0.0, 0.0));
          fit_frame[0].normalize();
          fit_frame[1] = crossProduct(normal, fit_frame[0]);
          fit_frame[2] = normal;
          IRL::Paraboloid paraboloid =
              IRL::Paraboloid(cell_centroid, fit_frame, 1.0e-3, -1.0e-3);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                              paraboloid);

          if (solver_distance.getDistance() == -DBL_MAX) {
            (*a_interface)(i, j, k) =
                IRL::Paraboloid(cell_centroid, fit_frame, 1.0e-3, -1.0e-3);
          } else {
            auto new_datum =
                IRL::Pt(paraboloid.getDatum() +
                        solver_distance.getDistance() * fit_frame[2]);
            paraboloid.setDatum(new_datum);
            (*a_interface)(i, j, k) = paraboloid;
          }
        }
      }
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

struct MOF2Functor {
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  typedef double Scalar;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // Variables
  int m_inputs, m_values;
  double m_cell_volume, m_vfrac, m_a, m_b;
  IRL::RectangularCuboid m_cell, m_cell_constraint;
  IRL::GeneralMoments3D<2> m_liquid_moments, m_gas_moments;
  IRL::ReferenceFrame m_frame;
  IRL::Pt m_datum, m_liquid_centroid, m_gas_centroid, m_cell_centroid;
  double m_length_scale, m_m0_scale, m_m1_scale_liquid, m_m1_scale_gas,
      m_m2_scale_liquid, m_m2_scale_gas, m_vfrac_constraint;

  // Constructor
  MOF2Functor(int inputs, int values, const IRL::RectangularCuboid& cell,
              const IRL::GeneralMoments3D<2>& liquid_moments,
              const IRL::GeneralMoments3D<2>& gas_moments,
              const IRL::RectangularCuboid& constraint_cell,
              const double constraint_vfrac)
      : m_inputs(inputs),
        m_values(values),
        m_cell(cell),
        m_cell_volume(cell.calculateVolume()),
        m_cell_centroid(cell.calculateCentroid()),
        m_liquid_moments(liquid_moments),
        m_gas_moments(gas_moments),
        m_vfrac(liquid_moments[0] / cell.calculateVolume()),
        m_vfrac_constraint(constraint_vfrac),
        m_cell_constraint(constraint_cell) {
    m_length_scale = std::cbrt(cell.calculateVolume());
    m_liquid_centroid = IRL::Pt(
        IRL::Pt(liquid_moments[1], liquid_moments[2], liquid_moments[3]) /
        liquid_moments[0]);
    m_gas_centroid =
        IRL::Pt(IRL::Pt(gas_moments[1], gas_moments[2], gas_moments[3]) /
                gas_moments[0]);
    m_datum =
        IRL::Pt((1.0 - m_vfrac) * m_liquid_centroid + m_vfrac * m_gas_centroid);
    RecenterMoments(&m_liquid_moments, m_liquid_centroid);
    RecenterMoments(&m_gas_moments, m_gas_centroid);
  }

  void setframe(const IRL::Normal& guess_direction,
                const IRL::Paraboloid& guess_paraboloid) {
    // Set scales
    m_m0_scale = 1.0 / m_cell_volume;
    // m_m1_scale_liquid = std::pow(m_cell_volume, -4.0 / 3.0);
    // m_m1_scale_gas = std::pow(m_cell_volume, -4.0 / 3.0);
    // m_m2_scale_liquid = std::pow(m_cell_volume, -5.0 / 3.0);
    // m_m2_scale_gas = std::pow(m_cell_volume, -5.0 / 3.0);
    m_m1_scale_liquid = std::pow(m_liquid_moments[0], -4.0 / 3.0);
    m_m1_scale_gas = std::pow(m_gas_moments[0], -4.0 / 3.0);
    m_m2_scale_liquid = std::pow(m_liquid_moments[0], -5.0 / 3.0);
    m_m2_scale_gas = std::pow(m_gas_moments[0], -5.0 / 3.0);

    // Try paraboloid guess and compute error
    const auto guess_coeffs = guess_paraboloid.getAlignedParaboloid();
    const auto guess_frame = guess_paraboloid.getReferenceFrame();
    m_frame = guess_frame;
    m_a = guess_coeffs.a();
    m_b = guess_coeffs.b();
    const double k1dx = 2.0 * m_a * m_length_scale;
    const double k2dx = 2.0 * m_b * m_length_scale;
    // If error is too large, revert to planar initial guess
    if (std::abs(k1dx) > 6.0 || std::abs(k2dx) > 6.0) {
      IRL::Normal normal = guess_frame[2];  // guess_direction;
      normal.normalize();
      m_frame[0] = crossProduct(guess_frame[1], normal);
      if (IRL::magnitude(m_frame[0]) < 1.0e-9) {
        m_frame[0] = crossProduct(normal, guess_frame[0]);
      }
      m_frame[0].normalize();
      m_frame[1] = crossProduct(normal, m_frame[0]);
      m_frame[2] = normal;
      m_a = 0.0;  // 1.0e-3;
      m_b = 0.0;  //-1.0e-3;
    }
  }

  const IRL::Paraboloid getparaboloid(const Eigen::VectorXd& x) const {
    const IRL::Pt datum = m_datum;
    IRL::ReferenceFrame newframe;
    IRL::Paraboloid paraboloid;
    const double vfrac_tol = 1.0e-2;
    if (m_vfrac < vfrac_tol || m_vfrac > 1.0 - vfrac_tol) {
      const auto rotation = IRL::UnitQuaternion(2.0 * M_PI * x(1), m_frame[1]) *
                            IRL::UnitQuaternion(2.0 * M_PI * x(0), m_frame[0]);
      newframe = rotation * m_frame;
      double delta = (x(2) - 0.999) / m_length_scale;
      if (std::abs(delta * m_length_scale) > 3.0)
        delta = std::copysign(3.0 / m_length_scale, delta);
      paraboloid = IRL::Paraboloid(datum, newframe, delta, delta);
    } else {
      const auto rotation = IRL::UnitQuaternion(2.0 * M_PI * x(4), m_frame[2]) *
                            IRL::UnitQuaternion(2.0 * M_PI * x(1), m_frame[1]) *
                            IRL::UnitQuaternion(2.0 * M_PI * x(0), m_frame[0]);
      newframe = rotation * m_frame;
      const double delta_a = (x(2) - 0.999) / m_length_scale;
      const double delta_b = (x(3) - 1.001) / m_length_scale;
      double coeff_a = m_a + delta_a;
      double coeff_b = m_b + delta_b;
      if (std::abs(coeff_a * m_length_scale) > 3.0)
        coeff_a = std::copysign(3.0 / m_length_scale, coeff_a);
      if (std::abs(coeff_b * m_length_scale) > 3.0)
        coeff_b = std::copysign(3.0 / m_length_scale, coeff_b);
      paraboloid = IRL::Paraboloid(datum, newframe, coeff_a, coeff_b);
    }
    IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
        solver_distance(m_cell_constraint, m_vfrac_constraint, 1.0e-13,
                        paraboloid);
    auto new_datum =
        IRL::Pt(datum + solver_distance.getDistance() * newframe[2]);
    paraboloid.setDatum(new_datum);
    return paraboloid;
  }

  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto paraboloid = this->getparaboloid(x);
    const auto sep_moments =
        IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::GeneralMoments3D<2>>>(
            m_cell, paraboloid);
    // const auto exact_sep_moments =
    //     IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
    //         m_cell, paraboloid);
    auto liquid_moments = sep_moments[0];
    auto gas_moments = sep_moments[1];
    // liquid_moments[0] = exact_sep_moments[0].volume();
    // liquid_moments[1] = exact_sep_moments[0].centroid()[0];
    // liquid_moments[2] = exact_sep_moments[0].centroid()[1];
    // liquid_moments[3] = exact_sep_moments[0].centroid()[2];
    // gas_moments[0] = exact_sep_moments[1].volume();
    // gas_moments[1] = exact_sep_moments[1].centroid()[0];
    // gas_moments[2] = exact_sep_moments[1].centroid()[1];
    // gas_moments[3] = exact_sep_moments[1].centroid()[2];
    RecenterMoments(&liquid_moments, m_liquid_centroid);
    RecenterMoments(&gas_moments, m_gas_centroid);
    fvec.setZero();
    const double k1dx =
        2.0 * (m_a + (x(2) - 0.999) / m_length_scale) * m_length_scale;
    const double k2dx =
        2.0 * (m_b + (x(3) - 1.001) / m_length_scale) * m_length_scale;
    const double maxkdx = 6.0;
    const double mu = 50.0;

    // Penalty to prevent kappa * dx > 1
    fvec(0) = mu * ((std::max(0.0, std::abs(k1dx) - maxkdx)) +
                    (std::max(0.0, std::abs(k2dx) - maxkdx)) +
                    (std::max(0.0, std::abs(1.0 - x(0)) - 0.05)) +
                    (std::max(0.0, std::abs(1.0 - x(1)) - 0.05)));
    // fvec(1) = m_m0_scale * (liquid_moments[0] - m_liquid_moments[0]);
    for (int i = 0; i < 3; ++i) {
      fvec(2 + i) =
          m_m1_scale_liquid * (liquid_moments[1 + i] - m_liquid_moments[1 + i]);
      fvec(5 + i) =
          m_m1_scale_gas * (gas_moments[1 + i] - m_gas_moments[1 + i]);
    }
    for (int i = 0; i < 6; ++i) {
      fvec(8 + i) =
          m_m2_scale_liquid * (liquid_moments[4 + i] - m_liquid_moments[4 + i]);
      fvec(14 + i) =
          m_m2_scale_gas * (gas_moments[4 + i] - m_gas_moments[4 + i]);
    }
    fvec(9) *= 2.0;
    fvec(10) *= 2.0;
    fvec(12) *= 2.0;
    fvec(15) *= 2.0;
    fvec(16) *= 2.0;
    fvec(18) *= 2.0;
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    this->errorvec(x, fvec);
    return 0;
  }

  // Get methods for users to determine function input and output dimensions
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

void PPICMOF2::getReconstruction(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::GeneralMoments3D<2>>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::Paraboloid>* a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  const BasicMesh& mesh = a_liquid_moments.getMesh();

  Jibben::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_W, a_interface, a_link_localized_paraboloid);

  const double cell_volume = mesh.cell_volume();
  const double cell_dx = std::pow(cell_volume, 1.0 / 3.0);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / cell_volume;
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
        } else {
          auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
              IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          const IRL::Pt cell_centroid =
              0.5 * (IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)) +
                     IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          auto liquid_moments = a_liquid_moments(i, j, k);
          auto gas_moments = a_gas_moments(i, j, k);
          IRL::Pt liquid_centroid =
              IRL::Pt(liquid_moments[1], liquid_moments[2], liquid_moments[3]);
          IRL::Pt gas_centroid =
              IRL::Pt(gas_moments[1], gas_moments[2], gas_moments[3]);
          liquid_centroid *= 1.0 / liquid_moments[0];
          gas_centroid *= 1.0 / gas_moments[0];

          ///////////////////////// EIGEN TEST
          auto centroid_line = IRL::Normal(gas_centroid - liquid_centroid);

          MOF2Functor myMOF2Functor(4, 20, cell, liquid_moments, gas_moments,
                                    cell, liquid_volume_fraction);
          myMOF2Functor.setframe(centroid_line, (*a_interface)(i, j, k));
          Eigen::NumericalDiff<MOF2Functor> NDMOF2Functor(myMOF2Functor);
          Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF2Functor>, double>
              MOF2LM(NDMOF2Functor);
          // MOF2LM.parameters.ftol = 1.0e-2;
          // MOF2LM.parameters.xtol = 1.0e-2;
          // MOF2LM.parameters.factor = 1.0;
          // MOF2LM.parameters.maxfev = 100;  // Max
          // iterations
          Eigen::VectorXd x(4);
          x.setZero();
          // MOF2LM.minimize(x);
          Eigen::LevenbergMarquardtSpace::Status status =
              MOF2LM.minimizeInit(x);
          do {
            status = MOF2LM.minimizeOneStep(x);
            // if (std::abs(x(2)) > 2.0) x(2) = std::copysign(2.0, x(2));
            // if (std::abs(x(3)) > 2.0) x(3) = std::copysign(2.0, x(3));
          } while (status == Eigen::LevenbergMarquardtSpace::Running);
          IRL::Paraboloid paraboloid = myMOF2Functor.getparaboloid(x);
          auto coeffs = paraboloid.getAlignedParaboloid();
          IRL::AlignedParaboloid new_coeffs = coeffs;
          if (std::abs(2.0 * new_coeffs.a() * cell_dx) > 2.0)
            new_coeffs.a() = std::copysign(1.0 / cell_dx, new_coeffs.a());
          if (std::abs(2.0 * new_coeffs.b() * cell_dx) > 2.0)
            new_coeffs.b() = std::copysign(1.0 / cell_dx, new_coeffs.b());
          paraboloid.setAlignedParaboloid(new_coeffs);
          const auto fit_frame = paraboloid.getReferenceFrame();
          ///////////////////////// EIGEN TEST
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                              paraboloid);
          auto new_datum =
              IRL::Pt(paraboloid.getDatum() +
                      solver_distance.getDistance() * fit_frame[2]);
          paraboloid.setDatum(new_datum);
          (*a_interface)(i, j, k) = paraboloid;
        }
      }
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void PPICSuperMOF2::getReconstruction(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::GeneralMoments3D<2>>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::Paraboloid>* a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  const BasicMesh& mesh = a_liquid_moments.getMesh();
  const bool second_pass = false;

  Jibben::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_W, a_interface, a_link_localized_paraboloid);

  const double cell_volume = mesh.cell_volume();
  const double cell_dx = mesh.dx();

  int nmixed_global = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / cell_volume;
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
        } else {
          nmixed_global++;
        }
      }
    }
  }

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  IRL::Paraboloid dummy_par;
  IRL::ByteBuffer dummy_buffer;
  dummy_buffer.resize(0);
  dummy_buffer.resetBufferPointer();
  IRL::serializeAndPack(dummy_par, &dummy_buffer);
  const int size_paraboloid = dummy_buffer.size();

  int nmixed_local = std::max(nmixed_global / size, 1);
  std::vector<int> proc_offset(size + 1);
  proc_offset[0] = 0;
  for (int r = 0; r < size; r++)
    proc_offset[r + 1] = proc_offset[r] + nmixed_local;
  proc_offset[size] = nmixed_global;
  for (int r = 1; r < size + 1; r++)
    proc_offset[r] = std::min(proc_offset[r], nmixed_global);
  nmixed_local = proc_offset[rank + 1] - proc_offset[rank];
  IRL::ByteBuffer interface_local, interface_global;
  interface_local.resize(nmixed_local * sizeof(IRL::Paraboloid));
  interface_global.resize(0);
  interface_local.resetBufferPointer();
  interface_global.resetBufferPointer();

  int count = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / cell_volume;
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          if (count >= proc_offset[rank] && count < proc_offset[rank + 1]) {
            const auto cell = IRL::RectangularCuboid::fromBoundingPts(
                IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
                IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
            const auto super_cell = IRL::RectangularCuboid::fromBoundingPts(
                IRL::Pt(mesh.x(i - 1), mesh.y(j - 1), mesh.z(k - 1)),
                IRL::Pt(mesh.x(i + 2), mesh.y(j + 2), mesh.z(k + 2)));
            const auto cell_centroid =
                IRL::Pt(mesh.xm(i), mesh.ym(j), mesh.zm(k));
            auto super_liquid_moments =
                IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);
            auto super_gas_moments =
                IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);
            for (int ii = -1; ii <= 1; ++ii) {
              for (int jj = -1; jj <= 1; ++jj) {
                for (int kk = -1; kk <= 1; ++kk) {
                  super_liquid_moments +=
                      a_liquid_moments(i + ii, j + jj, k + kk);
                  super_gas_moments += a_gas_moments(i + ii, j + jj, k + kk);
                }
              }
            }
            IRL::Pt liquid_centroid =
                IRL::Pt(super_liquid_moments[1], super_liquid_moments[2],
                        super_liquid_moments[3]);
            IRL::Pt gas_centroid =
                IRL::Pt(super_gas_moments[1], super_gas_moments[2],
                        super_gas_moments[3]);
            liquid_centroid *=
                1.0 / IRL::safelyEpsilon(super_liquid_moments[0]);
            gas_centroid *= 1.0 / IRL::safelyEpsilon(super_gas_moments[0]);

            ///////////////////////// EIGEN TEST
            auto centroid_line = IRL::Normal(gas_centroid - liquid_centroid);

            MOF2Functor myMOF2Functor(5, 20, super_cell, super_liquid_moments,
                                      super_gas_moments, cell,
                                      liquid_volume_fraction);
            myMOF2Functor.setframe(centroid_line, (*a_interface)(i, j, k));
            Eigen::NumericalDiff<MOF2Functor, Eigen::Forward> NDMOF2Functor(
                myMOF2Functor, 1.0e-12);
            Eigen::LevenbergMarquardt<
                Eigen::NumericalDiff<MOF2Functor, Eigen::Forward>, double>
                MOF2LM(NDMOF2Functor);
            MOF2LM.parameters.ftol = 1.0e-6;
            MOF2LM.parameters.xtol = 1.0e-6;
            // MOF2LM.parameters.factor = 1.0e6;
            // MOF2LM.parameters.epsfcn = std::sqrt(1.0e-6);
            MOF2LM.parameters.maxfev = 400;  // Max
            // iterations
            Eigen::VectorXd x(5);
            for (int ii = 0; ii < 5; ii++) x(ii) = 1.0;
            // MOF2LM.minimize(x);
            Eigen::VectorXd fvec(20);
            myMOF2Functor.errorvec(x, fvec);
            const double init_error_li = fvec.lpNorm<Eigen::Infinity>();
            double final_error_li = init_error_li;
            int it = 0, itmax = 1;
            do {
              for (int ii = 0; ii < 5; ii++) x(ii) = 1.0;
              Eigen::LevenbergMarquardtSpace::Status status =
                  MOF2LM.minimizeInit(x);
              do {
                status = MOF2LM.minimizeOneStep(x);
              } while (status == Eigen::LevenbergMarquardtSpace::Running);
              myMOF2Functor.errorvec(x, fvec);
              final_error_li = fvec.lpNorm<Eigen::Infinity>();
              MOF2LM.parameters.factor *= 10.0;
              it++;
            } while (init_error_li < final_error_li && it < itmax);

            IRL::Paraboloid paraboloid = myMOF2Functor.getparaboloid(x);
            if (0) {  // i == 42 && j == 31 && k == 16) {  // 0 &&
              //         // final_error_li > 0.01 && final_error_li >
              //         //  init_error_li) {
              // std::cout << "LINF = " << final_error_li << " from "
              //           << init_error_li << " at " << i << ", " << j << ", "
              //           << k << "VFRAC = "
              //           << super_liquid_moments[0] /
              //                  (volume_scale * cell_volume)
              //           << std::endl;
              // std::cout << "JIBBEN Paraboloid:\n"
              //           << (*a_interface)(i, j, k) << std::endl;
              // std::cout << "   MOF paraboloid:\n"
              //           << myMOF2Functor.getparaboloid(x) << std::endl;
              // std::cout << " fvec final =\n" << fvec << std::endl;

              // x(0) = 1.0;
              // x(1) = 1.0;
              // x(2) = 1.0;
              // x(3) = 1.0;
              // x(4) = 1.0;
              // myMOF2Functor.errorvec(x, fvec);
              // Eigen::Matrix<double, -1, -1> fjac;
              // fjac.resize(20, 5);
              // Eigen::Index df_ret = NDMOF2Functor.df(x, fjac);
              // std::cout << " fvec init =\n" << fvec << std::endl;
              // // std::cout << " fjac init =\n" << fjac << std::endl;
              // auto cell_moments =
              //     IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
              //         super_cell, IRL::Paraboloid::createAlwaysAbove());
              // RecenterMoments(&super_liquid_moments, cell_centroid);
              // RecenterMoments(&super_gas_moments, cell_centroid);
              // RecenterMoments(&cell_moments, cell_centroid);
              // std::cout << " liquid moments =\n"
              //           << super_liquid_moments * (1.0 / (27.0 *
              //           cell_volume))
              //           << std::endl;
              // std::cout << "   gas moments =\n"
              //           << super_gas_moments * (1.0 / (27.0 * cell_volume))
              //           << std::endl;
              // std::cout << "  cell moments =\n"
              //           << cell_moments * (1.0 / (27.0 * cell_volume))
              //           << std::endl;
              // // std::cout << "J^T fvec = " << fjac.transpose() * fvec
              // //           << std::endl;
              // // Eigen::Matrix<double, 5, 5> M =
              // //     fjac.transpose() * fjac +
              // //     100.0 * Eigen::Matrix<double, 5, 5>::Identity();
              // // std::cout << "J^T J + lambda I = " << M << std::endl;
              // // Eigen::VectorXd delta =
              // //     M.fullPivHouseholderQr().solve(fjac.transpose() * fvec);
              // // std::cout << "delta = " << delta << std::endl;
              // // exit(0);
            }
            auto coeffs = paraboloid.getAlignedParaboloid();
            double coeff_a = coeffs.a();
            double coeff_b = coeffs.b();
            if (std::abs(coeff_a * cell_dx) > 1.0)
              coeff_a = std::copysign(1.0 / cell_dx, coeff_a);
            if (std::abs(coeff_b * cell_dx) > 1.0)
              coeff_b = std::copysign(1.0 / cell_dx, coeff_b);
            paraboloid.setAlignedParaboloid(
                IRL::AlignedParaboloid({coeff_a, coeff_b}));
            const auto datum = paraboloid.getDatum();
            const auto fit_frame = paraboloid.getReferenceFrame();
            paraboloid = IRL::Paraboloid(datum, fit_frame, coeff_a, coeff_b);
            ///////////////////////// EIGEN TEST
            IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
                solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                                paraboloid);
            const auto test_volume =
                IRL::getVolumeMoments<IRL::Volume>(cell, paraboloid);
            if (test_volume < -1.0) {
              std::cout << "TEST VOLUME FAILED" << std::endl;
              exit(0);
            }
            auto new_datum =
                IRL::Pt(paraboloid.getDatum() +
                        solver_distance.getDistance() * fit_frame[2]);
            paraboloid.setDatum(new_datum);
            IRL::serializeAndPack(paraboloid, &interface_local);
            // (*a_interface)(i, j, k) = paraboloid;
          }
          count++;
        }
      }
    }
  }

  std::vector<int> proc_count(size);
  for (int r = 0; r < size; r++) {
    proc_count[r] = size_paraboloid * (proc_offset[r + 1] - proc_offset[r]);
    proc_offset[r] = size_paraboloid * proc_offset[r];
  }

  interface_global.resize(size_paraboloid * nmixed_global);
  MPI_Allgatherv(interface_local.data(), size_paraboloid * nmixed_local,
                 MPI_BYTE, interface_global.data(), proc_count.data(),
                 proc_offset.data(), MPI_BYTE, MPI_COMM_WORLD);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / cell_volume;
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          IRL::Paraboloid paraboloid;
          IRL::unpackAndStore(&paraboloid, &interface_global);
          (*a_interface)(i, j, k) = paraboloid;
        }
        //////////////////////////
        const auto coeffs = (*a_interface)(i, j, k).getAlignedParaboloid();
        if (std::abs(2.0 * coeffs.a() * cell_dx) > 2.0 ||
            std::abs(2.0 * coeffs.b() * cell_dx) > 2.0) {
          std::cout << "Error: coeffs = " << coeffs
                    << "; vfrac = " << liquid_volume_fraction << std::endl;
          exit(0);
        }
        //////////////////////
      }
    }
  }

  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);

  MPI_Barrier(MPI_COMM_WORLD);

  // if (second_pass) {
  //   Data<IRL::Paraboloid>
  //   corrected_interface(&mesh); for (int i =
  //   mesh.imin(); i <= mesh.imax(); ++i) {
  //     for (int j = mesh.jmin(); j <= mesh.jmax();
  //     ++j) {
  //       for (int k = mesh.kmin(); k <= mesh.kmax();
  //       ++k) {
  //         const double liquid_volume_fraction =
  //             a_liquid_moments(i, j, k)[0] /
  //             cell_volume;
  //         if (liquid_volume_fraction <
  //         IRL::global_constants::VF_LOW) {
  //           corrected_interface(i, j, k) =
  //           IRL::Paraboloid::createAlwaysBelow();
  //         } else if (liquid_volume_fraction >
  //         IRL::global_constants::VF_HIGH)
  //         {
  //           corrected_interface(i, j, k) =
  //           IRL::Paraboloid::createAlwaysAbove();
  //         } else {
  //           auto cell =
  //           IRL::RectangularCuboid::fromBoundingPts(
  //               IRL::Pt(mesh.x(i), mesh.y(j),
  //               mesh.z(k)), IRL::Pt(mesh.x(i + 1),
  //               mesh.y(j + 1), mesh.z(k + 1)));
  //           const auto datum = (*a_interface)(i, j,
  //           k).getDatum(); const auto ref_frame =
  //           (*a_interface)(i, j,
  //           k).getReferenceFrame(); const auto
  //           normal = ref_frame[2]; const double
  //           distance = normal * datum; const auto
  //           separator =
  //           IRL::PlanarSeparator::fromOnePlane(
  //               IRL::Plane(normal, distance));
  //           const auto polygon =
  //               IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
  //                   cell, separator, separator[0]);
  //           const IRL::Pt surface_centroid =
  //           polygon.calculateCentroid(); auto
  //           super_cell =
  //           IRL::RectangularCuboid::fromBoundingPts(
  //               IRL::Pt(surface_centroid[0] - 1.5 *
  //               mesh.dx(),
  //                       surface_centroid[1] - 1.5 *
  //                       mesh.dy(),
  //                       surface_centroid[2] - 1.5 *
  //                       mesh.dz()),
  //               IRL::Pt(surface_centroid[0] + 1.5 *
  //               mesh.dx(),
  //                       surface_centroid[1] + 1.5 *
  //                       mesh.dy(),
  //                       surface_centroid[2] + 1.5 *
  //                       mesh.dz()));
  //           const auto cell_centroid =
  //               IRL::Pt(mesh.xm(i), mesh.ym(j),
  //               mesh.zm(k));
  //           auto super_moments =
  //               IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
  //                   super_cell,
  //                   (*a_link_localized_paraboloid)(i,
  //                   j, k));
  //           const auto exact_super_moments =
  //               IRL::getVolumeMoments<IRL::VolumeMoments>(
  //                   super_cell,
  //                   (*a_link_localized_paraboloid)(i,
  //                   j, k));
  //           super_moments[0] =
  //           exact_super_moments.volume();
  //           super_moments[1] =
  //           exact_super_moments.centroid()[0];
  //           super_moments[2] =
  //           exact_super_moments.centroid()[1];
  //           super_moments[3] =
  //           exact_super_moments.centroid()[2];
  //           IRL::Pt liquid_centroid =
  //               IRL::Pt(super_moments[1],
  //               super_moments[2],
  //               super_moments[3]);
  //           IRL::Pt gas_centroid =
  //               27.0 * cell_volume * cell_centroid
  //               - liquid_centroid;
  //           liquid_centroid *= 1.0 /
  //           super_moments[0]; gas_centroid *= 1.0 /
  //           (27.0 * cell_volume -
  //           super_moments[0]);

  //           ///////////////////////// EIGEN TEST
  //           auto centroid_line =
  //           IRL::Normal(gas_centroid -
  //           liquid_centroid);

  //           MOF2Functor myMOF2Functor(4, 10,
  //           super_cell, super_moments);
  //           myMOF2Functor.setframe(centroid_line);
  //           Eigen::NumericalDiff<MOF2Functor>
  //           NDMOF2Functor(myMOF2Functor);
  //           Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF2Functor>,
  //           double>
  //               MOF2LM(NDMOF2Functor);
  //           // MOF2LM.parameters.ftol = 1.0e-5;
  //           // MOF2LM.parameters.xtol = 1.0e-5;
  //           // MOF2LM.parameters.factor = 1.0e3;
  //           MOF2LM.parameters.maxfev = 1000;  //
  //           Max iterations Eigen::VectorXd x(4);
  //           x.setZero(); MOF2LM.minimize(x);
  //           IRL::Paraboloid paraboloid =
  //           myMOF2Functor.getparaboloid(x); const
  //           auto fit_frame =
  //           paraboloid.getReferenceFrame();
  //           ///////////////////////// EIGEN TEST
  //           IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
  //               solver_distance(cell,
  //               liquid_volume_fraction, 1.0e-14,
  //                               paraboloid);
  //           if (solver_distance.getDistance() ==
  //           -DBL_MAX) {
  //             corrected_interface(i, j, k) =
  //                 IRL::Paraboloid(cell_centroid,
  //                 fit_frame, 1.0e-3, -1.0e-3);
  //           } else {
  //             auto new_datum =
  //                 IRL::Pt(paraboloid.getDatum() +
  //                         solver_distance.getDistance()
  //                         * fit_frame[2]);
  //             paraboloid.setDatum(new_datum);
  //             corrected_interface(i, j, k) =
  //             paraboloid;
  //           }
  //         }
  //       }
  //     }
  //   }
  //   for (int i = mesh.imin(); i <= mesh.imax();
  //   ++i) {
  //     for (int j = mesh.jmin(); j <= mesh.jmax();
  //     ++j) {
  //       for (int k = mesh.kmin(); k <= mesh.kmax();
  //       ++k) {
  //         const auto paraboloid =
  //         IRL::Paraboloid(corrected_interface(i, j,
  //         k));
  //         (*a_interface)(i, j, k) = paraboloid;
  //       }
  //     }
  //   }
  //   a_interface->updateBorder();
  //   correctInterfacePlaneBorders(a_interface);
  // }
}

void Jibben::getReconstruction(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::GeneralMoments3D<2>>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::Paraboloid>* a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  const BasicMesh& mesh = a_U.getMesh();

  Data<IRL::PlanarSeparator> interface(&mesh);
  updateReconstructionELVIRA(a_liquid_moments, &interface);
  updateReconstructionLVIRA(a_liquid_moments, 1, &interface);
  Data<IRL::Polygon> polygon(&mesh);
  updatePolygon(a_liquid_moments, interface, &polygon);
  polygon.updateBorder();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] -= mesh.lx();
        }
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] += mesh.lx();
        }
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] -= mesh.ly();
        }
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] += mesh.ly();
        }
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] -= mesh.lz();
        }
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] += mesh.lz();
        }
      }
    }
  }

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          // continue;
        } else {
          const IRL::Normal norm_poly = polygon(i, j, k).calculateNormal();
          const double poly_area = polygon(i, j, k).calculateVolume();
          const IRL::Pt pref = polygon(i, j, k).calculateCentroid();
          IRL::ReferenceFrame fit_frame;
          int largest_dir = 0;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[1]))
            largest_dir = 1;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[2]))
            largest_dir = 2;
          if (largest_dir == 0)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 1)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 0.0, 1.0));
          else
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(1.0, 0.0, 0.0));
          fit_frame[0].normalize();
          fit_frame[1] = crossProduct(norm_poly, fit_frame[0]);
          fit_frame[2] = norm_poly;
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
          IRL::Paraboloid paraboloid;

          double sum_vfrac = 0.0;
          for (int kk = -1; kk < 2; ++kk) {
            for (int jj = -1; jj < 2; ++jj) {
              for (int ii = -1; ii < 2; ++ii) {
                sum_vfrac += a_liquid_moments(i + ii, j + jj, k + kk)[0];
              }
            }
          }
          sum_vfrac /= mesh.cell_volume();

          auto sol_fit = fitParaboloidToPLICHeights(
              polygon, a_liquid_moments, pref, fit_frame, i, j, k, 1, 0.0);
          const double a = sol_fit[0], b = sol_fit[1], c = sol_fit[2],
                       d = sol_fit[3], e = sol_fit[4], f = sol_fit[5];
          const double theta = 0.5 * std::atan2(e, (IRL::safelyTiny(d - f)));
          const double cos_t = std::cos(theta);
          const double sin_t = std::sin(theta);
          const double A =
              -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
          const double B =
              -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);
          // Translation to coordinate system R' where aligned paraboloid valid
          //     Translation is R ' = {x' = x + u, y ' = y + v, z' = z + w }
          const double denominator = IRL::safelyTiny(4.0 * d * f - e * e);
          const double u = (2.0 * b * f - c * e) / denominator;
          const double v = -(b * e - 2.0 * d * c) / denominator;
          const double w =
              -(a + (-b * b * f + b * c * e - c * c * d) / denominator);

          IRL::UnitQuaternion rotation(theta, fit_frame[2]);
          IRL::Pt datum =
              pref - u * fit_frame[0] - v * fit_frame[1] - w * fit_frame[2];
          auto new_frame = rotation * fit_frame;
          const double max_curvature_dx = 2.0;
          double a_coeff = A;
          double b_coeff = B;
          if (std::sqrt(u * u + v * v + w * w) > 10.0 * mesh.dx() ||
              std::fabs(A) * mesh.dx() > max_curvature_dx ||
              std::fabs(B) * mesh.dx() > max_curvature_dx) {
            paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
          } else {
            if (fabs(a_coeff) < 1.0e-3) {
              a_coeff = std::copysign(1.0e-3, a_coeff);
            }
            if (fabs(b_coeff) < 1.0e-3) {
              b_coeff = std::copysign(1.0e-3, b_coeff);
            }
            paraboloid = IRL::Paraboloid(datum, new_frame, a_coeff, b_coeff);
          }

          // New approach
          // Eigen::Matrix2d W{{2.0 * d, e}, {e, 2.0 * f}};
          // Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigendecomp(W);
          // double a_coeff = -0.5 * eigendecomp.eigenvalues()[0];
          // double b_coeff = -0.5 * eigendecomp.eigenvalues()[1];
          // auto new_frame = fit_frame;
          // new_frame[0] =
          //     IRL::Normal(eigendecomp.eigenvectors().col(0)[0] * fit_frame[0]
          //     +
          //                 eigendecomp.eigenvectors().col(0)[1] *
          //                 fit_frame[1]);
          // new_frame[0].normalize();
          // new_frame[1] = IRL::crossProduct(fit_frame[2], new_frame[0]);
          // new_frame[1].normalize();
          // new_frame[2] =
          //     IRL::Normal(fit_frame[2] - b * fit_frame[0] - c *
          //     fit_frame[1]);
          // new_frame[2].normalize();
          // new_frame[0] = IRL::crossProduct(new_frame[1], new_frame[2]);
          // new_frame[0].normalize();
          // new_frame[1] = IRL::crossProduct(new_frame[2], new_frame[0]);
          // new_frame[1].normalize();
          // IRL::Pt datum = pref + a * fit_frame[2];
          // paraboloid = IRL::Paraboloid(datum, new_frame, a_coeff, b_coeff);

          // const double max_curvature_dx = 2.0;
          // if (std::fabs(a_coeff) * mesh.dx() > max_curvature_dx ||
          //     std::fabs(b_coeff) * mesh.dx() > max_curvature_dx) {
          //   paraboloid = IRL::Paraboloid(datum, new_frame, 1.0e-3, -1.0e-3);
          // }

          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                              paraboloid);

          auto new_datum =
              IRL::Pt(paraboloid.getDatum() +
                      solver_distance.getDistance() * new_frame[2]);
          paraboloid.setDatum(new_datum);
          (*a_interface)(i, j, k) = paraboloid;
        }
      }
    }
  }

  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void Centroid::getReconstruction(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::GeneralMoments3D<2>>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::Paraboloid>* a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  const BasicMesh& mesh = a_U.getMesh();

  Data<IRL::PlanarSeparator> interface(&mesh);
  updateReconstructionELVIRA(a_liquid_moments, &interface);
  updateReconstructionLVIRA(a_liquid_moments, 1, &interface);
  Data<IRL::Polygon> polygon(&mesh);
  updatePolygon(a_liquid_moments, interface, &polygon);
  polygon.updateBorder();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] -= mesh.lx();
        }
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] += mesh.lx();
        }
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] -= mesh.ly();
        }
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] += mesh.ly();
        }
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] -= mesh.lz();
        }
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] += mesh.lz();
        }
      }
    }
  }

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          // continue;
        } else {
          const IRL::Normal norm_poly = polygon(i, j, k).calculateNormal();
          const double poly_area = polygon(i, j, k).calculateVolume();
          const IRL::Pt pref = polygon(i, j, k).calculateCentroid();
          IRL::ReferenceFrame fit_frame;
          int largest_dir = 0;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[1]))
            largest_dir = 1;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[2]))
            largest_dir = 2;
          if (largest_dir == 0)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 1)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 0.0, 1.0));
          else
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(1.0, 0.0, 0.0));
          fit_frame[0].normalize();
          fit_frame[1] = crossProduct(norm_poly, fit_frame[0]);
          fit_frame[2] = norm_poly;
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
          IRL::Paraboloid paraboloid;

          double sum_vfrac = 0.0;
          for (int kk = -1; kk < 2; ++kk) {
            for (int jj = -1; jj < 2; ++jj) {
              for (int ii = -1; ii < 2; ++ii) {
                sum_vfrac += a_liquid_moments(i + ii, j + jj, k + kk)[0];
              }
            }
          }
          sum_vfrac /= mesh.cell_volume();

          if (std::fabs(sum_vfrac - 27.0 * 0.5) < 10.0 &&
              liquid_volume_fraction > 1.0e-6 &&
              liquid_volume_fraction < 1.0 - 1.0e-6) {
            auto sol_fit = fitParaboloidToCentroids(
                polygon, a_liquid_moments, pref, fit_frame, i, j, k, 1, 0.0);
            const double a = sol_fit[0], b = sol_fit[1], c = sol_fit[2],
                         d = sol_fit[3], e = sol_fit[4], f = sol_fit[5];
            const double theta = 0.5 * std::atan2(e, (IRL::safelyTiny(d - f)));
            const double cos_t = std::cos(theta);
            const double sin_t = std::sin(theta);
            const double A =
                -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
            const double B =
                -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);
            // Translation to coordinate system R' where aligned paraboloid
            // valid Translation is R' = {x' = x + u, y' = y + v, z' = z + w}
            const double denominator = IRL::safelyTiny(4.0 * d * f - e * e);
            const double u = (2.0 * b * f - c * e) / denominator;
            const double v = -(b * e - 2.0 * d * c) / denominator;
            const double w =
                -(a + (-b * b * f + b * c * e - c * c * d) / denominator);

            IRL::UnitQuaternion rotation(theta, fit_frame[2]);
            IRL::Pt datum =
                pref - u * fit_frame[0] - v * fit_frame[1] - w * fit_frame[2];
            auto new_frame = rotation * fit_frame;
            const double max_curvature_dx = 1.0;
            double a_coeff = A;
            double b_coeff = B;
            if (std::sqrt(u * u + v * v + w * w) > 10.0 * mesh.dx() ||
                std::fabs(A) * mesh.dx() > max_curvature_dx ||
                std::fabs(B) * mesh.dx() > max_curvature_dx) {
              paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
            } else {
              paraboloid = IRL::Paraboloid(datum, new_frame, a_coeff, b_coeff);
            }
          } else {
            paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
          }

          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                              paraboloid);
          if (solver_distance.getDistance() == -DBL_MAX) {
            (*a_interface)(i, j, k) =
                IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
          } else {
            auto new_datum =
                IRL::Pt(paraboloid.getDatum() +
                        solver_distance.getDistance() * fit_frame[2]);
            paraboloid.setDatum(new_datum);
            (*a_interface)(i, j, k) = paraboloid;
          }
        }
      }
    }
  }

  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void PLIC::getReconstruction(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::GeneralMoments3D<2>>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::Paraboloid>* a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  const BasicMesh& mesh = a_U.getMesh();

  Data<IRL::PlanarSeparator> interface(&mesh);
  updateReconstructionELVIRA(a_liquid_moments, &interface);
  updateReconstructionLVIRA(a_liquid_moments, 1, &interface);
  Data<IRL::Polygon> polygon(&mesh);
  updatePolygon(a_liquid_moments, interface, &polygon);
  polygon.updateBorder();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] -= mesh.lx();
        }
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] += mesh.lx();
        }
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] -= mesh.ly();
        }
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] += mesh.ly();
        }
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] -= mesh.lz();
        }
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] += mesh.lz();
        }
      }
    }
  }

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          // continue;
        } else {
          const IRL::Normal norm_poly = polygon(i, j, k).calculateNormal();
          const IRL::Pt pref = polygon(i, j, k).calculateCentroid();
          IRL::ReferenceFrame fit_frame;
          int largest_dir = 0;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[1]))
            largest_dir = 1;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[2]))
            largest_dir = 2;
          if (largest_dir == 0)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 1)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 0.0, 1.0));
          else
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(1.0, 0.0, 0.0));
          fit_frame[0].normalize();
          fit_frame[1] = crossProduct(norm_poly, fit_frame[0]);
          fit_frame[2] = norm_poly;
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
          IRL::Paraboloid paraboloid;

          paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);

          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                              paraboloid);

          if (solver_distance.getDistance() == -DBL_MAX) {
            (*a_interface)(i, j, k) =
                IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
            std::cout << "Distance solver failed: error = " << std::endl;
          } else {
            auto new_datum =
                IRL::Pt(paraboloid.getDatum() +
                        solver_distance.getDistance() * fit_frame[2]);
            paraboloid.setDatum(new_datum);
            (*a_interface)(i, j, k) = paraboloid;
          }
        }
      }
    }
  }

  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void correctInterfacePlaneBorders(Data<IRL::Paraboloid>* a_interface) {
  const BasicMesh& mesh = (*a_interface).getMesh();
  // Fix distances in reconstruction for periodic boundary

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[0] -= mesh.lx();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[0] += mesh.lx();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[1] -= mesh.ly();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[1] += mesh.ly();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[2] -= mesh.lz();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[2] += mesh.lz();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }
}
