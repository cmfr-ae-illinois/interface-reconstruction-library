// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include <mpi.h>
#include <limits>
#include <tuple>

#include "examples/variant_advector/reconstruction_types.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/interface_reconstruction_methods/constrained_optimization_behavior.h"
#include "irl/interface_reconstruction_methods/elvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/jibben_neighborhood.h"
#include "irl/interface_reconstruction_methods/plvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"
#include "irl/optimization/constrained_levenberg_marquardt.h"
#include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
#include "irl/paraboloid_reconstruction/hessian_paraboloid.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/planar_separator.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/Polynomials>
#include "examples/variant_advector/basic_mesh.h"
#include "examples/variant_advector/data.h"
#include "examples/variant_advector/vof_advection.h"

void getReconstruction(const std::string& a_reconstruction_method,
                       const Data<IRL::VolumeMoments>& a_liq_moments,
                       const Data<IRL::VolumeMoments>& a_gas_moments,
                       const double a_dt, const Data<double>& a_U,
                       const Data<double>& a_V, const Data<double>& a_W,
                       Data<IRL::SeparatorVariant>* a_interface) {
  if (a_reconstruction_method == "ELVIRA") {
    ELVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface);
  } else if (a_reconstruction_method == "LVIRA" ||
             a_reconstruction_method == "PLIC") {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface);
  } else if (a_reconstruction_method == "Jibben") {
    Jibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface);
  } else if (a_reconstruction_method == "PU") {
    PU::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                          a_interface);
  } else if (a_reconstruction_method == "MixedJibben") {
    MixedPLICJibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                       a_V, a_W, a_interface);
  } else if (a_reconstruction_method == "SlicesParabola") {
    SlicesParabola::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                      a_V, a_W, a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout
        << "Valid entries are: PLIC, Jibben, MixedJibben, SlicesParabola. \n";
    std::exit(-1);
  }
}

void ELVIRA::getReconstruction(const Data<IRL::VolumeMoments>& a_liq_moments,
                               const Data<IRL::VolumeMoments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V, const Data<double>& a_W,
                               Data<IRL::SeparatorVariant>* a_interface) {
  const BasicMesh& mesh = a_liq_moments.getMesh();

  IRL::ELVIRANeighborhood neighborhood;
  neighborhood.resize(27);
  IRL::RectangularCuboid cells[27];
  std::array<double, 27> cells_vfrac;
  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
            liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, liquid_volume_fraction - 0.5);
          (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }
        // Build surrounding stencil information for ELVIRA.
        for (int kk = k - 1; kk < k + 2; ++kk) {
          for (int jj = j - 1; jj < j + 2; ++jj) {
            for (int ii = i - 1; ii < i + 2; ++ii) {
              // Reversed order, bad for cache locality but thats okay..
              const int neigh_id =
                  (kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1);
              cells[neigh_id] = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                  IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
              cells_vfrac[neigh_id] =
                  a_liq_moments(ii, jj, kk).volume() / mesh.cell_volume();
              neighborhood.setMember(&cells[neigh_id], &cells_vfrac[neigh_id],
                                     ii - i, jj - j, kk - k);
            }
          }
        }
        // Now perform actual ELVIRA and obtain interface PlanarSeparator
        (*a_interface)(i, j, k) = reconstructionWithELVIRA3D(neighborhood);
      }
    }
  }
  // Update border with simple ghost-cell fill and correct distances for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

void LVIRA::getReconstruction(const Data<IRL::VolumeMoments>& a_liq_moments,
                              const Data<IRL::VolumeMoments>& a_gas_moments,
                              const double a_dt, const Data<double>& a_U,
                              const Data<double>& a_V, const Data<double>& a_W,
                              Data<IRL::SeparatorVariant>* a_interface,
                              const bool a_plic_already_built) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    ELVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface);
  }

  const BasicMesh& mesh = a_liq_moments.getMesh();

  IRL::LVIRANeighborhood<IRL::RectangularCuboid> neighborhood;
  neighborhood.resize(27);
  IRL::RectangularCuboid cells[27];
  std::array<double, 27> cells_vfrac;
  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
            liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, liquid_volume_fraction - 0.5);
          (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }
        // Build surrounding stencil information for LVIRA.
        for (int kk = k - 1; kk < k + 2; ++kk) {
          for (int jj = j - 1; jj < j + 2; ++jj) {
            for (int ii = i - 1; ii < i + 2; ++ii) {
              // Reversed order, bad for cache locality but thats okay..
              const int neigh_id =
                  (kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1);
              cells[neigh_id] = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                  IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
              cells_vfrac[neigh_id] =
                  a_liq_moments(ii, jj, kk).volume() / mesh.cell_volume();
              neighborhood.setMember(
                  static_cast<IRL::UnsignedIndex_t>(neigh_id), &cells[neigh_id],
                  &cells_vfrac[neigh_id]);
              // Trap center cell
              if (ii == i && jj == j && kk == k) {
                neighborhood.setCenterOfStencil(neigh_id);
              }
            }
          }
        }
        // Now perform actual LVIRA and obtain interface PlanarSeparator
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        (*a_interface)(i, j, k) =
            IRL::reconstructionWithLVIRA3D(neighborhood, planar_separator);
      }
    }
  }
  // Update border with simple ghost-cell fill and correct distances for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

void updatePtBorder(Data<IRL::Pt>* a_pts) {
  const BasicMesh& mesh = a_pts->getMesh();
  a_pts->updateBorder();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_pts)(i, j, k)[0] -= mesh.lx();
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_pts)(i, j, k)[0] += mesh.lx();
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_pts)(i, j, k)[1] -= mesh.ly();
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_pts)(i, j, k)[1] += mesh.ly();
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        (*a_pts)(i, j, k)[2] -= mesh.lz();
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        (*a_pts)(i, j, k)[2] += mesh.lz();
      }
    }
  }
}

void updatePolygonBorder(Data<IRL::Polygon>* a_polygons) {
  const BasicMesh& mesh = a_polygons->getMesh();
  a_polygons->updateBorder();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : (*a_polygons)(i, j, k)) {
          pt[0] -= mesh.lx();
        }
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : (*a_polygons)(i, j, k)) {
          pt[0] += mesh.lx();
        }
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : (*a_polygons)(i, j, k)) {
          pt[1] -= mesh.ly();
        }
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : (*a_polygons)(i, j, k)) {
          pt[1] += mesh.ly();
        }
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        for (auto& pt : (*a_polygons)(i, j, k)) {
          pt[2] -= mesh.lz();
        }
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        for (auto& pt : (*a_polygons)(i, j, k)) {
          pt[2] += mesh.lz();
        }
      }
    }
  }
}

const IRL::ReferenceFrame referenceFrameFromNormal(const IRL::Normal a_normal) {
  IRL::ReferenceFrame frame;
  int largest_dir = 0;
  if (std::fabs(a_normal[largest_dir]) < std::fabs(a_normal[1]))
    largest_dir = 1;
  if (std::fabs(a_normal[largest_dir]) < std::fabs(a_normal[2]))
    largest_dir = 2;
  if (largest_dir == 0)
    frame[0] = IRL::crossProduct(a_normal, IRL::Normal(0.0, 1.0, 0.0));
  else if (largest_dir == 1)
    frame[0] = IRL::crossProduct(a_normal, IRL::Normal(0.0, 0.0, 1.0));
  else
    frame[0] = IRL::crossProduct(a_normal, IRL::Normal(1.0, 0.0, 0.0));
  frame[0].normalize();
  frame[1] = IRL::crossProduct(a_normal, frame[0]);
  frame[2] = a_normal;
  return frame;
}

// This is modified from cut_polygon.tpp
void getIntersectionPts(const IRL::Polygon& a_polygon,
                        const IRL::Plane& a_cutting_plane,
                        IRL::StackVector<IRL::Pt, 2>* a_intersection_pts) {
  a_intersection_pts->resize(0);
  double distance = a_cutting_plane.signedDistanceToPoint(a_polygon[0]);
  for (int n = 0; n < a_polygon.getNumberOfVertices(); ++n) {
    const int next_id = (n + 1) % a_polygon.getNumberOfVertices();
    double next_distance =
        a_cutting_plane.signedDistanceToPoint(a_polygon[next_id]);
    if (distance * next_distance < 0.0) {
      a_intersection_pts->push_back(IRL::Pt::fromEdgeIntersection(
          a_polygon[n], distance, a_polygon[next_id], next_distance));
      if (a_intersection_pts->size() == 2) {
        break;
      }
    }
    distance = next_distance;
  }
}

// This is modified from general_polygon.tpp
const IRL::Normal calculatePolygonNormal(const IRL::Polygon& a_polygon) {
  const int nverts = a_polygon.getNumberOfVertices();
  if (nverts < 3) {
    return IRL::Normal(0, 0, 0);
  }
  for (int n = 0; n < nverts; ++n) {
    const IRL::Pt p0 = a_polygon[n];
    const IRL::Pt p1 = a_polygon[(n + 1) % nverts];
    const IRL::Pt p2 = a_polygon[(n + 2) % nverts];
    const IRL::Normal normal = IRL::crossProduct(p1 - p0, p2 - p0);
    const double normal_magnitude = IRL::magnitude(normal);
    if (normal_magnitude > 100.0 * DBL_EPSILON) {
      return normal / normal_magnitude;
    }
  }
  return IRL::Normal(0, 0, 0);
}

// Based on the paper by Chernova and Wijewickremab, JCAM 2013
const IRL::Pt closestPointOnParaboloid(const IRL::Paraboloid& a_paraboloid,
                                       const IRL::Pt& a_pt) {
  // Get paraboloid information
  const double a = a_paraboloid.getAlignedParaboloid().a();
  const double b = a_paraboloid.getAlignedParaboloid().b();
  const IRL::Pt& datum = a_paraboloid.getDatum();
  const IRL::ReferenceFrame& frame = a_paraboloid.getReferenceFrame();

  // Bring point in local frame of reference
  const IRL::Pt pt_tmp = a_pt - datum;
  IRL::Pt local_pt;
  for (int n = 0; n < 3; ++n) {
    local_pt[n] = frame[n] * pt_tmp;
  }

  // Compute quintic polynomial coeffs
  const double eps = std::numeric_limits<double>::epsilon();
  auto A = -16. * (a * a) * (b * b);
  auto B = 16. * a * b * (-b + a * (-1. + b * local_pt[2]));
  auto C =
      -4. * (a * a + 4. * a * b + b * b) + 16. * a * b * (a + b) * local_pt[2];
  auto D =
      4. *
      (a * a * (b * (local_pt[1] * local_pt[1]) + local_pt[2]) +
       b * (-1. + b * local_pt[2]) +
       a * (-1. + b * b * (local_pt[0] * local_pt[0]) + 4. * b * local_pt[2]));
  auto E =
      -1. +
      4. * a * b * (local_pt[0] * local_pt[0] + local_pt[1] * local_pt[1]) +
      4. * (a + b) * local_pt[2];
  auto F = a * (local_pt[0] * local_pt[0]) + b * (local_pt[1] * local_pt[1]) +
           local_pt[2];
  std::vector<double> t_vals(0);
  t_vals.reserve(5);

  if (std::abs(A) > eps) {
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    Eigen::Matrix<double, 6, 1> coeff(F, E, D, C, B, A);
    solver.compute(coeff);
    for (int i = 0; i < solver.roots().size(); ++i) {
      if (std::fabs(std::imag(solver.roots()[i])) < 1.0e4 * eps) {
        const double sol = std::real(solver.roots()[i]);
        t_vals.push_back(sol);
      }
    }
  } else if (std::abs(B) > eps) {
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    Eigen::Matrix<double, 5, 1> coeff(F, E, D, C, B);
    solver.compute(coeff);
    for (int i = 0; i < solver.roots().size(); ++i) {
      if (std::fabs(std::imag(solver.roots()[i])) < 1.0e4 * eps) {
        const double sol = std::real(solver.roots()[i]);
        t_vals.push_back(sol);
      }
    }
  } else if (std::abs(C) > eps) {
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    Eigen::Matrix<double, 4, 1> coeff(F, E, D, C);
    solver.compute(coeff);
    for (int i = 0; i < solver.roots().size(); ++i) {
      if (std::fabs(std::imag(solver.roots()[i])) < 1.0e4 * eps) {
        const double sol = std::real(solver.roots()[i]);
        t_vals.push_back(sol);
      }
    }

  } else if (std::abs(D) > eps) {
    double discriminant = E * E - 4.0 * D * F;
    if (discriminant > 0.0) {
      discriminant = std::sqrt(discriminant);
      const double q = -0.5 * (E + std::copysign(discriminant, E));
      const double safe_q = std::fabs(q) < eps ? std::copysign(eps, q) : q;
      const double sol1 = q / D;
      const double sol2 = F / safe_q;
      t_vals.push_back(sol1);
      t_vals.push_back(sol2);
    }
  } else if (std::abs(E) > eps) {
    const double sol = -F / E;
    t_vals.push_back(sol);
  }

  if (t_vals.size() > 0) {
    std::vector<std::pair<double, double>> t_and_dist(t_vals.size());
    for (int i = 0; i < t_vals.size(); i++) {
      const double t = t_vals[i];
      const double d = IRL::magnitude(
          IRL::Pt(local_pt[0] / (1. + 2. * a * t) - local_pt[0],
                  local_pt[1] / (1. + 2. * b * t) - local_pt[0], -t));
      t_and_dist[i] = std::make_pair(d, t);
    }
    std::sort(t_and_dist.begin(), t_and_dist.end());
    const double t = t_and_dist[0].second;
    auto closest_pt_local =
        IRL::Pt(local_pt[0] / (1. + 2. * a * t),
                local_pt[1] / (1. + 2. * b * t), local_pt[2] - t);
    closest_pt_local[2] = -a * closest_pt_local[0] * closest_pt_local[0] -
                          b * closest_pt_local[1] * closest_pt_local[1];
    auto closest_pt = IRL::Pt(0., 0., 0.);
    for (int d = 0; d < 3; ++d) {
      for (int n = 0; n < 3; ++n) {
        closest_pt[n] += frame[d][n] * closest_pt_local[d];
      }
    }
    closest_pt += datum;
    return closest_pt;
  }

  return a_pt;
}

const IRL::Paraboloid regenerateParaboloid(const IRL::Paraboloid& a_paraboloid,
                                           const IRL::Pt& a_pt) {
  const IRL::Pt closest_pt = closestPointOnParaboloid(a_paraboloid, a_pt);

  // Get paraboloid information
  const double a = a_paraboloid.getAlignedParaboloid().a();
  const double b = a_paraboloid.getAlignedParaboloid().b();
  const IRL::Pt& datum = a_paraboloid.getDatum();
  const IRL::ReferenceFrame& frame = a_paraboloid.getReferenceFrame();

  // Bring point in local frame of reference
  const IRL::Pt pt_tmp = closest_pt - datum;
  IRL::Pt local_pt;
  for (int n = 0; n < 3; ++n) {
    local_pt[n] = frame[n] * pt_tmp;
  }

  // Compute local derivatives
  const Eigen::Vector3d gradF(2. * a * local_pt[0], 2. * b * local_pt[1], 1.);
  Eigen::Matrix3d hessF = Eigen::Matrix3d::Zero();
  hessF(0, 0) = 2. * a;
  hessF(1, 1) = 2. * b;
  Eigen::Matrix3d adjHessF = Eigen::Matrix3d::Zero();
  adjHessF(2, 2) = 4. * a * b;
  auto new_normal = IRL::Normal(gradF(0), gradF(1), gradF(2));
  new_normal.normalize();

  // Based on Goldman 2005
  double H = gradF.transpose() * (hessF * gradF);
  H -= gradF.squaredNorm() * hessF.trace();
  H /= 2. * gradF.squaredNorm() * gradF.norm();
  double K = gradF.transpose() * (adjHessF * gradF);
  K /= gradF.squaredNorm() * gradF.squaredNorm();
  const double k1 = -H + std::sqrt(std::max(H * H - K, 0.));
  const double k2 = -H - std::sqrt(std::max(H * H - K, 0.));

  // Compute principal directions
  const double B = a - b, C = -gradF(1) * a, E = gradF(0) * b;
  const double U = 2. * gradF(0) * gradF(1) * a;
  const double V = 2. * (B - C * gradF(1) - E * gradF(0));
  const double W = -2. * gradF(1) * gradF(0) * b;
  const double delta = V * V - 4. * U * W;
  const double X1 = -V - std::sqrt(std::max(0., delta));
  auto T1 = IRL::Normal(X1, 2. * U, -X1 * gradF(0) - 2.0 * U * gradF(1));
  const double normT1 = IRL::magnitude(T1);
  IRL::ReferenceFrame new_local_frame;
  if (normT1 > DBL_EPSILON) {
    T1 *= 1.0 / normT1;
    new_local_frame[2] = new_normal;
    new_local_frame[0] = T1;
    new_local_frame[1] =
        IRL::crossProduct(new_local_frame[2], new_local_frame[0]);
  } else {  /// The point is umbilical
    new_local_frame = referenceFrameFromNormal(new_normal);
  }

  // Now move frame back to canonical frame and return paraboloid
  IRL::ReferenceFrame new_frame;
  for (int n = 0; n < 3; n++) {
    new_frame[n] = IRL::Normal(0., 0., 0.);
    for (int d = 0; d < 3; d++) {
      new_frame[n] += new_local_frame[n][d] * frame[d];
    }
    new_frame[n].normalize();
  }
  return IRL::Paraboloid(closest_pt, new_frame, 0.5 * k1, 0.5 * k2);
}

const IRL::Paraboloid generateLocalParaboloid(const IRL::Pt& a_pt,
                                              const Eigen::Vector3d& a_gradF,
                                              const Eigen::Matrix3d& a_hessF) {
  const Eigen::Matrix3d hessF = 0.5 * (a_hessF + a_hessF.transpose());

#if 0  // This uses the method of Che, Paul, Zhang, CAGD 2007 (-> Che2007 on
       // Zotero)
  // Compute adjugate hessian
  Eigen::Matrix3d adjHessF = Eigen::Matrix3d::Zero();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      adjHessF(i, j) =
          hessF((i + 1) % 3, (j + 1) % 3) * hessF((i + 2) % 3, (j + 2) % 3) -
          hessF((i + 1) % 3, (j + 2) % 3) * hessF((i + 2) % 3, (j + 1) % 3);
    }
  }

  // Based on Goldman 2005
  double H = a_gradF.transpose() * hessF * a_gradF;
  H -= a_gradF.squaredNorm() * hessF.trace();
  H /= 2. * a_gradF.squaredNorm() * a_gradF.norm();
  double K = a_gradF.transpose() * adjHessF * a_gradF;
  K /= a_gradF.squaredNorm() * a_gradF.squaredNorm();
  const double k1 = H + std::sqrt(std::max(H * H - K, 0.));
  const double k2 = H - std::sqrt(std::max(H * H - K, 0.));

  // Compute principal directions
  const double fx = a_gradF(0), fy = a_gradF(1), fz = a_gradF(2);
  const double fxx = hessF(0, 0), fxy = hessF(0, 1), fxz = hessF(0, 2);
  const double fyx = hessF(1, 0), fyy = hessF(1, 1), fyz = hessF(1, 2);
  const double fzx = hessF(2, 0), fzy = hessF(2, 1), fzz = hessF(2, 2);
  const double A = fy * fzx - fz * fyx;
  const double B = .5 * (fz * fxx - fx * fzx + fy * fzy - fz * fyy);
  const double C = .5 * (fy * fzz - fz * fyz + fx * fyx - fy * fxx);
  const double D = fz * fxy - fx * fzy;
  const double E = .5 * (fx * fyy - fy * fxy + fz * fxz - fx * fzz);
  const double F = fx * fyz - fy * fxz;
  const double U = A * fz * fz - 2. * C * fx * fz + F * fx * fx;
  const double V = 2. * (B * fz * fz - C * fy * fz - E * fx * fz + F * fx * fy);
  const double W = D * fz * fz - 2. * E * fy * fz + F * fy * fy;
  const double delta = V * V - 4. * U * W;
  const double sqrt_delta = std::sqrt(std::max(0., delta));
  const double X1 = -V + std::copysign(sqrt_delta, fz);
  auto T1 = IRL::Normal(X1 * fz, 2. * U * fz, -X1 * fx - 2. * U * fy);
  const double normT1 = IRL::magnitude(T1);

  auto new_normal = IRL::Normal(fx, fy, fz);
  new_normal.normalize();
  IRL::ReferenceFrame new_frame;
  if (normT1 > DBL_EPSILON) {
    T1 *= 1.0 / normT1;
    new_frame[2] = new_normal;
    new_frame[0] = T1;
    new_frame[1] = IRL::crossProduct(new_normal, T1);
  } else {  /// The point is umbilical
    new_frame = referenceFrameFromNormal(new_normal);
  }
  return IRL::Paraboloid(a_pt, new_frame, -0.5 * k1, -0.5 * k2);
#endif

  // This uses the method described in
  // https://www.geometrictools.com/Documentation/PrincipalCurvature.pdf
  const double inv_gradF_norm = 1. / a_gradF.norm();
  const IRL::Normal normal =
      IRL::Normal(a_gradF(0), a_gradF(1), a_gradF(2)) * inv_gradF_norm;
  IRL::ReferenceFrame frame = referenceFrameFromNormal(normal);
  Eigen::MatrixXd J(3, 2);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      J(i, j) = frame[j][i];
    }
  }
  const Eigen::Matrix2d A = (J.transpose() * hessF * J) * inv_gradF_norm;
  Eigen::EigenSolver<Eigen::Matrix2d> eigensolver(A);
  const double eval1 = eigensolver.eigenvalues()(0).real();
  const double eval2 = eigensolver.eigenvalues()(1).real();
  const Eigen::Vector2d evec1 =
      Eigen::Vector2d(eigensolver.eigenvectors()(0, 0).real(),
                      eigensolver.eigenvectors()(1, 0).real());
  const Eigen::Vector3d T1 = J * evec1;
  frame[0] = IRL::Normal(T1(0), T1(1), T1(2));
  frame[0].normalize();
  frame[1] = IRL::crossProduct(frame[2], frame[0]);
  return IRL::Paraboloid(a_pt, frame, 0.5 * eval1, 0.5 * eval2);

  ////////// TODO: compute eigenvalues and eigenvectors "by hand"
  // const double m11 = A(0, 0);
  // const double m22 = A(1, 1);
  // const double m12 = 0.5 * (A(1, 0) + A(0, 1));
  // // Compute eigenvalues and Givens rotation angle
  // const double tmp_sqrt = std::sqrt((m11 - m22) * (m11 - m22) + 4.0 * m12 *
  // m12); const double lambda1 = 0.5 * (m11 + m22 + tmp_sqrt); const double
  // lambda2 = 0.5 * (m11 + m22 - tmp_sqrt); const double theta = 0.5 *
  // std::atan2(2.0 * m12, m11 - m22);
  // // Compute principal directions (Darboux frame)
  // const IRL::UnitQuaternion givens_rotation(-theta, normal);
  // const auto darboux_frame = givens_rotation * frame;
  // return IRL::Paraboloid(a_pt, darboux_frame, 0.5 * lambda2, 0.5 * lambda1);
}

void Jibben::getReconstruction(const Data<IRL::VolumeMoments>& a_liq_moments,
                               const Data<IRL::VolumeMoments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V, const Data<double>& a_W,
                               Data<IRL::SeparatorVariant>* a_interface,
                               const bool a_plic_already_built,
                               Data<IRL::Pt>* a_centroids,
                               Data<double>* a_areas, Data<double>* a_errors) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface);
  }

  // Then, let's compute the PLIC polygons
  const BasicMesh& mesh = a_liq_moments.getMesh();

  Data<IRL::Polygon> polygon(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        polygon(i, j, k) = IRL::Polygon();
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
            liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        // Note: this std::get call is safe since all interfaces should be
        // of the type PlanarSeparator after calling the PLIC reconstruction
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, planar_separator, planar_separator[0]);
      }
    }
  }
  updatePolygonBorder(&polygon);

  if (a_centroids != nullptr && a_areas != nullptr) {
    for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
      for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
        for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
          const double liquid_volume_fraction =
              a_liq_moments(i, j, k).volume() / mesh.cell_volume();
          if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
              liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
            continue;
          }
          (*a_centroids)(i, j, k) = polygon(i, j, k).calculateCentroid();
          (*a_areas)(i, j, k) = polygon(i, j, k).calculateVolume();
        }
      }
    }
  }

  // Now let's compute the Jibben parabolic fit
  IRL::JibbenNeighborhood neighborhood;
  const int nlayers = 1;
  const int nstencil =
      (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);
  neighborhood.reserve(nstencil);
  neighborhood.setDelta(2.5 * mesh.dx());

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          // Fill neighborhood with polygons
          neighborhood.emptyNeighborhood();
          int count = 0;
          for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
            for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
              for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
                if (polygon(ii, jj, kk).getNumberOfVertices() > 0) {
                  neighborhood.addMember(polygon(ii, jj, kk));
                  if (i == ii && j == jj && k == kk) {
                    neighborhood.setCenterOfStencil(count);
                  }
                  count++;
                }
              }
            }
          }
          neighborhood.localize();
          // Now perform actual the Jibben fit and obtain interface
          (*a_interface)(i, j, k) =
              IRL::reconstructionWithJibben3D(neighborhood);

          // Match to volume fraction
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::setDistanceToMatchVolumeFraction(
              cell, liquid_volume_fraction, &(*a_interface)(i, j, k), 1.0e-14);
        }
      }
    }
  }

  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

const std::pair<double, Eigen::Vector3d> getPUAndGrad(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::SeparatorVariant>& a_interface,
    const Data<IRL::Pt>& a_centroid, const Data<double>& a_area,
    const int a_nlayers, const double a_delta, const int a_i, const int a_j,
    const int a_k, const IRL::Pt& a_pt) {
  const BasicMesh& mesh = a_liq_moments.getMesh();

  const double inv_delta = 1. / a_delta;
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d gradwxF_plus_wxgradF_sum = Eigen::Vector3d::Zero();
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();

  for (int k = a_k - a_nlayers; k <= a_k + a_nlayers; ++k) {
    for (int j = a_j - a_nlayers; j <= a_j + a_nlayers; ++j) {
      for (int i = a_i - a_nlayers; i <= a_i + a_nlayers; ++i) {
        const double vfrac =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        const double area_weight = a_area(i, j, k) / (mesh.dx() * mesh.dx());

        // Compute volume fraction weight
        double vfrac_weight = 1.0;
        const double limit_vfrac = 0.1;
        if (vfrac < 0.1) {
          vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * M_PI * vfrac);
        } else if (vfrac > 0.9) {
          vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * M_PI * (1.0 - vfrac));
        }

        // Compute distance weight
        const IRL::Pt x = a_pt - a_centroid(i, j, k);
        const double r = IRL::magnitude(x);
        const double rhat = r * inv_delta;
        if (rhat < 1.0) {  // Then weight is non zero;
          const double weight = area_weight * vfrac_weight * (1. + 4. * rhat) *
                                (1. - rhat) * (1. - rhat) * (1. - rhat) *
                                (1. - rhat);
          const Eigen::Vector3d grad_weight =
              Eigen::Vector3d(x[0], x[1], x[2]) * area_weight * vfrac_weight *
              20. * (rhat - 1.) * (rhat - 1.) * (rhat - 1.) * inv_delta *
              inv_delta;
          weight_sum += weight;
          grad_weight_sum += grad_weight;

          // Compute PU function
          if (const IRL::PlanarSeparator* separator =
                  std::get_if<IRL::PlanarSeparator>(
                      &(a_interface(i, j, k)))) {  // If plane
            if (separator->getNumberOfPlanes() > 0) {
              const IRL::Plane plane = (*separator)[0];
              // Compute plane normal
              const IRL::Normal n = plane.normal();
              const double F = n * x;
              const Eigen::Vector3d gradF = Eigen::Vector3d(n[0], n[1], n[2]);
              F_sum += weight * F;
              gradwxF_plus_wxgradF_sum += grad_weight * F + weight * gradF;
            }
          } else if (const IRL::Paraboloid* paraboloid =
                         std::get_if<IRL::Paraboloid>(
                             &(a_interface(i, j, k)))) {  // If paraboloid
            // Get paraboloid properties
            const IRL::ReferenceFrame& frame = paraboloid->getReferenceFrame();
            const double a = paraboloid->getAlignedParaboloid().a();
            const double b = paraboloid->getAlignedParaboloid().b();
            // Move sample point to local frame
            const IRL::Pt tmp = a_pt - paraboloid->getDatum();
            IRL::Pt xloc;
            for (int d = 0; d < 3; ++d) {
              xloc[d] = frame[d] * tmp;
            }
            const Eigen::Vector3d e0(frame[0][0], frame[0][1], frame[0][2]);
            const Eigen::Vector3d e1(frame[1][0], frame[1][1], frame[1][2]);
            const Eigen::Vector3d e2(frame[2][0], frame[2][1], frame[2][2]);
            const double dist_norm =
                1. / std::sqrt(1. + 4. * a * a * xloc[0] * xloc[0] +
                               4. * b * b * xloc[1] * xloc[1]);
            const Eigen::Vector3d grad_dist_norm =
                -4. * dist_norm * dist_norm * dist_norm *
                (a * a * e0 * xloc[0] + b * b * e1 * xloc[1]);
            const double F_alg =
                xloc[2] + a * xloc[0] * xloc[0] + b * xloc[1] * xloc[1];
            const Eigen::Vector3d grad_F_alg =
                e2 + 2. * (a * e0 * xloc[0] + b * e1 * xloc[1]);
            const double F = F_alg * dist_norm;
            const Eigen::Vector3d gradF =
                F_alg * grad_dist_norm + grad_F_alg * dist_norm;
            F_sum += weight * F;
            gradwxF_plus_wxgradF_sum += grad_weight * F + weight * gradF;
          }
        }
      }
    }
  }
  const double inv_weight_sum = 1.0 / weight_sum;
  const double PU_F = F_sum * inv_weight_sum;
  const Eigen::Vector3d PU_gradF =
      (gradwxF_plus_wxgradF_sum - F_sum * grad_weight_sum * inv_weight_sum) *
      inv_weight_sum;
  return std::make_pair(PU_F, PU_gradF);
}

const std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>
getPUAndGradAndHessian(const Data<IRL::VolumeMoments>& a_liq_moments,
                       const Data<IRL::SeparatorVariant>& a_interface,
                       const Data<IRL::Pt>& a_centroid,
                       const Data<double>& a_area, const int a_nlayers,
                       const double a_delta, const int a_i, const int a_j,
                       const int a_k, const IRL::Pt& a_pt) {
  const BasicMesh& mesh = a_liq_moments.getMesh();

  const double inv_delta = 1. / a_delta;
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hess_weight_sum = Eigen::Matrix3d::Zero();
  Eigen::Vector3d grad_product_sum = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hess_product_sum = Eigen::Matrix3d::Zero();

  for (int k = a_k - a_nlayers; k <= a_k + a_nlayers; ++k) {
    for (int j = a_j - a_nlayers; j <= a_j + a_nlayers; ++j) {
      for (int i = a_i - a_nlayers; i <= a_i + a_nlayers; ++i) {
        const double vfrac =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        const double area_weight = a_area(i, j, k) / (mesh.dx() * mesh.dx());

        // Compute volume fraction weight
        double vfrac_weight = 1.0;
        const double limit_vfrac = 0.1;
        if (vfrac < 0.1) {
          vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * M_PI * vfrac);
        } else if (vfrac > 0.9) {
          vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * M_PI * (1.0 - vfrac));
        }

        // Compute distance weight
        const IRL::Pt x = a_pt - a_centroid(i, j, k);
        const double r = IRL::magnitude(x);
        const double rhat = r * inv_delta;
        if (rhat < 1.0) {  // Then weight is non zero;
          const double weight = area_weight * vfrac_weight * (1. + 4. * rhat) *
                                (1. - rhat) * (1. - rhat) * (1. - rhat) *
                                (1. - rhat);
          const Eigen::Vector3d x_vector = Eigen::Vector3d(x[0], x[1], x[2]);
          const Eigen::Vector3d grad_weight =
              x_vector * area_weight * vfrac_weight * 20. * (rhat - 1.) *
              (rhat - 1.) * (rhat - 1.) * inv_delta * inv_delta;
          const double hess_factor = area_weight * vfrac_weight * 20. *
                                     (rhat - 1.) * (rhat - 1.) * inv_delta *
                                     inv_delta * inv_delta * inv_delta;
          Eigen::Matrix3d hess_weight = hess_factor * 3. * x_vector *
                                        x_vector.transpose() *
                                        (rhat > DBL_EPSILON ? 1.0 / rhat : 1.0);
          hess_weight += hess_factor * (a_delta * a_delta * (rhat - 1.)) *
                         Eigen::Matrix3d::Identity();

          // Increment weight sums
          weight_sum += weight;
          grad_weight_sum += grad_weight;
          hess_weight_sum += hess_weight;

          // Compute PU function
          if (const IRL::PlanarSeparator* separator =
                  std::get_if<IRL::PlanarSeparator>(
                      &(a_interface(i, j, k)))) {  // If plane
            if (separator->getNumberOfPlanes() > 0) {
              const IRL::Plane plane = (*separator)[0];
              // Compute plane normal
              const IRL::Normal n = plane.normal();
              const double F = n * x;
              const Eigen::Vector3d gradF = Eigen::Vector3d(n[0], n[1], n[2]);
              F_sum += weight * F;
              grad_product_sum += grad_weight * F + weight * gradF;
              hess_product_sum += hess_weight * F +
                                  grad_weight * gradF.transpose() +
                                  gradF * grad_weight.transpose();
            }
          } else if (const IRL::Paraboloid* paraboloid =
                         std::get_if<IRL::Paraboloid>(
                             &(a_interface(i, j, k)))) {  // If paraboloid
            // Get paraboloid properties
            const IRL::ReferenceFrame frame = paraboloid->getReferenceFrame();
            const double a = paraboloid->getAlignedParaboloid().a();
            const double b = paraboloid->getAlignedParaboloid().b();
            // Move sample point to local frame
            const IRL::Pt tmp = a_pt - paraboloid->getDatum();
            IRL::Pt xloc;
            for (int d = 0; d < 3; ++d) {
              xloc[d] = frame[d] * tmp;
            }
            const Eigen::Vector3d e0(frame[0][0], frame[0][1], frame[0][2]);
            const Eigen::Vector3d e1(frame[1][0], frame[1][1], frame[1][2]);
            const Eigen::Vector3d e2(frame[2][0], frame[2][1], frame[2][2]);
            // Compute Taubin distance norm term and grad and hessian
            const double dist_norm =
                1. / std::sqrt(1. + 4. * a * a * xloc[0] * xloc[0] +
                               4. * b * b * xloc[1] * xloc[1]);
            const Eigen::Vector3d grad_dist_norm =
                -4. * dist_norm * dist_norm * dist_norm *
                (a * a * e0 * xloc[0] + b * b * e1 * xloc[1]);
            const Eigen::Matrix3d hess_dist_norm =
                4. * dist_norm * dist_norm * dist_norm * dist_norm * dist_norm *
                (a * a *
                     (8. * a * a * xloc[0] * xloc[0] -
                      4. * b * b * xloc[1] * xloc[1] - 1.) *
                     e0 * e0.transpose() +
                 b * b *
                     (8. * b * b * xloc[1] * xloc[1] -
                      4. * a * a * xloc[0] * xloc[0] - 1.) *
                     e1 * e1.transpose() +
                 12. * a * a * b * b * xloc[0] * xloc[1] *
                     (e1 * e0.transpose() + e0 * e1.transpose()));
            // Compute algebraic distance term and grad and hessian
            const double F_alg =
                xloc[2] + a * xloc[0] * xloc[0] + b * xloc[1] * xloc[1];
            const Eigen::Vector3d grad_F_alg =
                e2 + 2. * (a * e0 * xloc[0] + b * e1 * xloc[1]);
            const Eigen::Matrix3d hess_F_alg =
                2. * (a * e0 * e0.transpose() + b * e1 * e1.transpose());
            // Compute signed distance and grad and hessian
            const double F = F_alg * dist_norm;
            const Eigen::Vector3d gradF =
                F_alg * grad_dist_norm + grad_F_alg * dist_norm;
            const Eigen::Matrix3d hessF =
                F_alg * hess_dist_norm +
                grad_F_alg * grad_dist_norm.transpose() +
                grad_dist_norm * grad_F_alg.transpose() +
                dist_norm * hess_F_alg;
            // Add to sums
            F_sum += weight * F;
            grad_product_sum += grad_weight * F + weight * gradF;
            hess_product_sum +=
                hess_weight * F + grad_weight * gradF.transpose() +
                gradF * grad_weight.transpose() + weight * hessF;
          }
        }
      }
    }
  }
  const double inv_weight_sum = 1.0 / weight_sum;
  const double PU_F = F_sum * inv_weight_sum;
  const Eigen::Vector3d PU_gradF =
      (grad_product_sum - F_sum * grad_weight_sum * inv_weight_sum) *
      inv_weight_sum;
  const Eigen::Matrix3d PU_hessF =
      (hess_product_sum + (grad_product_sum * grad_weight_sum.transpose() -
                           grad_weight_sum * grad_product_sum.transpose() -
                           F_sum * hess_weight_sum) *
                              inv_weight_sum) *
          inv_weight_sum -
      2. * (grad_product_sum - F_sum * grad_weight_sum * inv_weight_sum) *
          grad_weight_sum.transpose() * inv_weight_sum * inv_weight_sum;
  return std::make_tuple(PU_F, PU_gradF, PU_hessF);
}

const IRL::Pt projectCentroidOnPU(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::SeparatorVariant>& a_interface,
    const Data<IRL::Pt>& a_centroid, const Data<double>& a_area,
    const int a_nlayers, const double a_delta, const int a_i, const int a_j,
    const int a_k) {
  IRL::Pt projected_pt = a_centroid(a_i, a_j, a_k);
  const int itmax = 5;
  for (int i = 0; i < itmax; i++) {
    const auto F_and_gradF =
        getPUAndGrad(a_liq_moments, a_interface, a_centroid, a_area, a_nlayers,
                     a_delta, a_i, a_j, a_k, projected_pt);
    const double F = std::get<double>(F_and_gradF);
    if (F < a_delta * 1.e-6) {
      break;
    }
    const Eigen::Vector3d gradF = std::get<Eigen::Vector3d>(F_and_gradF);
    const double grad_norm_inv = 1.0 / gradF.squaredNorm();
    for (int d = 0; d < 3; d++) {
      projected_pt[d] -= F * gradF(d) * grad_norm_inv;
    }
  }
  return projected_pt;
}

void PU::getReconstruction(const Data<IRL::VolumeMoments>& a_liq_moments,
                           const Data<IRL::VolumeMoments>& a_gas_moments,
                           const double a_dt, const Data<double>& a_U,
                           const Data<double>& a_V, const Data<double>& a_W,
                           Data<IRL::SeparatorVariant>* a_interface,
                           const bool a_plic_already_built) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface);
  }

  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::SeparatorVariant> plic_reconstruction(&mesh);
  Data<IRL::SeparatorVariant> jibben_reconstruction(&mesh);
  // A element-wise copy is needed since std::memcpy is not safe with variants
  for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
        plic_reconstruction(i, j, k) = (*a_interface)(i, j, k);
        jibben_reconstruction(i, j, k) = (*a_interface)(i, j, k);
      }
    }
  }

  // Compute Jibben fit and plic centroids
  Data<IRL::Pt> interface_centroids(&mesh);
  Data<double> interface_areas(&mesh), jibben_errors(&mesh);
  Jibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                            &jibben_reconstruction, true, &interface_centroids,
                            &interface_areas, &jibben_errors);

  // Cleanup jibben
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          if (IRL::Paraboloid* paraboloid = std::get_if<IRL::Paraboloid>(
                  &jibben_reconstruction(i, j, k))) {
            const auto& aligned_paraboloid = paraboloid->getAlignedParaboloid();
            if (std::fabs(aligned_paraboloid.a()) > 1.0 / mesh.dx() ||
                std::fabs(aligned_paraboloid.b()) > 1.0 / mesh.dx()) {
              jibben_reconstruction(i, j, k) = plic_reconstruction(i, j, k);
            }
          }
        }
      }
    }
  }

  const int nlayers = 1;
  const double delta = 5.0 * mesh.dx();
  // const double jibben_error_threshold = 5.0e-3;

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          // Project local interface centroid on PU approximation
          const IRL::Pt pt_on_PU = projectCentroidOnPU(
              a_liq_moments, jibben_reconstruction, interface_centroids,
              interface_areas, nlayers, delta, i, j, k);

          // Compute local gradient and hessian of PU approximation
          const auto F_and_gradF_and_hessF = getPUAndGradAndHessian(
              a_liq_moments, jibben_reconstruction, interface_centroids,
              interface_areas, nlayers, delta, i, j, k, pt_on_PU);
          const Eigen::Vector3d gradF =
              std::get<Eigen::Vector3d>(F_and_gradF_and_hessF);
          const Eigen::Matrix3d hessF =
              std::get<Eigen::Matrix3d>(F_and_gradF_and_hessF);
          auto new_normal = IRL::Normal(gradF(0), gradF(1), gradF(2));
          new_normal.normalize();

          // Generate paraboloid and match volume fraction
          if (IRL::magnitude(new_normal) > 0.9) {
            // Generate paraboloid from gradient and hessian of implicit
            // surface
            auto paraboloid = generateLocalParaboloid(pt_on_PU, gradF, hessF);
            const double A = paraboloid.getAlignedParaboloid().a();
            const double B = paraboloid.getAlignedParaboloid().b();
            if (std::fabs(A) * mesh.dx() > 4.0 ||
                std::fabs(B) * mesh.dx() > 4.0) {
              continue;
            }
            // Translate paraboloid to match volume fraction
            const auto new_datum = paraboloid.getDatum();
            const auto new_frame = paraboloid.getReferenceFrame();
            const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
            const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                        mesh.z(k + 1));
            const auto cell = IRL::RectangularCuboid::fromBoundingPts(
                lower_cell_pt, upper_cell_pt);
            IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
                solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                                paraboloid);
            paraboloid.setDatum(IRL::Pt(
                new_datum + solver_distance.getDistance() * new_frame[2]));
            (*a_interface)(i, j, k) = paraboloid;
          }
        }
      }
    }
  }

  // Cleanup
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          double volume_supercell = 0.0;
          for (int kk = k - 1; kk <= k + 1; ++kk) {
            for (int jj = j - 1; jj <= j + 1; ++jj) {
              for (int ii = i - 1; ii <= i + 1; ++ii) {
                volume_supercell +=
                    a_liq_moments(ii, jj, kk).volume() / mesh.cell_volume();
              }
            }
          }
          if (liquid_volume_fraction < 1.0e-2 && volume_supercell < 1.0e-1) {
            (*a_interface)(i, j, k) = plic_reconstruction(i, j, k);
          }
        }
      }
    }
  }
}

void MixedPLICJibben::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface) {
  // First, we need to build the plic and copy it
  LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                           a_interface);
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::SeparatorVariant> plic_reconstruction(&mesh);
  // A element-wise copy is needed since std::memcpy is not safe with variants
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        plic_reconstruction(i, j, k) = (*a_interface)(i, j, k);
      }
    }
  }

  // Second, we build the Jibben reconstruction
  Jibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                            a_interface, true);

  // Choose between PLIC and Jibben
  const double vfrac_threshold = 1.0e-4;
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction < vfrac_threshold ||
            liquid_volume_fraction > 1.0 - vfrac_threshold) {
          (*a_interface)(i, j, k) = plic_reconstruction(i, j, k);
        }
      }
    }
  }

  // Update border with simple ghost-cell fill and correct distances for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

void SlicesParabola::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface, const bool a_plic_already_built) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface);
  }

  // Then, let's compute the PLIC polygons
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::Polygon> polygon(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
            liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        // Note: this std::get call is safe since all interfaces should be
        // of the type PlanarSeparator after calling the PLIC reconstruction
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, planar_separator, planar_separator[0]);
      }
    }
  }
  updatePolygonBorder(&polygon);

  // Now let's fit parabolas based on the PLIC polygons
  const int nsamples_per_segment = 10;
  const int nslices = 10;
  const int nlayers = 1;
  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);
  IRL::StackVector<IRL::Pt, 2> intersections;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (polygon(i, j, k).getNumberOfVertices() > 2) {
          // Compute local frame of reference based on polygon
          const IRL::Normal polygon_normal =
              calculatePolygonNormal(polygon(i, j, k));
          const IRL::ReferenceFrame polygon_frame =
              referenceFrameFromNormal(polygon_normal);
          const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();
          // Build list of polygons in 5x5x5 stencil
          polygon_vfrac_list.resize(0);
          for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
            for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
              for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
                if (polygon(ii, jj, kk).getNumberOfVertices() > 2) {
                  const double vfrac =
                      a_liq_moments(ii, jj, kk).volume() / mesh.cell_volume();
                  polygon_vfrac_list.push_back(
                      std::make_pair(polygon(ii, jj, kk), vfrac));
                }
              }
            }
          }
          // Initialize approximate tensor of curvature
          Eigen::Matrix3d Mtilde = Eigen::Matrix3d::Zero();
          // Loop over slices and compute local directional curvature
          // contribution to approximate tensor of curvature
          for (int s = 0; s < nslices; s++) {
            // Compute slice weight
            const double weight_n = 1.0 / static_cast<double>(nslices);
            // Compute rotation angle
            const double theta_n = M_PI * weight_n * static_cast<double>(s);
            // Compute rotated tangent
            const IRL::UnitQuaternion rotation(theta_n, polygon_frame[2]);
            const auto local_frame = rotation * polygon_frame;
            // Create plane with local normal and rotated tangent
            const auto slicing_plane =
                IRL::Plane(local_frame[1], local_frame[1] * polygon_centroid);
            // Compute intersections of plane with PLIC polygons and the
            // resulting directional curvature
            Eigen::Matrix2d Q = Eigen::Matrix2d::Zero();
            Eigen::Vector2d r = Eigen::Vector2d::Zero();
            for (int p = 0; p < polygon_vfrac_list.size(); p++) {
              getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                                 &intersections);
              // If PLIC does intersection with plane
              if (intersections.size() == 2) {
                const IRL::Pt pt0 = intersections[0] - polygon_centroid;
                const IRL::Pt pt1 = intersections[1] - polygon_centroid;
                // Compute weight
                const double distance =
                    IRL::magnitude(IRL::Pt(0.5 * (pt0 + pt1)));
                const double distance_ndim = distance / 2.5 * mesh.dx();
                const double distance_weight =
                    distance_ndim >= 1.0
                        ? 0.0
                        : (1.0 + 4.0 * distance_ndim) *
                              std::pow(1.0 - distance_ndim, 4.0);
                const double vfrac = polygon_vfrac_list[p].second;
                double vfrac_weight = 1.0;
                const double limit_vfrac = 0.1;
                if (vfrac < 0.1) {
                  vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * M_PI * vfrac);
                } else if (vfrac > 0.9) {
                  vfrac_weight =
                      0.5 - 0.5 * std::cos(10.0 * M_PI * (1.0 - vfrac));
                }
                const double weight = vfrac_weight * distance_weight /
                                      static_cast<double>(nsamples_per_segment);
                // Loop over segment samples
                const double nsample_norm =
                    1.0 / static_cast<double>(nsamples_per_segment - 1);
                for (int pp = 0; pp < nsamples_per_segment; pp++) {
                  const IRL::Pt pt = pt0 + (pt1 - pt0) *
                                               static_cast<double>(pp) *
                                               nsample_norm;
                  const double x = pt * local_frame[0];
                  const double y = pt * local_frame[2];
                  const double x2 = x * x;
                  Q(0, 0) += weight * x2 * x2;
                  Q(0, 1) += weight * x2;
                  Q(1, 0) += weight * x2;
                  Q(1, 1) += weight;
                  r(0) += weight * x2 * y;
                  r(1) += weight * y;
                }
              }
            }
            const auto coeffs = Q.fullPivHouseholderQr().solve(r);
            const double curvature_n = -2.0 * coeffs(0);
            // Add contribution to tensor of curvature
            const Eigen::Vector3d Ttheta_n(local_frame[0][0], local_frame[0][1],
                                           local_frame[0][2]);
            Mtilde += weight_n * curvature_n * Ttheta_n * Ttheta_n.transpose();
          }
          // Compute Householder matrix
          IRL::Normal W_plus = IRL::Normal(1, 0, 0) + polygon_frame[2];
          IRL::Normal W_minus = IRL::Normal(1, 0, 0) - polygon_frame[2];
          IRL::Normal W_max =
              (IRL::squaredMagnitude(W_plus) > IRL::squaredMagnitude(W_plus))
                  ? W_plus
                  : W_minus;
          W_max.normalize();
          const Eigen::Vector3d W(W_max[0], W_max[1], W_max[2]);
          const Eigen::Matrix3d Q =
              Eigen::Matrix3d::Identity() - 2.0 * W * W.transpose();
          // Compute minor matrix
          const Eigen::Matrix3d QMQ = Q.transpose() * Mtilde * Q;
          const double m11 = QMQ(1, 1);
          const double m12 = 0.5 * (QMQ(1, 2) + QMQ(2, 1));
          const double m22 = QMQ(2, 2);
          // Compute eigenvalues and Givens rotation angle
          const double tmp_sqrt =
              std::sqrt((m11 - m22) * (m11 - m22) + 4.0 * m12 * m12);
          const double lambda1 = 0.5 * (m11 + m22 + tmp_sqrt);
          const double lambda2 = 0.5 * (m11 + m22 - tmp_sqrt);
          const double theta = 0.5 * std::atan2(2.0 * m12, m11 - m22);
          // Compute principal directions (Darboux frame)
          const IRL::UnitQuaternion givens_rotation(theta, polygon_frame[2]);
          const auto darboux_frame = givens_rotation * polygon_frame;
          // Compute principal curvatures
          const double k1 = 3.0 * lambda1 - lambda2;
          const double k2 = 3.0 * lambda2 - lambda1;
          auto paraboloid = IRL::Paraboloid(polygon_centroid, darboux_frame,
                                            0.5 * k1, 0.5 * k2);
          // Translate paraboloid to match volume fraction
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                              paraboloid);
          paraboloid.setDatum(
              IRL::Pt(polygon_centroid +
                      solver_distance.getDistance() * darboux_frame[2]));
          (*a_interface)(i, j, k) = paraboloid;
        }
      }
    }
  }

  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

void recenterMoments(IRL::VolumeMoments* moments, const IRL::Pt& center) {
  for (int i = 0; i < 3; i++) {
    (*moments).centroid()[i] - (*moments).volume() * center[i];
  }
}

void correctInterfaceBorders(Data<IRL::SeparatorVariant>* a_interface) {
  const BasicMesh& mesh = (*a_interface).getMesh();
  // Fix distances in reconstruction for periodic boundary
  Data<IRL::Pt> shift(&mesh);
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        shift(i, j, k) = IRL::Pt(0.0, 0.0, 0.0);
      }
    }
  }

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        shift(i, j, k)[0] -= mesh.lx();
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        shift(i, j, k)[0] += mesh.lx();
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        shift(i, j, k)[1] -= mesh.ly();
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        shift(i, j, k)[1] += mesh.ly();
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        shift(i, j, k)[2] -= mesh.lz();
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        shift(i, j, k)[2] += mesh.lz();
      }
    }
  }

  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        if (auto interface_ptr =
                std::get_if<IRL::PlanarSeparator>(&(*a_interface)(i, j, k))) {
          for (auto& plane : *interface_ptr) {
            plane.distance() += plane.normal() * shift(i, j, k);
          }
        } else if (auto interface_ptr =
                       std::get_if<IRL::Paraboloid>(&(*a_interface)(i, j, k))) {
          const IRL::Pt datum = (*interface_ptr).getDatum() + shift(i, j, k);
          (*interface_ptr).setDatum(datum);
        } else if (auto interface_ptr =
                       std::get_if<IRL::Cylinder>(&(*a_interface)(i, j, k))) {
          const IRL::Pt datum = (*interface_ptr).getDatum() + shift(i, j, k);
          (*interface_ptr).setDatum(datum);
        }
      }
    }
  }
}
