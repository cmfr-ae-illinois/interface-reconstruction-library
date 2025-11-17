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

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "examples/implicit_surface_reconstruction/data.h"
#include "examples/implicit_surface_reconstruction/reconstruction_types.h"
// #include "examples/variant_advector/vof_advection.h"

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
  } else if (a_reconstruction_method == "Taubin") {
    Taubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface);
  } else if (a_reconstruction_method == "SlicesTaubin") {
    SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                    a_V, a_W, a_interface);
  } else if (a_reconstruction_method == "PLICAligned") {
    PLICAligned::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                   a_W, a_interface);
  } else if (a_reconstruction_method == "SlicesParticle") {
    SlicesParticle::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                      a_V, a_W, a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: PLIC, Jibben, PU, MixedJibben, "
                 "SlicesParabola, Taubin. \n";
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
            // Generate paraboloid from gradient and hessian of surface
            auto paraboloid =
                IRL::Paraboloid::fromDerivatives(pt_on_PU, gradF, hessF);
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
                const double distance_ndim = distance / (2.5 * mesh.dx());
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
              (IRL::squaredMagnitude(W_plus) > IRL::squaredMagnitude(W_minus))
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

void Taubin::getReconstruction(const Data<IRL::VolumeMoments>& a_liq_moments,
                               const Data<IRL::VolumeMoments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V, const Data<double>& a_W,
                               Data<IRL::SeparatorVariant>* a_interface,
                               const bool a_plic_already_built) {
  using namespace IRL;

  if (!a_plic_already_built) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface);
  }

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

        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));

        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, planar_separator, planar_separator[0]);
      }
    }
  }
  updatePolygonBorder(&polygon);

  const int nsamples_per_segment = 10;
  const int nslices = 101;
  const int nlayers = 2;

  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();

        if (polygon(i, j, k).getNumberOfVertices() <= 2) continue;

        // Local frame from PLIC polygon
        const IRL::Normal polygon_normal =
            calculatePolygonNormal(polygon(i, j, k));
        IRL::ReferenceFrame polygon_frame =
            referenceFrameFromNormal(polygon_normal);
        const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();

        // Gather stencil polygons and their volume fractions
        polygon_vfrac_list.clear();
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              if (polygon(ii, jj, kk).getNumberOfVertices() > 2) {
                const double vfrac_n =
                    a_liq_moments(ii, jj, kk).volume() / mesh.cell_volume();
                polygon_vfrac_list.emplace_back(polygon(ii, jj, kk), vfrac_n);
              }
            }
          }
        }

        //  accumulators for least squares
        Eigen::Matrix3d XtX = Eigen::Matrix3d::Zero();
        Eigen::Vector3d Xty = Eigen::Vector3d::Zero();

        // Sweep angles over [0, pi)
        for (int s = 0; s < nslices; ++s) {
          const double w_s =
              static_cast<double>(s) / static_cast<double>(nslices);
          const double theta_s = -M_PI + (2.0 * M_PI) * w_s;

          const IRL::UnitQuaternion rotation(theta_s, polygon_frame[2]);
          const IRL::ReferenceFrame local_frame = rotation * polygon_frame;

          const IRL::Plane slicing_plane(local_frame[1],
                                         local_frame[1] * polygon_centroid);

          // ---- Weighted linear LS circle fit on the slice ----
          // Model: x^2 + y^2 + A x + B y + C = 0
          // Normal eqs H * [A,B,C]^T = rhs
          double Sx = 0, Sy = 0, Sw = 0;
          double Sxx = 0, Syy = 0, Sxy = 0;
          double Sz = 0, Sxz = 0, Syz = 0;

          for (int p = 0; p < static_cast<int>(polygon_vfrac_list.size());
               ++p) {
            IRL::StackVector<IRL::Pt, 2> intersections;
            getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                               &intersections);
            if (intersections.size() != 2) continue;

            // Segment endpoints relative to centroid
            const IRL::Pt pt0 = intersections[0] - polygon_centroid;
            const IRL::Pt pt1 = intersections[1] - polygon_centroid;

            const double distance = IRL::magnitude(IRL::Pt(0.5 * (pt0 + pt1)));
            const double distance_ndim = distance / (2.5 * mesh.dx());
            const double distance_weight =
                (distance_ndim >= 1.0) ? 0.0
                                       : (1.0 + 4.0 * distance_ndim) *
                                             std::pow(1.0 - distance_ndim, 4.0);

            const double vfrac = polygon_vfrac_list[p].second;
            double vfrac_weight = 1.0;
            if (vfrac < 0.1) {
              vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * M_PI * vfrac);
            } else if (vfrac > 0.9) {
              vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * M_PI * (1.0 - vfrac));
            }

            IRL::Normal nloc =
                calculatePolygonNormal(polygon_vfrac_list[p].first);
            double n_dot = polygon_normal[0] * nloc[0] +
                           polygon_normal[1] * nloc[1] +
                           polygon_normal[2] * nloc[2];
            double normal_weight = 1.;
            if (n_dot <= 0) normal_weight = 0;
            normal_weight = std::max(0.0, n_dot);

            const double seg_weight = vfrac_weight * distance_weight /
                                      static_cast<double>(nsamples_per_segment);

            const double nsample_norm =
                1.0 / static_cast<double>(nsamples_per_segment - 1);
            for (int pp = 0; pp < nsamples_per_segment; ++pp) {
              const IRL::Pt pt =
                  pt0 + (pt1 - pt0) * (static_cast<double>(pp) * nsample_norm);
              const double x = pt * local_frame[0];  // in-plane tangent
              const double y = pt * local_frame[2];  // height along normal
              const double z = x * x + y * y;

              const double w = seg_weight;

              Sx += w * x;
              Sy += w * y;
              Sw += w;

              Sxx += w * x * x;
              Syy += w * y * y;
              Sxy += w * x * y;

              Sz += w * z;
              Sxz += w * x * z;
              Syz += w * y * z;
            }
          }

          Eigen::Matrix3d H;
          H.setZero();
          H(0, 0) = Sxx;
          H(0, 1) = Sxy;
          H(0, 2) = Sx;
          H(1, 0) = Sxy;
          H(1, 1) = Syy;
          H(1, 2) = Sy;
          H(2, 0) = Sx;
          H(2, 1) = Sy;
          H(2, 2) = Sw;

          Eigen::Vector3d rhs;
          rhs << -Sxz, -Syz, -Sz;

          double k_theta = 0.0;
          bool slice_has_k = false;

          if (H.norm() > 0.0) {
            // Tiny ridge for robustness
            const double ridge = 1e-14;
            H(0, 0) += ridge;
            H(1, 1) += ridge;
            H(2, 2) += ridge;

            Eigen::LDLT<Eigen::Matrix3d> ldlt3(H);
            if (ldlt3.info() == Eigen::Success) {
              const Eigen::Vector3d abc = ldlt3.solve(rhs);
              const double A = abc(0), B = abc(1), C = abc(2);

              // Center and radius in local slice frame
              const double xc = -0.5 * A;
              const double yc = -0.5 * B;
              const double R2 = std::max(0.0, xc * xc + yc * yc - C);
              const double R = std::sqrt(R2);

              if (R > 0.0 && std::isfinite(R)) {
                const double sign = (yc >= 0.0) ? -1.0 : +1.0;
                k_theta = sign * (1.0 / R);
                slice_has_k = std::isfinite(k_theta);
              }
            }
          }

          if (!slice_has_k) continue;

          const double x0 = 1.0;
          const double x1 = std::cos(2.0 * theta_s);
          const double x2 = std::sin(2.0 * theta_s);
          const double w_row = 1.0;  // could use a slice confidence if desired

          XtX(0, 0) += w_row * x0 * x0;
          XtX(0, 1) += w_row * x0 * x1;
          XtX(0, 2) += w_row * x0 * x2;
          XtX(1, 0) += w_row * x1 * x0;
          XtX(1, 1) += w_row * x1 * x1;
          XtX(1, 2) += w_row * x1 * x2;
          XtX(2, 0) += w_row * x2 * x0;
          XtX(2, 1) += w_row * x2 * x1;
          XtX(2, 2) += w_row * x2 * x2;

          Xty(0) += w_row * x0 * k_theta;
          Xty(1) += w_row * x1 * k_theta;
          Xty(2) += w_row * x2 * k_theta;
        }  // end sweep over theta

        // Solve normal equations for [alpha, beta, gamma]
        bool fit_ok = false;
        double k1 = 0.0, k2 = 0.0, phi = 0.0;
        if (XtX.norm() > 0.0) {
          const double ridge = 1e-14;
          XtX(0, 0) += ridge;
          XtX(1, 1) += ridge;
          XtX(2, 2) += ridge;

          Eigen::LDLT<Eigen::Matrix3d> ldlt(XtX);
          if (ldlt.info() == Eigen::Success) {
            const Eigen::Vector3d abg = ldlt.solve(Xty);
            const double alpha = abg(0);
            const double beta = abg(1);
            const double gamma = abg(2);

            const double Rmag = std::sqrt(beta * beta + gamma * gamma);
            k1 = alpha + Rmag;  // max principal curvature
            k2 = alpha - Rmag;  // min principal curvature
            phi = 0.5 * std::atan2(gamma, beta);

            fit_ok =
                std::isfinite(k1) && std::isfinite(k2) && std::isfinite(phi);
          }
        }

        if (!fit_ok) {
          // Fallback: keep planar if fit fails
          std::cout << "Taubin circle fit failed!" << std::endl;
          (*a_interface)(i, j, k) =
              std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
          continue;
        }

        // Rotate polygon_frame around its normal by phi to get Darboux frame
        const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
        const IRL::ReferenceFrame darboux_frame = rotate_phi * polygon_frame;

        // Build paraboloid with coefficients a=0.5*k1, b=0.5*k2 (curvature =
        // 2a, 2b)
        IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 * k1,
                                   0.5 * k2);

        // std::cout << "a = " << 0.5 * k1 << " , b = " << 0.5 * k2 <<
        // std::endl;

        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);

        IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
            solver_distance(cell, liquid_volume_fraction, 1.0e-14, paraboloid);

        paraboloid.setDatum(
            IRL::Pt(polygon_centroid +
                    solver_distance.getDistance() * darboux_frame[2]));

        (*a_interface)(i, j, k) = paraboloid;
      }
    }
  }

  // 5) Border handling
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// volume fraction weight
double getVfracWeight(const double& a_vfrac) {
  const double limit_vfrac = 0.1;
  if (a_vfrac < limit_vfrac) {
    return 0.5 - 0.5 * std::cos(M_PI * a_vfrac / limit_vfrac);
  } else if (a_vfrac > (1.0 - limit_vfrac)) {
    return 0.5 - 0.5 * std::cos(M_PI * (1.0 - a_vfrac) / limit_vfrac);
  } else {
    return 1.0;
  }
}

// distance weight
double getDistanceWeight(const IRL::Pt& a_pref, const IRL::Pt& a_ploc,
                         const double& h) {
  IRL::Pt d_vec = (a_ploc - a_pref);
  double distance = std::sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] +
                              d_vec[2] * d_vec[2]) /
                    h;
  // if (distance < 2.5) {
  //   return (1.0 + 4.0 * distance / 2.5) * std::pow(1.0 - distance
  //   / 2.5, 4.0);
  // } else {
  //   return 0.0;
  // }
  return 1.0;
}

// normal weight
double getNormalWeight(const IRL::Normal& a_nref, const IRL::Normal& a_nloc) {
  double n_dot =
      a_nref[0] * a_nloc[0] + a_nref[1] * a_nloc[1] + a_nref[2] * a_nloc[2];
  return std::max(0.0, n_dot);
}

struct TaubinCircleData {
  double k;                // signed curvbature
  double R;                // circle radius
  Eigen::Vector2d center;  // circle center
};

// taubin circle coefficient
TaubinCircleData getTaubinData(
    const std::vector<std::pair<IRL::Pt, IRL::Pt>>& a_end_points,
    const int& a_nsamples, const std::vector<double>& a_vfrac_list,
    const std::vector<IRL::Normal>& a_normal_list, const double& a_h,
    const IRL::Normal& a_target_normal, const IRL::Pt& a_target_centroid,
    const IRL::ReferenceFrame& a_local_frame, const IRL::Pt& a_local_origin) {
  TaubinCircleData tcd;

  // rotation matrix
  Eigen::Vector3d e1(a_local_frame[0][0], a_local_frame[0][1],
                     a_local_frame[0][2]);
  Eigen::Vector3d e2(a_local_frame[2][0], a_local_frame[2][1],
                     a_local_frame[2][2]);
  Eigen::Vector3d e3(a_local_frame[1][0], a_local_frame[1][1],
                     a_local_frame[1][2]);
  Eigen::Matrix3d R;
  R.col(0) = e1;
  R.col(1) = e2;
  R.col(2) = e3;
  Eigen::Vector3d o(a_local_origin[0], a_local_origin[1], a_local_origin[2]);

  // sampling points in local frame
  std::vector<Eigen::Vector3d> points;
  std::vector<double> vfw, dw, nw;  // weights
  for (int i = 0; i < a_end_points.size(); i++) {
    Eigen::Vector3d x0(a_end_points[i].first[0], a_end_points[i].first[1],
                       a_end_points[i].first[2]);
    Eigen::Vector3d x1(a_end_points[i].second[0], a_end_points[i].second[1],
                       a_end_points[i].second[2]);
    // end points in local frame
    Eigen::Vector3d x0_local = R.transpose() * (x0 - o);
    Eigen::Vector3d x1_local = R.transpose() * (x1 - o);

    // sampling and finding weights
    for (int j = 0; j < a_nsamples; j++) {
      double t = static_cast<double>(j) / (static_cast<double>(a_nsamples) - 1);
      Eigen::Vector3d pt = x0_local * (1.0 - t) + x1_local * t;
      IRL::Pt p(pt[0], pt[1], pt[2]);
      points.push_back(pt);
      double pt_vfw = getVfracWeight(a_vfrac_list[i]);
      double pt_dw = getDistanceWeight(a_target_centroid, p, a_h);
      double pt_nw = getNormalWeight(a_target_normal, a_normal_list[i]);
      vfw.push_back(pt_vfw);
      dw.push_back(pt_dw);
      nw.push_back(pt_nw);
    }
  }

  // moment matrix
  int n = points.size();
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    double w = vfw[i] * dw[i] * nw[i];  // applying weights
    w = nw[i];
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();  // weighted outer product
  }

  // constraint matrix
  Eigen::Matrix4d C;
  C.setZero();
  C(0, 0) = 4.0 * M(0, 3);
  C(0, 1) = 2.0 * M(1, 3);
  C(0, 2) = 2.0 * M(2, 3);
  C(1, 0) = C(0, 1);
  C(1, 1) = n;
  C(2, 0) = C(0, 2);
  C(2, 2) = n;

  // solving the generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // extracting smallest positive eigenvalue and its eigenvector
  std::vector<std::pair<double, int>> positive_eigs;
  for (int i = 0; i < eigenvalues.size(); i++) {
    double real_part = eigenvalues[i].real();
    if (real_part > 0) {
      positive_eigs.emplace_back(real_part, i);
    }
  }
  std::sort(positive_eigs.begin(), positive_eigs.end());
  double A = eigenvectors.col(positive_eigs[0].second)[0].real();
  double B = eigenvectors.col(positive_eigs[0].second)[1].real();
  double C_ = eigenvectors.col(positive_eigs[0].second)[2].real();
  double D = eigenvectors.col(positive_eigs[0].second)[3].real();

  // scaling parameters
  double constraint = 4.0 * A * A * M(0, 3) + 4.0 * A * B * M(1, 3) +
                      4.0 * A * C_ * M(2, 3) + B * B * static_cast<double>(n) +
                      C_ * C_ * static_cast<double>(n);
  double scale_factor = 1.0 / std::sqrt(constraint) * std::sqrt(80.0);
  A *= scale_factor;
  B *= scale_factor;
  C_ *= scale_factor;
  D *= scale_factor;

  // radius of circle fit
  double radius = std::sqrt((B * B + C_ * C_ - 4.0 * A * D) /
                            IRL::safelyEpsilon(4.0 * A * A));

  // y-coordinate of center
  double yc = -C_ / IRL::safelyEpsilon(2.0 * A);
  double xc = -B / IRL::safelyEpsilon(2.0 * A);
  double sign = (yc >= 0.0) ? -1.0 : +1.0;

  tcd.k = sign * (1.0 / IRL::safelyEpsilon(radius));
  tcd.R = radius;
  tcd.center = Eigen::Vector2d(xc, yc);

  return tcd;
}

// taubin circle coefficient
double getTaubinCurvature(
    const std::vector<std::pair<IRL::Pt, IRL::Pt>>& a_end_points,
    const int& a_nsamples, const std::vector<double>& a_vfrac_list,
    const std::vector<IRL::Normal>& a_normal_list, const double& a_h,
    const IRL::Normal& a_target_normal, const IRL::Pt& a_target_centroid,
    const IRL::ReferenceFrame& a_local_frame, const IRL::Pt& a_local_origin) {
  // rotation matrix
  Eigen::Vector3d e1(a_local_frame[0][0], a_local_frame[0][1],
                     a_local_frame[0][2]);
  Eigen::Vector3d e2(a_local_frame[2][0], a_local_frame[2][1],
                     a_local_frame[2][2]);
  Eigen::Vector3d e3(a_local_frame[1][0], a_local_frame[1][1],
                     a_local_frame[1][2]);
  Eigen::Matrix3d R;
  R.col(0) = e1;
  R.col(1) = e2;
  R.col(2) = e3;
  Eigen::Vector3d o(a_local_origin[0], a_local_origin[1], a_local_origin[2]);

  // sampling points in local frame
  std::vector<Eigen::Vector3d> points;
  std::vector<double> vfw, dw, nw;  // weights
  for (int i = 0; i < a_end_points.size(); i++) {
    Eigen::Vector3d x0(a_end_points[i].first[0], a_end_points[i].first[1],
                       a_end_points[i].first[2]);
    Eigen::Vector3d x1(a_end_points[i].second[0], a_end_points[i].second[1],
                       a_end_points[i].second[2]);
    // end points in local frame
    Eigen::Vector3d x0_local = R.transpose() * (x0 - o);
    Eigen::Vector3d x1_local = R.transpose() * (x1 - o);

    // sampling and finding weights
    for (int j = 0; j < a_nsamples; j++) {
      double t = static_cast<double>(j) / (static_cast<double>(a_nsamples) - 1);
      Eigen::Vector3d pt = x0_local * (1.0 - t) + x1_local * t;
      IRL::Pt p(pt[0], pt[1], pt[2]);
      points.push_back(pt);
      double pt_vfw = getVfracWeight(a_vfrac_list[i]);
      double pt_dw = getDistanceWeight(a_target_centroid, p, a_h);
      double pt_nw = getNormalWeight(a_target_normal, a_normal_list[i]);
      vfw.push_back(pt_vfw);
      dw.push_back(pt_dw);
      nw.push_back(pt_nw);
    }
  }

  // moment matrix
  int n = points.size();
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    double w = vfw[i] * dw[i] * nw[i];  // applying weights
    w = nw[i];
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();  // weighted outer product
  }

  // constraint matrix
  Eigen::Matrix4d C;
  C.setZero();
  C(0, 0) = 4.0 * M(0, 3);
  C(0, 1) = 2.0 * M(1, 3);
  C(0, 2) = 2.0 * M(2, 3);
  C(1, 0) = C(0, 1);
  C(1, 1) = n;
  C(2, 0) = C(0, 2);
  C(2, 2) = n;

  // solving the generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // extracting smallest positive eigenvalue and its eigenvector
  std::vector<std::pair<double, int>> positive_eigs;
  for (int i = 0; i < eigenvalues.size(); i++) {
    double real_part = eigenvalues[i].real();
    if (real_part > 0) {
      positive_eigs.emplace_back(real_part, i);
    }
  }
  std::sort(positive_eigs.begin(), positive_eigs.end());
  double A = eigenvectors.col(positive_eigs[0].second)[0].real();
  double B = eigenvectors.col(positive_eigs[0].second)[1].real();
  double C_ = eigenvectors.col(positive_eigs[0].second)[2].real();
  double D = eigenvectors.col(positive_eigs[0].second)[3].real();

  // scaling parameters
  double constraint = 4.0 * A * A * M(0, 3) + 4.0 * A * B * M(1, 3) +
                      4.0 * A * C_ * M(2, 3) + B * B * static_cast<double>(n) +
                      C_ * C_ * static_cast<double>(n);
  double scale_factor = 1.0 / std::sqrt(constraint) * std::sqrt(80.0);
  A *= scale_factor;
  B *= scale_factor;
  C_ *= scale_factor;
  D *= scale_factor;

  // radius of circle fit
  double radius = std::sqrt((B * B + C_ * C_ - 4.0 * A * D) / (4.0 * A * A));

  // y-coordinate of center
  double yc = -C_ / (2.0 * A);
  double sign = (yc >= 0.0) ? -1.0 : +1.0;

  return sign * (1.0 / radius);  // return signed curvature
}

// estimating reference frame
Eigen::Vector2d estimateTaubinNormal(const Eigen::Vector2d& circle_center,
                                     const double& R) {
  Eigen::Vector2d plic_centroid(0.,
                                0.);  // plic centroid is origin in local frame
  Eigen::Vector2d dir = plic_centroid - circle_center;
  dir.normalize();
  Eigen::Vector2d datum = circle_center + R * dir;
  Eigen::Vector2d normal = datum - circle_center;
  normal.normalize();
  return normal;
}

void SlicesTaubin::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface, const bool a_plic_already_built) {
  // plic
  if (!a_plic_already_built) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface);
  }

  // clipped polygon from planes
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::Polygon> polygon(&mesh);
  Data<double> vfrac(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        vfrac(i, j, k) = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vfrac(i, j, k) < IRL::global_constants::VF_LOW ||
            vfrac(i, j, k) > IRL::global_constants::VF_HIGH) {
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, planar_separator, planar_separator[0]);
      }
    }
  }
  updatePolygonBorder(&polygon);

  // slicing params
  const int nsamples_per_segment = 10;
  const int nslices = 18;
  const int nlayers = 1;

  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);
  IRL::StackVector<IRL::Pt, 2> intersections;

  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        if (polygon(i, j, k).getNumberOfVertices() > 2) {
          // target cell local frame
          const IRL::Normal polygon_normal =
              calculatePolygonNormal(polygon(i, j, k));
          IRL::ReferenceFrame polygon_frame =
              referenceFrameFromNormal(polygon_normal);
          const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();

          // building stencil
          polygon_vfrac_list.resize(0);
          for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
            for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
              for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
                if (polygon(ii, jj, kk).getNumberOfVertices() > 2) {
                  polygon_vfrac_list.push_back(
                      std::make_pair(polygon(ii, jj, kk), vfrac(ii, jj, kk)));
                }
              }
            }
          }

          // accumulator for principal curvature least squares
          Eigen::Matrix3d AtA = Eigen::Matrix3d::Zero();
          Eigen::Vector3d Atb = Eigen::Vector3d::Zero();

          // sweeping over slices from [0, pi)
          Eigen::Vector3d n_global_sum = Eigen::Vector3d::Zero();
          for (int s = 0; s < nslices; s++) {
            // plane rotation angle
            double theta_s =
                M_PI * static_cast<double>(s) / static_cast<double>(nslices);

            // rotating target polygon frame about normal
            const IRL::UnitQuaternion rotation(theta_s, polygon_frame[2]);
            const auto local_frame = rotation * polygon_frame;

            // slicing plane
            const IRL::Plane slicing_plane(local_frame[1],
                                           local_frame[1] * polygon_centroid);

            // intersection of slicing plane with polygons in stencil
            std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
            std::vector<double> vfrac_list;
            std::vector<IRL::Normal> normal_list;
            for (int p = 0; p < static_cast<int>(polygon_vfrac_list.size());
                 p++) {
              getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                                 &intersections);
              if (intersections.size() != 2) continue;
              IRL::Pt start_point = intersections[0];
              IRL::Pt end_point = intersections[1];
              IRL::Normal neighbor_normal =
                  calculatePolygonNormal(polygon_vfrac_list[p].first);
              end_points_list.push_back({start_point, end_point});
              vfrac_list.push_back(polygon_vfrac_list[p].second);
              normal_list.push_back(neighbor_normal);
            }
            // taubin circle data
            TaubinCircleData tcd =
                getTaubinData(end_points_list, nsamples_per_segment, vfrac_list,
                              normal_list, mesh.dx(), polygon_normal,
                              polygon_centroid, local_frame, polygon_centroid);

            // local normal
            Eigen::Vector2d n_taubin = estimateTaubinNormal(tcd.center, tcd.R);
            Eigen::Vector3d n_loc(n_taubin[0], n_taubin[1], 0.0);
            // if (i == 6 && j == 8 && k == 10) {
            //   std::cout << "-------Slice " << s << std::endl;
            //   std::cout << "Taubin Normal = " << n_taubin.transpose()
            //             << std::endl;
            // }

            // convertint to global frame
            Eigen::Vector3d e1(local_frame[0][0], local_frame[0][1],
                               local_frame[0][2]);
            Eigen::Vector3d e2(local_frame[2][0], local_frame[2][1],
                               local_frame[2][2]);
            Eigen::Vector3d e3(local_frame[1][0], local_frame[1][1],
                               local_frame[1][2]);
            Eigen::Matrix3d R;
            R.col(0) = e1;
            R.col(1) = e2;
            R.col(2) = e3;
            Eigen::Vector3d n_global = R * n_loc;
            n_global_sum += n_global;

            // directional curvature
            double k_theta = tcd.k;
            const double x0 = 1.0;
            const double x1 = std::cos(2.0 * theta_s);
            const double x2 = std::sin(2.0 * theta_s);
            const double w_row = 1.0;  // slicing weight (if any)

            AtA(0, 0) += w_row * x0 * x0;
            AtA(0, 1) += w_row * x0 * x1;
            AtA(0, 2) += w_row * x0 * x2;
            AtA(1, 0) += w_row * x1 * x0;
            AtA(1, 1) += w_row * x1 * x1;
            AtA(1, 2) += w_row * x1 * x2;
            AtA(2, 0) += w_row * x2 * x0;
            AtA(2, 1) += w_row * x2 * x1;
            AtA(2, 2) += w_row * x2 * x2;

            Atb(0) += w_row * x0 * k_theta;
            Atb(1) += w_row * x1 * k_theta;
            Atb(2) += w_row * x2 * k_theta;

          }  // end slicing
          IRL::Normal n_global_avg(n_global_sum[0], n_global_sum[1],
                                   n_global_sum[2]);
          n_global_avg.normalize();
          // n_global_avg = IRL::Normal(polygon_centroid[0] - 0.35,
          //                            polygon_centroid[1] - 0.35,
          //                            polygon_centroid[2] - 0.35);
          // n_global_avg.normalize();

          // solving least squares for principal curvatures
          bool fit_ok = false;
          double k1 = 0.0, k2 = 0.0, phi = 0.0;
          Eigen::LDLT<Eigen::Matrix3d> ldlt(AtA);
          if (ldlt.info() == Eigen::Success) {
            const Eigen::Vector3d abg = ldlt.solve(Atb);
            const double alpha = abg(0);
            const double beta = abg(1);
            const double gamma = abg(2);

            const double Rmag = std::sqrt(beta * beta + gamma * gamma);
            k1 = alpha + Rmag;  // max principal curvature
            k2 = alpha - Rmag;  // min principal curvature
            phi = 0.5 * std::atan2(gamma, beta);

            fit_ok =
                std::isfinite(k1) && std::isfinite(k2) && std::isfinite(phi);
          }

          if (!fit_ok) {
            // Fallback: keep planar if fit fails
            std::cout << "Taubin circle fit failed!" << std::endl;
            (*a_interface)(i, j, k) =
                std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
            continue;
          }

          //  Darboux frame
          IRL::ReferenceFrame temp_frame;
          temp_frame[0] = IRL::crossProduct(polygon_frame[1], n_global_avg);
          temp_frame[0].normalize();
          temp_frame[2] = n_global_avg;
          temp_frame[1] = IRL::crossProduct(temp_frame[2], temp_frame[0]);

          polygon_frame = temp_frame;

          // accumulator for principal curvature least squares
          AtA = Eigen::Matrix3d::Zero();
          Atb = Eigen::Vector3d::Zero();

          // sweeping over slices from [0, pi)
          for (int s = 0; s < nslices; s++) {
            // plane rotation angle
            double theta_s =
                M_PI * static_cast<double>(s) / static_cast<double>(nslices);

            // rotating target polygon frame about normal
            const IRL::UnitQuaternion rotation(theta_s, polygon_frame[2]);
            const auto local_frame = rotation * polygon_frame;

            // slicing plane
            const IRL::Plane slicing_plane(local_frame[1],
                                           local_frame[1] * polygon_centroid);

            // intersection of slicing plane with polygons in stencil
            std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
            std::vector<double> vfrac_list;
            std::vector<IRL::Normal> normal_list;
            for (int p = 0; p < static_cast<int>(polygon_vfrac_list.size());
                 p++) {
              getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                                 &intersections);
              if (intersections.size() != 2) continue;
              IRL::Pt start_point = intersections[0];
              IRL::Pt end_point = intersections[1];
              IRL::Normal neighbor_normal =
                  calculatePolygonNormal(polygon_vfrac_list[p].first);
              end_points_list.push_back({start_point, end_point});
              vfrac_list.push_back(polygon_vfrac_list[p].second);
              normal_list.push_back(neighbor_normal);
            }
            // taubin circle data
            TaubinCircleData tcd =
                getTaubinData(end_points_list, nsamples_per_segment, vfrac_list,
                              normal_list, mesh.dx(), polygon_normal,
                              polygon_centroid, local_frame, polygon_centroid);

            // local normal
            Eigen::Vector2d n_taubin = estimateTaubinNormal(tcd.center, tcd.R);
            Eigen::Vector3d n_loc(n_taubin[0], n_taubin[1], 0.0);

            // convertint to global frame
            Eigen::Vector3d e1(local_frame[0][0], local_frame[0][1],
                               local_frame[0][2]);
            Eigen::Vector3d e2(local_frame[2][0], local_frame[2][1],
                               local_frame[2][2]);
            Eigen::Vector3d e3(local_frame[1][0], local_frame[1][1],
                               local_frame[1][2]);
            Eigen::Matrix3d R;
            R.col(0) = e1;
            R.col(1) = e2;
            R.col(2) = e3;
            Eigen::Vector3d n_global = R * n_loc;
            n_global_sum += n_global;

            // directional curvature
            double k_theta = tcd.k;
            const double x0 = 1.0;
            const double x1 = std::cos(2.0 * theta_s);
            const double x2 = std::sin(2.0 * theta_s);
            const double w_row = 1.0;  // slicing weight (if any)

            AtA(0, 0) += w_row * x0 * x0;
            AtA(0, 1) += w_row * x0 * x1;
            AtA(0, 2) += w_row * x0 * x2;
            AtA(1, 0) += w_row * x1 * x0;
            AtA(1, 1) += w_row * x1 * x1;
            AtA(1, 2) += w_row * x1 * x2;
            AtA(2, 0) += w_row * x2 * x0;
            AtA(2, 1) += w_row * x2 * x1;
            AtA(2, 2) += w_row * x2 * x2;

            Atb(0) += w_row * x0 * k_theta;
            Atb(1) += w_row * x1 * k_theta;
            Atb(2) += w_row * x2 * k_theta;

          }  // end slicing

          // solving least squares for principal curvatures
          fit_ok = false;
          k1 = 0.0, k2 = 0.0, phi = 0.0;
          Eigen::LDLT<Eigen::Matrix3d> ldlt2(AtA);
          if (ldlt2.info() == Eigen::Success) {
            const Eigen::Vector3d abg = ldlt2.solve(Atb);
            const double alpha = abg(0);
            const double beta = abg(1);
            const double gamma = abg(2);

            const double Rmag = std::sqrt(beta * beta + gamma * gamma);
            k1 = alpha + Rmag;  // max principal curvature
            k2 = alpha - Rmag;  // min principal curvature
            phi = 0.5 * std::atan2(gamma, beta);

            fit_ok =
                std::isfinite(k1) && std::isfinite(k2) && std::isfinite(phi);
          }

          if (!fit_ok) {
            // Fallback: keep planar if fit fails
            std::cout << "Taubin circle fit failed!" << std::endl;
            (*a_interface)(i, j, k) =
                std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
            continue;
          }

          // const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
          // const IRL::ReferenceFrame darboux_frame = rotate_phi *
          // polygon_frame;

          const IRL::UnitQuaternion rotate_phi(phi, temp_frame[2]);
          const IRL::ReferenceFrame darboux_frame = rotate_phi * temp_frame;
          // new paraboloid
          IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 * k1,
                                     0.5 * k2);

          // exact normal
          // IRL::Normal sphere_normal(polygon_centroid[0] - 0.35,
          //                           polygon_centroid[1] - 0.35,
          //                           polygon_centroid[2] - 0.35);
          // sphere_normal.normalize();
          // IRL::ReferenceFrame sphere_frame =
          //     referenceFrameFromNormal(sphere_normal);

          // translate paraboloid to match volume fraction
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, vfrac(i, j, k), 1.0e-14, paraboloid);
          paraboloid.setDatum(
              IRL::Pt(polygon_centroid +
                      solver_distance.getDistance() * darboux_frame[2]));
          (*a_interface)(i, j, k) = paraboloid;
        }
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// void SlicesTaubin::getReconstruction(
//     const Data<IRL::VolumeMoments>& a_liq_moments,
//     const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
//     const Data<double>& a_U, const Data<double>& a_V, const Data<double>&
//     a_W, Data<IRL::SeparatorVariant>* a_interface, const bool
//     a_plic_already_built) {
//   // plic
//   if (!a_plic_already_built) {
//     LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
//     a_W,
//                              a_interface);
//   }

//   // clipped polygon from planes
//   const BasicMesh& mesh = a_liq_moments.getMesh();
//   Data<IRL::Polygon> polygon(&mesh);
//   Data<double> vfrac(&mesh);
//   for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
//       for (int i = mesh.imin(); i <= mesh.imax(); i++) {
//         vfrac(i, j, k) = a_liq_moments(i, j, k).volume() /
//         mesh.cell_volume(); if (vfrac(i, j, k) <
//         IRL::global_constants::VF_LOW ||
//             vfrac(i, j, k) > IRL::global_constants::VF_HIGH) {
//           continue;
//         }
//         auto cell = IRL::RectangularCuboid::fromBoundingPts(
//             IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
//             IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
//         const auto planar_separator =
//             std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
//         polygon(i, j, k) =
//         IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
//             cell, planar_separator, planar_separator[0]);
//       }
//     }
//   }
//   updatePolygonBorder(&polygon);

//   // outputting params
//   // std::string path = "/home/parinht2/Desktop/Curvature
//   sampling/paraview/";
//   // std::string csv_path = "/home/parinht2/Desktop/Curvature sampling/";
//   // std::string filename = "slicing_taubin_normals.vtk";
//   // std::string filepath = path + filename;
//   // std::vector<IRL::Pt> vtk_centroid_list;
//   // std::vector<IRL::Normal> vtk_normal_list;
//   // std::ofstream csvfile(csv_path + "sphere_normals.csv");
//   // csvfile << "x0,y0,z0,nx,ny,nz\n";

//   // slicing params
//   const int nsamples_per_segment = 10;
//   const int nslices = 18;
//   const int nlayers = 1;

//   std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
//   polygon_vfrac_list.reserve(125);
//   IRL::StackVector<IRL::Pt, 2> intersections;

//   for (int i = mesh.imin(); i <= mesh.imax(); i++) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
//       for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
//         if (polygon(i, j, k).getNumberOfVertices() > 2) {
//           // target cell local frame
//           const IRL::Normal polygon_normal =
//               calculatePolygonNormal(polygon(i, j, k));
//           const IRL::ReferenceFrame polygon_frame =
//               referenceFrameFromNormal(polygon_normal);
//           const IRL::Pt polygon_centroid = polygon(i, j,
//           k).calculateCentroid();

//           // building stencil
//           polygon_vfrac_list.resize(0);
//           for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
//             for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
//               for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
//                 if (polygon(ii, jj, kk).getNumberOfVertices() > 2) {
//                   polygon_vfrac_list.push_back(
//                       std::make_pair(polygon(ii, jj, kk), vfrac(ii, jj,
//                       kk)));
//                 }
//               }
//             }
//           }

//           // accumulator for principal curvature least squares
//           Eigen::Matrix3d AtA = Eigen::Matrix3d::Zero();
//           Eigen::Vector3d Atb = Eigen::Vector3d::Zero();

//           // sweeping over slices from [0, pi)
//           for (int s = 0; s < nslices; s++) {
//             // plane rotation angle
//             double theta_s =
//                 M_PI * static_cast<double>(s) / static_cast<double>(nslices);

//             // rotating target polygon frame about normal
//             const IRL::UnitQuaternion rotation(theta_s, polygon_frame[2]);
//             const auto local_frame = rotation * polygon_frame;

//             // slicing plane
//             const IRL::Plane slicing_plane(local_frame[1],
//                                            local_frame[1] *
//                                            polygon_centroid);

//             // intersection of slicing plane with polygons in stencil
//             std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
//             std::vector<double> vfrac_list;
//             std::vector<IRL::Normal> normal_list;
//             for (int p = 0; p < static_cast<int>(polygon_vfrac_list.size());
//                  p++) {
//               getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
//                                  &intersections);
//               if (intersections.size() != 2) continue;
//               IRL::Pt start_point = intersections[0];
//               IRL::Pt end_point = intersections[1];
//               IRL::Normal neighbor_normal =
//                   calculatePolygonNormal(polygon_vfrac_list[p].first);
//               end_points_list.push_back({start_point, end_point});
//               vfrac_list.push_back(polygon_vfrac_list[p].second);
//               normal_list.push_back(neighbor_normal);
//             }
//             // signed curvature using taubin fit
//             double k_theta = getTaubinCurvature(
//                 end_points_list, nsamples_per_segment, vfrac_list,
//                 normal_list, mesh.dx(), polygon_normal, polygon_centroid,
//                 local_frame, polygon_centroid);
//             const double x0 = 1.0;
//             const double x1 = std::cos(2.0 * theta_s);
//             const double x2 = std::sin(2.0 * theta_s);
//             const double w_row = 1.0;  // slicing weight (if any)

//             AtA(0, 0) += w_row * x0 * x0;
//             AtA(0, 1) += w_row * x0 * x1;
//             AtA(0, 2) += w_row * x0 * x2;
//             AtA(1, 0) += w_row * x1 * x0;
//             AtA(1, 1) += w_row * x1 * x1;
//             AtA(1, 2) += w_row * x1 * x2;
//             AtA(2, 0) += w_row * x2 * x0;
//             AtA(2, 1) += w_row * x2 * x1;
//             AtA(2, 2) += w_row * x2 * x2;

//             Atb(0) += w_row * x0 * k_theta;
//             Atb(1) += w_row * x1 * k_theta;
//             Atb(2) += w_row * x2 * k_theta;

//           }  // end slicing
//           // solving least squares for principal curvatures
//           bool fit_ok = false;
//           double k1 = 0.0, k2 = 0.0, phi = 0.0;
//           Eigen::LDLT<Eigen::Matrix3d> ldlt(AtA);
//           if (ldlt.info() == Eigen::Success) {
//             const Eigen::Vector3d abg = ldlt.solve(Atb);
//             const double alpha = abg(0);
//             const double beta = abg(1);
//             const double gamma = abg(2);

//             const double Rmag = std::sqrt(beta * beta + gamma * gamma);
//             k1 = alpha + Rmag;  // max principal curvature
//             k2 = alpha - Rmag;  // min principal curvature
//             phi = 0.5 * std::atan2(gamma, beta);

//             fit_ok =
//                 std::isfinite(k1) && std::isfinite(k2) && std::isfinite(phi);
//           }

//           if (!fit_ok) {
//             // Fallback: keep planar if fit fails
//             std::cout << "Taubin circle fit failed!" << std::endl;
//             (*a_interface)(i, j, k) =
//                 std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
//             continue;
//           }

//           //  Darboux frame
//           const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
//           const IRL::ReferenceFrame darboux_frame = rotate_phi *
//           polygon_frame;

//           // // exact normal
//           // IRL::Normal sphere_normal(polygon_centroid[0] - 0.35,
//           //                           polygon_centroid[1] - 0.35,
//           //                           polygon_centroid[2] - 0.35);
//           // sphere_normal.normalize();
//           // IRL::ReferenceFrame sphere_frame =
//           //     referenceFrameFromNormal(sphere_normal);
//           // IRL::Paraboloid paraboloid(polygon_centroid, sphere_frame, 0.5 *
//           // k1,
//           //                            0.5 * k2);
//           // vtk_centroid_list.push_back(polygon_centroid);
//           // vtk_normal_list.push_back(sphere_normal);

//           // // normals to csv
//           // csvfile << polygon_centroid[0] << "," << polygon_centroid[1] <<
//           ","
//           //         << polygon_centroid[2] << "," << sphere_normal[0] << ","
//           //         << sphere_normal[1] << "," << sphere_normal[2] << "\n";

//           // new paraboloid
//           IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 *
//           k1,
//                                      0.5 * k2);
//           // vtk_centroid_list.push_back(polygon_centroid);
//           // vtk_normal_list.push_back(darboux_frame[2]);

//           // translate paraboloid to match volume fraction
//           const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
//           const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
//                                       mesh.z(k + 1));
//           const auto cell = IRL::RectangularCuboid::fromBoundingPts(
//               lower_cell_pt, upper_cell_pt);
//           IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
//               solver_distance(cell, vfrac(i, j, k), 1.0e-14, paraboloid);
//           paraboloid.setDatum(
//               IRL::Pt(polygon_centroid +
//                       solver_distance.getDistance() * darboux_frame[2]));
//           (*a_interface)(i, j, k) = paraboloid;
//         }
//       }
//     }
//   }
//   // writeVectorsVTK(filepath, vtk_centroid_list, vtk_normal_list);
//   // writeScatterVTK(vtk_centroid_list, path +
//   "slicing_taubin_centroids.vtk");
//   // csvfile.close();

//   a_interface->updateBorder();
//   correctInterfaceBorders(a_interface);
// }

void PLICAligned::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface, const bool a_plic_already_built) {
  // build plic
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface);
  }

  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::SeparatorVariant> plic_reconstruction(&mesh);
  Data<IRL::SeparatorVariant> jibben_reconstruction(&mesh);
  for (int i = mesh.imino(); i <= mesh.imaxo(); i++) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); j++) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); k++) {
        plic_reconstruction(i, j, k) = (*a_interface)(i, j, k);
        jibben_reconstruction(i, j, k) = (*a_interface)(i, j, k);
      }
    }
  }

  using VolumeMomentsAndSurface =
      IRL::AddSurfaceOutput<IRL::VolumeMoments,
                            IRL::ParaboloidParametrizedSurfaceOutput>;
  Data<IRL::Pt> interface_centroids(&mesh);
  Data<double> interface_areas(&mesh), jibben_errors(&mesh);
  // Jibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
  // a_W,
  //                           &jibben_reconstruction, true,
  //                           &interface_centroids, &interface_areas,
  //                           &jibben_errors);
  SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                  a_W, &jibben_reconstruction, true);
  Data<IRL::Normal> average_normals(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          if (IRL::Paraboloid* paraboloid = std::get_if<IRL::Paraboloid>(
                  &jibben_reconstruction(i, j, k))) {
            // finding average normals
            auto cell = IRL::RectangularCuboid::fromBoundingPts(
                IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
                IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
            auto surface = IRL::getVolumeMoments<VolumeMomentsAndSurface>(
                               cell, *paraboloid)
                               .getSurface();
            average_normals(i, j, k) = surface.getAverageNormalNonAligned();
            const auto plane = std::get_if<IRL::PlanarSeparator>(
                &(plic_reconstruction(i, j, k)));
            (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
                IRL::Plane(average_normals(i, j, k), (*plane)[0].distance()));
            IRL::setDistanceToMatchVolumeFraction(cell, liquid_volume_fraction,
                                                  &(*a_interface)(i, j, k),
                                                  1.0e-14);
          }
        }
      }
    }
  }
  SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                  a_W, a_interface, true);
}

// functions for particle method --------------------------------------

bool isTargetSegment(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                     double tol = 1e-12) {
  // collinearity test
  double cross = (p2[0] - p1[0]) * (-p1[1]) - (p2[1] - p1[1]) * (-p1[1]);
  if (std::fabs(cross) > tol) return false;

  // bounding box test
  if (std::min(p1[0], p2[0]) - tol <= 0. &&
      0. <= std::max(p1[0], p2[0]) + tol &&
      std::min(p1[1], p2[1]) - tol <= 0. &&
      0. <= std::max(p1[1], p2[1]) + tol) {
    return true;
  }
  return false;
}

// computing particle positions
std::vector<Eigen::Vector2d> ComputeParticlePositions(const int& N,
                                                      const Eigen::Vector2d& p,
                                                      const double& phi,
                                                      const double& theta,
                                                      const double& hp) {
  // N: Number of particles (odd)
  // p: coordinate of central particle (origin of local frame)
  // phi: orientation angle (angle between tangent of p and x-axis)
  // theta: bending angle (turning angle between chords of the circle)
  // hp: distance between particles (or chord length)

  std::vector<Eigen::Vector2d> particle_positions(N, Eigen::Vector2d::Zero());

  const int c = (N - 1) / 2;  // central particle index

  // parameterized coordinates of all other particles on the arc
  for (int i = 0; i < N; i++) {
    if (i > c) {
      for (int j = 1; j <= i - c; j++) {
        particle_positions[i] +=
            hp * Eigen::Vector2d(
                     std::cos(phi + (static_cast<double>(j) - 0.5) * theta),
                     std::sin(phi + (static_cast<double>(j) - 0.5) * theta));
      }
      particle_positions[i] = p + particle_positions[i];
    } else if (i < c) {
      for (int j = 1; j <= c - i; j++) {
        particle_positions[i] +=
            hp * Eigen::Vector2d(
                     std::cos(phi - (static_cast<double>(j) - 0.5) * theta),
                     std::sin(phi - (static_cast<double>(j) - 0.5) * theta));
      }
      particle_positions[i] = p - particle_positions[i];
    } else {
      particle_positions[i] = p;
    }
  }

  return particle_positions;
}

// compute single particle force
Eigen::Vector2d ComputeParticleForce(
    const Eigen::Vector2d& x,
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>&
        line_seg_endpoints,
    const double& eta) {
  // x: position of particle
  // line_seg_endpoints: endpoints {a,b} of line segments that are cloest to the
  // particle

  Eigen::Vector2d particle_force = Eigen::Vector2d::Zero();

  // Computing closest distance to all line segments in the vicinity of the
  // particle
  for (int i = 0; i < line_seg_endpoints.size(); i++) {
    Eigen::Vector2d a = line_seg_endpoints[i].first;
    Eigen::Vector2d b = line_seg_endpoints[i].second;

    // finding t using projection
    Eigen::Vector2d ab = b - a;
    Eigen::Vector2d ax = x - a;
    double t = ax.dot(ab) / ab.squaredNorm();
    double t_clamped = std::max(0.0, std::min(1.0, t));

    // finding closet point on the line segment to the point
    Eigen::Vector2d y = a + t_clamped * ab;

    // Finding xy distance and keeping minimum value of "force"
    Eigen::Vector2d xy = y - x;

    if (i == 0) {
      particle_force = xy;
    } else {
      if (xy.norm() < particle_force.norm()) {
        particle_force = xy;
      }
    }
  }
  return (eta * particle_force);
}

// initializing particle positions
std::vector<Eigen::Vector2d> InitializeParticlePositions(
    const std::pair<Eigen::Vector2d, Eigen::Vector2d>& target_endpoints,
    const double& hp, const int& N) {
  // target_endpoints: end points of the target interface where curvature is to
  // be estimated hp: spacing between particles along line segment N: number of
  // particles (odd)

  std::vector<Eigen::Vector2d> initial_particle_positions(
      N, Eigen::Vector2d::Zero());

  // line segment end points
  Eigen::Vector2d a = target_endpoints.first;
  Eigen::Vector2d b = target_endpoints.second;

  // unit vector along the line segment
  Eigen::Vector2d unit_ab = (b - a) / (b - a).norm();

  // central particle at midpoint of line segment
  initial_particle_positions[(N - 1) / 2] = (a + b) / 2.0;

  // other particles are spaced by hp on either side of the central particle
  // along the line segment
  for (int i = 1; i <= (N - 1) / 2; i++) {
    initial_particle_positions[(N - 1) / 2 + i] =
        initial_particle_positions[(N - 1) / 2] + hp * unit_ab * i;
    initial_particle_positions[(N - 1) / 2 - i] =
        initial_particle_positions[(N - 1) / 2] - hp * unit_ab * i;
  }

  return initial_particle_positions;
}

// particle force projections
double ComputeParticleForceProjection(
    const int& N, const double& phi, const double& theta, const double& hp,
    const bool& iswrtPhi, const std::vector<Eigen::Vector2d> particle_forces) {
  // iswrtPhi: "true" will compute derivative wrt phi else wrt theta

  int c = (N - 1) / 2;  // central particle index

  std::vector<Eigen::Vector2d> position_derivative(N, Eigen::Vector2d::Zero());

  position_derivative[c] = Eigen::Vector2d(0.0, 0.0);  // central particle

  if (iswrtPhi == true) {
    for (int i = 1; i <= c; i++) {
      // i > c
      position_derivative[c + i] =
          position_derivative[c + (i - 1)] +
          hp * Eigen::Vector2d(std::cos(phi +
                                        (static_cast<double>(i) -
                                         static_cast<double>(c) - 0.5) *
                                            theta +
                                        M_PI / 2.0),
                               std::sin(phi +
                                        (static_cast<double>(i) -
                                         static_cast<double>(c) - 0.5) *
                                            theta +
                                        M_PI / 2.0));
      // i < c
      position_derivative[c - i] =
          position_derivative[c - (i - 1)] -
          hp * Eigen::Vector2d(std::cos(phi -
                                        (static_cast<double>(c) -
                                         static_cast<double>(i) - 0.5) *
                                            theta +
                                        M_PI / 2.0),
                               std::sin(phi -
                                        (static_cast<double>(c) -
                                         static_cast<double>(i) - 0.5) *
                                            theta +
                                        M_PI / 2.0));
    }
  } else {
    for (int i = 1; i <= c; i++) {
      // i > c
      position_derivative[c + i] =
          position_derivative[c + (i - 1)] +
          hp * (static_cast<double>(i) - static_cast<double>(c) - 0.5) *
              Eigen::Vector2d(std::cos(phi +
                                       (static_cast<double>(i) -
                                        static_cast<double>(c) - 0.5) *
                                           theta +
                                       M_PI / 2.0),
                              std::sin(phi +
                                       (static_cast<double>(i) -
                                        static_cast<double>(c) - 0.5) *
                                           theta +
                                       M_PI / 2.0));
      // i < c
      position_derivative[c - i] =
          position_derivative[c - (i - 1)] -
          hp * (static_cast<double>(c) - static_cast<double>(i) - 0.5) *
              Eigen::Vector2d(std::cos(phi -
                                       (static_cast<double>(c) -
                                        static_cast<double>(i) - 0.5) *
                                           theta +
                                       M_PI / 2.0),
                              std::sin(phi -
                                       (static_cast<double>(c) -
                                        static_cast<double>(i) - 0.5) *
                                           theta +
                                       M_PI / 2.0));
    }
  }

  // projecting force on derivative
  double num = 0.0, denom = 0.0;
  for (int i = 0; i < N; i++) {
    num += particle_forces[i].dot(position_derivative[i]);
    denom += position_derivative[i].dot(position_derivative[i]);
  }

  return (num / denom);
}

// finding circle center
Eigen::Vector2d findCircleCenter(const std::vector<Eigen::Vector2d>& points) {
  // input is 3 non-collinear points on circle
  Eigen::Vector2d A = points[0], B = points[1], C = points[2];

  // midpoint of AB and BC
  Eigen::Vector2d midAB((A.x() + B.x()) / 2.0, (A.y() + B.y()) / 2.0);
  Eigen::Vector2d midBC((B.x() + C.x()) / 2.0, (B.y() + C.y()) / 2.0);

  // slopes
  double dxAB = B.x() - A.x();
  double dyAB = B.y() - A.y();
  double dxBC = C.x() - B.x();
  double dyBC = C.y() - B.y();

  // collinearity check (area of triangle)
  double area = A.x() * (B.y() - C.y()) + B.x() * (C.y() - A.y()) +
                C.x() * (A.y() - B.y());
  if (std::abs(area) < 1e-10) {
    // throw std::runtime_error("Points are collinear. Circle is undefined");
    std::cout << "Points are collinear. Circle is undefined" << std::endl;
  }

  // slopes for perpendicular bisectors
  double slopePerpAB = 0.0, slopePerpBC = 0.0;

  // handling vertical line segments
  bool verticalAB = (std::abs(dxAB) < 1e-10);
  bool verticalBC = (std::abs(dxBC) < 1e-10);
  if (!verticalAB) slopePerpAB = -dxAB / dyAB;
  if (!verticalBC) slopePerpBC = -dxBC / dyBC;

  // solving for center of circle
  double cx, cy;
  if (verticalAB) {
    // AB is vertical → perpendicular bisector is horizontal
    cy = midAB.y();
    cx = slopePerpBC * (cy - midBC.y()) + midBC.x();
  } else if (verticalBC) {
    // BC is vertical → perpendicular bisector is horizontal
    cy = midBC.y();
    cx = slopePerpAB * (cy - midAB.y()) + midAB.x();
  } else {
    // Solve intersection of two lines
    cx = (slopePerpAB * midAB.x() - slopePerpBC * midBC.x() + midBC.y() -
          midAB.y()) /
         (slopePerpAB - slopePerpBC);
    cy = slopePerpAB * (cx - midAB.x()) + midAB.y();
  }

  return Eigen::Vector2d(cx, cy);
}

// finding curvature
double getParticleMethodCurvature(
    const std::pair<Eigen::Vector2d, Eigen::Vector2d>& target_endpoints,
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& endpoints,
    const int& N, const double& Hp, const double& h, const double& eta, int s) {
  double theta, phi;
  std::vector<Eigen::Vector2d> particle_positions(N, Eigen::Vector2d::Zero()),
      particle_forces(N, Eigen::Vector2d::Zero()),
      particle_positions_prev(N, Eigen::Vector2d::Zero()),
      particle_forces_prev(N, Eigen::Vector2d::Zero()),
      particle_positions_s(N, Eigen::Vector2d::Zero()),
      particle_positions_ss(N, Eigen::Vector2d::Zero()),
      particle_forces_s(N, Eigen::Vector2d::Zero()),
      particle_forces_ss(N, Eigen::Vector2d::Zero());

  // particle spacing
  double hp = Hp * h / (static_cast<double>(N) - 1.0);

  // initializing particle positions
  particle_positions = InitializeParticlePositions(target_endpoints, hp, N);

  // initialize forces
  for (int i = 0; i < N; i++) {
    particle_forces[i] =
        ComputeParticleForce(particle_positions[i], endpoints, eta);
  }

  // initialize orientation and bending angle
  theta = 0.0;
  Eigen::Vector2d ab_star = target_endpoints.second - target_endpoints.first;
  phi = std::atan2(ab_star[1], ab_star[0]);
  if (phi < 0.0) {
    phi += 2.0 * M_PI;
  }

  // iteration parameters
  int max_iter = 100;
  double tol = 1e-3;
  int iter = 0;
  double residual = 1.0;

  // index of central particle
  int c = (N - 1) / 2;

  // update positions and forces
  while (std::abs(residual) > tol) {
    iter++;

    // prev iter
    particle_positions_prev = particle_positions;
    particle_forces_prev = particle_forces;

    // step 1: correct central particle position using force
    particle_positions[c] += particle_forces[c];

    // step 1: change in position for other particles
    particle_positions_s =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 1: subtracting change of position from forces
    for (int i = 0; i < N; i++) {
      particle_forces_s[i] =
          particle_forces_prev[i] -
          (particle_positions_s[i] - particle_positions_prev[i]);
    }

    // step 2: correct phi by projection of force
    phi += ComputeParticleForceProjection(N, phi, theta, hp, true,
                                          particle_forces_s);

    // step 2: change in position for other particles
    particle_positions_ss =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 2: subtracting change of position from forces
    for (int i = 0; i < N; i++) {
      particle_forces_ss[i] = particle_forces_s[i] - (particle_positions_ss[i] -
                                                      particle_positions_s[i]);
    }

    // step 3: correct theta by projection of force
    theta -= ComputeParticleForceProjection(N, phi, theta, hp, false,
                                            particle_forces_ss);

    // step 3: update particle positions
    particle_positions =
        ComputeParticlePositions(N, particle_positions[c], phi, theta, hp);

    // step 3: update particle forces
    for (int i = 0; i < N; i++) {
      particle_forces[i] =
          ComputeParticleForce(particle_positions[i], endpoints, eta);
    }

    // residual: change in position for all particles (max value among all
    // particles)
    residual = 0.0;
    for (int i = 0; i < N; i++) {
      residual =
          std::max(
              residual,
              (particle_positions[i] - particle_positions_prev[i]).norm()) /
          (eta * h);
    }

    if (iter == max_iter) {
      // std::cout << "Maximum iterations reached for particle method. Residual
      // = "
      //           << residual << std::endl;
      break;
    }
  }

  // std::cout << "phi: " << phi * 180.0 / M_PI
  //           << " theta: " << theta * 180.0 / M_PI << std::endl;
  // std::cout << "residual: " << residual << " iter: " << iter << std::endl;

  // outputting final positions and forces to csv
  // std::string csv_path =
  //     "/home/parinht2/Desktop/Curvature "
  //     "sampling/particle_method_csv/";
  // std::string csv_filename =
  //     csv_path + "final_particle_data_" + std::to_string(s) + ".csv";
  // std::ofstream csvfile(csv_filename);
  // csvfile << "particle_id,x,y,Fx,Fy\n";
  // for (int p = 0; p < particle_positions.size(); p++) {
  //   csvfile << p + 1 << "," << particle_positions[p][0] << ","
  //           << particle_positions[p][1] << "," << particle_forces[p][0] <<
  //           ","
  //           << particle_forces[p][1] << "\n";
  // }
  // csvfile.close();

  // circle center for sign
  Eigen::Vector2d xc = findCircleCenter(particle_positions);
  double yc = xc.y();
  double sign = (yc >= 0.0) ? -1.0 : +1.0;

  return sign * std::abs((2 * std::sin(theta / 2.0) / hp));  // curvature
}

// reconstruction using particle method on each slice
void SlicesParticle::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface, const bool a_plic_already_built) {
  // plic
  if (!a_plic_already_built) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface);
  }

  // clipped polygon from planes
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::Polygon> polygon(&mesh);
  Data<double> vfrac(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        vfrac(i, j, k) = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vfrac(i, j, k) < IRL::global_constants::VF_LOW ||
            vfrac(i, j, k) > IRL::global_constants::VF_HIGH) {
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, planar_separator, planar_separator[0]);
      }
    }
  }
  updatePolygonBorder(&polygon);

  // slicing params
  const int nslices = 18;
  const int nlayers = 2;

  // particle method params
  const int N = 7;
  const double Hp = 5.0;
  const double eta = 1.0;
  const double h = mesh.dx();
  const double hp = Hp * h / (static_cast<double>(N) - 1.0);

  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);
  IRL::StackVector<IRL::Pt, 2> intersections;

  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        if (polygon(i, j, k).getNumberOfVertices() > 2) {
          // if (i == 6 && j == 8 && k == 10) {  // REMOVE AFTER DONE
          // target cell local frame
          const IRL::Normal polygon_normal =
              calculatePolygonNormal(polygon(i, j, k));
          const IRL::ReferenceFrame polygon_frame =
              referenceFrameFromNormal(polygon_normal);
          const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();

          // building stencil
          polygon_vfrac_list.resize(0);
          for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
            for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
              for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
                if (polygon(ii, jj, kk).getNumberOfVertices() > 2) {
                  polygon_vfrac_list.push_back(
                      std::make_pair(polygon(ii, jj, kk), vfrac(ii, jj, kk)));
                }
              }
            }
          }

          // accumulator for principal curvature least squares
          Eigen::Matrix3d AtA = Eigen::Matrix3d::Zero();
          Eigen::Vector3d Atb = Eigen::Vector3d::Zero();

          // sweeping over slices from [0, pi)
          for (int s = 0; s < nslices; s++) {
            // std::cout << "slice = " << s << " ------------------------"
            //           << std::endl;
            // plane rotation angle
            double theta_s =
                M_PI * static_cast<double>(s) / static_cast<double>(nslices);

            // rotating target polygon frame about normal
            const IRL::UnitQuaternion rotation(theta_s, polygon_frame[2]);
            const auto local_frame = rotation * polygon_frame;

            // rotation matrix for global to local frame
            Eigen::Vector3d e1(local_frame[0][0], local_frame[0][1],
                               local_frame[0][2]);
            Eigen::Vector3d e2(local_frame[2][0], local_frame[2][1],
                               local_frame[2][2]);
            Eigen::Vector3d e3(local_frame[1][0], local_frame[1][1],
                               local_frame[1][2]);
            Eigen::Matrix3d R;
            R.col(0) = e1;
            R.col(1) = e2;
            R.col(2) = e3;

            // local frame origin
            Eigen::Vector3d o(polygon_centroid[0], polygon_centroid[1],
                              polygon_centroid[2]);

            // slicing plane
            const IRL::Plane slicing_plane(local_frame[1],
                                           local_frame[1] * polygon_centroid);

            // intersection of slicing plane with polygons in stencil
            std::vector<std::pair<IRL::Pt, IRL::Pt>> endpoints_list;
            for (int p = 0; p < static_cast<int>(polygon_vfrac_list.size());
                 p++) {
              getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                                 &intersections);
              if (intersections.size() != 2) continue;
              IRL::Pt start_point = intersections[0];
              IRL::Pt end_point = intersections[1];
              endpoints_list.push_back({start_point, end_point});
            }

            // converting endpoints to local frame (2D)
            std::pair<Eigen::Vector2d, Eigen::Vector2d> target_endpoints = {
                Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero()};
            std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>
                endpoints_list_local;
            int target_index = 0;
            for (int e = 0; e < endpoints_list.size(); e++) {
              Eigen::Vector3d x0(endpoints_list[e].first[0],
                                 endpoints_list[e].first[1],
                                 endpoints_list[e].first[2]);
              Eigen::Vector3d x1(endpoints_list[e].second[0],
                                 endpoints_list[e].second[1],
                                 endpoints_list[e].second[2]);
              Eigen::Vector3d x0_local = R.transpose() * (x0 - o);
              Eigen::Vector3d x1_local = R.transpose() * (x1 - o);
              endpoints_list_local.push_back(
                  {Eigen::Vector2d(x0_local[0], x0_local[1]),
                   Eigen::Vector2d(x1_local[0], x1_local[1])});
              // target line segment
              if (isTargetSegment(Eigen::Vector2d(x0_local[0], x0_local[1]),
                                  Eigen::Vector2d(x1_local[0], x1_local[1]))) {
                target_index = e;
                target_endpoints = {Eigen::Vector2d(x0_local[0], x0_local[1]),
                                    Eigen::Vector2d(x1_local[0], x1_local[1])};
              }
            }

            // finding directional curvature
            double k_theta = getParticleMethodCurvature(
                target_endpoints, endpoints_list_local, N, Hp, h, eta, s);
            // k_theta = std::abs(k_theta);
            // std::cout << "Curvature = " << k_theta << std::endl;

            // OUTPUTTING
            // std::cout << "s = " << s << " k = " << k_theta << std::endl;
            // std::string csv_path =
            //     "/home/parinht2/Desktop/Curvature "
            //     "sampling/particle_method_csv/";
            // std::string csv_filename = csv_path + "initial_particle_data_"
            // +
            //                            std::to_string(s) + ".csv";
            // std::ofstream csvfile(csv_filename);
            // csvfile << "particle_id,x,y,Fx,Fy\n";
            // std::vector<Eigen::Vector2d> initial_particle_positions =
            //     InitializeParticlePositions(target_endpoints, hp, N);
            // std::vector<Eigen::Vector2d> initial_particle_forces(
            //     initial_particle_positions.size(),
            //     Eigen::Vector2d::Zero());
            // for (int p = 0; p < initial_particle_positions.size(); p++) {
            //   initial_particle_forces[p] = ComputeParticleForce(
            //       initial_particle_positions[p], endpoints_list_local,
            //       eta);
            //   csvfile << p + 1 << "," << initial_particle_positions[p][0]
            //           << "," << initial_particle_positions[p][1] << ","
            //           << initial_particle_forces[p][0] << ","
            //           << initial_particle_forces[p][1] << "\n";
            // }
            // csvfile.close();

            const double x0 = 1.0;
            const double x1 = std::cos(2.0 * theta_s);
            const double x2 = std::sin(2.0 * theta_s);
            const double w_row = 1.0;  // slicing weight (if any)

            AtA(0, 0) += w_row * x0 * x0;
            AtA(0, 1) += w_row * x0 * x1;
            AtA(0, 2) += w_row * x0 * x2;
            AtA(1, 0) += w_row * x1 * x0;
            AtA(1, 1) += w_row * x1 * x1;
            AtA(1, 2) += w_row * x1 * x2;
            AtA(2, 0) += w_row * x2 * x0;
            AtA(2, 1) += w_row * x2 * x1;
            AtA(2, 2) += w_row * x2 * x2;

            Atb(0) += w_row * x0 * k_theta;
            Atb(1) += w_row * x1 * k_theta;
            Atb(2) += w_row * x2 * k_theta;

          }  // end slicing

          // solving least squares for principal curvatures
          bool fit_ok = false;
          double k1 = 0.0, k2 = 0.0, phi = 0.0;
          Eigen::LDLT<Eigen::Matrix3d> ldlt(AtA);
          if (ldlt.info() == Eigen::Success) {
            const Eigen::Vector3d abg = ldlt.solve(Atb);
            const double alpha = abg(0);
            const double beta = abg(1);
            const double gamma = abg(2);

            const double Rmag = std::sqrt(beta * beta + gamma * gamma);
            k1 = alpha + Rmag;  // max principal curvature
            k2 = alpha - Rmag;  // min principal curvature
            phi = 0.5 * std::atan2(gamma, beta);

            fit_ok =
                std::isfinite(k1) && std::isfinite(k2) && std::isfinite(phi);
          }

          if (!fit_ok) {
            // Fallback: keep planar if fit fails
            std::cout << "Unable to obtain principal curvatures" << std::endl;
            (*a_interface)(i, j, k) =
                std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
            continue;
          }

          //  Darboux frame
          const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
          const IRL::ReferenceFrame darboux_frame = rotate_phi * polygon_frame;

          // new paraboloid
          IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 * k1,
                                     0.5 * k2);
          // std::cout << "a = " << 0.5 * k1 << "  b = " << 0.5 * k2 <<
          // std::endl;

          // translate paraboloid to match volume fraction
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, vfrac(i, j, k), 1.0e-14, paraboloid);
          paraboloid.setDatum(
              IRL::Pt(polygon_centroid +
                      solver_distance.getDistance() * darboux_frame[2]));
          (*a_interface)(i, j, k) = paraboloid;
          // }
        }
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// ------------------------------------------------------------------------------------

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