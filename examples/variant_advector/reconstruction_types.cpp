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
#include <map>
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
#include "irl/quadratic_reconstruction/gauss_legendre_integrator.h"

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
                       Data<IRL::SeparatorVariant>* a_interface,
                       std::vector<InterfaceScalarField>* a_scalar_fields) {
  if (a_reconstruction_method == "ELVIRA") {
    ELVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "LVIRA" ||
             a_reconstruction_method == "PLIC") {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "Jibben") {
    Jibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "Jibben2") {
    Jibben2::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                               a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "JibbenCubic") {
    JibbenCubic::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                   a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "JibbenM") {
    JibbenM::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                               a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "PU") {
    PU::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                          a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "MixedJibben") {
    MixedPLICJibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                       a_V, a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "SlicesParabola") {
    SlicesParabola::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                      a_V, a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "SlicesTaubin") {
    SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                    a_V, a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "SlicesTaubin2") {
    SlicesTaubin2::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                     a_V, a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "SlicesTaubin3") {
    SlicesTaubin3::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                     a_V, a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "SlicesTaubinLS") {
    SlicesTaubinLS::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                      a_V, a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "SlicesTaubinS") {
    SlicesTaubinS::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                     a_V, a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "PLICalign") {
    PLICalign::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                 a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "PLICalign2") {
    PLICalign2::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                  a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "MossoSwartz") {
    MossoSwartz::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                   a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "Hybrid") {
    Hybrid::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "Hybrid2") {
    Hybrid2::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                               a_W, a_interface, a_scalar_fields);
  } else if (a_reconstruction_method == "SlicesParticle") {
    SlicesParticle::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                      a_V, a_W, a_interface, a_scalar_fields);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: PLIC, Jibben, PU, MixedJibben, "
                 "SlicesParabola, SlicesTaubin, SlicesParticle. \n";
    std::exit(-1);
  }
}

void ELVIRA::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields) {
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

void LVIRA::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    ELVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface, a_scalar_fields);
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

void Jibben::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built, Data<IRL::Pt>* a_centroids,
    Data<double>* a_areas, Data<double>* a_errors) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
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

// ============ Jibben based on squared volume ===============
void Jibben2::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built, Data<IRL::Pt>* a_centroids,
    Data<double>* a_areas, Data<double>* a_errors) {
  // First, we need to build the plic
  LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                           a_interface, a_scalar_fields);

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

  // Now let's compute the Jibben parabolic fit
  IRL::JibbenNeighborhood neighborhood;
  const int nlayers = 1;
  const int nstencil =
      (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);
  const double delta = 2.5 * mesh.dx();
  neighborhood.reserve(nstencil);
  neighborhood.setDelta(delta);

  InterfaceScalarField normal_error("normal_metric", &mesh);
  InterfaceScalarField normal_eigen_error("normal_eigen_metric", &mesh);
  InterfaceScalarField normal_std_error("normal_std_metric", &mesh);

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

          // jibben fit
          IRL::Jibben_3D jibben_solver;
          (*a_interface)(i, j, k) = jibben_solver.solve2(&neighborhood, delta);

          // error metrics
          // double err = jibben_solver.getVolumeErrorSquared(mesh.dx());
          double normal_err = jibben_solver.getNormalMetric();
          double normal_eigen_err = jibben_solver.getNormalEigenMetric();
          double normal_std_err = jibben_solver.getNormalVarianceMetric();

          // Match to volume fraction
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::setDistanceToMatchVolumeFraction(
              cell, liquid_volume_fraction, &(*a_interface)(i, j, k), 1.0e-14);

          if (const auto ptr =
                  std::get_if<IRL::Paraboloid>(&(*a_interface)(i, j, k))) {
            normal_error.paraboloid_scalar_data(i, j, k) = normal_err;
            normal_eigen_error.paraboloid_scalar_data(i, j, k) =
                normal_eigen_err;
            normal_std_error.paraboloid_scalar_data(i, j, k) = normal_std_err;
            // squared_vol_error.paraboloid_scalar_data(i, j, k) = err;
          } else if (const auto ptr = std::get_if<IRL::PlanarSeparator>(
                         &(*a_interface)(i, j, k))) {
            normal_error.polygon_scalar_data(i, j, k) = normal_err;
            normal_eigen_error.polygon_scalar_data(i, j, k) = normal_eigen_err;
            normal_std_error.polygon_scalar_data(i, j, k) = normal_std_err;
            // squared_vol_error.polygon_scalar_data(i, j, k) = err;
          }
        }
      }
    }
  }
  a_scalar_fields->push_back(normal_error);
  a_scalar_fields->push_back(normal_eigen_error);
  a_scalar_fields->push_back(normal_std_error);

  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// ============ Jibben based on cubic fit ===============
void JibbenCubic::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built, Data<IRL::Pt>* a_centroids,
    Data<double>* a_areas, Data<double>* a_errors) {
  // First, we need to build the plic
  LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                           a_interface, a_scalar_fields);

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

          // jibben fit
          IRL::Jibben_3D jibben_solver;
          (*a_interface)(i, j, k) = jibben_solver.solveCubic(&neighborhood);

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

// ============ Modified jibben to output metrics ===============
void JibbenM::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields) {
  // First, we need to build the plic
  LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                           a_interface, a_scalar_fields);

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

  // Now let's compute the Jibben parabolic fit
  IRL::JibbenNeighborhood neighborhood;
  const int nlayers = 1;
  const int nstencil =
      (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);
  neighborhood.reserve(nstencil);
  neighborhood.setDelta(2.5 * mesh.dx());

  // InterfaceScalarField squared_vol_error("squared_volume_error", &mesh);
  InterfaceScalarField normal_error("normal_metric", &mesh);
  InterfaceScalarField normal_eigen_error("normal_eigen_metric", &mesh);
  InterfaceScalarField normal_std_error("normal_std_metric", &mesh);

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

          // jibben fit
          IRL::Jibben_3D jibben_solver;
          (*a_interface)(i, j, k) = jibben_solver.solve2(&neighborhood);
          // IRL::Jibben_3D jibben_solver;
          // (*a_interface)(i, j, k) = jibben_solver.solve2(&neighborhood);

          // error metrics
          // double err = jibben_solver.getVolumeErrorSquared(mesh.dx());
          double normal_err = jibben_solver.getNormalMetric();
          double normal_eigen_err = jibben_solver.getNormalEigenMetric();
          double normal_std_err = jibben_solver.getNormalVarianceMetric();

          // Match to volume fraction
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::setDistanceToMatchVolumeFraction(
              cell, liquid_volume_fraction, &(*a_interface)(i, j, k), 1.0e-14);

          if (const auto ptr =
                  std::get_if<IRL::Paraboloid>(&(*a_interface)(i, j, k))) {
            normal_error.paraboloid_scalar_data(i, j, k) = normal_err;
            normal_eigen_error.paraboloid_scalar_data(i, j, k) =
                normal_eigen_err;
            normal_std_error.paraboloid_scalar_data(i, j, k) = normal_std_err;
            // squared_vol_error.paraboloid_scalar_data(i, j, k) = err;
          } else if (const auto ptr = std::get_if<IRL::PlanarSeparator>(
                         &(*a_interface)(i, j, k))) {
            normal_error.polygon_scalar_data(i, j, k) = normal_err;
            normal_eigen_error.polygon_scalar_data(i, j, k) = normal_eigen_err;
            normal_std_error.polygon_scalar_data(i, j, k) = normal_std_err;
            // squared_vol_error.polygon_scalar_data(i, j, k) = err;
          }
        }
      }
    }
  }
  // a_scalar_fields->push_back(squared_vol_error);
  a_scalar_fields->push_back(normal_error);
  a_scalar_fields->push_back(normal_eigen_error);
  a_scalar_fields->push_back(normal_std_error);

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
  const int itmax = 50;
  for (int i = 0; i < itmax; i++) {
    const auto F_and_gradF =
        getPUAndGrad(a_liq_moments, a_interface, a_centroid, a_area, a_nlayers,
                     a_delta, a_i, a_j, a_k, projected_pt);
    const double F = std::get<double>(F_and_gradF);
    if (F < a_delta * 1.e-6) {
      break;
    }
    // if ((i + 1) == itmax)
    //   std::cout << "projection incomplete F = " << F << std::endl;
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
                           std::vector<InterfaceScalarField>* a_scalar_fields,
                           const bool a_plic_already_built) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
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
                            &jibben_reconstruction, a_scalar_fields, true,
                            &interface_centroids, &interface_areas,
                            &jibben_errors);

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

  const int nlayers = 3;
  // const double delta = 5.0 * mesh.dx();
  const double delta = 2.5 * mesh.dx();
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
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields) {
  // First, we need to build the plic and copy it
  LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                           a_interface, a_scalar_fields);
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
                            a_interface, a_scalar_fields, true);

  // Choose between PLIC and Jibben
  const double vfrac_threshold = 1.0e-4;
  const double kdx_threshold = 4.0;
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        // vfrac check
        if (liquid_volume_fraction < vfrac_threshold ||
            liquid_volume_fraction > 1.0 - vfrac_threshold) {
          (*a_interface)(i, j, k) = plic_reconstruction(i, j, k);
          continue;
        }
        // curvature check
        if (IRL::Paraboloid* paraboloid =
                std::get_if<IRL::Paraboloid>(&(*a_interface)(i, j, k))) {
          const auto& aligned = paraboloid->getAlignedParaboloid();
          const double a = aligned.a();
          const double b = aligned.b();
          if (std::abs(a) * mesh.dx() > kdx_threshold ||
              std::abs(b) * mesh.dx() > kdx_threshold) {
            (*a_interface)(i, j, k) = plic_reconstruction(i, j, k);
          }
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
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
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

// Helping functions for Taubin fit -----------------------------------------

// vf weight
double getVfWeight(const double& a_vfrac) {
  const double limit_vfrac = 0.05;
  if (a_vfrac < limit_vfrac) {
    return 0.5 - 0.5 * std::cos(M_PI * a_vfrac / limit_vfrac);
  } else if (a_vfrac > (1.0 - limit_vfrac)) {
    return 0.5 - 0.5 * std::cos(M_PI * (1.0 - a_vfrac) / limit_vfrac);
  } else {
    return 1.0;
  }
}

// normal weight
double getNormalWeight(const IRL::Normal& a_nref, const IRL::Normal& a_nloc) {
  const double n_dot =
      a_nref[0] * a_nloc[0] + a_nref[1] * a_nloc[1] + a_nref[2] * a_nloc[2];
  return std::max(0.0, n_dot);
}

// distance weight
double getDistanceWeight(const IRL::Pt& a_pref, const IRL::Pt& a_ploc,
                         const double& h) {
  const double distance = IRL::magnitude(a_ploc - a_pref);
  const double distance_ndim = distance / (2.5 * h);
  const double distance_weight =
      distance_ndim >= 1.0
          ? 0.0
          : (1.0 + 4.0 * distance_ndim) * std::pow(1.0 - distance_ndim, 4.0);
  return distance_weight;
}

// total weight
double getTotalWeight(const double& a_vfrac, const IRL::Normal& a_nref,
                      const IRL::Normal& a_nloc, const IRL::Pt& a_pref,
                      const IRL::Pt& a_ploc, const double& h) {
  double vf_weight = getVfWeight(a_vfrac);
  double normal_weight = getNormalWeight(a_nref, a_nloc);
  double distance_weight = getDistanceWeight(a_pref, a_ploc, h);
  return vf_weight * normal_weight * distance_weight;
}

// sampling points in local frame
void sampleLocalPoints(
    const IRL::ReferenceFrame& local_frame, const IRL::Pt local_datum,
    const std::vector<std::pair<IRL::Pt, IRL::Pt>>& end_points_list,
    const int& num_samples_per_segment,
    std::vector<Eigen::Vector2d>* points_local_frame) {
  // rotation matrix and origin
  Eigen::Vector3d e1(local_frame[0][0], local_frame[0][1], local_frame[0][2]);
  Eigen::Vector3d e2(local_frame[2][0], local_frame[2][1], local_frame[2][2]);
  Eigen::Vector3d e3(local_frame[1][0], local_frame[1][1], local_frame[1][2]);
  Eigen::Matrix3d R;
  R.col(0) = e1;
  R.col(1) = e2;
  R.col(2) = e3;
  Eigen::Vector3d o(local_datum[0], local_datum[1], local_datum[2]);

  // sampling points
  for (int i = 0; i < end_points_list.size(); i++) {
    // starting and ending points
    Eigen::Vector3d start(end_points_list[i].first[0],
                          end_points_list[i].first[1],
                          end_points_list[i].first[2]);
    Eigen::Vector3d end(end_points_list[i].second[0],
                        end_points_list[i].second[1],
                        end_points_list[i].second[2]);

    // start and end in local rotated frame
    Eigen::Vector3d start_local = R.transpose() * (start - o);
    Eigen::Vector3d end_local = R.transpose() * (end - o);

    // sampling points
    for (int j = 0; j < num_samples_per_segment; j++) {
      double t = static_cast<double>(j) /
                 (static_cast<double>(num_samples_per_segment) - 1);
      Eigen::Vector3d pt_local = start_local * (1.0 - t) + end_local * t;
      Eigen::Vector2d pt_local_2d(pt_local[0], pt_local[1]);
      points_local_frame->push_back(pt_local_2d);
    }
  }
}

// taubin circle fit using points to get curvature
double getSignedTaubinCurvature(const std::vector<Eigen::Vector2d>& points,
                                const std::vector<double>& weights) {
  // number of points
  const int n = points.size();  // in local frame
  if (n < 3)
    return std::numeric_limits<double>::quiet_NaN();  // not enough points

  // moment matrix
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    double w = weights[i];
    // w = 1.0;
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
  C(1, 1) = static_cast<double>(n);
  C(2, 0) = C(0, 2);
  C(2, 2) = static_cast<double>(n);

  // collinearity check for safety (again)
  const double Mxx = M(1, 1) / static_cast<double>(n);
  const double Myy = M(2, 2) / static_cast<double>(n);
  const double Mxy = M(1, 2) / static_cast<double>(n);
  const double cov_det = Mxx * Myy - Mxy * Mxy;
  if (std::abs(cov_det) < 1e-14) {
    return 0.0;  // zero curvature if collinear
  }

  // solving generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);

  // eigen values/vectors
  const auto& evals = ges.eigenvalues();
  const auto& evecs = ges.eigenvectors();

  // need to find eigenvector with smallest non-negative real eigenvalue
  int bestIndex = -1;
  double bestLambda = std::numeric_limits<double>::infinity();
  for (int i = 0; i < 4; i++) {
    const auto lam = evals(i);
    if (std::abs(lam.imag()) > 1e-9) continue;
    const double lambdaReal = lam.real();
    if (lambdaReal <= 0.0) continue;
    if (lambdaReal < bestLambda) {
      bestLambda = lambdaReal;
      bestIndex = i;
    }
  }
  if (bestIndex < 0) {  // no eigenvalue found
    return std::numeric_limits<double>::quiet_NaN();
  }

  // extracting eigenvector components
  Eigen::Vector4cd v_c = evecs.col(bestIndex);
  Eigen::Vector4d a = v_c.real();  // imaginary parts should be tiny
  double A = a(0);
  double B = a(1);
  double Cc = a(2);
  double D = a(3);

  if (std::abs(A) < 1e-14) {
    return 0.0;  // very small A is infinite radius again
  }

  // finding radius of curvature
  const double num = B * B + Cc * Cc - 4.0 * A * D;
  if (num <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double R = 0.5 * std::sqrt(num) / std::abs(A);

  // signed curvature
  double yc_local = -Cc / (2.0 * A);
  double sign = (yc_local >= 0.0) ? -1.0 : +1.0;
  return sign * (1.0 / R);
}

// solving for principal directions and curvatures
std::tuple<double, double, double, bool> getPrincipalCurvatures(
    const std::vector<double>& theta, const std::vector<double>& k_theta) {
  const std::size_t m = theta.size();
  if (m < 3 || m != k_theta.size()) {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(), false};
  }

  // building least squares system
  Eigen::MatrixXd X(m, 3);
  Eigen::VectorXd y(m);

  for (std::size_t s = 0; s < m; ++s) {
    const double th = theta[s];
    const double ks = k_theta[s];

    X(s, 0) = 1.0;
    X(s, 1) = std::cos(2.0 * th);
    X(s, 2) = std::sin(2.0 * th);
    y(s) = ks;
  }

  // solving minimization problem
  Eigen::Vector3d a = X.colPivHouseholderQr().solve(y);
  const double res_norm = (X * a - y).norm();
  if (std::isnan(res_norm)) {  // NaN check
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(), false};
  }

  const double a0 = a(0);
  const double a1 = a(1);
  const double a2 = a(2);

  const double Delta = std::sqrt(a1 * a1 + a2 * a2);
  double k1, k2, phi;

  const double eps = 1e-14;
  if (Delta < eps) {
    k1 = a0;
    k2 = a0;
    phi = 0.0;  // arbitrary
    return {k1, k2, phi, true};
  }

  k1 = a0 + Delta;
  k2 = a0 - Delta;
  phi = 0.5 * std::atan2(a2, a1);

  return {k1, k2, phi, true};
}

// Taubin's approximation to gradient weighted algrabric fit
void SlicesTaubin::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // PLIC
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  // clipped PLIC polygons
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::Polygon> polygon(&mesh);
  Data<double> vf(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        vf(i, j, k) = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vf(i, j, k) < IRL::global_constants::VF_LOW ||
            vf(i, j, k) > IRL::global_constants::VF_HIGH) {
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

  // slicing and circle fit params
  const int nsamples_per_segment = 10;
  const int nslices = 18;
  const int nlayers = 2;

  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);
  IRL::StackVector<IRL::Pt, 2> intersections;

  // computing circle fit for each mixed cell (or each cell with polygon)
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        if (polygon(i, j, k).getNumberOfVertices() <= 2) continue;

        // polygon local frame
        const IRL::Normal polygon_normal =
            calculatePolygonNormal(polygon(i, j, k));
        IRL::ReferenceFrame polygon_frame =
            referenceFrameFromNormal(polygon_normal);
        const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();

        // stencil polygon and volume fraction data
        polygon_vfrac_list.clear();
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              if (polygon(ii, jj, kk).getNumberOfVertices() <= 2) continue;
              polygon_vfrac_list.emplace_back(polygon(ii, jj, kk),
                                              vf(ii, jj, kk));
            }
          }
        }
        const int num_polygons = polygon_vfrac_list.size();

        // revert back to plane if there are no polygons in the neighboring
        // stencil. Circle fit will fail if all points are collinear
        if (num_polygons < 2) continue;  // back to LVIRA

        // slicing about local normal [0, pi)
        std::vector<double> theta_list, k_theta_list;
        for (int s = 0; s < nslices; s++) {
          // rotation angle
          double theta =
              M_PI * static_cast<double>(s) / static_cast<double>(nslices);

          // rotating local polygon frame
          const IRL::UnitQuaternion rotation(theta, polygon_frame[2]);
          const auto rotated_polygon_frame = rotation * polygon_frame;

          // slicing plane
          const IRL::Plane slicing_plane(
              rotated_polygon_frame[1],
              rotated_polygon_frame[1] * polygon_centroid);

          // slicing polygon stencial with plane and collecting info
          std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
          std::vector<double> weights;
          for (int p = 0; p < num_polygons; p++) {
            getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                               &intersections);
            if (intersections.size() != 2) continue;
            // end points for sampling
            IRL::Pt start_point = intersections[0];
            IRL::Pt end_point = intersections[1];
            end_points_list.push_back({start_point, end_point});
            // computing weights
            IRL::Normal neighbor_normal =
                calculatePolygonNormal(polygon_vfrac_list[p].first);
            IRL::Normal neighbor_centroid =
                (polygon_vfrac_list[p].first).calculateCentroid();
            double neighbor_weight =
                getTotalWeight(polygon_vfrac_list[p].second,
                               rotated_polygon_frame[2], neighbor_normal,
                               polygon_centroid, neighbor_centroid, mesh.dx());
            weights.insert(weights.end(), nsamples_per_segment,
                           neighbor_weight);  // all points sampled on a
                                              // segment have the same weight
          }

          // getting points in local frame
          const int num_points = end_points_list.size() * nsamples_per_segment;
          std::vector<Eigen::Vector2d> points_rotated_frame;
          points_rotated_frame.reserve(num_points);
          sampleLocalPoints(rotated_polygon_frame, polygon_centroid,
                            end_points_list, nsamples_per_segment,
                            &points_rotated_frame);

          // taubin circle fit using points in local frame for curvature
          double k_theta =
              getSignedTaubinCurvature(points_rotated_frame, weights);
          if (std::isfinite(k_theta)) {
            theta_list.push_back(theta);
            k_theta_list.push_back(k_theta);
          }
        }
        // principal curvatures and directions
        auto [k1, k2, phi, ok] =
            getPrincipalCurvatures(theta_list, k_theta_list);
        if (!ok) continue;  // plane

        // Darboux frame
        const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
        const IRL::ReferenceFrame darboux_frame = rotate_phi * polygon_frame;

        // new paraboloid
        IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 * k1,
                                   0.5 * k2);

        // volume fraction matching
        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
            solver_distance(cell, vf(i, j, k), 1.0e-14, paraboloid);
        paraboloid.setDatum(
            IRL::Pt(polygon_centroid +
                    solver_distance.getDistance() * darboux_frame[2]));
        (*a_interface)(i, j, k) = paraboloid;
      }
    }
  }
}

// taubin fit using PU to get normal for cutting plane
void SlicesTaubin2::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // PLIC
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  // clipped PLIC polygons
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::Polygon> polygon(&mesh);
  Data<IRL::Pt> polygon_centroids(&mesh);
  Data<double> polygon_areas(&mesh);
  Data<double> vf(&mesh);
  Data<IRL::SeparatorVariant> plic_reconstruction(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        plic_reconstruction(i, j, k) = (*a_interface)(i, j, k);
        vf(i, j, k) = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vf(i, j, k) < IRL::global_constants::VF_LOW ||
            vf(i, j, k) > IRL::global_constants::VF_HIGH) {
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, planar_separator, planar_separator[0]);
        polygon_centroids(i, j, k) = polygon(i, j, k).calculateCentroid();
        polygon_areas(i, j, k) = std::abs(polygon(i, j, k).calculateVolume());
      }
    }
  }
  updatePolygonBorder(&polygon);

  // storing normals obtained from PU for the orientation of the cutting plane
  Data<IRL::Normal> pu_normals(&mesh);
  const int pu_layers = 2;
  const double pu_delta = 2.5 * mesh.dx();
  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        if (vf(i, j, k) < IRL::global_constants::VF_LOW ||
            vf(i, j, k) > IRL::global_constants::VF_HIGH)
          continue;
        const auto F_and_gradF =
            getPUAndGrad(a_liq_moments, plic_reconstruction, polygon_centroids,
                         polygon_areas, pu_layers, pu_delta, i, j, k,
                         polygon_centroids(i, j, k));
        Eigen::Vector3d pu_normal = F_and_gradF.second;
        pu_normals(i, j, k) =
            IRL::Normal(pu_normal[0], pu_normal[1], pu_normal[2]);
      }
    }
  }

  // slicing and circle fit params
  const int nsamples_per_segment = 10;
  const int nslices = 18;
  const int nlayers = 2;

  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);
  IRL::StackVector<IRL::Pt, 2> intersections;

  // computing circle fit for each mixed cell (or each cell with polygon)
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        if (polygon(i, j, k).getNumberOfVertices() <= 2) continue;

        // polygon local frame
        // const IRL::Normal polygon_normal =
        //     calculatePolygonNormal(polygon(i, j, k));
        const IRL::Normal polygon_normal = pu_normals(i, j, k);
        IRL::ReferenceFrame polygon_frame =
            referenceFrameFromNormal(polygon_normal);
        const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();

        // stencil polygon and volume fraction data
        polygon_vfrac_list.clear();
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              if (polygon(ii, jj, kk).getNumberOfVertices() <= 2) continue;
              polygon_vfrac_list.emplace_back(polygon(ii, jj, kk),
                                              vf(ii, jj, kk));
            }
          }
        }
        const int num_polygons = polygon_vfrac_list.size();

        // revert back to plane if there are no polygons in the neighboring
        // stencil. Circle fit will fail if all points are collinear
        if (num_polygons < 2) continue;  // back to LVIRA

        // slicing about local normal [0, pi)
        std::vector<double> theta_list, k_theta_list;
        for (int s = 0; s < nslices; s++) {
          // rotation angle
          double theta =
              M_PI * static_cast<double>(s) / static_cast<double>(nslices);

          // rotating local polygon frame
          const IRL::UnitQuaternion rotation(theta, polygon_frame[2]);
          const auto rotated_polygon_frame = rotation * polygon_frame;

          // slicing plane
          const IRL::Plane slicing_plane(
              rotated_polygon_frame[1],
              rotated_polygon_frame[1] * polygon_centroid);

          // slicing polygon stencial with plane and collecting info
          std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
          std::vector<double> weights;
          for (int p = 0; p < num_polygons; p++) {
            getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                               &intersections);
            if (intersections.size() != 2) continue;
            // end points for sampling
            IRL::Pt start_point = intersections[0];
            IRL::Pt end_point = intersections[1];
            end_points_list.push_back({start_point, end_point});
            // computing weights
            IRL::Normal neighbor_normal =
                calculatePolygonNormal(polygon_vfrac_list[p].first);
            IRL::Normal neighbor_centroid =
                (polygon_vfrac_list[p].first).calculateCentroid();
            double neighbor_weight =
                getTotalWeight(polygon_vfrac_list[p].second,
                               rotated_polygon_frame[2], neighbor_normal,
                               polygon_centroid, neighbor_centroid, mesh.dx());
            weights.insert(weights.end(), nsamples_per_segment,
                           neighbor_weight);  // all points sampled on a
                                              // segment have the same weight
          }

          // getting points in local frame
          const int num_points = end_points_list.size() * nsamples_per_segment;
          std::vector<Eigen::Vector2d> points_rotated_frame;
          points_rotated_frame.reserve(num_points);
          sampleLocalPoints(rotated_polygon_frame, polygon_centroid,
                            end_points_list, nsamples_per_segment,
                            &points_rotated_frame);

          // taubin circle fit using points in local frame for curvature
          double k_theta =
              getSignedTaubinCurvature(points_rotated_frame, weights);
          if (std::isfinite(k_theta)) {
            theta_list.push_back(theta);
            k_theta_list.push_back(k_theta);
          }
        }
        // principal curvatures and directions
        auto [k1, k2, phi, ok] =
            getPrincipalCurvatures(theta_list, k_theta_list);
        if (!ok) continue;  // plane

        // Darboux frame
        const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
        const IRL::ReferenceFrame darboux_frame = rotate_phi * polygon_frame;

        // new paraboloid
        IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 * k1,
                                   0.5 * k2);

        // volume fraction matching
        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
            solver_distance(cell, vf(i, j, k), 1.0e-14, paraboloid);
        paraboloid.setDatum(
            IRL::Pt(polygon_centroid +
                    solver_distance.getDistance() * darboux_frame[2]));
        (*a_interface)(i, j, k) = paraboloid;
      }
    }
  }
}

void SlicesTaubin3::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // PLIC
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  // clipped PLIC polygons
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::Polygon> polygon(&mesh);
  Data<IRL::Pt> polygon_centroids(&mesh);
  Data<double> polygon_areas(&mesh);
  Data<double> vf(&mesh);
  Data<IRL::SeparatorVariant> plic_reconstruction(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        plic_reconstruction(i, j, k) = (*a_interface)(i, j, k);
        vf(i, j, k) = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vf(i, j, k) < IRL::global_constants::VF_LOW ||
            vf(i, j, k) > IRL::global_constants::VF_HIGH) {
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, planar_separator, planar_separator[0]);
        polygon_centroids(i, j, k) = polygon(i, j, k).calculateCentroid();
        polygon_areas(i, j, k) = std::abs(polygon(i, j, k).calculateVolume());
      }
    }
  }
  updatePolygonBorder(&polygon);

  // getting normal based on eigenvector
  Data<IRL::Normal> new_normals(&mesh);
  const int layers = 2;
  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        if (vf(i, j, k) < IRL::global_constants::VF_LOW ||
            vf(i, j, k) > IRL::global_constants::VF_HIGH)
          continue;
        Eigen::Vector3d o(polygon_centroids(i, j, k)[0],
                          polygon_centroids(i, j, k)[1],
                          polygon_centroids(i, j, k)[2]);
        IRL::ReferenceFrame local_frame = referenceFrameFromNormal(
            polygon(i, j, k).getPlaneOfExistence().normal());
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
        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
        IRL::UnsignedIndex_t count = 0;
        for (int kk = k - layers; kk <= k + layers; ++kk) {
          for (int jj = j - layers; jj <= j + layers; ++jj) {
            for (int ii = i - layers; ii <= i + layers; ++ii) {
              if (polygon(ii, jj, kk).getNumberOfVertices() <= 2) continue;
              count++;
              IRL::Normal normal =
                  polygon(ii, jj, kk).getPlaneOfExistence().normal();
              Eigen::Vector3d n(normal[0], normal[1], normal[2]);
              Eigen::Vector3d n_local = R.transpose() * n;
              C += n_local * n_local.transpose();
            }
          }
        }
        if (count > 0) {
          C /= static_cast<double>(count);
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
          // eigenvector corresponding to largest eigenvalue
          Eigen::Vector3d evec_loc = es.eigenvectors().col(2);
          // converting to global frame
          Eigen::Vector3d evec_glob = R * evec_loc;
          // checking orientation with plic normal otherwise flipping
          IRL::Normal plic_normal =
              polygon(i, j, k).getPlaneOfExistence().normal();
          new_normals(i, j, k) =
              IRL::Normal(evec_glob[0], evec_glob[1], evec_glob[2]);
          if (new_normals(i, j, k) * plic_normal < 0.0) {
            new_normals(i, j, k) =
                IRL::Normal(-evec_glob[0], -evec_glob[1], -evec_glob[2]);
          }
        } else {
          new_normals(i, j, k) = IRL::Normal(
              local_frame[2][0], local_frame[2][1], local_frame[2][2]);
        }
      }
    }
  }

  // slicing and circle fit params
  const int nsamples_per_segment = 10;
  const int nslices = 18;
  const int nlayers = 2;

  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);
  IRL::StackVector<IRL::Pt, 2> intersections;

  // computing circle fit for each mixed cell (or each cell with polygon)
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        if (polygon(i, j, k).getNumberOfVertices() <= 2) continue;

        // polygon local frame
        // const IRL::Normal polygon_normal =
        //     calculatePolygonNormal(polygon(i, j, k));
        const IRL::Normal polygon_normal = new_normals(i, j, k);
        IRL::ReferenceFrame polygon_frame =
            referenceFrameFromNormal(polygon_normal);
        const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();

        // stencil polygon and volume fraction data
        polygon_vfrac_list.clear();
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              if (polygon(ii, jj, kk).getNumberOfVertices() <= 2) continue;
              polygon_vfrac_list.emplace_back(polygon(ii, jj, kk),
                                              vf(ii, jj, kk));
            }
          }
        }
        const int num_polygons = polygon_vfrac_list.size();

        // revert back to plane if there are no polygons in the neighboring
        // stencil. Circle fit will fail if all points are collinear
        if (num_polygons < 2) continue;  // back to LVIRA

        // slicing about local normal [0, pi)
        std::vector<double> theta_list, k_theta_list;
        for (int s = 0; s < nslices; s++) {
          // rotation angle
          double theta =
              M_PI * static_cast<double>(s) / static_cast<double>(nslices);

          // rotating local polygon frame
          const IRL::UnitQuaternion rotation(theta, polygon_frame[2]);
          const auto rotated_polygon_frame = rotation * polygon_frame;

          // slicing plane
          const IRL::Plane slicing_plane(
              rotated_polygon_frame[1],
              rotated_polygon_frame[1] * polygon_centroid);

          // slicing polygon stencial with plane and collecting info
          std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
          std::vector<double> weights;
          for (int p = 0; p < num_polygons; p++) {
            getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                               &intersections);
            if (intersections.size() != 2) continue;
            // end points for sampling
            IRL::Pt start_point = intersections[0];
            IRL::Pt end_point = intersections[1];
            end_points_list.push_back({start_point, end_point});
            // computing weights
            IRL::Normal neighbor_normal =
                calculatePolygonNormal(polygon_vfrac_list[p].first);
            IRL::Normal neighbor_centroid =
                (polygon_vfrac_list[p].first).calculateCentroid();
            double neighbor_weight =
                getTotalWeight(polygon_vfrac_list[p].second,
                               rotated_polygon_frame[2], neighbor_normal,
                               polygon_centroid, neighbor_centroid, mesh.dx());
            weights.insert(weights.end(), nsamples_per_segment,
                           neighbor_weight);  // all points sampled on a
                                              // segment have the same weight
          }

          // getting points in local frame
          const int num_points = end_points_list.size() * nsamples_per_segment;
          std::vector<Eigen::Vector2d> points_rotated_frame;
          points_rotated_frame.reserve(num_points);
          sampleLocalPoints(rotated_polygon_frame, polygon_centroid,
                            end_points_list, nsamples_per_segment,
                            &points_rotated_frame);

          // taubin circle fit using points in local frame for curvature
          double k_theta =
              getSignedTaubinCurvature(points_rotated_frame, weights);
          if (std::isfinite(k_theta)) {
            theta_list.push_back(theta);
            k_theta_list.push_back(k_theta);
          }
        }
        // principal curvatures and directions
        auto [k1, k2, phi, ok] =
            getPrincipalCurvatures(theta_list, k_theta_list);
        if (!ok) continue;  // plane

        // Darboux frame
        const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
        const IRL::ReferenceFrame darboux_frame = rotate_phi * polygon_frame;

        // new paraboloid
        IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 * k1,
                                   0.5 * k2);

        // volume fraction matching
        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
            solver_distance(cell, vf(i, j, k), 1.0e-14, paraboloid);
        paraboloid.setDatum(
            IRL::Pt(polygon_centroid +
                    solver_distance.getDistance() * darboux_frame[2]));
        (*a_interface)(i, j, k) = paraboloid;
      }
    }
  }
}

// ================ taubin line segment fit =================

// moment matrix terms
std::vector<double> getTaubinMomentTerms(const double& x0, const double& y0,
                                         const double& dx, const double& dy) {
  std::vector<double> terms(10, 0.0);

  // length of line segment
  double Le = std::sqrt(dx * dx + dy * dy);
  // if (Le < 1e-14) return terms;  // degenerate segment

  // Mzz
  terms[0] = (dx * dx * dx * dx) / 5. + (2. * (dx * dx) * (dy * dy)) / 5. +
             (dy * dy * dy * dy) / 5. + dx * dx * dx * x0 +
             dx * (dy * dy) * x0 + 2. * (dx * dx) * (x0 * x0) +
             (2. * (dy * dy) * (x0 * x0)) / 3. + 2. * dx * (x0 * x0 * x0) +
             x0 * x0 * x0 * x0 + dx * dx * dy * y0 + dy * dy * dy * y0 +
             (8. * dx * dy * x0 * y0) / 3. + 2. * dy * (x0 * x0) * y0 +
             (2. * (dx * dx) * (y0 * y0)) / 3. + 2. * (dy * dy) * (y0 * y0) +
             2. * dx * x0 * (y0 * y0) + 2. * (x0 * x0) * (y0 * y0) +
             2. * dy * (y0 * y0 * y0) + y0 * y0 * y0 * y0;

  // Mxz
  terms[1] = (dx * dx * dx) / 4. + (dx * (dy * dy)) / 4. + dx * dx * x0 +
             (dy * dy * x0) / 3. + (3. * dx * (x0 * x0)) / 2. + x0 * x0 * x0 +
             (2. * dx * dy * y0) / 3. + dy * x0 * y0 + (dx * (y0 * y0)) / 2. +
             x0 * (y0 * y0);

  // Myz
  terms[2] = (dx * dx * dy) / 4. + (dy * dy * dy) / 4. +
             (2. * dx * dy * x0) / 3. + (dy * (x0 * x0)) / 2. +
             (dx * dx * y0) / 3. + dy * dy * y0 + dx * x0 * y0 + x0 * x0 * y0 +
             (3. * dy * (y0 * y0)) / 2. + y0 * y0 * y0;

  // Mz
  terms[3] =
      (dx * dx) / 3. + (dy * dy) / 3. + dx * x0 + x0 * x0 + dy * y0 + y0 * y0;

  // Mxx
  terms[4] = (dx * dx) / 3. + dx * x0 + x0 * x0;

  // Mxy
  terms[5] = (dx * dy) / 3. + (dy * x0) / 2. + (dx * y0) / 2. + x0 * y0;

  // Mx
  terms[6] = dx / 2. + x0;

  // Myy
  terms[7] = (dy * dy) / 3. + dy * y0 + y0 * y0;

  // My
  terms[8] = dy / 2. + y0;

  // Le
  terms[9] = 1.0;

  // for (auto& term : terms) {
  //   term *= Le;
  // }

  return terms;
}

// moment matrix terms using GL quadrature
std::vector<double> getTaubinMomentTermsGL(const double& x0, const double& y0,
                                           const double& dx, const double& dy) {
  std::vector<double> terms(10, 0.0);

  // quadrature params
  const int nGL = 3;
  const int a = 0.0;
  const int b = 1.0;

  // functions for integration
  auto f_x = [&](const double& t) { return x0 + dx * t; };
  auto f_y = [&](const double& t) { return y0 + dy * t; };
  auto f_z = [&](const double& t) {
    double x = f_x(t);
    double y = f_y(t);
    return x * x + y * y;
  };

  // first order moments
  auto f_Mx = [&](const double& t) { return f_x(t); };
  auto f_My = [&](const double& t) { return f_y(t); };
  auto f_Mz = [&](const double& t) { return f_z(t); };

  // second order moments
  auto f_Mxx = [&](const double& t) {
    double x = f_x(t);
    return x * x;
  };
  auto f_Myy = [&](const double& t) {
    double y = f_y(t);
    return y * y;
  };
  auto f_Mzz = [&](const double& t) {
    double z = f_z(t);
    return z * z;
  };
  auto f_Mxy = [&](const double& t) {
    double x = f_x(t);
    double y = f_y(t);
    return x * y;
  };
  auto f_Mxz = [&](const double& t) {
    double x = f_x(t);
    double z = f_z(t);
    return x * z;
  };
  auto f_Myz = [&](const double& t) {
    double y = f_y(t);
    double z = f_z(t);
    return y * z;
  };

  // GL integrator
  IRL::GaussLegendreIntegrator<double, 1> integrator(nGL);

  // Mzz
  terms[0] = integrator.integrate(f_Mzz, a, b);

  // Mxz
  terms[1] = integrator.integrate(f_Mxz, a, b);

  // Myz
  terms[2] = integrator.integrate(f_Myz, a, b);

  // Mz
  terms[3] = integrator.integrate(f_Mz, a, b);

  // Mxx
  terms[4] = integrator.integrate(f_Mxx, a, b);

  // Mxy
  terms[5] = integrator.integrate(f_Mxy, a, b);

  // Mx
  terms[6] = integrator.integrate(f_Mx, a, b);

  // Myy
  terms[7] = integrator.integrate(f_Myy, a, b);

  // My
  terms[8] = integrator.integrate(f_My, a, b);

  // Le
  terms[9] = 1.0;

  return terms;
}

// building taubin moment and contraint matrix
void getTaubinMatrices(
    const std::vector<std::pair<IRL::Pt, IRL::Pt>>& end_points,
    const std::vector<double>& weights, Eigen::Matrix4d* M,
    Eigen::Matrix4d* C) {
  M->setZero();                     // moment matrix
  C->setZero();                     // constraint matrix
  const int n = end_points.size();  // number of line segments

  // building moment matrix
  for (std::size_t i = 0; i < n; i++) {
    // line segment end points
    const IRL::Pt x0 = end_points[i].first;
    const IRL::Pt x1 = end_points[i].second;

    const double dx = x1[0] - x0[0];
    const double dy = x1[1] - x0[1];

    // moment terms
    std::vector<double> terms = getTaubinMomentTerms(
        x0[0], x0[1], dx,
        dy);  // getTaubinMomentTermsGL for eval using Gauss-Legendre
    // 0 -> Mzz, 1 -> Mxz, 2 -> Myz, 3 -> Mz, 4 -> Mxx, 5 -> Mxy,
    // 6 -> Mx, 7 -> Myy, 8 -> My, 9 -> Le

    // weight
    const double w = weights[i];

    // updating moment matrix
    // M->operator()(0, 0) += w * terms[0];  // Mzz
    // M->operator()(0, 1) += w * terms[1];  // Mxz
    // M->operator()(0, 2) += w * terms[2];  // Myz
    // M->operator()(0, 3) += w * terms[3];  // Mz
    // M->operator()(1, 0) += w * terms[1];  // Mxz
    // M->operator()(1, 1) += w * terms[4];  // Mxx
    // M->operator()(1, 2) += w * terms[5];  // Mxy
    // M->operator()(1, 3) += w * terms[6];  // Mx
    // M->operator()(2, 0) += w * terms[2];  // Myz
    // M->operator()(2, 1) += w * terms[5];  // Mxy
    // M->operator()(2, 2) += w * terms[7];  // Myy
    // M->operator()(2, 3) += w * terms[8];  // My
    // M->operator()(3, 0) += w * terms[3];  // Mz
    // M->operator()(3, 1) += w * terms[6];  // Mx
    // M->operator()(3, 2) += w * terms[8];  // My
    // M->operator()(3, 3) += w * terms[9];  // Le

    double Mzz = terms[0], Mxz = terms[1], Myz = terms[2], Mz = terms[3],
           Mxx = terms[4], Mxy = terms[5], Mx = terms[6], Myy = terms[7],
           My = terms[8];

    M->operator()(0, 0) += Mzz;
    M->operator()(0, 1) += Mxz;
    M->operator()(0, 2) += Myz;
    M->operator()(0, 3) += Mz;
    M->operator()(1, 0) += Mxz;
    M->operator()(1, 1) += Mxx;
    M->operator()(1, 2) += Mxy;
    M->operator()(1, 3) += Mx;
    M->operator()(2, 0) += Myz;
    M->operator()(2, 1) += Mxy;
    M->operator()(2, 2) += Myy;
    M->operator()(2, 3) += My;
    M->operator()(3, 0) += Mz;
    M->operator()(3, 1) += Mx;
    M->operator()(3, 2) += My;
    M->operator()(3, 3) += terms[9];
  }

  // building constraint matrix
  C->operator()(0, 0) = 4.0 * M->operator()(0, 3);
  C->operator()(0, 1) = 2.0 * M->operator()(1, 3);
  C->operator()(0, 2) = 2.0 * M->operator()(2, 3);
  C->operator()(1, 0) = C->operator()(0, 1);
  C->operator()(1, 1) = M->operator()(3, 3);
  C->operator()(2, 0) = C->operator()(0, 2);
  C->operator()(2, 2) = M->operator()(3, 3);
}

// signed curvature
double getSignedTaubinCurvatureLS(
    const std::vector<std::pair<IRL::Pt, IRL::Pt>>& end_points,
    const std::vector<double>& weights) {
  // collinearity check
  const int n = end_points.size();
  if (n < 2) return 0.0;  // zero curvature if only one line segment

  // moment and constraint matrices
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  Eigen::Matrix4d C = Eigen::Matrix4d::Zero();
  getTaubinMatrices(end_points, weights, &M, &C);

  // solving generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);

  // eigen values/vectors
  const auto& evals = ges.eigenvalues();
  const auto& evecs = ges.eigenvectors();

  // need to find eigenvector with smallest non-negative real eigenvalue
  int bestIndex = -1;
  double bestLambda = std::numeric_limits<double>::infinity();
  for (int i = 0; i < 4; i++) {
    const auto lam = evals(i);
    if (std::abs(lam.imag()) > 1e-9) continue;
    const double lambdaReal = lam.real();
    if (lambdaReal <= 0.0) continue;
    if (lambdaReal < bestLambda) {
      bestLambda = lambdaReal;
      bestIndex = i;
    }
  }
  if (bestIndex < 0) {  // no eigenvalue found
    return std::numeric_limits<double>::quiet_NaN();
  }

  // extracting eigenvector components
  Eigen::Vector4cd v_c = evecs.col(bestIndex);
  Eigen::Vector4d a = v_c.real();  // imaginary parts should be tiny
  double A = a(0);
  double B = a(1);
  double Cc = a(2);
  double D = a(3);

  if (std::abs(A) < 1e-14) {
    return 0.0;  // very small A is infinite radius again
  }

  // finding radius of curvature
  const double num = B * B + Cc * Cc - 4.0 * A * D;
  if (num <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double R = 0.5 * std::sqrt(num) / std::abs(A);

  // signed curvature
  double yc_local = -Cc / (2.0 * A);
  double sign = (yc_local >= 0.0) ? -1.0 : +1.0;
  return sign * (1.0 / R);
}

// to local
IRL::Pt toLocal(const IRL::ReferenceFrame& local_frame,
                const IRL::Pt& local_datum, const IRL::Pt& global_pt) {
  // rotation matrix and origin
  Eigen::Vector3d e1(local_frame[0][0], local_frame[0][1], local_frame[0][2]);
  Eigen::Vector3d e2(local_frame[2][0], local_frame[2][1], local_frame[2][2]);
  Eigen::Vector3d e3(local_frame[1][0], local_frame[1][1], local_frame[1][2]);
  Eigen::Matrix3d R;
  R.col(0) = e1;
  R.col(1) = e2;
  R.col(2) = e3;
  Eigen::Vector3d o(local_datum[0], local_datum[1], local_datum[2]);

  // point in local rotated frame
  Eigen::Vector3d gpt(global_pt[0], global_pt[1], global_pt[2]);
  Eigen::Vector3d pt_local = R.transpose() * (gpt - o);
  return IRL::Pt(pt_local[0], pt_local[1], pt_local[2]);
}

void SlicesTaubinLS::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // PLIC
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  // clipped PLIC polygons
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::Polygon> polygon(&mesh);
  Data<double> vf(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        vf(i, j, k) = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vf(i, j, k) < IRL::global_constants::VF_LOW ||
            vf(i, j, k) > IRL::global_constants::VF_HIGH) {
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

  // slicing and circle fit params
  const int nslices = 18;
  const int nlayers = 2;

  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);
  IRL::StackVector<IRL::Pt, 2> intersections;

  // computing circle fit for each mixed cell (or each cell with polygon)
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        if (polygon(i, j, k).getNumberOfVertices() <= 2) continue;

        // polygon local frame
        const IRL::Normal polygon_normal =
            calculatePolygonNormal(polygon(i, j, k));
        IRL::ReferenceFrame polygon_frame =
            referenceFrameFromNormal(polygon_normal);
        const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();

        // stencil polygon and volume fraction data
        polygon_vfrac_list.clear();
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              if (polygon(ii, jj, kk).getNumberOfVertices() <= 2) continue;
              polygon_vfrac_list.emplace_back(polygon(ii, jj, kk),
                                              vf(ii, jj, kk));
            }
          }
        }
        const int num_polygons = polygon_vfrac_list.size();

        // revert back to plane if there are no polygons in the neighboring
        // stencil. Circle fit will fail if all points are collinear
        if (num_polygons < 2) continue;  // back to LVIRA

        // accumulator for prinical curvature and direction least squares

        // slicing about local normal [0, pi)
        std::vector<double> theta_list, k_theta_list;
        for (int s = 0; s < nslices; s++) {
          // rotation angle
          double theta =
              M_PI * static_cast<double>(s) / static_cast<double>(nslices);

          // rotating local polygon frame
          const IRL::UnitQuaternion rotation(theta, polygon_frame[2]);
          const auto rotated_polygon_frame = rotation * polygon_frame;

          // slicing plane
          const IRL::Plane slicing_plane(
              rotated_polygon_frame[1],
              rotated_polygon_frame[1] * polygon_centroid);

          // slicing polygon stencial with plane and collecting info
          std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
          std::vector<double> weights;
          for (int p = 0; p < num_polygons; p++) {
            getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                               &intersections);
            if (intersections.size() != 2) continue;
            // end points for sampling
            IRL::Pt start_point = intersections[0];
            IRL::Pt end_point = intersections[1];
            end_points_list.push_back({start_point, end_point});
            // computing weights
            IRL::Normal neighbor_normal =
                calculatePolygonNormal(polygon_vfrac_list[p].first);
            IRL::Normal neighbor_centroid =
                (polygon_vfrac_list[p].first).calculateCentroid();
            double neighbor_weight =
                getTotalWeight(polygon_vfrac_list[p].second,
                               rotated_polygon_frame[2], neighbor_normal,
                               polygon_centroid, neighbor_centroid, mesh.dx());
            weights.push_back(neighbor_weight);
          }

          // end points in local frame
          std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_local_frame;
          for (const auto& ep : end_points_list) {
            IRL::Pt p1_local =
                toLocal(rotated_polygon_frame, polygon_centroid, ep.first);
            IRL::Pt p2_local =
                toLocal(rotated_polygon_frame, polygon_centroid, ep.second);
            end_points_local_frame.push_back({p1_local, p2_local});
          }

          // taubin circle fit using points in local frame for curvature
          double k_theta =
              getSignedTaubinCurvatureLS(end_points_local_frame, weights);
          if (std::isfinite(k_theta)) {
            theta_list.push_back(theta);
            k_theta_list.push_back(k_theta);
          }
        }
        // principal curvatures and directions
        auto [k1, k2, phi, ok] =
            getPrincipalCurvatures(theta_list, k_theta_list);
        if (!ok) continue;  // plane

        // Darboux frame
        const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
        const IRL::ReferenceFrame darboux_frame = rotate_phi * polygon_frame;

        // new paraboloid
        IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 * k1,
                                   0.5 * k2);

        // volume fraction matching
        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
            solver_distance(cell, vf(i, j, k), 1.0e-14, paraboloid);
        paraboloid.setDatum(
            IRL::Pt(polygon_centroid +
                    solver_distance.getDistance() * darboux_frame[2]));
        (*a_interface)(i, j, k) = paraboloid;
      }
    }
  }
}

// ====================== Taubin fit using sphere =======================

// monomial integral for barycentric coordinates
double monomialIntegral(const int& i, const int& j, const int& k,
                        const double& area) {
  auto factorial = [](int n) {
    if (n < 0) throw std::invalid_argument("negative n");
    double result = 1.0;
    for (int i = 2; i <= n; ++i) result *= i;
    return result;
  };
  return 2.0 * area * factorial(i) * factorial(j) * factorial(k) /
         factorial(i + j + k + 2);
}

using Poly = std::map<std::tuple<int, int, int>, double>;  // (i, j, k) -> coeff

Poly addPoly(const Poly& a, const Poly& b) {
  Poly out = a;
  for (auto& [e, c] : b) out[e] += c;
  return out;
}

Poly mulPoly(const Poly& a, const Poly& b, int maxDeg = 4) {
  Poly out;
  for (auto& [ea, ca] : a) {
    auto [ia, ja, ka] = ea;
    for (auto& [eb, cb] : b) {
      auto [ib, jb, kb] = eb;
      int i = ia + ib, j = ja + jb, k = ka + kb;
      if (i + j + k > maxDeg) continue;
      out[{i, j, k}] += ca * cb;
    }
  }
  return out;
}

double integratePoly(const Poly& p, const double& A) {
  double I = 0.0;
  for (auto& [e, c] : p) {
    auto [i, j, k] = e;
    I += c * monomialIntegral(i, j, k, A);
  }
  return I;
}

// barycentric coordinate formulation
Poly makeX(const double& x1, const double& x2, const double& x3) {
  Poly x;
  x[{1, 0, 0}] = x1;
  x[{0, 1, 0}] = x2;
  x[{0, 0, 1}] = x3;
  return x;
}

void getSphereFitMatrices(
    Eigen::Matrix<double, 5, 5>* M, Eigen::Matrix<double, 5, 5>* C,
    const std::vector<std::pair<IRL::Polygon, double>>& polygons) {
  // loop over polygons
  for (const auto& [poly, w] : polygons) {
    const int n_simplices = poly.getNumberOfSimplicesInDecomposition();
    double A_polygon = 0.0;

    double Mx = 0.0, My = 0.0, Mz = 0.0;
    double Mxx = 0.0, Myy = 0.0, Mzz = 0.0;
    double Mxy = 0.0, Mxz = 0.0, Myz = 0.0;
    double Ms = 0.0;
    double Mxs = 0.0, Mys = 0.0, Mzs = 0.0;
    double Mss = 0.0;

    // loop over simplices
    for (int t = 0; t < n_simplices; t++) {
      const auto& triangle = poly.getSimplexFromDecomposition(t);

      double A = triangle.calculateVolume();  // area of triangle
      A_polygon += A;

      // linear barycentric polynomials
      Poly x = makeX(triangle[0][0], triangle[1][0], triangle[2][0]);
      Poly y = makeX(triangle[0][1], triangle[1][1], triangle[2][1]);
      Poly z = makeX(triangle[0][2], triangle[1][2], triangle[2][2]);
      Poly s = addPoly(mulPoly(x, x), addPoly(mulPoly(y, y), mulPoly(z, z)));

      Mx += integratePoly(x, A);
      My += integratePoly(y, A);
      Mz += integratePoly(z, A);

      Mxx += integratePoly(mulPoly(x, x), A);
      Myy += integratePoly(mulPoly(y, y), A);
      Mzz += integratePoly(mulPoly(z, z), A);
      Mxy += integratePoly(mulPoly(x, y), A);
      Mxz += integratePoly(mulPoly(x, z), A);
      Myz += integratePoly(mulPoly(y, z), A);

      Ms += integratePoly(s, A);
      Mxs += integratePoly(mulPoly(x, s), A);
      Mys += integratePoly(mulPoly(y, s), A);
      Mzs += integratePoly(mulPoly(z, s), A);
      Mss += integratePoly(mulPoly(s, s), A);
    }
    // polygon moment matrix contribution
    (*M)(0, 0) += Mss / A_polygon * w;
    (*M)(0, 1) += Mxs / A_polygon * w;
    (*M)(0, 2) += Mys / A_polygon * w;
    (*M)(0, 3) += Mzs / A_polygon * w;
    (*M)(0, 4) += Ms / A_polygon * w;
    (*M)(1, 0) += Mxs / A_polygon * w;
    (*M)(1, 1) += Mxx / A_polygon * w;
    (*M)(1, 2) += Mxy / A_polygon * w;
    (*M)(1, 3) += Mxz / A_polygon * w;
    (*M)(1, 4) += Mx / A_polygon * w;
    (*M)(2, 0) += Mys / A_polygon * w;
    (*M)(2, 1) += Mxy / A_polygon * w;
    (*M)(2, 2) += Myy / A_polygon * w;
    (*M)(2, 3) += Myz / A_polygon * w;
    (*M)(2, 4) += My / A_polygon * w;
    (*M)(3, 0) += Mzs / A_polygon * w;
    (*M)(3, 1) += Mxz / A_polygon * w;
    (*M)(3, 2) += Myz / A_polygon * w;
    (*M)(3, 3) += Mzz / A_polygon * w;
    (*M)(3, 4) += Mz / A_polygon * w;
    (*M)(4, 0) += Ms / A_polygon * w;
    (*M)(4, 1) += Mx / A_polygon * w;
    (*M)(4, 2) += My / A_polygon * w;
    (*M)(4, 3) += Mz / A_polygon * w;
    (*M)(4, 4) += 1.0 * w;
  }
  // constraint matrix
  (*C)(0, 0) = 4.0 * ((*M)(0, 4));
  (*C)(0, 1) = 2.0 * ((*M)(1, 4));
  (*C)(0, 2) = 2.0 * ((*M)(2, 4));
  (*C)(0, 3) = 2.0 * ((*M)(3, 4));
  (*C)(1, 0) = 2.0 * ((*M)(1, 4));
  (*C)(2, 0) = 2.0 * ((*M)(2, 4));
  (*C)(3, 0) = 2.0 * ((*M)(3, 4));
  (*C)(1, 1) = (*M)(4, 4);
  (*C)(2, 2) = (*M)(4, 4);
  (*C)(3, 3) = (*M)(4, 4);
}

std::pair<IRL::Pt, double> getSphereParams(
    const std::vector<std::pair<IRL::Polygon, double>>& polygons) {  //,
  // const std::vector<double>& weights) {
  const int n = polygons.size();
  if (n < 2) return {IRL::Pt(0.0, 0.0, 0.0), 0.0};  // only one polygon segment

  // moment and constraint matrices
  Eigen::Matrix<double, 5, 5> M = Eigen::Matrix<double, 5, 5>::Zero();
  Eigen::Matrix<double, 5, 5> C = Eigen::Matrix<double, 5, 5>::Zero();
  getSphereFitMatrices(&M, &C, polygons);

  // solving generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix<double, 5, 5>> ges;
  ges.compute(M, C);

  // eigen values/vectors
  const auto& evals = ges.eigenvalues();
  const auto& evecs = ges.eigenvectors();

  // need to find eigenvector with smallest non-negative real eigenvalue
  int bestIndex = -1;
  double bestLambda = std::numeric_limits<double>::infinity();
  for (int i = 0; i < 5; i++) {
    const auto lam = evals(i);
    if (std::abs(lam.imag()) > 1e-9) continue;
    const double lambdaReal = lam.real();
    if (lambdaReal <= 0.0) continue;
    if (lambdaReal < bestLambda) {
      bestLambda = lambdaReal;
      bestIndex = i;
    }
  }
  if (bestIndex < 0) {  // no eigenvalue found
    return std::numeric_limits<std::pair<IRL::Pt, double>>::quiet_NaN();
  }

  // extracting eigenvector components
  Eigen::Matrix<std::complex<double>, 5, 1> v_c = evecs.col(bestIndex);
  Eigen::Matrix<double, 5, 1> a = v_c.real();
  double A = a(0);
  double B = a(1);
  double Cc = a(2);
  double D = a(3);
  double E = a(4);

  if (std::abs(A) < 1e-14) {
    return {IRL::Pt(0.0, 0.0, 0.0), 0.0};  // infinite radius
  }

  // finding radius of curvature
  const double num = B * B + Cc * Cc + D * D - 4.0 * A * E;
  if (num <= 0.0) {
    return std::numeric_limits<std::pair<IRL::Pt, double>>::quiet_NaN();
  }
  const double R = 0.5 * std::sqrt(num) / std::abs(A);
  IRL::Pt xc(-a(1) / (2.0 * A), -a(2) / (2.0 * A), -a(3) / (2.0 * A));

  return {xc, R};
}

// normal from sphere
IRL::Normal getNormalFromSphere(const IRL::Pt& center, const IRL::Pt& p,
                                const IRL::Normal& n_plic) {
  IRL::Pt v = p - center;
  IRL::Normal n(v[0], v[1], v[2]);
  n.normalize();

  // changing orientation based on plic
  if (n * n_plic < 0) {
    return -n;
  } else {
    return n;
  }
}

void localizePolygons(std::vector<std::pair<IRL::Polygon, double>>* polygons,
                      const int& center_index) {
  const IRL::Polygon& center_polygon = (*polygons)[center_index].first;
  const IRL::Pt datum = center_polygon.calculateCentroid();
  const IRL::ReferenceFrame frame = IRL::ReferenceFrame::fromNormal(
      center_polygon.getPlaneOfExistence().normal());

  for (auto& [polygon, weight] : *polygons) {
    // moving vertices
    for (auto& pt : polygon) {
      const IRL::Pt tmp_pt = pt - datum;
      for (int d = 0; d < 3; d++) {
        pt[d] = frame[d] * tmp_pt;
      }
    }
    if (polygon.getNumberOfVertices() > 0) {
      // aligning normals
      const auto& plane = polygon.getPlaneOfExistence();
      IRL::Normal normal = plane.normal();
      const IRL::Normal tmp_normal = plane.normal();
      for (int d = 0; d < 3; ++d) {
        normal[d] = frame[d] * tmp_normal;
      }
      polygon.setPlaneOfExistence(IRL::Plane(normal, normal * polygon[0]));
    }
  }
}

// Using sphere fit for reference frame
void SlicesTaubinS::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // PLIC
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  // clipped PLIC polygons
  // int count = 0;
  const BasicMesh& mesh = a_liq_moments.getMesh();
  Data<IRL::Polygon> polygon(&mesh);
  Data<double> vf(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        vf(i, j, k) = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vf(i, j, k) < IRL::global_constants::VF_LOW ||
            vf(i, j, k) > IRL::global_constants::VF_HIGH) {
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

  // storing normals
  Data<IRL::Normal> aligned_normals(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        if (vf(i, j, k) < IRL::global_constants::VF_LOW ||
            vf(i, j, k) > IRL::global_constants::VF_HIGH)
          continue;
        std::vector<std::pair<IRL::Polygon, double>> polygons;
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        IRL::Normal n_plic = planar_separator[0].normal();
        IRL::Pt centroid = polygon(i, j, k).calculateCentroid();
        const int nlayers = 2;
        int center_index = 0;
        int index_count = 0;
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              if (polygon(ii, jj, kk).getNumberOfVertices() <= 2) continue;
              const auto planar_separator_n =
                  std::get<IRL::PlanarSeparator>((*a_interface)(ii, jj, kk));
              IRL::Normal n_plic_n = planar_separator_n[0].normal();
              IRL::Pt centroid_n = polygon(ii, jj, kk).calculateCentroid();
              double weight = getTotalWeight(vf(ii, jj, kk), n_plic, n_plic_n,
                                             centroid, centroid_n, mesh.dx());
              // weight = 1.0;
              polygons.push_back(std::make_pair(polygon(ii, jj, kk), weight));
              if (kk == k && jj == j && ii == i) center_index = index_count;
              index_count++;
            }
          }
        }
        localizePolygons(&polygons, center_index);

        // sphere center and radius
        std::pair<IRL::Pt, double> sphere_params = getSphereParams(polygons);

        // liq barycenter
        IRL::VolumeMoments liq_moment = a_liq_moments(i, j, k);
        liq_moment.normalizeByVolume();
        IRL::Pt liq_centroid = liq_moment.centroid();

        // sphere center to global frame
        IRL::ReferenceFrame frame = IRL::ReferenceFrame::fromNormal(
            polygon(i, j, k).getPlaneOfExistence().normal());
        const IRL::Pt x0 = polygon(i, j, k).calculateCentroid();
        IRL::Pt xc_global;
        for (int e = 0; e < 3; e++) {
          for (int d = 0; d < 3; d++) {
            xc_global[e] += frame[d][e] * sphere_params.first[d];
          }
        }
        xc_global += x0;

        // normal from sphere
        aligned_normals(i, j, k) =
            getNormalFromSphere(xc_global, liq_centroid, n_plic);
      }
    }
  }

  // updating plics based on aligned normals
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        double vf = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vf < IRL::global_constants::VF_LOW ||
            vf > IRL::global_constants::VF_HIGH)
          continue;
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        IRL::Pt cell_center(mesh.xm(i), mesh.ym(j), mesh.zm(k));
        IRL::Plane plane(aligned_normals(i, j, k),
                         aligned_normals(i, j, k) * cell_center);
        IRL::PlanarSeparator updated_ps =
            IRL::PlanarSeparator::fromOnePlane(plane);
        IRL::setDistanceToMatchVolumeFraction(cell, vf, &updated_ps);
        // (*a_interface)(i, j, k) = updated_ps;
        // updating polygon
        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, updated_ps, updated_ps[0]);
      }
    }
  }
  updatePolygonBorder(&polygon);

  // slicing and circle fit params
  const int nsamples_per_segment = 10;
  const int nslices = 18;
  const int nlayers = 2;

  std::vector<std::pair<IRL::Polygon, double>> polygon_vfrac_list;
  polygon_vfrac_list.reserve(125);
  IRL::StackVector<IRL::Pt, 2> intersections;

  // computing circle fit for each mixed cell (or each cell with polygon)
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        if (polygon(i, j, k).getNumberOfVertices() <= 2) continue;

        // polygon local frame
        // const IRL::Normal polygon_normal =
        //     calculatePolygonNormal(polygon(i, j, k));
        const IRL::Normal polygon_normal = aligned_normals(i, j, k);
        IRL::ReferenceFrame polygon_frame =
            referenceFrameFromNormal(polygon_normal);
        const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();

        // stencil polygon and volume fraction data
        polygon_vfrac_list.clear();
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              if (polygon(ii, jj, kk).getNumberOfVertices() <= 2) continue;
              polygon_vfrac_list.emplace_back(polygon(ii, jj, kk),
                                              vf(ii, jj, kk));
            }
          }
        }
        const int num_polygons = polygon_vfrac_list.size();

        // revert back to plane if there are no polygons in the neighboring
        // stencil. Circle fit will fail if all points are collinear
        if (num_polygons < 2) continue;  // back to LVIRA

        // accumulator for prinical curvature and direction least squares

        // slicing about local normal [0, pi)
        std::vector<double> theta_list, k_theta_list;
        for (int s = 0; s < nslices; s++) {
          // rotation angle
          double theta =
              M_PI * static_cast<double>(s) / static_cast<double>(nslices);

          // rotating local polygon frame
          const IRL::UnitQuaternion rotation(theta, polygon_frame[2]);
          const auto rotated_polygon_frame = rotation * polygon_frame;

          // slicing plane
          const IRL::Plane slicing_plane(
              rotated_polygon_frame[1],
              rotated_polygon_frame[1] * polygon_centroid);

          // slicing polygon stencial with plane and collecting info
          std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
          std::vector<double> weights;
          for (int p = 0; p < num_polygons; p++) {
            getIntersectionPts(polygon_vfrac_list[p].first, slicing_plane,
                               &intersections);
            if (intersections.size() != 2) continue;
            // end points for sampling
            IRL::Pt start_point = intersections[0];
            IRL::Pt end_point = intersections[1];
            end_points_list.push_back({start_point, end_point});
            // computing weights
            IRL::Normal neighbor_normal =
                calculatePolygonNormal(polygon_vfrac_list[p].first);
            IRL::Normal neighbor_centroid =
                (polygon_vfrac_list[p].first).calculateCentroid();
            double neighbor_weight =
                getTotalWeight(polygon_vfrac_list[p].second,
                               rotated_polygon_frame[2], neighbor_normal,
                               polygon_centroid, neighbor_centroid, mesh.dx());
            weights.insert(weights.end(), nsamples_per_segment,
                           neighbor_weight);  // all points sampled on a
                                              // segment have the same weight
          }

          // getting points in local frame
          const int num_points = end_points_list.size() * nsamples_per_segment;
          std::vector<Eigen::Vector2d> points_rotated_frame;
          points_rotated_frame.reserve(num_points);
          sampleLocalPoints(rotated_polygon_frame, polygon_centroid,
                            end_points_list, nsamples_per_segment,
                            &points_rotated_frame);

          // taubin circle fit using points in local frame for curvature
          double k_theta =
              getSignedTaubinCurvature(points_rotated_frame, weights);
          if (std::isfinite(k_theta)) {
            theta_list.push_back(theta);
            k_theta_list.push_back(k_theta);
          }
        }
        // principal curvatures and directions
        auto [k1, k2, phi, ok] =
            getPrincipalCurvatures(theta_list, k_theta_list);
        if (!ok) continue;  // plane

        // Darboux frame
        const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
        const IRL::ReferenceFrame darboux_frame = rotate_phi * polygon_frame;

        // new paraboloid
        IRL::Paraboloid paraboloid(polygon_centroid, darboux_frame, 0.5 * k1,
                                   0.5 * k2);

        // volume fraction matching
        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
            solver_distance(cell, vf(i, j, k), 1.0e-14, paraboloid);
        paraboloid.setDatum(
            IRL::Pt(polygon_centroid +
                    solver_distance.getDistance() * darboux_frame[2]));
        (*a_interface)(i, j, k) = paraboloid;
      }
    }
  }
  // a_interface->updateBorder();
  // correctInterfaceBorders(a_interface);
}

// ====================== Relaligning PLIC normals =======================
// normal of projecting point on paraboloid
IRL::Normal getNormalAtProjectedPoint(const IRL::Pt& a_point,
                                      const IRL::Paraboloid& a_paraboloid,
                                      const int max_iter = 100,
                                      const double tol = 1.0e-10) {
  // reference frame and datum
  const IRL::ReferenceFrame frame = a_paraboloid.getReferenceFrame();
  const IRL::Pt datum = a_paraboloid.getDatum();

  // coefficients
  const double a = a_paraboloid.getAlignedParaboloid().a();
  const double b = a_paraboloid.getAlignedParaboloid().b();

  auto F = [&](const double& x, const double& y, const double& z) -> double {
    return z + a * x * x + b * y * y;
  };

  auto gradF = [&](const double& x, const double& y,
                   const double& z) -> Eigen::Vector3d {
    return Eigen::Vector3d(2.0 * a * x, 2.0 * b * y, 1.0);
  };

  // point in local frame
  IRL::Pt local_pt = toLocal(frame, datum, a_point);
  Eigen::Vector3d x_proj(local_pt[0], local_pt[1], local_pt[2]);

  // newton iterations
  for (int i = 0; i < max_iter; i++) {
    const double f = F(x_proj[0], x_proj[1], x_proj[2]);
    const Eigen::Vector3d g = gradF(x_proj[0], x_proj[1], x_proj[2]);
    const double g_norm2 = g.squaredNorm();
    if (g_norm2 < 1.0e-14) break;
    x_proj -= (f / g_norm2) * g;
    if (std::abs(f) < tol) break;
    if (i == max_iter - 1) {
      std::cout << "Max iterations reached. Projection incomplete. "
                << "f = " << std::abs(f) << std::endl;
    }
  }

  // normal at projected point
  Eigen::Vector3d normal_vec = gradF(x_proj[0], x_proj[1], x_proj[2]);
  normal_vec.normalize();

  // bring it back to global frame
  Eigen::Vector3d e1(frame[0][0], frame[0][1], frame[0][2]);
  Eigen::Vector3d e2(frame[2][0], frame[2][1], frame[2][2]);
  Eigen::Vector3d e3(frame[1][0], frame[1][1], frame[1][2]);
  Eigen::Matrix3d R;
  R.col(0) = e1;
  R.col(1) = e2;
  R.col(2) = e3;
  Eigen::Vector3d normal_global = R * normal_vec;

  return IRL::Normal(normal_global[0], normal_global[1], normal_global[2]);
}

// ===================== normal averaging =============================
void PLICalign::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  const BasicMesh& mesh = a_liq_moments.getMesh();

  // planar interface that will be updated iteratively
  Data<IRL::SeparatorVariant> plic_interface(&mesh);
  plic_interface = *a_interface;

  using VolumeMomentsAndSurface =
      IRL::AddSurfaceOutput<IRL::VolumeMoments,
                            IRL::ParaboloidParametrizedSurfaceOutput>;

  const int n_iters = 100;  // iterations of realigning normals

  for (int it = 0; it < n_iters; ++it) {
    // taubin paraboloid fit
    Data<IRL::SeparatorVariant> taubin_interface(&mesh);
    taubin_interface = plic_interface;
    SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                    a_V, a_W, &taubin_interface,
                                    a_scalar_fields, true);

    // realigning plic normals using taubin paraboloid
    for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
        for (int i = mesh.imin(); i <= mesh.imax(); i++) {
          double vf = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
          if (vf < IRL::global_constants::VF_LOW ||
              vf > IRL::global_constants::VF_HIGH)
            continue;
          auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
              IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));

          const auto* parab =
              std::get_if<IRL::Paraboloid>(&(taubin_interface(i, j, k)));
          if (!parab) continue;
          auto* planar =
              std::get_if<IRL::PlanarSeparator>(&(plic_interface(i, j, k)));
          if (!planar) continue;

          // parameterized paraboloid surface
          auto surface =
              IRL::getVolumeMoments<VolumeMomentsAndSurface>(cell, *parab)
                  .getSurface();

          // average normal on paraboloid surface
          IRL::Normal n_avg = surface.getAverageNormalNonAligned();

          IRL::Pt cell_center(mesh.xm(i), mesh.ym(j), mesh.zm(k));
          IRL::Plane updated_plane(n_avg, n_avg * cell_center);

          IRL::PlanarSeparator updated_ps =
              IRL::PlanarSeparator::fromOnePlane(updated_plane);

          IRL::setDistanceToMatchVolumeFraction(cell, vf, &updated_ps);

          plic_interface(i, j, k) = updated_ps;
        }
      }
    }
  }

  // taubin paraboloid with new plics
  SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                  a_W, &plic_interface, a_scalar_fields, true);

  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        (*a_interface)(i, j, k) = plic_interface(i, j, k);
      }
    }
  }
}

void PLICalign2::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  const BasicMesh& mesh = a_liq_moments.getMesh();

  // planar interface that will be updated iteratively
  Data<IRL::SeparatorVariant> plic_interface(&mesh);
  plic_interface = *a_interface;

  // taubin fit
  Data<IRL::SeparatorVariant> taubin_interface(&mesh);
  taubin_interface = *a_interface;
  SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                  a_W, &taubin_interface, a_scalar_fields,
                                  true);

  using VolumeMomentsAndSurface =
      IRL::AddSurfaceOutput<IRL::VolumeMoments,
                            IRL::ParaboloidParametrizedSurfaceOutput>;

  Data<std::pair<IRL::Normal, int>> normal_contribution(&mesh);

  const int nlayers = 2;

  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        double vf = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vf < IRL::global_constants::VF_LOW ||
            vf > IRL::global_constants::VF_HIGH)
          continue;
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              double vf_neighbor =
                  a_liq_moments(ii, jj, kk).volume() / mesh.cell_volume();
              if (vf_neighbor < IRL::global_constants::VF_LOW ||
                  vf_neighbor > IRL::global_constants::VF_HIGH)
                continue;

              const auto* parab =
                  std::get_if<IRL::Paraboloid>(&(taubin_interface(ii, jj, kk)));
              if (!parab) continue;

              auto cell = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                  IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));

              // parameterized paraboloid surface
              auto surface =
                  IRL::getVolumeMoments<VolumeMomentsAndSurface>(cell, *parab)
                      .getSurface();

              // average normal on paraboloid surface
              IRL::Normal n_avg = surface.getAverageNormalNonAligned();

              normal_contribution(ii, jj, kk).first += n_avg;
              normal_contribution(ii, jj, kk).second += 1;
            }
          }
        }
      }
    }
  }

  // updating plic normals
  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        double vf = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (vf < IRL::global_constants::VF_LOW ||
            vf > IRL::global_constants::VF_HIGH)
          continue;
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        IRL::Normal n_avg =
            normal_contribution(i, j, k).first /
            static_cast<double>(normal_contribution(i, j, k).second);
        n_avg.normalize();
        IRL::Pt cell_center(mesh.xm(i), mesh.ym(j), mesh.zm(k));
        IRL::Plane plane(n_avg, n_avg * cell_center);
        IRL::PlanarSeparator updated_ps =
            IRL::PlanarSeparator::fromOnePlane(plane);
        IRL::setDistanceToMatchVolumeFraction(cell, vf, &updated_ps);
        taubin_interface(i, j, k) = updated_ps;
      }
    }
  }
  // final taubin paraboloid with new plics
  SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                  a_W, &taubin_interface, a_scalar_fields,
                                  true);

  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        (*a_interface)(i, j, k) = taubin_interface(i, j, k);
      }
    }
  }
}

// ====================== Mosso Swartz =======================
// function to rotate a vector about an axis
Eigen::Vector3d rotate_about_axis(const Eigen::Vector3d& v,
                                  const Eigen::Vector3d& axis,
                                  const double& theta) {
  return v + std::sin(theta) * axis.cross(v) +
         (1.0 - std::cos(theta)) * axis.cross(axis.cross(v));
}

void MossoSwartz::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  const BasicMesh& mesh = a_liq_moments.getMesh();

  const int max_iter = 4;

  for (int iter = 0; iter < max_iter; ++iter) {
    Data<IRL::SeparatorVariant> updated_interfaces(&mesh);
    updated_interfaces = *a_interface;

    for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
        for (int i = mesh.imin(); i <= mesh.imax(); i++) {
          double vf = a_liq_moments(i, j, k).volume() / mesh.cell_volume();
          if (vf < IRL::global_constants::VF_LOW ||
              vf > IRL::global_constants::VF_HIGH)
            continue;
          auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
              IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          const auto planar_separator =
              std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
          IRL::Polygon target_polygon =
              IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                  cell, planar_separator, planar_separator[0]);
          IRL::Pt targ_centroid = target_polygon.calculateCentroid();
          IRL::Normal targ_normal = planar_separator[0].normal();
          Eigen::Vector3d target_centroid(targ_centroid[0], targ_centroid[1],
                                          targ_centroid[2]);
          Eigen::Vector3d target_normal(targ_normal[0], targ_normal[1],
                                        targ_normal[2]);
          std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>
              stencil_quantities;  // centroids and normals

          // gathering stencil quantities
          for (int kk = k - 1; kk <= k + 1; ++kk) {
            for (int jj = j - 1; jj <= j + 1; ++jj) {
              for (int ii = i - 1; ii <= i + 1; ++ii) {
                if (i == ii && j == jj && k == kk) continue;
                double vf_neighbor =
                    a_liq_moments(ii, jj, kk).volume() / mesh.cell_volume();
                if (vf_neighbor < IRL::global_constants::VF_LOW ||
                    vf_neighbor > IRL::global_constants::VF_HIGH)
                  continue;
                auto neighbor_cell = IRL::RectangularCuboid::fromBoundingPts(
                    IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                    IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
                const auto neighbor_planar_separator =
                    std::get<IRL::PlanarSeparator>((*a_interface)(ii, jj, kk));
                IRL::Polygon neighbor_polygon =
                    IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                        neighbor_cell, neighbor_planar_separator,
                        neighbor_planar_separator[0]);
                IRL::Pt neighbor_centroid =
                    neighbor_polygon.calculateCentroid();
                IRL::Normal neighbor_normal =
                    neighbor_planar_separator[0].normal();
                stencil_quantities.push_back(
                    {Eigen::Vector3d(neighbor_centroid[0], neighbor_centroid[1],
                                     neighbor_centroid[2]),
                     Eigen::Vector3d(neighbor_normal[0], neighbor_normal[1],
                                     neighbor_normal[2])});
              }
            }
          }
          // Mosso swartz update
          std::vector<Eigen::Vector3d> rotated_normals;
          Eigen::Vector3d rotated_normal_sum = Eigen::Vector3d::Zero();
          for (const auto& sq : stencil_quantities) {
            Eigen::Vector3d d = sq.first - target_centroid;
            if (d.norm() < 1e-12) continue;
            d.normalize();
            Eigen::Vector3d axis = d.cross(target_normal);
            if (axis.norm() < 1e-12) continue;
            axis.normalize();
            double angle = M_PI / 2.0;
            Eigen::Vector3d rotated_normal = rotate_about_axis(d, axis, angle);
            rotated_normal.normalize();
            if (rotated_normal.dot(target_normal) < 0.0) {
              rotated_normal = -rotated_normal;
            }
            rotated_normal_sum += rotated_normal;
            rotated_normals.push_back(rotated_normal);
          }
          IRL::Normal new_normal(rotated_normal_sum[0], rotated_normal_sum[1],
                                 rotated_normal_sum[2]);
          if (rotated_normals.empty()) continue;
          new_normal /= rotated_normals.size();
          new_normal.normalize();

          // plane with new normal
          IRL::Pt cell_center(mesh.xm(i), mesh.ym(j), mesh.zm(k));
          IRL::Plane plane(new_normal, new_normal * cell_center);
          IRL::PlanarSeparator updated_ps =
              IRL::PlanarSeparator::fromOnePlane(plane);
          IRL::setDistanceToMatchVolumeFraction(cell, vf, &updated_ps);
          updated_interfaces(i, j, k) = updated_ps;
        }
      }
    }
    for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
        for (int i = mesh.imin(); i <= mesh.imax(); i++) {
          (*a_interface)(i, j, k) = updated_interfaces(i, j, k);
        }
      }
    }
  }
}

// ====================== Hybrid Jibben and PU =======================
void Hybrid::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  const BasicMesh& mesh = a_liq_moments.getMesh();

  // plic for jibben and PU
  Data<IRL::SeparatorVariant> jibben_interface(&a_liq_moments.getMesh());
  Data<IRL::SeparatorVariant> pu_interface(&a_liq_moments.getMesh());
  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        jibben_interface(i, j, k) = (*a_interface)(i, j, k);
        pu_interface(i, j, k) = (*a_interface)(i, j, k);
      }
    }
  }

  // PU reconstruction
  PU::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                        &pu_interface, a_scalar_fields, true);

  // plic polygons
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
        const auto planar_separator =
            std::get<IRL::PlanarSeparator>((*a_interface)(i, j, k));
        polygon(i, j, k) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            cell, planar_separator, planar_separator[0]);
      }
    }
  }
  updatePolygonBorder(&polygon);

  // jibben fit
  IRL::JibbenNeighborhood neighborhood;
  const int nlayers = 1;
  const int nstencil =
      (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);
  neighborhood.reserve(nstencil);
  neighborhood.setDelta(2.5 * mesh.dx());

  InterfaceScalarField interface_type_field("interface_type", &mesh);

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

          // jibben fit
          IRL::Jibben_3D jibben_solver;
          jibben_interface(i, j, k) = jibben_solver.solve(&neighborhood);

          // check for squared volume error
          // double volume_error_sq =
          //     jibben_solver.getVolumeErrorSquared(mesh.dx());

          // normal error
          // double normal_error = jibben_solver.getNormalMetric();
          double normal_error = jibben_solver.getNormalEigenMetric();
          // double normal_error = jibben_solver.getNormalVarianceMetric();

          // double volume_error = jibben_solver.getVolumeError(mesh.dx());

          // if (volume_error_sq > 0.01) {
          // if (volume_error_sq > 0.05) {
          if (normal_error > 0.2) {
            // if (volume_error > 0.05) {
            (*a_interface)(i, j, k) = pu_interface(i, j, k);
            interface_type_field.paraboloid_scalar_data(i, j, k) = 1.0;
          } else {
            const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
            const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                        mesh.z(k + 1));
            auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                                upper_cell_pt);
            IRL::setDistanceToMatchVolumeFraction(cell, liquid_volume_fraction,
                                                  &(jibben_interface)(i, j, k),
                                                  1.0e-14);
            (*a_interface)(i, j, k) = jibben_interface(i, j, k);
            interface_type_field.paraboloid_scalar_data(i, j, k) = 0.0;
          }
        }
      }
    }
  }

  a_scalar_fields->push_back(interface_type_field);

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// ====================== Hybrid Circle Fit and PU =======================
void Hybrid2::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // building planar interface
  if (a_plic_already_built == false) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
  }

  const BasicMesh& mesh = a_liq_moments.getMesh();

  // plic for jibben and PU
  Data<IRL::SeparatorVariant> plic_reconstruction(&a_liq_moments.getMesh());
  Data<IRL::SeparatorVariant> taubin_reconstruction(&a_liq_moments.getMesh());
  for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int i = mesh.imin(); i <= mesh.imax(); i++) {
        plic_reconstruction(i, j, k) = (*a_interface)(i, j, k);
        taubin_reconstruction(i, j, k) = (*a_interface)(i, j, k);
      }
    }
  }

  // circle fit reconstruction
  SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                  a_W, &taubin_reconstruction, a_scalar_fields,
                                  true);

  // plic centroids and areas
  Data<IRL::Pt> interface_centroids(&mesh);
  Data<double> interface_areas(&mesh);
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        interface_centroids(i, j, k) = IRL::Pt();
        interface_areas(i, j, k) = 0.0;
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
        IRL::Polygon polygon =
            IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                cell, planar_separator, planar_separator[0]);
        interface_centroids(i, j, k) = polygon.calculateCentroid();
        interface_areas(i, j, k) = polygon.calculateVolume();
      }
    }
  }

  // pu params
  const int nlayers = 3;
  const double delta = 2.5 * mesh.dx();

  // pu fit on taubin interface
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
            liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          continue;
        }
        // projecting centroid on PU approximation
        const IRL::Pt pt_on_PU = projectCentroidOnPU(
            a_liq_moments, taubin_reconstruction, interface_centroids,
            interface_areas, nlayers, delta, i, j, k);

        // gradient and hessian on PU approximation
        const auto F_and_gradF_and_hessF = getPUAndGradAndHessian(
            a_liq_moments, taubin_reconstruction, interface_centroids,
            interface_areas, nlayers, delta, i, j, k, pt_on_PU);
        const Eigen::Vector3d gradF =
            std::get<Eigen::Vector3d>(F_and_gradF_and_hessF);
        const Eigen::Matrix3d hessF =
            std::get<Eigen::Matrix3d>(F_and_gradF_and_hessF);
        auto new_normal = IRL::Normal(gradF(0), gradF(1), gradF(2));
        new_normal.normalize();

        // generating paraboloid from PU using gradient and hessian
        // if (IRL::magnitude(new_normal) > 0.9) {
        //   auto paraboloid =
        //       IRL::Paraboloid::fromDerivatives(pt_on_PU, gradF, hessF);
        //   const double A = paraboloid.getAlignedParaboloid().a();
        //   const double B = paraboloid.getAlignedParaboloid().b();
        //   // if (std::fabs(A) * mesh.dx() > 4.0 ||
        //   //     std::fabs(B) * mesh.dx() > 4.0) {
        //   //   continue;
        //   // }
        //   // Translate paraboloid to match volume fraction
        //   const auto new_datum = paraboloid.getDatum();
        //   const auto new_frame = paraboloid.getReferenceFrame();
        //   const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        //   const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
        //                               mesh.z(k + 1));
        //   const auto cell = IRL::RectangularCuboid::fromBoundingPts(
        //       lower_cell_pt, upper_cell_pt);
        //   IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
        //       solver_distance(cell, liquid_volume_fraction, 1.0e-14,
        //                       paraboloid);
        //   paraboloid.setDatum(IRL::Pt(
        //       new_datum + solver_distance.getDistance() * new_frame[2]));
        //   (*a_interface)(i, j, k) = paraboloid;
        // }

        // realigning planar interface with new normal
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        IRL::Pt cell_center(mesh.xm(i), mesh.ym(j), mesh.zm(k));
        IRL::Plane plane(new_normal, new_normal * cell_center);
        IRL::PlanarSeparator updated_ps =
            IRL::PlanarSeparator::fromOnePlane(plane);
        IRL::setDistanceToMatchVolumeFraction(cell, liquid_volume_fraction,
                                              &updated_ps, 1.0e-14);
        // plic_reconstruction(i, j, k) = updated_ps;
        (*a_interface)(i, j, k) = updated_ps;
      }
    }
  }

  // circle fit reconstruction with new plics
  SlicesTaubin::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V,
                                  a_W, a_interface, a_scalar_fields, true);
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
  // line_seg_endpoints: endpoints {a,b} of line segments that are cloest to
  // the particle

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
  // target_endpoints: end points of the target interface where curvature is
  // to be estimated hp: spacing between particles along line segment N:
  // number of particles (odd)

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
      // std::cout << "Maximum iterations reached for particle method.
      // Residual = "
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
    Data<IRL::SeparatorVariant>* a_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields,
    const bool a_plic_already_built) {
  // plic
  if (!a_plic_already_built) {
    LVIRA::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                             a_interface, a_scalar_fields);
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
