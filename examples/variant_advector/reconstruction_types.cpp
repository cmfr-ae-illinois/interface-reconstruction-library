// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include <mpi.h>

#include "examples/variant_advector/reconstruction_types.h"

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
#include "examples/variant_advector/basic_mesh.h"
#include "examples/variant_advector/data.h"
#include "examples/variant_advector/vof_advection.h"

void getReconstruction(
    const std::string& a_reconstruction_method,
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::SeparatorVariant>* a_interface) {
  if (a_reconstruction_method == "PLIC") {
    PLIC::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                            a_interface, a_link_localized_interface);
  } else if (a_reconstruction_method == "Jibben") {
    Jibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                              a_interface, a_link_localized_interface);
  } else if (a_reconstruction_method == "Mixed") {
    MixedPLICJibben::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U,
                                       a_V, a_W, a_interface,
                                       a_link_localized_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: PLIC, Jibben, Mixed. \n";
    std::exit(-1);
  }
}

void PLIC::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface) {
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

void Jibben::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface,
    const bool a_plic_already_built) {
  // First, we need to build the plic
  if (a_plic_already_built == false) {
    PLIC::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                            a_interface, a_link_localized_interface);
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

  // Now let's compute the Jibben parabolic fit
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          // Compute local frame of reference based on polygon
          const IRL::Normal polygon_normal = polygon(i, j, k).calculateNormal();
          const IRL::ReferenceFrame polygon_frame =
              referenceFrameFromNormal(polygon_normal);
          const IRL::Pt polygon_centroid = polygon(i, j, k).calculateCentroid();
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
          // Compute number of existing polygons in 5x5x5 stencil
          const int nneigh = 1;
          int num_polygons = 0;
          for (int kk = k - nneigh; kk <= k + nneigh; ++kk) {
            for (int jj = j - nneigh; jj <= j + nneigh; ++jj) {
              for (int ii = i - nneigh; ii <= i + nneigh; ++ii) {
                if (polygon(ii, jj, kk).getNumberOfVertices() > 0) {
                  num_polygons++;
                }
              }
            }
          }
          // Allocate least-squares system (using Eigen)
          Eigen::MatrixXd Amat = Eigen::MatrixXd::Zero(num_polygons, 6);
          Eigen::VectorXd bvec = Eigen::VectorXd::Zero(num_polygons);
          num_polygons = 0;
          for (int kk = k - nneigh; kk <= k + nneigh; ++kk) {
            for (int jj = j - nneigh; jj <= j + nneigh; ++jj) {
              for (int ii = i - nneigh; ii <= i + nneigh; ++ii) {
                // Skip if no polygon
                const IRL::UnsignedIndex_t num_vertices =
                    polygon(ii, jj, kk).getNumberOfVertices();
                if (num_vertices == 0) {
                  continue;
                }
                // Local polygon normal and centroid
                IRL::Pt local_polygon_centroid =
                    polygon(ii, jj, kk).calculateCentroid() - polygon_centroid;
                IRL::Normal local_polygon_normal =
                    polygon(ii, jj, kk).calculateNormal();
                // Ignore polygons oriented more than 90 degrees compared to
                // central polygon
                if (polygon_frame[2] * local_polygon_normal <= 0.0) {
                  continue;
                }
                // Move centroid and normal to local frame
                const IRL::Pt tmp_c = local_polygon_centroid;
                const IRL::Normal tmp_n = local_polygon_normal;
                for (int d = 0; d < 3; ++d) {
                  local_polygon_centroid[d] = polygon_frame[d] * tmp_c;
                  local_polygon_normal[d] = polygon_frame[d] * tmp_n;
                }
                // Compute polygon plane coefficients
                Eigen::VectorXd plane_coeffs(3);
                plane_coeffs
                    << -(local_polygon_centroid * local_polygon_normal),
                    local_polygon_normal[0], local_polygon_normal[1];
                plane_coeffs /= -local_polygon_normal[2];
                // Compute distance and volume fraction weighting
                const double distance = IRL::magnitude(local_polygon_centroid);
                const double distance_ndim = distance / 2.5 * mesh.dx();
                const double distance_weight =
                    distance_ndim >= 1.0
                        ? 0.0
                        : (1.0 + 4.0 * distance_ndim) *
                              std::pow(1.0 - distance_ndim, 4.0);
                const double vfrac =
                    a_liq_moments(ii, jj, kk).volume() / mesh.cell_volume();
                double vfrac_weight = 1.0;
                const double limit_vfrac = 0.1;
                if (vfrac < 0.1) {
                  vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * M_PI * vfrac);
                } else if (vfrac > 0.9) {
                  vfrac_weight =
                      0.5 - 0.5 * std::cos(10.0 * M_PI * (1.0 - vfrac));
                }
                const double weight = vfrac_weight * distance_weight;
                // Compute momonial integrals and feed to LS system
                Eigen::VectorXd integrals = Eigen::VectorXd::Zero(6);
                double b_dot_sum = 0.0;
                for (int v = 0; v < num_vertices; ++v) {
                  const int vn = (v + 1) % num_vertices;
                  IRL::Pt vert1 = polygon(ii, jj, kk)[v] - polygon_centroid;
                  IRL::Pt vert2 = polygon(ii, jj, kk)[vn] - polygon_centroid;
                  IRL::Pt tmp_pt1 = vert1, tmp_pt2 = vert2;
                  for (int d = 0; d < 3; ++d) {
                    vert1[d] = polygon_frame[d] * tmp_pt1;
                    vert2[d] = polygon_frame[d] * tmp_pt2;
                  }
                  const double xv = vert1[0], yv = vert1[1];
                  const double xvn = vert2[0], yvn = vert2[1];
                  Eigen::VectorXd integral_to_add(6);
                  integral_to_add << (xv * yvn - xvn * yv) / 2.0,
                      (xv + xvn) * (xv * yvn - xvn * yv) / 6.0,
                      (yv + yvn) * (xv * yvn - xvn * yv) / 6.0,
                      (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0,
                      (yvn - yv) *
                          (3.0 * xv * xv * yv + xv * xv * yvn +
                           2.0 * xv * xvn * yv + 2.0 * xv * xvn * yvn +
                           xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
                          24.0,
                      (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;
                  integrals += integral_to_add;
                }
                if (weight > 0.0) {
                  Amat.row(num_polygons) = weight * integrals;
                  bvec(num_polygons) =
                      weight * integrals.head(3).dot(plane_coeffs);
                }
                num_polygons++;
              }
            }
          }
          // Unconstrained LS solution
          const Eigen::VectorXd sol =
              Amat.completeOrthogonalDecomposition().pseudoInverse() * bvec;
          const double a = sol[0], b = sol[1], c = sol[2], d = sol[3],
                       e = sol[4], f = sol[5];
          const double theta = 0.5 * std::atan2(e, (IRL::safelyTiny(d - f)));
          const double cos_t = std::cos(theta);
          const double sin_t = std::sin(theta);
          const double A =
              -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
          const double B =
              -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);
          // Translation to coordinate system R' where aligned paraboloid valid
          // Translation is R ' = {x' = x + u, y ' = y + v, z' = z + w }
          const double denom_inv = 1.0 / IRL::safelyTiny(4.0 * d * f - e * e);
          const double u = (2.0 * b * f - c * e) * denom_inv;
          const double v = -(b * e - 2.0 * d * c) * denom_inv;
          const double w =
              -(a + (-b * b * f + b * c * e - c * c * d) * denom_inv);
          const IRL::Pt paraboloid_datum =
              polygon_centroid - u * polygon_frame[0] - v * polygon_frame[1] -
              w * polygon_frame[2];
          const IRL::UnitQuaternion rotation(theta, polygon_frame[2]);
          const auto paraboloid_frame = rotation * polygon_frame;
          auto paraboloid =
              IRL::Paraboloid(paraboloid_datum, paraboloid_frame, A, B);

          // Translate paraboloid to match volume fraction
          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, liquid_volume_fraction, 1.0e-14,
                              paraboloid);
          paraboloid.setDatum(
              IRL::Pt(paraboloid_datum +
                      solver_distance.getDistance() * paraboloid_frame[2]));
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

void MixedPLICJibben::getReconstruction(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::SeparatorVariant>* a_interface,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface) {
  // First, we need to build the plic and copy it
  PLIC::getReconstruction(a_liq_moments, a_gas_moments, a_dt, a_U, a_V, a_W,
                          a_interface, a_link_localized_interface);
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
                            a_interface, a_link_localized_interface, true);

  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        if (j < int(mesh.jmax() / 2)) {
          (*a_interface)(i, j, k) = plic_reconstruction(i, j, k);
        }
      }
    }
  }
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
