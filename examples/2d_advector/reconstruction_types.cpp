// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include <mpi.h>

#include "examples/2d_advector/reconstruction_types.h"

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
#include "examples/2d_advector/basic_mesh.h"
#include "examples/2d_advector/data.h"
#include "examples/2d_advector/vof_advection.h"

void getReconstruction(const std::string& a_reconstruction_method,
                       const Data<IRL2D::Moments>& a_liquid_moments,
                       const Data<IRL2D::Moments>& a_gas_moments,
                       const double a_dt, const Data<double>& a_U,
                       const Data<double>& a_V,
                       Data<IRL2D::Parabola>* a_interface) {
  // if (a_reconstruction_method == "PLIC") {
  //   PLIC::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
  //                           a_W, a_interface, a_localized_paraboloid_link);
  // } else {
  std::cout << "Unknown reconstruction method of : " << a_reconstruction_method
            << '\n';
  std::cout << "Valid entries are: PLIC. \n";
  // std::exit(-1);
  // }
}

void RecenterMoments(IRL2D::Moments* moments, const IRL2D::Vec& center) {
  moments->m2() += -outer_product((*moments).m1(), center) -
                   outer_product(center, (*moments).m1()) +
                   (*moments).m0() * outer_product(center, center);
  moments->m1() -= (*moments).m0() * center;
}

// void PLIC::getReconstruction(
//     const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
//     const Data<IRL::GeneralMoments3D<2>>& a_gas_moments, const double a_dt,
//     const Data<double>& a_U, const Data<double>& a_V, const Data<double>&
//     a_W, Data<IRL::Paraboloid>* a_interface,
//     Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid)
//     {
//   const BasicMesh& mesh = a_U.getMesh();

//   Data<IRL::PlanarSeparator> interface(&mesh);
//   updateReconstructionELVIRA(a_liquid_moments, &interface);
//   updateReconstructionLVIRA(a_liquid_moments, 1, &interface);
//   Data<IRL::Polygon> polygon(&mesh);
//   updatePolygon(a_liquid_moments, interface, &polygon);
//   polygon.updateBorder();

//   // x- boundary
//   for (int i = mesh.imino(); i < mesh.imin(); ++i) {
//     for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
//       for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
//         for (auto& pt : polygon(i, j, k)) {
//           pt[0] -= mesh.lx();
//         }
//       }
//     }
//   }

//   // x+ boundary
//   for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
//     for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
//       for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
//         for (auto& pt : polygon(i, j, k)) {
//           pt[0] += mesh.lx();
//         }
//       }
//     }
//   }

//   // y- boundary
//   for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
//     for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
//       for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
//         for (auto& pt : polygon(i, j, k)) {
//           pt[1] -= mesh.ly();
//         }
//       }
//     }
//   }

//   // y+ boundary
//   for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
//     for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
//       for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
//         for (auto& pt : polygon(i, j, k)) {
//           pt[1] += mesh.ly();
//         }
//       }
//     }
//   }

//   // z- boundary
//   for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
//     for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
//       for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
//         for (auto& pt : polygon(i, j, k)) {
//           pt[2] -= mesh.lz();
//         }
//       }
//     }
//   }

//   // z+ boundary
//   for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
//     for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
//       for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
//         for (auto& pt : polygon(i, j, k)) {
//           pt[2] += mesh.lz();
//         }
//       }
//     }
//   }

//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
//         const double liquid_volume_fraction =
//             a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
//         if (liquid_volume_fraction < IRL::global_constants::VF_LOW) {
//           (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
//         } else if (liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
//           (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
//           // continue;
//         } else {
//           const IRL::Normal norm_poly = polygon(i, j, k).calculateNormal();
//           const IRL::Pt pref = polygon(i, j, k).calculateCentroid();
//           IRL::ReferenceFrame fit_frame;
//           int largest_dir = 0;
//           if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[1]))
//             largest_dir = 1;
//           if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[2]))
//             largest_dir = 2;
//           if (largest_dir == 0)
//             fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 1.0,
//             0.0));
//           else if (largest_dir == 1)
//             fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0,
//             0.0, 1.0));
//           else
//             fit_frame[0] = crossProduct(norm_poly, IRL::Normal(1.0, 0.0,
//             0.0));
//           fit_frame[0].normalize();
//           fit_frame[1] = crossProduct(norm_poly, fit_frame[0]);
//           fit_frame[2] = norm_poly;
//           const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
//           const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
//                                       mesh.z(k + 1));
//           const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
//           IRL::Paraboloid paraboloid;

//           paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);

//           auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
//                                                               upper_cell_pt);
//           IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
//               solver_distance(cell, liquid_volume_fraction, 1.0e-14,
//                               paraboloid);

//           if (solver_distance.getDistance() == -DBL_MAX) {
//             (*a_interface)(i, j, k) =
//                 IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
//             std::cout << "Distance solver failed: error = " << std::endl;
//           } else {
//             auto new_datum =
//                 IRL::Pt(paraboloid.getDatum() +
//                         solver_distance.getDistance() * fit_frame[2]);
//             paraboloid.setDatum(new_datum);
//             (*a_interface)(i, j, k) = paraboloid;
//           }
//         }
//       }
//     }
//   }

//   // Update border with simple ghost-cell fill and correct datum for
//   // assumed periodic boundary
//   a_interface->updateBorder();
//   correctInterfacePlaneBorders(a_interface);
// }

void correctInterfaceBorders(Data<IRL2D::Parabola>* a_interface) {
  const BasicMesh& mesh = (*a_interface).getMesh();
  // Fix distances in reconstruction for periodic boundary

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      (*a_interface)(i, j).datum()[0] -= mesh.lx();
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      (*a_interface)(i, j).datum()[0] += mesh.lx();
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      (*a_interface)(i, j).datum()[1] -= mesh.ly();
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      (*a_interface)(i, j).datum()[1] += mesh.ly();
    }
  }
}
