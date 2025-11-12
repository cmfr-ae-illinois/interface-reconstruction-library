// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include "examples/2d_advector/reconstruction_types.h"

#include "irl/interface_reconstruction_methods/constrained_optimization_behavior.h"
#include "irl/interface_reconstruction_methods/elvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/lvira_neighborhood.h"
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

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <nlopt.hpp>
#include <queue>
#include <set>

void getReconstruction(const std::string& a_reconstruction_method,
                       const Data<IRL2D::Moments>& a_liquid_moments,
                       const Data<IRL2D::Moments>& a_gas_moments,
                       const double a_dt, const Data<double>& a_U,
                       const Data<double>& a_V,
                       Data<IRL2D::Parabola>* a_interface) {
  if (a_reconstruction_method == "ELVIRA") {
    ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "LVIRA") {
    LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                             a_interface);
  } else if (a_reconstruction_method == "Jibben") {
    Jibben::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "Particle") {
    Particle::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                                a_interface);
  } else if (a_reconstruction_method == "JPHybrid") {
    JPHybrid::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                                a_interface);
  } else if (a_reconstruction_method == "LVIRAQ") {
    LVIRAQ::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "LMIRAQ") {
    LMIRAQ::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "MOF1") {
    MOF1::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);
  } else if (a_reconstruction_method == "MOF2") {
    MOF2::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);
  } else if (a_reconstruction_method == "MOF2AL") {
    MOF2AL::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "MOF2ALUnit") {
    MOF2ALUnit::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U,
                                  a_V, a_interface);
  } else if (a_reconstruction_method == "NLOPT") {
    NLOPT::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                             a_interface);
  } else if (a_reconstruction_method == "Pratt") {
    Pratt::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                             a_interface);
  } else if (a_reconstruction_method == "Taubin") {
    Taubin::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "LMA") {
    LMA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);
  } else if (a_reconstruction_method == "Pratt2") {
    Pratt2::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "Pratt3") {
    Pratt3::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "Cubic") {
    Cubic::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                             a_interface);
  } else if (a_reconstruction_method == "PU") {
    PU::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                          a_interface);
  } else if (a_reconstruction_method == "PUPLIC") {
    PUPLIC::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: ELVIRA, LVIRA, Jibben, Particle, LMIRAQ, "
                 "LVIRAQ, MOF1, MOF2, MOF2AL, NLOPT, Pratt, Taubin, LMA. \n";
    std::exit(-1);
  }
}

void RecenterMoments(IRL2D::Moments* moments, const IRL2D::Vec& center) {
  moments->m2() += -outer_product((*moments).m1(), center) -
                   outer_product(center, (*moments).m1()) +
                   (*moments).m0() * outer_product(center, center);
  moments->m1() -= (*moments).m0() * center;
}

void ELVIRA::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  const BasicMesh& mesh = a_U.getMesh();
  IRL::ELVIRANeighborhood neighborhood;
  neighborhood.resize(9);
  std::array<IRL::RectangularCuboid, 9> cells;
  std::array<double, 9> liquid_volume_fraction;

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double vfrac = a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (vfrac < IRL::global_constants::VF_LOW) {
        (*a_interface)(i, j).markAsAlwaysBelow();
      } else if (vfrac > IRL::global_constants::VF_HIGH) {
        (*a_interface)(i, j).markAsAlwaysAbove();
      } else {
        for (int jj = -1; jj <= 1; ++jj) {
          for (int ii = -1; ii <= 1; ++ii) {
            const double x_shift = static_cast<double>(ii);
            const double y_shift = static_cast<double>(jj);
            const IRL::UnsignedIndex_t linear_index = (jj + 1) * 3 + (ii + 1);
            cells[linear_index] = IRL::RectangularCuboid::fromBoundingPts(
                IRL::Pt(-0.5 + x_shift, -0.5 + y_shift, -0.5),
                IRL::Pt(0.5 + x_shift, 0.5 + y_shift, 0.5));
            liquid_volume_fraction[linear_index] =
                a_liquid_moments(i + ii, j + jj).m0() / mesh.cell_volume();
            neighborhood.setMember(&cells[linear_index],
                                   &liquid_volume_fraction[linear_index], ii,
                                   jj);
          }
        }
        const auto planar_separator =
            IRL::reconstructionWithELVIRA2D(neighborhood);
        const auto normal = planar_separator[0].normal();
        const auto frame =
            IRL2D::ReferenceFrame(IRL2D::Vec(normal[1], -normal[0]),
                                  IRL2D::Vec(normal[0], normal[1]));
        const auto datum = IRL2D::Vec(mesh.xm(i), mesh.ym(j));
        const auto x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
        const auto x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
        const auto cell = IRL2D::RectangleFromBounds(x0, x1);
        const auto parabola = IRL2D::Parabola(datum, frame, 0.0);
        (*a_interface)(i, j) =
            IRL2D::MatchToVolumeFraction(cell, parabola, vfrac);
      }
    }
  }
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

void LVIRA::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                              const Data<IRL2D::Moments>& a_gas_moments,
                              const double a_dt, const Data<double>& a_U,
                              const Data<double>& a_V,
                              Data<IRL2D::Parabola>* a_interface) {
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);

  const BasicMesh& mesh = a_U.getMesh();
  IRL::LVIRANeighborhood<IRL::RectangularCuboid> neighborhood;
  neighborhood.resize(9);
  std::array<IRL::RectangularCuboid, 9> cells;
  std::array<double, 9> liquid_volume_fraction;

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double vfrac = a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (vfrac < IRL::global_constants::VF_LOW) {
        (*a_interface)(i, j).markAsAlwaysBelow();
      } else if (vfrac > IRL::global_constants::VF_HIGH) {
        (*a_interface)(i, j).markAsAlwaysAbove();
      } else {
        IRL::UnsignedIndex_t ndata = 0;
        for (int jj = -1; jj <= 1; ++jj) {
          for (int ii = -1; ii <= 1; ++ii) {
            if (ii == 0 && jj == 0) {
              neighborhood.setCenterOfStencil(ndata);
            }
            const double x_shift = static_cast<double>(ii);
            const double y_shift = static_cast<double>(jj);
            const IRL::UnsignedIndex_t linear_index = (jj + 1) * 3 + (ii + 1);
            cells[ndata] = IRL::RectangularCuboid::fromBoundingPts(
                IRL::Pt(-0.5 + x_shift, -0.5 + y_shift, -0.5),
                IRL::Pt(0.5 + x_shift, 0.5 + y_shift, 0.5));
            liquid_volume_fraction[ndata] =
                a_liquid_moments(i + ii, j + jj).m0() / mesh.cell_volume();
            neighborhood.setMember(ndata, &cells[ndata],
                                   &liquid_volume_fraction[ndata]);
            ndata++;
          }
        }
        const auto guess_datum = (*a_interface)(i, j).datum();
        const auto guess_frame = (*a_interface)(i, j).frame();
        const auto planar_guess = IRL::PlanarSeparator::fromOnePlane(
            IRL::Plane(IRL::Normal(guess_frame[1][0], guess_frame[1][1], 0.0),
                       guess_datum.magnitude()));
        const auto planar_separator =
            IRL::reconstructionWithLVIRA2D(neighborhood, planar_guess);
        const auto normal = planar_separator[0].normal();
        const auto frame =
            IRL2D::ReferenceFrame(IRL2D::Vec(normal[1], -normal[0]),
                                  IRL2D::Vec(normal[0], normal[1]));
        const auto datum = IRL2D::Vec(mesh.xm(i), mesh.ym(j));
        const auto x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
        const auto x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
        const auto cell = IRL2D::RectangleFromBounds(x0, x1);
        const auto parabola = IRL2D::Parabola(datum, frame, 0.0);
        (*a_interface)(i, j) =
            IRL2D::MatchToVolumeFraction(cell, parabola, vfrac);
      }
    }
  }
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// void Jibben::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
//                                const Data<IRL2D::Moments>& a_gas_moments,
//                                const double a_dt, const Data<double>& a_U,
//                                const Data<double>& a_V,
//                                Data<IRL2D::Parabola>* a_interface) {

//   ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
//                             a_interface);

//   const BasicMesh& mesh = a_U.getMesh();

//   // flag for mixed cells
//   Data<int> band(&mesh);
//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       band(i, j) = 0;
//       const double liquid_volume_fraction =
//           (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
//       if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
//           liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
//         band(i, j) = 1;
//       }
//     }
//   }

//   // Jibben parabola estimation
//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       if (band(i, j) == 1) {
//         IRL2D::Parabola target_interface = (*a_interface)(i,j);
//         IRL2D::BezierList target_cell = IRL2D::RectangleFromBounds(
//           IRL2D::Vec(mesh.x(i), mesh.y(j)),
//           IRL2D::Vec(mesh.x(i+1), mesh.y(j+1))
//         );
//         IRL2D::Vec n_target = target_interface.frame()[1];
//         std::vector<IRL2D::Parabola> all_interfaces;
//         std::vector<IRL2D::BezierList> all_cells;
//         for (int jj = -2; jj <= 2; jj++){
//           for (int ii = -2; ii <= 2; ii++){
//             if (band(ii + i, jj + j) == 1){
//               IRL2D::Vec n_neighbor = (*a_interface)(ii + i, jj +
//               j).frame()[1]; if (n_neighbor * n_target <= -0.5){
//                 continue; // filtering based on normals
//               }
//               all_interfaces.push_back((*a_interface)(ii + i, jj + j));
//               all_cells.push_back(IRL2D::RectangleFromBounds(
//                 IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
//                 IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))
//               ));
//             }
//           }
//         }
//         IRL2D::Parabola parabolaJibben = IRL2D::getParabolaJibben(
//           target_interface, target_cell, all_interfaces, all_cells
//         );
//         double vfrac = (a_liquid_moments)(i,j).m0() /
//         IRL2D::ComputeArea(target_cell);
//         (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(target_cell,
//         parabolaJibben, vfrac);
//       }
//     }
//   }

//   a_interface->updateBorder();
//   correctInterfaceBorders(a_interface);
// }

double f_parabola(double x, double A, double B, double C) {
  return A * x * x + B * x + C;
}

double fp_parabola(double x, double A, double B) { return 2.0 * A * x + B; }

double fpp_parabola(double A) { return 2.0 * A; }

struct ClosestParabolaPointFunctor {
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // variables
  double A, B, C;

  // constructor
  ClosestParabolaPointFunctor(const double& A_, const double& B_,
                              const double& C_)
      : A(A_), B(B_), C(C_) {}

  int inputs() const { return 1; }
  int values() const { return 1; }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    double xp = x(0);
    double yp = A * xp * xp + B * xp + C;
    fvec(0) = std::sqrt(xp * xp + yp * yp);
    return 0;
  }
};

void Jibben::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);

  // storing PLIC
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();

  // local <-> global
  auto globalToLocal = [&](const IRL2D::Vec& p,
                           const IRL2D::Parabola& interface) -> IRL2D::Vec {
    IRL2D::Vec t = interface.frame()[0], n = interface.frame()[1];
    IRL2D::Vec ploc = {(p.x() - interface.datum().x()) * t.x() +
                           (p.y() - interface.datum().y()) * t.y(),
                       (p.x() - interface.datum().x()) * n.x() +
                           (p.y() - interface.datum().y()) * n.y()};
    return ploc;
  };
  auto localToGlobal = [&](const IRL2D::Vec& ploc,
                           const IRL2D::Parabola& interface) -> IRL2D::Vec {
    IRL2D::Vec t = interface.frame()[0], n = interface.frame()[1];
    IRL2D::Vec p = {
        interface.datum().x() + t.x() * ploc.x() + n.x() * ploc.y(),
        interface.datum().y() + t.y() * ploc.x() + n.y() * ploc.y()};
    return p;
  };
  auto vectorGlobalToLocal =
      [&](const IRL2D::Vec& v, const IRL2D::Parabola& interface) -> IRL2D::Vec {
    IRL2D::Vec t = interface.frame()[0];
    IRL2D::Vec n = interface.frame()[1];
    IRL2D::Vec vloc = {v.x() * t.x() + v.y() * t.y(),
                       v.x() * n.x() + v.y() * n.y()};
    return vloc;
  };
  auto vectorLocalToGlobal =
      [&](const IRL2D::Vec& vloc,
          const IRL2D::Parabola& interface) -> IRL2D::Vec {
    IRL2D::Vec t = interface.frame()[0], n = interface.frame()[1];
    IRL2D::Vec v = {t.x() * vloc.x() + n.x() * vloc.y(),
                    t.y() * vloc.x() + n.y() * vloc.y()};
    return v;
  };

  // flag for mixed cells
  Data<int> band(&mesh);
  Data<double> vf(&mesh);
  Data<IRL2D::Vec> plic_center(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
        vf(i, j) = liquid_volume_fraction;
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(cell, plic(i, j), true);
        plic_center(i, j) =
            0.5 * (clipped_plic[0].first + clipped_plic[1].first);
      }
    }
  }

  // Jibben parabola estimation
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (band(i, j) == 1) {
        IRL2D::Parabola target_interface = plic(i, j);
        IRL2D::BezierList target_cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        IRL2D::Vec pref_local =
            globalToLocal(plic_center(i, j), target_interface);
        IRL2D::Vec n_target = target_interface.frame()[1];
        IRL2D::Vec nref_local = vectorGlobalToLocal(n_target, target_interface);
        std::vector<double> weights;
        std::vector<IRL2D::NeighborInfo> neighbors;
        for (int jj = -2; jj <= 2; jj++) {
          for (int ii = -2; ii <= 2; ii++) {
            if (band(ii + i, jj + j) == 1) {
              IRL2D::Vec n_neighbor = plic(ii + i, jj + j).frame()[1];
              if (n_neighbor * n_target <= -0.5) {
                continue;
              }
              IRL2D::Vec n_local =
                  vectorGlobalToLocal(n_neighbor, target_interface);
              IRL2D::Vec ploc_local =
                  globalToLocal(plic_center(i + ii, j + jj), target_interface);
              double vf_weight = IRL2D::getVfracWeight(vf(ii + i, jj + j));
              double d_weight =
                  IRL2D::getDistanceWeight(pref_local, ploc_local, mesh.dx());
              double n_weight = IRL2D::getNormalWeight(nref_local, n_local);
              weights.push_back(vf_weight * d_weight * n_weight);
              neighbors.push_back(
                  {plic(ii + i, jj + j),
                   IRL2D::RectangleFromBounds(
                       IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
                       IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))),
                   ii + i, jj + j, vf(ii + i, jj + j)});
            }
          }
        }

        // coefficients of parabola in local frame
        IRL2D::Parabola parabolaJibben;
        std::vector<double> jibbenCoeffs = IRL2D::getJibbenCoeffs(
            target_interface, target_cell, neighbors, weights);
        double A = jibbenCoeffs[0], B = jibbenCoeffs[1], C = jibbenCoeffs[2];

        // datum
        ClosestParabolaPointFunctor functor(A, B, C);
        Eigen::NumericalDiff<ClosestParabolaPointFunctor> numDiff(functor);
        Eigen::LevenbergMarquardt<
            Eigen::NumericalDiff<ClosestParabolaPointFunctor>, double>
            lm(numDiff);
        Eigen::VectorXd x(1);
        x(0) = 0.0;
        lm.parameters.maxfev = 1000;
        lm.parameters.xtol = 1e-12;
        lm.minimize(x);
        double x_star = x(0);
        double y_star = f_parabola(x_star, A, B, C);
        parabolaJibben.datum() =
            localToGlobal({x_star, y_star}, target_interface);

        // curvature
        double fp = fp_parabola(x_star, A, B);
        double fpp = fpp_parabola(A);
        double curvature = -fpp / std::pow(1.0 + fp * fp, 1.5);
        parabolaJibben.coeff() = 0.5 * curvature;

        // reference frame
        IRL2D::Vec normal_local = {-fp, 1.0};
        normal_local.normalize();
        IRL2D::Vec normal = vectorLocalToGlobal(normal_local, target_interface);
        if ((normal * target_interface.frame()[1]) < 0) {
          normal *= -1.0;
          parabolaJibben.coeff() *= -1.0;
        }
        IRL2D::Vec tangent = {normal.y(), -normal.x()};
        parabolaJibben.frame() = {tangent, normal};

        // reverting black to plane is curvature is too large
        const double maxkdx = 4.0;
        const double length_scale = std::sqrt(ComputeArea(target_cell));
        const double kdx = 2.0 * parabolaJibben.coeff() * length_scale;
        if (std::abs(kdx) > maxkdx) {
          parabolaJibben = plic(i, j);
        }

        // IRL2D::Parabola parabolaJibben = IRL2D::getParabolaJibben(
        // target_interface, target_cell, neighbors, i, j);
        // parabolaJibben.datum() = IRL2D::ComputeMoments(target_cell).m1() /
        // IRL2D::ComputeMoments(target_cell).m0();

        // vf matching
        (*a_interface)(i, j) =
            IRL2D::MatchToVolumeFraction(target_cell, parabolaJibben, vf(i, j));
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// void Jibben::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
//                               const Data<IRL2D::Moments>& a_gas_moments,
//                               const double a_dt, const Data<double>& a_U,
//                               const Data<double>& a_V,
//                               Data<IRL2D::Parabola>* a_interface) {

//   LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
//                             a_interface);

//   const BasicMesh& mesh = a_U.getMesh();

//   // band for mixed cells
//   Data<int> band(&mesh);
//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       band(i, j) = 0;
//       const double liquid_volume_fraction =
//           (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
//       if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
//           liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
//         band(i, j) = 1;
//       }
//     }
//   }

//   auto getConnectedInterfaceNeighbors = [&](int i, int j, const IRL2D::Vec&
//   n_target) {
//     std::vector<std::pair<int, int>> connected;
//     std::queue<std::pair<int, int>> q;
//     std::set<std::pair<int, int>> visited;

//     q.push({i, j});
//     visited.insert({i, j});

//     while (!q.empty()) {
//       auto [ci, cj] = q.front();
//       q.pop();
//       connected.emplace_back(ci, cj);

//       for (int dj = -1; dj <= 1; ++dj) {
//         for (int di = -1; di <= 1; ++di) {
//           if (di == 0 && dj == 0) continue;

//           int ni = ci + di;
//           int nj = cj + dj;
//           if (visited.count({ni, nj})) continue;
//           if (ni < mesh.imin() || ni > mesh.imax() || nj < mesh.jmin() || nj
//           > mesh.jmax()) continue; if (std::abs(ni - i) > 2 || std::abs(nj -
//           j) > 2) continue;

//           if (band(ni, nj) == 1) {
//             IRL2D::Vec n_neighbor = (*a_interface)(ni, nj).frame()[1];
//             if (n_neighbor * n_target > -0.5) {
//               visited.insert({ni, nj});
//               q.push({ni, nj});
//             }
//           }
//         }
//       }
//     }

//     return connected;
//   };

//   // Jibben parabola estimation
//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       if (band(i, j) == 1) {
//         IRL2D::Parabola target_interface = (*a_interface)(i,j);
//         IRL2D::BezierList target_cell = IRL2D::RectangleFromBounds(
//           IRL2D::Vec(mesh.x(i), mesh.y(j)),
//           IRL2D::Vec(mesh.x(i+1), mesh.y(j+1))
//         );
//         IRL2D::Vec n_target = target_interface.frame()[1];
//         std::vector<IRL2D::Parabola> all_interfaces;
//         std::vector<IRL2D::BezierList> all_cells;
//         auto neighbors = getConnectedInterfaceNeighbors(i, j, n_target);
//         for (const auto& [ii, jj] : neighbors) {
//           all_interfaces.push_back((*a_interface)(ii, jj));
//           all_cells.push_back(IRL2D::RectangleFromBounds(
//               IRL2D::Vec(mesh.x(ii), mesh.y(jj)),
//               IRL2D::Vec(mesh.x(ii + 1), mesh.y(jj + 1))));
//         }
//         IRL2D::Parabola parabolaJibben = IRL2D::getParabolaJibben(
//           target_interface, target_cell, all_interfaces, all_cells
//         );
//         double vfrac = (a_liquid_moments)(i,j).m0() /
//         IRL2D::ComputeArea(target_cell);
//         (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(target_cell,
//         parabolaJibben, vfrac);
//       }
//     }
//   }

//   a_interface->updateBorder();
//   correctInterfaceBorders(a_interface);
// }

void Particle::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                                 const Data<IRL2D::Moments>& a_gas_moments,
                                 const double a_dt, const Data<double>& a_U,
                                 const Data<double>& a_V,
                                 Data<IRL2D::Parabola>* a_interface) {
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);

  // storing PLIC
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();

  // band for mixed cells
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }

  int N = 5;
  double Hp = 3.0;
  double h = mesh.dx();
  double eta = 0.5;

  // estimation curvature using particles
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (band(i, j) == 1) {
        IRL2D::Parabola target_interface = plic(i, j);
        IRL2D::BezierList target_cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        IRL2D::Vec n_target = target_interface.frame()[1];
        std::vector<IRL2D::Parabola> all_interfaces;
        std::vector<IRL2D::BezierList> all_cells;
        for (int jj = -2; jj <= 2; jj++) {
          for (int ii = -2; ii <= 2; ii++) {
            if (band(ii + i, jj + j) == 1) {
              IRL2D::Vec n_neighbor = plic(ii + i, jj + j).frame()[1];
              if (n_neighbor * n_target <= -0.5) {
                continue;  // filtering based on normals
              }
              all_interfaces.push_back(plic(ii + i, jj + j));
              all_cells.push_back(IRL2D::RectangleFromBounds(
                  IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
                  IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))));
            }
          }
        }
        IRL2D::Parabola parabolaParticle = target_interface;
        double coefficient =
            0.5 * IRL2D::getCurvature(target_interface, target_cell,
                                      all_interfaces, all_cells, N, Hp, h, eta);
        // std::cout << "Coefficient = " << coefficient << std::endl;
        parabolaParticle.coeff() = coefficient;
        double vfrac =
            (a_liquid_moments)(i, j).m0() / IRL2D::ComputeArea(target_cell);
        (*a_interface)(i, j) =
            IRL2D::MatchToVolumeFraction(target_cell, parabolaParticle, vfrac);

        // curvature check
        const double maxkdx = 4.0;
        const double length_scale = std::sqrt(ComputeArea(target_cell));
        const double kdx = 2.0 * parabolaParticle.coeff() * length_scale;
        if (std::abs(kdx) > maxkdx) {
          std::cout << "Warning: Curvature too large.\n";
          // parabolaJibben.coeff() = 0.0;
        }
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// void Particle::getReconstruction(const Data<IRL2D::Moments>&
// a_liquid_moments,
//                                   const Data<IRL2D::Moments>& a_gas_moments,
//                                   const double a_dt, const Data<double>& a_U,
//                                   const Data<double>& a_V,
//                                   Data<IRL2D::Parabola>* a_interface) {

//   LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
//                             a_interface);

//   const BasicMesh& mesh = a_U.getMesh();
//   Data<int> band(&mesh);

//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       band(i, j) = 0;
//       const double vfrac = (a_liquid_moments)(i, j).m0() /
//       mesh.cell_volume(); if (vfrac >= IRL::global_constants::VF_LOW && vfrac
//       <= IRL::global_constants::VF_HIGH) {
//         band(i, j) = 1;
//       }
//     }
//   }

//   auto getConnectedInterfaceNeighbors = [&](int i, int j, const IRL2D::Vec&
//   n_target) {
//     std::vector<std::pair<int, int>> connected;
//     std::queue<std::pair<int, int>> q;
//     std::set<std::pair<int, int>> visited;

//     q.push({i, j});
//     visited.insert({i, j});

//     while (!q.empty()) {
//       auto [ci, cj] = q.front();
//       q.pop();
//       connected.emplace_back(ci, cj);

//       for (int dj = -1; dj <= 1; ++dj) {
//         for (int di = -1; di <= 1; ++di) {
//           if (di == 0 && dj == 0) continue;

//           int ni = ci + di;
//           int nj = cj + dj;
//           if (visited.count({ni, nj})) continue;
//           if (ni < mesh.imin() || ni > mesh.imax() || nj < mesh.jmin() || nj
//           > mesh.jmax()) continue; if (std::abs(ni - i) > 2 || std::abs(nj -
//           j) > 2) continue;

//           if (band(ni, nj) == 1) {
//             IRL2D::Vec n_neighbor = (*a_interface)(ni, nj).frame()[1];
//             if (n_neighbor * n_target > -0.5) {
//               visited.insert({ni, nj});
//               q.push({ni, nj});
//             }
//           }
//         }
//       }
//     }

//     return connected;
//   };

//   int N = 5;
//   double Hp = 3.0;
//   double h = mesh.dx();
//   double eta = 1.0;

//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       if (band(i, j) == 1) {
//         IRL2D::Parabola target_interface = (*a_interface)(i, j);
//         IRL2D::BezierList target_cell = IRL2D::RectangleFromBounds(
//             IRL2D::Vec(mesh.x(i), mesh.y(j)),
//             IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
//         IRL2D::Vec n_target = target_interface.frame()[1];

//         std::vector<IRL2D::Parabola> all_interfaces;
//         std::vector<IRL2D::BezierList> all_cells;

//         auto neighbors = getConnectedInterfaceNeighbors(i, j, n_target);
//         for (const auto& [ii, jj] : neighbors) {
//           all_interfaces.push_back((*a_interface)(ii, jj));
//           all_cells.push_back(IRL2D::RectangleFromBounds(
//               IRL2D::Vec(mesh.x(ii), mesh.y(jj)),
//               IRL2D::Vec(mesh.x(ii + 1), mesh.y(jj + 1))));
//         }

//         IRL2D::Parabola parabolaParticle = target_interface;
//         double curvature = IRL2D::getCurvature(target_interface, target_cell,
//                                               all_interfaces, all_cells, N,
//                                               Hp, h, eta);
//         parabolaParticle.coeff() = 0.5 * curvature;

//         double vfrac = (a_liquid_moments)(i, j).m0() /
//         IRL2D::ComputeArea(target_cell);
//         (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(target_cell,
//         parabolaParticle, vfrac);
//       }
//     }
//   }

//   a_interface->updateBorder();
//   correctInterfaceBorders(a_interface);
// }

void JPHybrid::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                                 const Data<IRL2D::Moments>& a_gas_moments,
                                 const double a_dt, const Data<double>& a_U,
                                 const Data<double>& a_V,
                                 Data<IRL2D::Parabola>* a_interface) {
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);

  // storing PLIC
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();

  // band for mixed cells
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }

  int N = 5;
  double Hp = 3.0;
  double h = mesh.dx();
  double eta = 0.5;

  // estimation curvature using particles
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (band(i, j) == 1) {
        IRL2D::Parabola target_interface = plic(i, j);
        IRL2D::BezierList target_cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        IRL2D::Vec n_target = target_interface.frame()[1];
        std::vector<IRL2D::Parabola> all_interfaces;
        std::vector<IRL2D::BezierList> all_cells;
        for (int jj = -2; jj <= 2; jj++) {
          for (int ii = -2; ii <= 2; ii++) {
            if (band(ii + i, jj + j) == 1) {
              IRL2D::Vec n_neighbor = plic(ii + i, jj + j).frame()[1];
              if (n_neighbor * n_target <= -0.5) {
                continue;  // filtering based on normals
              }
              all_interfaces.push_back(plic(ii + i, jj + j));
              all_cells.push_back(IRL2D::RectangleFromBounds(
                  IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
                  IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))));
            }
          }
        }
        IRL2D::Parabola parabolaParticle = target_interface;
        double coefficient =
            0.5 * IRL2D::getCurvature(target_interface, target_cell,
                                      all_interfaces, all_cells, N, Hp, h, eta);
        // std::cout << "Coefficient = " << coefficient << std::endl;
        parabolaParticle.coeff() = coefficient;
        double vfrac =
            (a_liquid_moments)(i, j).m0() / IRL2D::ComputeArea(target_cell);
        (*a_interface)(i, j) =
            IRL2D::MatchToVolumeFraction(target_cell, parabolaParticle, vfrac);

        // curvature check
        const double maxkdx = 4.0;
        const double length_scale = std::sqrt(ComputeArea(target_cell));
        const double kdx = 2.0 * parabolaParticle.coeff() * length_scale;
        if (std::abs(kdx) > maxkdx) {
          std::cout << "Warning: Curvature too large.\n";
          // parabolaJibben.coeff() = 0.0;
        }
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

struct LVIRAQFunctor {
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
  std::array<std::array<double, 3>, 3> m_vfracs;
  std::array<std::array<IRL2D::BezierList, 3>, 3> m_cells;
  IRL2D::Vec m_datum;
  IRL2D::ReferenceFrame m_frame;
  double m_coeff;
  double m_length_scale;

  // Constructor
  LVIRAQFunctor(int inputs, int values,
                const std::array<std::array<IRL2D::BezierList, 3>, 3>& cells,
                std::array<std::array<double, 3>, 3>& vfracs)
      : m_inputs(inputs), m_values(values), m_cells(cells), m_vfracs(vfracs) {
    m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cells[1][1]));
  }

  void setframe(const IRL2D::Parabola& guess_paraboloid) {
    m_coeff = guess_paraboloid.coeff();
    m_frame = guess_paraboloid.frame();
    m_datum = IRL2D::ComputeMoments(m_cells[1][1]).m1() /
              IRL2D::ComputeArea(m_cells[1][1]);
    // If error is too large, revert to planar initial guess
    const double maxkdx = 4.0;
    const double kdx = 2.0 * m_coeff * m_length_scale;
    if (std::abs(kdx) > maxkdx) {
      m_coeff = 0.0;
    }
  }

  const IRL2D::Parabola getparabola(const Eigen::VectorXd& x) const {
    const auto rotation = IRL2D::ReferenceFrame(x(0) * M_PI);
    const auto new_frame =
        IRL2D::ReferenceFrame(rotation * m_frame[0], rotation * m_frame[1]);
    const double new_coeff = m_coeff + x(1) / m_length_scale;
    const auto parabola = IRL2D::Parabola(m_datum, new_frame, new_coeff);
    return IRL2D::MatchToVolumeFraction(m_cells[1][1], parabola,
                                        m_vfracs[1][1]);
  }

  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto parabola = this->getparabola(x);
    int count = 0;
    for (int jj = 0; jj < 3; ++jj) {
      for (int ii = 0; ii < 3; ++ii) {
        fvec(count++) =
            IRL2D::ComputeVFrac(m_cells[ii][jj], parabola) - m_vfracs[ii][jj];
      }
    }
    // Penalty to prevent kappa * dx > 6
    const double mu = 50.0;
    const double kdx = 2.0 * m_coeff * m_length_scale + x(1);
    const double maxkdx = 4.0;
    fvec(count++) = mu * std::max(0.0, std::abs(kdx) - maxkdx);
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    this->errorvec(x, fvec);
    return 0;
  }

  // Get methods for users to determine function input and output dimensions
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

void LVIRAQ::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  // First guess with ELVIRA
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);

  // Now fit parabola LVIRA-style
  const BasicMesh& mesh = a_U.getMesh();

#ifdef USE_MPI
  const double cell_volume = mesh.cell_volume();
  int nmixed_global = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / cell_volume;
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        nmixed_global++;
      }
    }
  }

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  IRL2D::Parabola dummy_par;
  IRL::ByteBuffer dummy_buffer;
  dummy_buffer.resize(0);
  dummy_buffer.resetBufferPointer();
  IRL::serializeAndPack(dummy_par, &dummy_buffer);
  const int size_parabola = dummy_buffer.size();

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
  interface_local.resize(nmixed_local * sizeof(IRL2D::Parabola));
  interface_global.resize(0);
  interface_local.resetBufferPointer();
  interface_global.resetBufferPointer();

  int count = 0;
#endif

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double vfrac = a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (vfrac < IRL::global_constants::VF_LOW) {
        (*a_interface)(i, j).markAsAlwaysBelow();
      } else if (vfrac > IRL::global_constants::VF_HIGH) {
        (*a_interface)(i, j).markAsAlwaysAbove();
      } else {
#ifdef USE_MPI
        if (count >= proc_offset[rank] && count < proc_offset[rank + 1]) {
#endif
          // Fill stencil of moments
          std::array<std::array<double, 3>, 3> vfracs;
          std::array<std::array<IRL2D::BezierList, 3>, 3> cells;
          for (int jj = 0; jj < 3; ++jj) {
            for (int ii = 0; ii < 3; ++ii) {
              vfracs[ii][jj] = a_liquid_moments(i + ii - 1, j + jj - 1).m0() /
                               mesh.cell_volume();
              const auto x0 =
                  IRL2D::Vec(mesh.x(i + ii - 1), mesh.y(j + jj - 1));
              const auto x1 = IRL2D::Vec(mesh.x(i + ii), mesh.y(j + jj));
              cells[ii][jj] = IRL2D::RectangleFromBounds(x0, x1);
            }
          }

          // Create functor for LM minimization
          LVIRAQFunctor myLVIRAQFunctor(2, 10, cells, vfracs);
          myLVIRAQFunctor.setframe((*a_interface)(i, j));
          Eigen::NumericalDiff<LVIRAQFunctor> NDLVIRAQFunctor(myLVIRAQFunctor);
          Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LVIRAQFunctor>, double>
              LVIRAQ_LM(NDLVIRAQFunctor);
          // LVIRAQ_LM.parameters.ftol = 1.0e-8;
          // LVIRAQ_LM.parameters.xtol = 1.0e-8;
          // LVIRAQ_LM.parameters.factor = 1.0;
          // LVIRAQ_LM.parameters.maxfev = 1000;  // Max
          Eigen::VectorXd x(2);
          x.setZero();
          Eigen::LevenbergMarquardtSpace::Status status =
              LVIRAQ_LM.minimizeInit(x);
          do {
            status = LVIRAQ_LM.minimizeOneStep(x);
          } while (status == Eigen::LevenbergMarquardtSpace::Running);
          const auto parabola = IRL2D::MatchToVolumeFractionBisection(
              cells[1][1], myLVIRAQFunctor.getparabola(x), vfracs[1][1], 500);
#ifdef USE_MPI
          IRL::serializeAndPack(parabola, &interface_local);
        }
        count++;
#else
        (*a_interface)(i, j) = parabola;
#endif
      }
    }
  }

#ifdef USE_MPI
  std::vector<int> proc_count(size);
  for (int r = 0; r < size; r++) {
    proc_count[r] = size_parabola * (proc_offset[r + 1] - proc_offset[r]);
    proc_offset[r] = size_parabola * proc_offset[r];
  }

  interface_global.resize(size_parabola * nmixed_global);
  MPI_Allgatherv(interface_local.data(), size_parabola * nmixed_local, MPI_BYTE,
                 interface_global.data(), proc_count.data(), proc_offset.data(),
                 MPI_BYTE, MPI_COMM_WORLD);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / cell_volume;
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        IRL2D::Parabola parabola;
        IRL::unpackAndStore(&parabola, &interface_global);
        (*a_interface)(i, j) = parabola;
      }
    }
  }
#endif

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// Functor for MOF2 with least squares
struct LMIRAQFunctor {
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
  std::array<std::array<IRL2D::Moments, 3>, 3> m_liq_moments;
  std::array<std::array<IRL2D::Moments, 3>, 3> m_gas_moments;
  std::array<std::array<double, 3>, 3> m_vfracs;
  std::array<std::array<IRL2D::Vec, 3>, 3> m_liq_centroids;
  std::array<std::array<IRL2D::Vec, 3>, 3> m_gas_centroids;
  std::array<std::array<IRL2D::Mat, 3>, 3> m_liq_m2s;
  std::array<std::array<IRL2D::Mat, 3>, 3> m_gas_m2s;
  std::array<std::array<IRL2D::BezierList, 3>, 3> m_cells;
  IRL2D::Vec m_datum;
  IRL2D::ReferenceFrame m_frame;
  double m_coeff;
  double m_length_scale;

  // Constructor
  LMIRAQFunctor(int inputs, int values,
                const std::array<std::array<IRL2D::BezierList, 3>, 3>& cells,
                const std::array<std::array<IRL2D::Moments, 3>, 3>& liq_moments,
                const std::array<std::array<IRL2D::Moments, 3>, 3>& gas_moments)
      : m_inputs(inputs),
        m_values(values),
        m_cells(cells),
        m_liq_moments(liq_moments),
        m_gas_moments(gas_moments) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        m_vfracs[i][j] =
            liq_moments[i][j].m0() / IRL2D::ComputeArea(cells[i][j]);
        m_liq_centroids[i][j] = liq_moments[i][j].m1() / liq_moments[i][j].m0();
        m_gas_centroids[i][j] = gas_moments[i][j].m1() / gas_moments[i][j].m0();
        RecenterMoments(&m_liq_moments[i][j], m_liq_centroids[i][j]);
        RecenterMoments(&m_gas_moments[i][j], m_gas_centroids[i][j]);
        m_liq_m2s[i][j] = m_liq_moments[i][j].m2();
        m_gas_m2s[i][j] = m_gas_moments[i][j].m2();
      }
    }
    m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cells[1][1]));
  }

  void setframe(const IRL2D::Parabola& guess_parabola) {
    m_coeff = guess_parabola.coeff();
    m_frame = guess_parabola.frame();
    m_datum = guess_parabola.datum();

    // check for curvature
    const double maxkdx = 4.0;
    const double kdx = 2.0 * m_coeff * m_length_scale;
    if (std::abs(kdx) > maxkdx) {
      m_coeff = 0.0;  // plane
    }
  }

  const IRL2D::Parabola getparabola(const Eigen::VectorXd& x) const {
    const auto rotation = IRL2D::ReferenceFrame(x(0) * M_PI);
    const auto new_frame =
        IRL2D::ReferenceFrame(rotation * m_frame[0], rotation * m_frame[1]);
    const double new_coeff = m_coeff + x(1) / m_length_scale;
    const auto parabola = IRL2D::Parabola(m_datum, new_frame, new_coeff);
    return IRL2D::MatchToVolumeFraction(m_cells[1][1], parabola,
                                        m_vfracs[1][1]);
  }

  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto parabola = this->getparabola(x);
    int count = 0;

    for (int jj = 0; jj < 3; ++jj) {
      for (int ii = 0; ii < 3; ++ii) {
        const auto& cell = m_cells[ii][jj];
        const double area = IRL2D::ComputeArea(cell);

        // liq and gas moments
        IRL2D::Moments lm = IRL2D::ComputeMoments(cell, parabola);
        IRL2D::Moments gm = IRL2D::ComputeMoments(cell) - lm;

        // Centroids
        IRL2D::Vec lc = lm.m1() / lm.m0();
        IRL2D::Vec gc = gm.m1() / gm.m0();

        // Recentering moments
        RecenterMoments(&lm, m_liq_centroids[ii][jj]);
        RecenterMoments(&gm, m_gas_centroids[ii][jj]);

        // second moments
        IRL2D::Mat lm2 = lm.m2();
        IRL2D::Mat gm2 = gm.m2();

        const double denom = std::pow(m_length_scale, 2.0);

        // Volume fraction
        fvec(count++) = (lm.m0() / area) - m_vfracs[ii][jj];

        // Centroids
        fvec(count++) = (lc[0] - m_liq_centroids[ii][jj][0]) / m_length_scale;
        fvec(count++) = (lc[1] - m_liq_centroids[ii][jj][1]) / m_length_scale;
        fvec(count++) = (gc[0] - m_gas_centroids[ii][jj][0]) / m_length_scale;
        fvec(count++) = (gc[1] - m_gas_centroids[ii][jj][1]) / m_length_scale;

        // second moment
        fvec(count++) =
            (lm2[0][0] - m_liq_m2s[ii][jj][0][0]) / (lm.m0() * denom);
        fvec(count++) =
            (lm2[1][0] - m_liq_m2s[ii][jj][1][0]) / (lm.m0() * denom);
        fvec(count++) =
            (lm2[1][1] - m_liq_m2s[ii][jj][1][1]) / (lm.m0() * denom);
        fvec(count++) =
            (gm2[0][0] - m_gas_m2s[ii][jj][0][0]) / (gm.m0() * denom);
        fvec(count++) =
            (gm2[1][0] - m_gas_m2s[ii][jj][1][0]) / (gm.m0() * denom);
        fvec(count++) =
            (gm2[1][1] - m_gas_m2s[ii][jj][1][1]) / (gm.m0() * denom);
      }
    }

    // Curvature penalty
    const double mu = 50.0;
    const double kdx = 2.0 * m_coeff * m_length_scale + x(1);
    const double maxkdx = 4.0;
    fvec(count++) = mu * std::max(0.0, std::abs(kdx) - maxkdx);
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    this->errorvec(x, fvec);
    return 0;
  }

  // Get methods for users to determine function input and output dimensions
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

void LMIRAQ::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  // PLIC for curvature estimation
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);

  const BasicMesh& mesh = a_U.getMesh();

  // finding mixed cells
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }

  // parameters for curvature estimation using particles (initial guess)
  int N = 5;
  double Hp = 3.0;
  double h = mesh.dx();
  double eta = 0.5;

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (band(i, j) == 1) {
        // current cell
        IRL2D::Vec x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
        IRL2D::Vec x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
        IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(x0, x1);

        // initial guess using particle
        IRL2D::Parabola target_interface = (*a_interface)(i, j);
        IRL2D::Vec n_target = target_interface.frame()[1];
        std::vector<IRL2D::Parabola> all_interfaces;
        std::vector<IRL2D::BezierList> all_cells;
        for (int jj = -2; jj <= 2; jj++) {
          for (int ii = -2; ii <= 2; ii++) {
            if (band(ii + i, jj + j) == 1) {
              IRL2D::Vec n_neighbor = (*a_interface)(ii + i, jj + j).frame()[1];
              if (n_neighbor * n_target <= 0) {
                continue;
              }
              all_interfaces.push_back((*a_interface)(ii + i, jj + j));
              all_cells.push_back(IRL2D::RectangleFromBounds(
                  IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
                  IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))));
            }
          }
        }
        IRL2D::Parabola parabola_guess = target_interface;
        parabola_guess.coeff() =
            0.5 * IRL2D::getCurvature(target_interface, rectangle,
                                      all_interfaces, all_cells, N, Hp, h, eta);

        // filling stencil for functor inputs
        std::array<std::array<IRL2D::BezierList, 3>, 3> cells;
        std::array<std::array<IRL2D::Moments, 3>, 3> liq_moments;
        std::array<std::array<IRL2D::Moments, 3>, 3> gas_moments;
        for (int jj = 0; jj < 3; ++jj) {
          for (int ii = 0; ii < 3; ++ii) {
            const auto x0 = IRL2D::Vec(mesh.x(i + ii - 1), mesh.y(j + jj - 1));
            const auto x1 = IRL2D::Vec(mesh.x(i + ii), mesh.y(j + jj));
            cells[ii][jj] = IRL2D::RectangleFromBounds(x0, x1);
            liq_moments[ii][jj] = a_liquid_moments(i + ii - 1, j + jj - 1);
            gas_moments[ii][jj] = a_gas_moments(i + ii - 1, j + jj - 1);
          }
        }

        // for LM solver
        int num_inputs = 2;
        int num_eq = 100;

        Eigen::VectorXd x(num_inputs);
        x.setZero();

        LMIRAQFunctor myLMIRAQFunctor(num_inputs, num_eq, cells, liq_moments,
                                      gas_moments);
        myLMIRAQFunctor.setframe(parabola_guess);
        Eigen::NumericalDiff<LMIRAQFunctor> numericalDiffMyFunctor(
            myLMIRAQFunctor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LMIRAQFunctor>, double>
            LM(numericalDiffMyFunctor);
        LM.minimize(x);

        (*a_interface)(i, j) = myLMIRAQFunctor.getparabola(x);
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// Functor for MOF
struct MOF1Functor {
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // variables
  const IRL2D::BezierList& m_cell;
  const IRL2D::Vec m_liq_centroid_star;
  const IRL2D::Vec m_gas_centroid_star;
  const double m_liq_f_star;
  IRL2D::Vec m_datum;
  IRL2D::ReferenceFrame m_frame;
  double m_coeff;
  double m_length_scale;

  MOF1Functor(const IRL2D::BezierList& cell,
              const IRL2D::Vec& liq_centroid_star,
              const IRL2D::Vec& gas_centroid_star, const double liq_f_star)
      : m_cell(cell),
        m_liq_centroid_star(liq_centroid_star),
        m_gas_centroid_star(gas_centroid_star),
        m_liq_f_star(liq_f_star) {
    m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cell));
  }

  void setframe(const IRL2D::Parabola& guess_plane) {
    m_coeff = 0.0;
    m_frame = guess_plane.frame();
    m_datum = IRL2D::ComputeMoments(m_cell).m1() / IRL2D::ComputeArea(m_cell);
    // m_datum = guess_plane.datum();
  }

  const IRL2D::Parabola getplane(const Eigen::VectorXd& x) const {
    const auto rotation = IRL2D::ReferenceFrame(x(0));
    const auto new_frame =
        IRL2D::ReferenceFrame(rotation * m_frame[0], rotation * m_frame[1]);
    const auto plane = IRL2D::Parabola(m_datum, new_frame, 0.0);
    return IRL2D::MatchToVolumeFraction(m_cell, plane, m_liq_f_star);
  }

  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto plane = this->getplane(x);
    IRL2D::Moments liq_moments = IRL2D::ComputeMoments(m_cell, plane);
    IRL2D::Moments gas_moments = IRL2D::ComputeMoments(m_cell) - liq_moments;
    IRL2D::Vec liq_centroid_h = liq_moments.m1() / liq_moments.m0();
    IRL2D::Vec gas_centroid_h = gas_moments.m1() / gas_moments.m0();
    fvec(0) = (m_liq_centroid_star[0] - liq_centroid_h[0]) / m_length_scale;
    fvec(1) = (m_liq_centroid_star[1] - liq_centroid_h[1]) / m_length_scale;
    fvec(2) = (m_gas_centroid_star[0] - gas_centroid_h[0]) / m_length_scale;
    fvec(3) = (m_gas_centroid_star[1] - gas_centroid_h[1]) / m_length_scale;
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    this->errorvec(x, fvec);
    return 0;
  }

  int inputs() const { return 1; }
  int values() const { return 4; }
};

void MOF1::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                             const Data<IRL2D::Moments>& a_gas_moments,
                             const double a_dt, const Data<double>& a_U,
                             const Data<double>& a_V,
                             Data<IRL2D::Parabola>* a_interface) {
  // initial guess
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);

  const BasicMesh& mesh = a_U.getMesh();

#ifdef USE_MPI
  const double cell_volume = mesh.cell_volume();
  int nmixed_global = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / cell_volume;
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        nmixed_global++;
      }
    }
  }

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  IRL2D::Parabola dummy_par;
  IRL::ByteBuffer dummy_buffer;
  dummy_buffer.resize(0);
  dummy_buffer.resetBufferPointer();
  IRL::serializeAndPack(dummy_par, &dummy_buffer);
  const int size_parabola = dummy_buffer.size();

  int nmixed_local = std::max(nmixed_global / size, 1);
  std::vector<int> proc_offset(size + 1);
  proc_offset[0] = 0;
  for (int r = 0; r < size; r++) {
    proc_offset[r + 1] = proc_offset[r] + nmixed_local;
  }
  proc_offset[size] = nmixed_global;
  for (int r = 1; r < size + 1; r++) {
    proc_offset[r] = std::min(proc_offset[r], nmixed_global);
  }
  nmixed_local = proc_offset[rank + 1] - proc_offset[rank];
  IRL::ByteBuffer interface_local, interface_global;
  interface_local.resize(nmixed_local * sizeof(IRL2D::Parabola));
  interface_global.resize(0);
  interface_local.resetBufferPointer();
  interface_global.resetBufferPointer();

  int count = 0;
#endif

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
#ifdef USE_MPI
        if (count >= proc_offset[rank] && proc_offset[rank + 1]) {
#endif

          IRL2D::Vec x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
          IRL2D::Vec x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
          IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(x0, x1);

          // matching volume fraction and centroid
          IRL2D::Vec liq_centroid_star =
              a_liquid_moments(i, j).m1() / a_liquid_moments(i, j).m0();
          IRL2D::Vec gas_centroid_star =
              a_gas_moments(i, j).m1() / a_gas_moments(i, j).m0();
          MOF1Functor myMOFFunctor(rectangle, liq_centroid_star,
                                   gas_centroid_star, liquid_volume_fraction);
          myMOFFunctor.setframe((*a_interface)(i, j));
          Eigen::NumericalDiff<MOF1Functor> numericalDiffMyFunctor(
              myMOFFunctor);
          Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF1Functor>, double>
              lm(numericalDiffMyFunctor);
          Eigen::VectorXd x(1);
          x.setZero();
          Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(x);
          do {
            status = lm.minimizeOneStep(x);
          } while (status == Eigen::LevenbergMarquardtSpace::Running);
          const auto plane = IRL2D::MatchToVolumeFractionBisection(
              rectangle, myMOFFunctor.getplane(x), liquid_volume_fraction, 500);
#ifdef USE_MPI
          IRL::serializeAndPack(plane, &interface_local);
        }
        count++;
#else
        (*a_interface)(i, j) = plane;
#endif
      }
    }
  }

#ifdef USE_MPI
  std::vector<int> proc_count(size);
  for (int r = 0; r < size; r++) {
    proc_count[r] = size_parabola * (proc_offset[r + 1] - proc_offset[r]);
    proc_offset[r] = size_parabola * proc_offset[r];
  }

  interface_global.resize(size_parabola * nmixed_global);
  MPI_Allgatherv(interface_local.data(), size_parabola * nmixed_local, MPI_BYTE,
                 interface_global.data(), proc_count.data(), proc_offset.data(),
                 MPI_BYTE, MPI_COMM_WORLD);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / cell_volume;
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        IRL2D::Parabola parabola;
        IRL::unpackAndStore(&parabola, &interface_global);
        (*a_interface)(i, j) = parabola;
      }
    }
  }

#endif

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// Functor for MOF2
struct MOF2Functor {
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // variables
  const int m_inputs, m_values;
  const IRL2D::BezierList& m_cell;
  IRL2D::Moments m_liq_moments;
  IRL2D::Moments m_gas_moments;
  double m_liq_f_star;
  IRL2D::Vec m_liq_centroid_star;
  IRL2D::Vec m_gas_centroid_star;
  IRL2D::Mat m_liq_M2_star;
  IRL2D::Mat m_gas_M2_star;
  IRL2D::Vec m_datum;
  IRL2D::ReferenceFrame m_frame;
  double m_coeff;
  double m_length_scale;

  // constructor
  MOF2Functor(int inputs, int values, const IRL2D::BezierList& cell,
              const IRL2D::Moments& liq_moments,
              const IRL2D::Moments& gas_moments)
      : m_inputs(inputs),
        m_values(values),
        m_cell(cell),
        m_liq_moments(liq_moments),
        m_gas_moments(gas_moments),
        m_liq_f_star(liq_moments.m0() / IRL2D::ComputeArea(cell)),
        m_liq_centroid_star(liq_moments.m1() / liq_moments.m0()),
        m_gas_centroid_star(gas_moments.m1() / gas_moments.m0()) {
    // m_datum = ( (1-m_liq_f_star) * m_liq_centroid_star + m_liq_f_star *
    // m_gas_centroid_star );
    m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cell));
    RecenterMoments(&m_liq_moments, m_liq_centroid_star);
    RecenterMoments(&m_gas_moments, m_gas_centroid_star);
    m_liq_M2_star = m_liq_moments.m2();
    m_gas_M2_star = m_gas_moments.m2();
  }

  void setframe(const IRL2D::Parabola& guess_parabola) {
    m_coeff = guess_parabola.coeff();
    m_frame = guess_parabola.frame();
    m_datum = guess_parabola.datum();

    // check for curvature
    const double maxkdx = 4.0;
    const double kdx = 2.0 * m_coeff * m_length_scale;
    if (std::abs(kdx) > maxkdx) {
      m_coeff = 0.0;  // plane
    }
  }

  // x(0): theta    x(1): alpha
  const IRL2D::Parabola getparabola(const Eigen::VectorXd& x) const {
    const auto rotation = IRL2D::ReferenceFrame(x(0));
    const auto new_frame =
        IRL2D::ReferenceFrame(rotation * m_frame[0], rotation * m_frame[1]);
    const double new_coeff = m_coeff + x(1) / m_length_scale;
    const auto parabola = IRL2D::Parabola(m_datum, new_frame, new_coeff);
    return IRL2D::MatchToVolumeFraction(m_cell, parabola, m_liq_f_star);
  }

  mutable int iteration = 0;  // counting iterations for convergence
  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto parabola = this->getparabola(x);
    IRL2D::Moments lm = IRL2D::ComputeMoments(m_cell, parabola);
    IRL2D::Moments gm = IRL2D::ComputeMoments(m_cell) - lm;
    IRL2D::Vec lc = lm.m1() / lm.m0();
    IRL2D::Vec gc = gm.m1() / gm.m0();
    RecenterMoments(&lm, m_liq_centroid_star);
    RecenterMoments(&gm, m_gas_centroid_star);
    IRL2D::Mat lm2 = lm.m2();
    IRL2D::Mat gm2 = gm.m2();
    fvec(0) = (m_liq_centroid_star[0] - lc[0]) / m_length_scale;
    fvec(1) = (m_liq_centroid_star[1] - lc[1]) / m_length_scale;
    fvec(2) = (m_gas_centroid_star[0] - gc[0]) / m_length_scale;
    fvec(3) = (m_gas_centroid_star[1] - gc[1]) / m_length_scale;
    fvec(4) = (m_liq_M2_star[0][0] - lm2[0][0]) /
              (lm.m0() * std::pow(m_length_scale, 2.0));
    fvec(5) = (m_liq_M2_star[1][0] - lm2[1][0]) /
              (lm.m0() * std::pow(m_length_scale, 2.0));
    fvec(6) = (m_liq_M2_star[1][1] - lm2[1][1]) /
              (lm.m0() * std::pow(m_length_scale, 2.0));
    fvec(7) = (m_gas_M2_star[0][0] - gm2[0][0]) /
              (gm.m0() * std::pow(m_length_scale, 2.0));
    fvec(8) = (m_gas_M2_star[1][0] - gm2[1][0]) /
              (gm.m0() * std::pow(m_length_scale, 2.0));
    fvec(9) = (m_gas_M2_star[1][1] - gm2[1][1]) /
              (gm.m0() * std::pow(m_length_scale, 2.0));

    // scaling based on Shashkov's paper
    //  fvec(0) = std::pow((liq_centroid_h[0] - m_liq_centroid_star[0]), 2.0) /
    //  m_liq_moments.m0(); fvec(1) = std::pow((liq_centroid_h[1] -
    //  m_liq_centroid_star[1]), 2.0) / m_liq_moments.m0(); fvec(2) =
    //  std::pow((gas_centroid_h[0] - m_gas_centroid_star[0]), 2.0) /
    //  m_gas_moments.m0(); fvec(3) = std::pow((gas_centroid_h[1] -
    //  m_gas_centroid_star[1]), 2.0) / m_gas_moments.m0(); fvec(4) =
    //  std::pow((liq_M2_h[0][0] - m_liq_M2_star[0][0]), 2.0) /
    //  (std::pow(m_liq_M2_star[0][0],2.0)) ; fvec(5) = std::pow((liq_M2_h[1][1]
    //  - m_liq_M2_star[1][1]), 2.0) / (std::pow(m_liq_M2_star[1][1],2.0));
    //  fvec(6) = std::pow((liq_M2_h[1][0] - m_liq_M2_star[1][0]), 2.0) /
    //  (std::pow(m_liq_M2_star[0][0],2.0) + std::pow(m_liq_M2_star[1][1],2.0));
    //  fvec(7) = std::pow((gas_M2_h[0][0] - m_gas_M2_star[0][0]), 2.0) /
    //  (std::pow(m_gas_M2_star[0][0],2.0)); fvec(8) = std::pow((gas_M2_h[1][1]
    //  - m_gas_M2_star[1][1]), 2.0) / (std::pow(m_gas_M2_star[1][1],2.0));
    //  fvec(9) = std::pow((gas_M2_h[1][0] - m_gas_M2_star[1][0]), 2.0) /
    //  (std::pow(m_gas_M2_star[0][0],2.0) + std::pow(m_gas_M2_star[1][1],2.0));

    // Penalty to prevent kappa * dx > 4
    const double mu = 50.0;
    const double kdx = 2.0 * m_coeff * m_length_scale + x(1);
    const double maxkdx = 4.0;
    fvec(10) = mu * std::max(0.0, std::abs(kdx) - maxkdx);

    // std::cout << "Iteration: " << iteration << ", ||fvec|| = " << fvec.norm()
    // <<std::endl;

    // if (iteration == 0){
    //   std::cout << fvec.transpose() << std::endl;
    //   std::cout << fvec.norm() << std::endl;
    // }

    // iteration++;
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    this->errorvec(x, fvec);
    return 0;
  }

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

void MOF2::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                             const Data<IRL2D::Moments>& a_gas_moments,
                             const double a_dt, const Data<double>& a_U,
                             const Data<double>& a_V,
                             Data<IRL2D::Parabola>* a_interface) {
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);

  const BasicMesh& mesh = a_U.getMesh();

  // band for mixed cells
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }

  int N = 5;
  double Hp = 3.0;
  double h = mesh.dx();
  double eta = 0.5;

  // estimation curvature using particles
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (band(i, j) == 1) {
        IRL2D::Parabola target_interface = (*a_interface)(i, j);
        IRL2D::BezierList target_cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        IRL2D::Vec n_target = target_interface.frame()[1];
        std::vector<IRL2D::Parabola> all_interfaces;
        std::vector<IRL2D::BezierList> all_cells;
        for (int jj = -2; jj <= 2; jj++) {
          for (int ii = -2; ii <= 2; ii++) {
            if (band(ii + i, jj + j) == 1) {
              IRL2D::Vec n_neighbor = (*a_interface)(ii + i, jj + j).frame()[1];
              if (n_neighbor * n_target <= 0) {
                continue;  // filtering based on normals
              }
              all_interfaces.push_back((*a_interface)(ii + i, jj + j));
              all_cells.push_back(IRL2D::RectangleFromBounds(
                  IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
                  IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))));
            }
          }
        }
        IRL2D::Parabola parabolaParticle = target_interface;
        double coefficient =
            0.5 * IRL2D::getCurvature(target_interface, target_cell,
                                      all_interfaces, all_cells, N, Hp, h, eta);
        // std::cout << "Coefficient = " << coefficient << std::endl;
        parabolaParticle.coeff() = coefficient;

        int num_inputs = 2;
        int num_eq = 11;

        MOF2Functor myMOF2Functor(num_inputs, num_eq, target_cell,
                                  (a_liquid_moments)(i, j),
                                  (a_gas_moments)(i, j));
        myMOF2Functor.setframe(parabolaParticle);
        Eigen::NumericalDiff<MOF2Functor> numericalDiffMyFunctor(myMOF2Functor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF2Functor>, double> LM(
            numericalDiffMyFunctor);
        Eigen::VectorXd x(2);
        x.setZero();
        // lm.parameters.ftol = 1e-14;
        // lm.parameters.xtol = 1e-14;
        //  lm.parameters.maxfev = 500;
        LM.minimize(x);
        (*a_interface)(i, j) = myMOF2Functor.getparabola(x);
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// Functor for MOF2
// struct MOF2AugmentedLagrangianFunctor{
//   typedef double Scalar;
//   typedef Eigen::VectorXd InputType;
//   typedef Eigen::VectorXd ValueType;
//   typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
//   enum{
//     InputsAtCompileTime = Eigen::Dynamic,
//     ValuesAtCompileTime = Eigen::Dynamic
//   };

//   // variables
//   const int m_inputs, m_values;
//   const IRL2D::BezierList& m_cell;
//   IRL2D::Moments m_liq_moments;
//   IRL2D::Moments m_gas_moments;
//   double m_liq_f_star;
//   IRL2D::Vec m_liq_centroid_star;
//   IRL2D::Vec m_gas_centroid_star;
//   IRL2D::Mat m_liq_M2_star;
//   IRL2D::Mat m_gas_M2_star;
//   IRL2D::Vec m_datum;
//   IRL2D::ReferenceFrame m_frame;
//   IRL2D::Vec m_cell_centroid;
//   double m_coeff;
//   double m_length_scale;
//   std::vector<double> m_lambda;
//   std::vector<double> m_mu;

//   // constructor
//   MOF2AugmentedLagrangianFunctor(const int inputs, const int values, const
//   IRL2D::BezierList& cell,
//                                  const IRL2D::Moments& liq_moments, const
//                                  IRL2D::Moments& gas_moments, const
//                                  std::vector<double>& lambda, const
//                                  std::vector<double>& mu)
//     : m_inputs(inputs),
//       m_values(values),
//       m_cell(cell),
//       m_liq_moments(liq_moments),
//       m_gas_moments(gas_moments),
//       m_liq_f_star(liq_moments.m0() / IRL2D::ComputeArea(cell)),
//       m_liq_centroid_star(liq_moments.m1() / liq_moments.m0()),
//       m_gas_centroid_star(gas_moments.m1() / gas_moments.m0()),
//       m_lambda(lambda),
//       m_mu(mu) {
//       m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cell));
//       RecenterMoments(&m_liq_moments, m_liq_centroid_star);
//       RecenterMoments(&m_gas_moments, m_gas_centroid_star);
//       m_liq_M2_star = m_liq_moments.m2();
//       m_gas_M2_star = m_gas_moments.m2();
//       m_cell_centroid = IRL2D::ComputeMoments(m_cell).m1() /
//       IRL2D::ComputeArea(m_cell);
//   }

//   void setframe(const IRL2D::Parabola& guess_parabola){
//     m_coeff = guess_parabola.coeff();
//     m_frame = guess_parabola.frame();
//     m_datum = guess_parabola.datum();

//     // check for curvature
//     const double maxkdx = 4.0;
//     const double kdx = 2.0 * m_coeff * m_length_scale;
//     if (std::abs(kdx) > maxkdx){
//       m_coeff = 0.0; // plane
//     }
//   }

//   // x(0): theta    x(1): alpha   x(2): datum[0]    x(3): datum[1]
//   const IRL2D::Parabola getparabola(const Eigen::VectorXd& x) const {
//     const auto rotation = IRL2D::ReferenceFrame(x(0));
//     const auto new_frame = IRL2D::ReferenceFrame(rotation * m_frame[0],
//     rotation * m_frame[1]); const double new_coeff = m_coeff + x(1) /
//     m_length_scale; const auto new_datum = IRL2D::Vec(m_datum[0] +
//     x(2)*m_length_scale , m_datum[1] + x(3)*m_length_scale); const auto
//     parabola = IRL2D::Parabola(new_datum, new_frame, new_coeff); return
//     parabola;
//   }

//   double getVolfrac (const Eigen::VectorXd& x) const{
//     const auto parabola = this->getparabola(x);
//     const auto moments = IRL2D::ComputeMoments(m_cell, parabola);
//     return (moments.m0() / IRL2D::ComputeArea(m_cell));
//   }

//   // std::vector<double> getIneqFunction(const Eigen::VectorXd& x) const{
//   //   std::vector<double> IneqFunctions;
//   //   const double dx = m_length_scale;
//   //   const double x0 = (m_cell[0].first)[0], x1 = (m_cell[1].first)[0],
//   //                y0 = (m_cell[0].first)[1], y1 = (m_cell[2].first)[1];
//   //   std::vector<double> bounds = { (x0 - (m_datum[0] + x(2))) ,
//   ((m_datum[0] + x(2)) - x1) ,
//   //                                  (y0 - (m_datum[1] + x(3))) ,
//   ((m_datum[1] + x(3)) - y1) };
//   //   for (int bound = 0; bound < bounds.size(); bound++){
//   //     IneqFunctions[bound] = std::max(0.0 , bounds[bound]);
//   //   }
//   //   return IneqFunctions;
//   // }

//   mutable int iteration = 0; // counting iterations for convergence
//   void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
//     const auto parabola = this->getparabola(x);
//     IRL2D::Moments liq_mom = IRL2D::ComputeMoments(m_cell, parabola);
//     IRL2D::Moments gas_mom = IRL2D::ComputeMoments(m_cell) - liq_mom;
//     IRL2D::Vec liq_centroid_h = liq_mom.m1() / liq_mom.m0();
//     IRL2D::Vec gas_centroid_h = gas_mom.m1() / gas_mom.m0();
//     RecenterMoments(&liq_mom, m_liq_centroid_star);
//     RecenterMoments(&gas_mom, m_gas_centroid_star);
//     IRL2D::Mat liq_M2_h = liq_mom.m2();
//     IRL2D::Mat gas_M2_h = gas_mom.m2();

//     // centroids
//     fvec(0) = (m_liq_centroid_star[0] - liq_centroid_h[0]) / m_length_scale;
//     fvec(1) = (m_liq_centroid_star[1] - liq_centroid_h[1]) / m_length_scale;
//     fvec(2) = (m_gas_centroid_star[0] - gas_centroid_h[0]) / m_length_scale;
//     fvec(3) = (m_gas_centroid_star[1] - gas_centroid_h[1]) / m_length_scale;

//     // second moments
//     fvec(4) = (m_liq_M2_star[0][0] - liq_M2_h[0][0]) /
//     (liq_mom.m0()*std::pow(m_length_scale, 2.0)); fvec(5) =
//     (m_liq_M2_star[1][0] - liq_M2_h[1][0]) /
//     (liq_mom.m0()*std::pow(m_length_scale, 2.0)); fvec(6) =
//     (m_liq_M2_star[1][1] - liq_M2_h[1][1]) /
//     (liq_mom.m0()*std::pow(m_length_scale, 2.0)); fvec(7) =
//     (m_gas_M2_star[0][0] - gas_M2_h[0][0]) /
//     (gas_mom.m0()*std::pow(m_length_scale, 2.0)); fvec(8) =
//     (m_gas_M2_star[1][0] - gas_M2_h[1][0]) /
//     (gas_mom.m0()*std::pow(m_length_scale, 2.0)); fvec(9) =
//     (m_gas_M2_star[1][1] - gas_M2_h[1][1]) /
//     (gas_mom.m0()*std::pow(m_length_scale, 2.0));

//     // volume fraction constraint
//     double liq_f_h = liq_mom.m0() / IRL2D::ComputeArea(m_cell);
//     double f_constraint = liq_f_h - m_liq_f_star;
//     fvec(10) = std::sqrt(m_mu[0]) * f_constraint + 0.5 / (std::sqrt(m_mu[0]))
//     * m_lambda[0];

//     // curvature
//     const double mu_kdx = 50.0;
//     const double kdx = 2.0 * m_coeff * m_length_scale + x(1);
//     const double maxkdx = 4.0;
//     fvec(11) = mu_kdx * std::max(0.0, std::abs(kdx) - maxkdx);

//     // std::cout << iteration++ << ". residual = [" << fvec(0) << " " <<
//     fvec(1) << " " << fvec(2) << " " << fvec(3) << " " << fvec(4) << " "
//     //           << fvec(5) << " " << fvec(6) << " " << fvec(7) << " " <<
//     fvec(8) << " " << fvec(9) << " " << fvec(10) << "] "
//     //           << "mu = " << m_mu[0] << " lambda = " << m_lambda[0] << "
//     norm = " << fvec.norm() << std::endl;
//     // confining the datum within the cell
//     // const double dx = m_length_scale;
//     // const double x0 = (m_cell[0].first)[0], x1 = (m_cell[1].first)[0],
//     //              y0 = (m_cell[0].first)[1], y1 = (m_cell[2].first)[1];
//     //std::vector<double> bounds = { -(dx/2.0 + x(2)) , (x(2) - dx/2.0) ,
//     -(dx/2.0 + x(3)) , (x(3) - dx/2.0) };
//     // std::vector<double> bounds = { (x0 - (m_datum[0] + x(2))) ,
//     ((m_datum[0] + x(2)) - x1) ,
//     //                                (y0 - (m_datum[1] + x(3))) ,
//     ((m_datum[1] + x(3)) - y1) };
//     // for (int bound = 0; bound < bounds.size(); bound++){
//     //   fvec(10 + (bound + 1)) = std::sqrt(m_mu[bound + 1]) * std::max(0.0 ,
//     bounds[bound]) +
//     //                            0.5 / (std::sqrt(m_mu[bound + 1])) *
//     m_lambda[bound + 1];
//     // }

//     // std::vector<double> IneqFunctions = this->getIneqFunction(x);
//     // for (int bound = 0; bound < IneqFunctions.size(); bound++){
//     //   fvec(10 + (bound + 1)) = std::sqrt(m_mu[bound + 1]) *
//     IneqFunctions[bound] +
//     //                            0.5 / (std::sqrt(m_mu[bound + 1])) *
//     m_lambda[bound + 1];
//     // }

//     // // curvature inequality
//     // const double kdx = 2.0 * m_coeff * m_length_scale + x(1);
//     // const double maxkdx = 4.0;
//     // fvec(15) = std::sqrt(m_mu[5]) * std::max(0.0, std::abs(kdx) - maxkdx);
//     +
//     //            0.5 / (std::sqrt(m_mu[5])) * m_lambda[5];

//   }

//   int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
//     this->errorvec(x, fvec);
//     return 0;
//   }

//   int inputs() const { return m_inputs; }
//   int values() const { return m_values; }

// };

struct MOF2AugmentedLagrangianFunctor {
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  const int m_inputs, m_values;
  const IRL2D::BezierList& m_cell;
  IRL2D::Moments m_liq_moments;
  IRL2D::Moments m_gas_moments;
  double m_liq_f_star;
  IRL2D::Vec m_liq_centroid_star;
  IRL2D::Vec m_gas_centroid_star;
  IRL2D::Mat m_liq_M2_star;
  IRL2D::Mat m_gas_M2_star;
  IRL2D::Vec m_datum;
  IRL2D::ReferenceFrame m_frame;
  IRL2D::Vec m_cell_centroid;
  double m_coeff;
  double m_length_scale;
  std::vector<double> m_lambda;
  std::vector<double> m_mu;
  IRL2D::Vec m_cell_min, m_cell_max;

  MOF2AugmentedLagrangianFunctor(const int inputs, const int values,
                                 const IRL2D::BezierList& cell,
                                 const IRL2D::Moments& liq_moments,
                                 const IRL2D::Moments& gas_moments,
                                 const std::vector<double>& lambda,
                                 const std::vector<double>& mu)
      : m_inputs(inputs),
        m_values(values),
        m_cell(cell),
        m_liq_moments(liq_moments),
        m_gas_moments(gas_moments),
        m_liq_f_star(liq_moments.m0() / IRL2D::ComputeArea(cell)),
        m_liq_centroid_star(liq_moments.m1() / liq_moments.m0()),
        m_gas_centroid_star(gas_moments.m1() / gas_moments.m0()),
        m_lambda(lambda),
        m_mu(mu) {
    m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cell));
    RecenterMoments(&m_liq_moments, m_liq_centroid_star);
    RecenterMoments(&m_gas_moments, m_gas_centroid_star);
    m_liq_M2_star = m_liq_moments.m2();
    m_gas_M2_star = m_gas_moments.m2();
    m_cell_centroid =
        IRL2D::ComputeMoments(m_cell).m1() / IRL2D::ComputeArea(m_cell);

    m_cell_min = IRL2D::Vec(m_cell[0].first.x(), m_cell[0].first.y());
    m_cell_max = IRL2D::Vec(m_cell[2].first.x(), m_cell[2].first.y());
  }

  void setframe(const IRL2D::Parabola& guess_parabola) {
    m_coeff = guess_parabola.coeff();
    m_frame = guess_parabola.frame();
    m_datum = guess_parabola.datum();
    const double maxkdx = 4.0;
    const double kdx = 2.0 * m_coeff * m_length_scale;
    if (std::abs(kdx) > maxkdx) {
      m_coeff = 0.0;
    }
  }

  const IRL2D::Parabola getparabola(const Eigen::VectorXd& x) const {
    auto rotation = IRL2D::ReferenceFrame(x(0));
    auto new_frame =
        IRL2D::ReferenceFrame(rotation * m_frame[0], rotation * m_frame[1]);
    double new_coeff = m_coeff + x(1) / m_length_scale;
    auto new_datum = IRL2D::Vec(m_datum[0] + x(2) * m_length_scale,
                                m_datum[1] + x(3) * m_length_scale);
    return IRL2D::Parabola(new_datum, new_frame, new_coeff);
  }

  double getVolfrac(const Eigen::VectorXd& x) const {
    return IRL2D::ComputeMoments(m_cell, this->getparabola(x)).m0() /
           IRL2D::ComputeArea(m_cell);
  }

  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    auto parabola = this->getparabola(x);
    auto liq_mom = IRL2D::ComputeMoments(m_cell, parabola);
    auto gas_mom = IRL2D::ComputeMoments(m_cell) - liq_mom;
    auto liq_centroid_h = liq_mom.m1() / liq_mom.m0();
    auto gas_centroid_h = gas_mom.m1() / gas_mom.m0();
    RecenterMoments(&liq_mom, m_liq_centroid_star);
    RecenterMoments(&gas_mom, m_gas_centroid_star);
    auto liq_M2_h = liq_mom.m2();
    auto gas_M2_h = gas_mom.m2();

    fvec(0) = (m_liq_centroid_star[0] - liq_centroid_h[0]) / m_length_scale;
    fvec(1) = (m_liq_centroid_star[1] - liq_centroid_h[1]) / m_length_scale;
    fvec(2) = (m_gas_centroid_star[0] - gas_centroid_h[0]) / m_length_scale;
    fvec(3) = (m_gas_centroid_star[1] - gas_centroid_h[1]) / m_length_scale;

    fvec(4) = (m_liq_M2_star[0][0] - liq_M2_h[0][0]) /
              (liq_mom.m0() * std::pow(m_length_scale, 2.0));
    fvec(5) = (m_liq_M2_star[1][0] - liq_M2_h[1][0]) /
              (liq_mom.m0() * std::pow(m_length_scale, 2.0));
    fvec(6) = (m_liq_M2_star[1][1] - liq_M2_h[1][1]) /
              (liq_mom.m0() * std::pow(m_length_scale, 2.0));
    fvec(7) = (m_gas_M2_star[0][0] - gas_M2_h[0][0]) /
              (gas_mom.m0() * std::pow(m_length_scale, 2.0));
    fvec(8) = (m_gas_M2_star[1][0] - gas_M2_h[1][0]) /
              (gas_mom.m0() * std::pow(m_length_scale, 2.0));
    fvec(9) = (m_gas_M2_star[1][1] - gas_M2_h[1][1]) /
              (gas_mom.m0() * std::pow(m_length_scale, 2.0));

    double liq_f_h = liq_mom.m0() / IRL2D::ComputeArea(m_cell);
    double f_constraint = liq_f_h - m_liq_f_star;
    fvec(10) = std::sqrt(m_mu[0]) * f_constraint +
               0.5 / std::sqrt(m_mu[0]) * m_lambda[0];

    const double mu_kdx = 50.0;
    const double kdx = 2.0 * m_coeff * m_length_scale + x(1);
    const double maxkdx = 4.0;
    fvec(11) = mu_kdx * std::max(0.0, std::abs(kdx) - maxkdx);

    auto datum = this->getparabola(x).datum();
    fvec(12) = std::sqrt(m_mu[1]) * std::max(0.0, m_cell_min.x() - datum.x()) +
               0.5 / std::sqrt(m_mu[1]) * m_lambda[1];
    fvec(13) = std::sqrt(m_mu[2]) * std::max(0.0, datum.x() - m_cell_max.x()) +
               0.5 / std::sqrt(m_mu[2]) * m_lambda[2];
    fvec(14) = std::sqrt(m_mu[3]) * std::max(0.0, m_cell_min.y() - datum.y()) +
               0.5 / std::sqrt(m_mu[3]) * m_lambda[3];
    fvec(15) = std::sqrt(m_mu[4]) * std::max(0.0, datum.y() - m_cell_max.y()) +
               0.5 / std::sqrt(m_mu[4]) * m_lambda[4];
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    this->errorvec(x, fvec);
    return 0;
  }

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

// void MOF2AL::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
//                             const Data<IRL2D::Moments>& a_gas_moments,
//                             const double a_dt, const Data<double>& a_U,
//                             const Data<double>& a_V,
//                             Data<IRL2D::Parabola>* a_interface){

//   // PLIC for curvature estimation
//   ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
//                             a_interface);

//   const BasicMesh& mesh = a_U.getMesh();

//   // finding mixed cells
//   Data<int> band(&mesh);
//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       band(i, j) = 0;
//       const double liquid_volume_fraction =
//           (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
//       if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
//           liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
//         band(i, j) = 1;
//       }
//     }
//   }

//   int count_noconv = 0; // cells that do not converge

//   // parameters for curvature estimation using particles
//   int N = 5; double Hp = 3.0; double h = mesh.dx(); double eta = 0.5;

//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       const double liquid_volume_fraction = a_liquid_moments(i, j).m0() /
//       mesh.cell_volume(); if (liquid_volume_fraction >=
//       IRL::global_constants::VF_LOW &&
//           liquid_volume_fraction <= IRL::global_constants::VF_HIGH){ // && i
//           == 23 && j == 51) {

//         // current cell
//         IRL2D::Vec x0 = IRL2D::Vec(mesh.x(i) , mesh.y(j));
//         IRL2D::Vec x1 = IRL2D::Vec(mesh.x(i+1) , mesh.y(j+1));
//         IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(x0, x1);

//         // Using Jibben for parabola fit
//         IRL2D::Parabola target_interface = (*a_interface)(i,j);
//         IRL2D::Vec n_target = target_interface.frame()[1];
//         std::vector<IRL2D::Parabola> all_interfaces;
//         std::vector<IRL2D::BezierList> all_cells;
//         for (int jj = -2; jj <= 2; jj++){
//           for (int ii = -2; ii <= 2; ii++){
//             if (band(ii + i, jj + j) == 1){
//               IRL2D::Vec n_neighbor = (*a_interface)(ii + i, jj +
//               j).frame()[1]; if (n_neighbor * n_target <= 0){
//                 continue; // filtering based on normals
//               }
//               all_interfaces.push_back((*a_interface)(ii + i, jj + j));
//               all_cells.push_back(IRL2D::RectangleFromBounds(
//                 IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
//                 IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))
//               ));
//             }
//           }
//         }
//         IRL2D::Parabola parabola_guess = target_interface;
//         parabola_guess.coeff() = 0.5 * IRL2D::getCurvature(target_interface,
//         rectangle,
//                                        all_interfaces, all_cells, N, Hp, h,
//                                        eta);

//         // moments
//         const auto liq_moments = a_liquid_moments(i,j);
//         const auto gas_moments = a_gas_moments(i,j);

//         // parameters for augmented lagrangian
//         std::vector<double> lambda(1, 0.0); // CHANGE 1 -> 6 with
//         inequalities std::vector<double> mu(1, 1.0); double mu_max = 10000.0;
//         // 100 before double tol = 1e-8; // -12 before double
//         f_constraint_current = 1; double f_constraint_prev; int iter = 0; int
//         max_iter = 1000; // 1000 before

//         // initial guess for first iteration
//         // auto guess_interface = (*a_interface)(i,j);
//         auto guess_interface = parabola_guess;
//         IRL2D::Parabola parabola;

//         // for LM solver
//         int num_inputs = 4;
//         int num_eq = 12;

//         Eigen::VectorXd x(num_inputs);
//         x.setZero();

//         // residual vector
//         Eigen::VectorXd residuals(num_eq);
//         residuals.setZero();

//         while (std::abs(f_constraint_current) > tol){

//           iter++;
//           // std::cout << "-----------------Augmented Lagrangian Iteration: "
//           << iter << " (" << i << "," << j<< ")-------------------" <<
//           std::endl;
//           // minimization variables
//           // Eigen::VectorXd x(num_inputs);
//           // x.setZero();

//           // setting up LM solver
//           MOF2AugmentedLagrangianFunctor myMOFFunctor(num_inputs, num_eq,
//           rectangle,
//                                                       liq_moments,
//                                                       gas_moments, lambda,
//                                                       mu);
//           myMOFFunctor.setframe(guess_interface);
//           Eigen::NumericalDiff<MOF2AugmentedLagrangianFunctor>
//           numericalDiffMyFunctor(myMOFFunctor);
//           Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF2AugmentedLagrangianFunctor>,
//           double> lm(numericalDiffMyFunctor);

//           //LM solver parameters
//           // lm.parameters.ftol = 1e-12;
//           // lm.parameters.xtol = 1e-12;

//           auto x_prev = x; // storing prev iter for AL multiplier updates
//           lm.minimize(x); // minimization

//           // compute volume fraction with new and prev x
//           f_constraint_current = myMOFFunctor.getVolfrac(x) -
//           liquid_volume_fraction; f_constraint_prev =
//           myMOFFunctor.getVolfrac(x_prev) - liquid_volume_fraction;

//           // update Lagrangian multiplier for constraint
//           lambda[0] += 2.0 * mu[0] * f_constraint_current;

//           // update penalty parameter for constraint
//           if (std::abs(f_constraint_current) < 0.25 *
//           std::abs(f_constraint_prev)){
//             //mu[0] = mu[0];
//             mu[0] = std::min(mu[0], mu_max);
//           } else {
//             //mu[0] *= 2.0;
//             mu[0] = std::min(2*mu[0], mu_max);
//           }

//           // updating parameters for inequalities
//           // std::vector<double> IneqFunctions =
//           myMOFFunctor.getIneqFunction(x);
//           // for (int ineq = 0 ; ineq < IneqFunctions.size(); ineq++){
//           //   // Lagrangian multiplier update
//           //   lambda[ineq + 1] += 2.0 * mu[ineq + 1] * IneqFunctions[ineq];

//           //   // penalty term update
//           // }

//           // updating interface guess for next AL iteration
//           parabola = myMOFFunctor.getparabola(x); // new parabola
//           //guess_interface = parabola;

//           // break if it reaches max iter
//           if (iter == max_iter){
//             count_noconv++;
//             std::cout << "------------------------------------------------"
//             << std::endl;
//             // cell location
//             std::cout << "(i,j) = (" << i << " , " << j << ")" << std::endl;

//             // calculating residuals
//             myMOFFunctor.errorvec(x, residuals);
//             std::cout << "residual = [" << residuals(0) << " " <<
//             residuals(1) << " " << residuals(2) << " " << residuals(3) << " "
//             << residuals(4) << " "
//                       << residuals(5) << " " << residuals(6) << " " <<
//                       residuals(7) << " " << residuals(8) << " " <<
//                       residuals(9) << " " << residuals(10) << "] "
//                       << " norm = " << residuals.norm() << std::endl;
//             std::cout << "Volfrac residual = " <<
//             std::abs(f_constraint_current) << std::endl;

//             break;
//           }

//         }
//         //std::cout << "AL iterations: " << iter - 1000 << std::endl;
//         //(*a_interface)(i,j) = parabola; // final interface

//         (*a_interface)(i, j) = parabola;
//       }
//     }
//   }
//   std::cout << "Cells not converged = " << count_noconv << std::endl;

//   a_interface->updateBorder();
//   correctInterfaceBorders(a_interface);

// }

void MOF2AL::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  // PLIC for curvature estimation
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);

  const BasicMesh& mesh = a_U.getMesh();

  // finding mixed cells
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }

  int count_noconv = 0;  // cells that do not converge

  // parameters for curvature estimation using particles
  int N = 5;
  double Hp = 3.0;
  double h = mesh.dx();
  double eta = 0.5;

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <=
              IRL::global_constants::VF_HIGH) {  // && i == 23 && j == 51) {

        // current cell
        IRL2D::Vec x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
        IRL2D::Vec x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
        IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(x0, x1);

        // Using Jibben for parabola fit
        IRL2D::Parabola target_interface = (*a_interface)(i, j);
        IRL2D::Vec n_target = target_interface.frame()[1];
        std::vector<IRL2D::Parabola> all_interfaces;
        std::vector<IRL2D::BezierList> all_cells;
        for (int jj = -2; jj <= 2; jj++) {
          for (int ii = -2; ii <= 2; ii++) {
            if (band(ii + i, jj + j) == 1) {
              IRL2D::Vec n_neighbor = (*a_interface)(ii + i, jj + j).frame()[1];
              if (n_neighbor * n_target <= -0.5) {
                continue;  // filtering based on normals
              }
              all_interfaces.push_back((*a_interface)(ii + i, jj + j));
              all_cells.push_back(IRL2D::RectangleFromBounds(
                  IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
                  IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))));
            }
          }
        }
        IRL2D::Parabola parabola_guess = target_interface;
        parabola_guess.coeff() =
            0.5 * IRL2D::getCurvature(target_interface, rectangle,
                                      all_interfaces, all_cells, N, Hp, h, eta);

        // moments
        const auto liq_moments = a_liquid_moments(i, j);
        const auto gas_moments = a_gas_moments(i, j);

        // parameters for augmented lagrangian
        std::vector<double> lambda(5, 0.0);  // CHANGE 1 -> 6 with inequalities
        std::vector<double> mu(5, 1.0);
        double mu_max = 10000.0;  // 100 before
        double tol = 1e-8;        // -12 before
        double f_constraint_current = 1;
        double f_constraint_prev;
        int iter = 0;
        int max_iter = 1000;  // 1000 before

        // initial guess for first iteration
        // auto guess_interface = (*a_interface)(i,j);
        auto guess_interface = parabola_guess;
        IRL2D::Parabola parabola;

        // for LM solver
        int num_inputs = 4;
        int num_eq = 16;

        Eigen::VectorXd x(num_inputs);
        x.setZero();

        // residual vector
        Eigen::VectorXd residuals(num_eq);
        residuals.setZero();

        while (std::abs(f_constraint_current) > tol) {
          iter++;
          // std::cout << "-----------------Augmented Lagrangian Iteration: " <<
          // iter << " (" << i << "," << j<< ")-------------------" <<
          // std::endl; minimization variables Eigen::VectorXd x(num_inputs);
          // x.setZero();

          // setting up LM solver
          MOF2AugmentedLagrangianFunctor myMOFFunctor(num_inputs, num_eq,
                                                      rectangle, liq_moments,
                                                      gas_moments, lambda, mu);
          myMOFFunctor.setframe(guess_interface);
          Eigen::NumericalDiff<MOF2AugmentedLagrangianFunctor>
              numericalDiffMyFunctor(myMOFFunctor);
          Eigen::LevenbergMarquardt<
              Eigen::NumericalDiff<MOF2AugmentedLagrangianFunctor>, double>
              lm(numericalDiffMyFunctor);

          // LM solver parameters
          //  lm.parameters.ftol = 1e-12;
          //  lm.parameters.xtol = 1e-12;

          auto x_prev = x;  // storing prev iter for AL multiplier updates
          lm.minimize(x);   // minimization

          // compute volume fraction with new and prev x
          f_constraint_current =
              myMOFFunctor.getVolfrac(x) - liquid_volume_fraction;
          f_constraint_prev =
              myMOFFunctor.getVolfrac(x_prev) - liquid_volume_fraction;

          // update Lagrangian multiplier for constraint
          lambda[0] += 2.0 * mu[0] * f_constraint_current;

          // update penalty parameter for constraint
          if (std::abs(f_constraint_current) <
              0.25 * std::abs(f_constraint_prev)) {
            mu[0] = std::min(mu[0], mu_max);
          } else {
            mu[0] = std::min(2 * mu[0], mu_max);
          }

          IRL2D::Vec datum = myMOFFunctor.getparabola(x).datum();
          std::vector<double> datum_ineq = {
              myMOFFunctor.m_cell_min.x() - datum.x(),
              datum.x() - myMOFFunctor.m_cell_max.x(),
              myMOFFunctor.m_cell_min.y() - datum.y(),
              datum.y() - myMOFFunctor.m_cell_max.y()};
          for (int k = 0; k < 4; ++k) {
            if (datum_ineq[k] > 0) {
              lambda[k + 1] += 2.0 * mu[k + 1] * datum_ineq[k];
              mu[k + 1] = std::min(2 * mu[k + 1], mu_max);
            }
          }

          parabola = myMOFFunctor.getparabola(x);  // new parabola

          // break if it reaches max iter
          if (iter == max_iter) {
            count_noconv++;
            std::cout << "------------------------------------------------"
                      << std::endl;
            // cell location
            std::cout << "(i,j) = (" << i << " , " << j << ")" << std::endl;

            // calculating residuals
            myMOFFunctor.errorvec(x, residuals);
            std::cout << "residual = [" << residuals(0) << " " << residuals(1)
                      << " " << residuals(2) << " " << residuals(3) << " "
                      << residuals(4) << " " << residuals(5) << " "
                      << residuals(6) << " " << residuals(7) << " "
                      << residuals(8) << " " << residuals(9) << " "
                      << residuals(10) << "] "
                      << " norm = " << residuals.norm() << std::endl;
            std::cout << "Volfrac residual = " << std::abs(f_constraint_current)
                      << std::endl;

            break;
          }
        }
        // std::cout << "AL iterations: " << iter - 1000 << std::endl;
        //(*a_interface)(i,j) = parabola; // final interface

        (*a_interface)(i, j) = parabola;
      }
    }
  }
  std::cout << "Cells not converged = " << count_noconv << std::endl;

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// using unit cell
// -----------------------------------------------------------------------------------------

void MOF2ALUnit::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                                   const Data<IRL2D::Moments>& a_gas_moments,
                                   const double a_dt, const Data<double>& a_U,
                                   const Data<double>& a_V,
                                   Data<IRL2D::Parabola>* a_interface) {
  // initial guess
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);

  const BasicMesh& mesh = a_U.getMesh();

  // #ifdef USE_MPI
  //   const double cell_volume = mesh.cell_volume();
  //   int nmixed_global = 0;
  //   for (int i = mesh.imin(); i <= mesh.imax(); ++i){
  //     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j){
  //       const double liquid_volume_fraction = a_liquid_moments(i,j).m0() /
  //       cell_volume; if (liquid_volume_fraction >=
  //       IRL::global_constants::VF_LOW &&
  //           liquid_volume_fraction <= IRL::global_constants::VF_HIGH){
  //         nmixed_global++;
  //       }
  //     }
  //   }

  //   int rank, size;
  //   MPI_Comm_size(MPI_COMM_WORLD, &size);
  //   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //   IRL2D::Parabola dummy_par;
  //   IRL::ByteBuffer dummy_buffer;
  //   dummy_buffer.resize(0);
  //   dummy_buffer.resetBufferPointer();
  //   IRL::serializeAndPack(dummy_par , &dummy_buffer);
  //   const int size_parabola = dummy_buffer.size();

  //   int nmixed_local = std::max(nmixed_global / size , 1);
  //   std::vector<int> proc_offset(size + 1);
  //   proc_offset[0] = 0;
  //   for (int r = 0; r < size; r++){
  //     proc_offset[r + 1] = proc_offset[r] + nmixed_local;
  //   }
  //   proc_offset[size] = nmixed_global;
  //   for (int r = 1; r < size + 1; r++){
  //     proc_offset[r] = std::min(proc_offset[r], nmixed_global);
  //   }
  //   nmixed_local = proc_offset[rank + 1] - proc_offset[rank];
  //   IRL::ByteBuffer interface_local, interface_global;
  //   interface_local.resize(nmixed_local * sizeof(IRL2D::Parabola));
  //   interface_global.resize(0);
  //   interface_local.resetBufferPointer();
  //   interface_global.resetBufferPointer();

  //   int count = 0;
  // #endif

  // counting how many cells did not converge
  int count_noconv = 0;

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <=
              IRL::global_constants::VF_HIGH) {  // && i == 23 && j == 51) {

        // #ifdef USE_MPI
        //       if (count >= proc_offset[rank] && proc_offset[rank + 1]) {
        // #endif

        // current cell
        IRL2D::Vec x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
        IRL2D::Vec x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
        IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(x0, x1);

        // moments
        const auto liq_moments = a_liquid_moments(i, j);
        const auto gas_moments = a_gas_moments(i, j);

        // parameters for augmented lagrangian
        std::vector<double> lambda(1, 0.0);  // CHANGE 1 -> 6 with inequalities
        std::vector<double> mu(1, 1.0);
        double mu_max = 10000.0;  // 100 before
        double tol = 1e-12;       // -12 before
        double f_constraint_current = 1;
        double f_constraint_prev;
        int iter = 0;
        int max_iter = 10000;  // 1000 before

        // initial guess for first iteration
        auto guess_interface = (*a_interface)(i, j);
        IRL2D::Parabola parabola;

        // transform interface to unit cell

        // for LM solver
        int num_inputs = 4;
        int num_eq = 12;

        Eigen::VectorXd x(num_inputs);
        x.setZero();

        // residual vector
        Eigen::VectorXd residuals(num_eq);
        residuals.setZero();

        // transfored quantities
        IRL2D::BezierList UnitRectangle =
            IRL2D::ComputeTransformedCell(rectangle, true);
        const auto unit_liq_moments = IRL2D::ComputeTransformedCellMoments(
            rectangle, (*a_interface)(i, j), true);
        const auto unit_gas_moments =
            IRL2D::ComputeMoments(UnitRectangle) - unit_liq_moments;
        IRL2D::Parabola unit_guess_interface =
            IRL2D::ComputeTransformedParabola(rectangle, guess_interface, true);

        while (std::abs(f_constraint_current) > tol) {
          iter++;

          // std::cout << "-----------------Augmented Lagrangian Iteration: " <<
          // iter << "-------------------" << std::endl;

          // minimization variables
          // Eigen::VectorXd x(num_inputs);
          // x.setZero();

          // setting up LM solver
          MOF2AugmentedLagrangianFunctor myMOFFunctor(
              num_inputs, num_eq, UnitRectangle, unit_liq_moments,
              unit_gas_moments, lambda, mu);
          myMOFFunctor.setframe(unit_guess_interface);  // guess for interface
          Eigen::NumericalDiff<MOF2AugmentedLagrangianFunctor>
              numericalDiffMyFunctor(myMOFFunctor);
          Eigen::LevenbergMarquardt<
              Eigen::NumericalDiff<MOF2AugmentedLagrangianFunctor>, double>
              lm(numericalDiffMyFunctor);

          // LM solver parameters
          lm.parameters.ftol = 1e-12;
          lm.parameters.xtol = 1e-12;

          auto x_prev = x;  // storing prev iter for AL multiplier updates
          lm.minimize(x);   // minimization

          // compute volume fraction with new and prev x
          f_constraint_current =
              myMOFFunctor.getVolfrac(x) - liquid_volume_fraction;
          f_constraint_prev =
              myMOFFunctor.getVolfrac(x_prev) - liquid_volume_fraction;

          // update Lagrangian multiplier for constraint
          lambda[0] += 2.0 * mu[0] * f_constraint_current;

          // update penalty parameter for constraint
          if (std::abs(f_constraint_current) <
              0.25 * std::abs(f_constraint_prev)) {
            // mu[0] = mu[0];
            mu[0] = std::min(mu[0], mu_max);
          } else {
            // mu[0] *= 2.0;
            mu[0] = std::min(2 * mu[0], mu_max);
          }

          // updating parameters for inequalities
          // std::vector<double> IneqFunctions =
          // myMOFFunctor.getIneqFunction(x); for (int ineq = 0 ; ineq <
          // IneqFunctions.size(); ineq++){
          //   // Lagrangian multiplier update
          //   lambda[ineq + 1] += 2.0 * mu[ineq + 1] * IneqFunctions[ineq];

          //   // penalty term update
          // }

          // updating interface guess for next AL iteration
          parabola = myMOFFunctor.getparabola(x);  // new parabola
          // guess_interface = parabola;

          // break if it reaches max iter
          if (iter == max_iter) {
            count_noconv++;
            std::cout << "------------------------------------------------"
                      << std::endl;
            // cell location
            std::cout << "(i,j) = (" << i << " , " << j << ")" << std::endl;

            // calculating residuals
            myMOFFunctor.errorvec(x, residuals);
            std::cout << "residual = [" << residuals(0) << " " << residuals(1)
                      << " " << residuals(2) << " " << residuals(3) << " "
                      << residuals(4) << " " << residuals(5) << " "
                      << residuals(6) << " " << residuals(7) << " "
                      << residuals(8) << " " << residuals(9) << " "
                      << residuals(10) << "] "
                      << " norm = " << residuals.norm() << std::endl;
            std::cout << "Volfrac residual = " << std::abs(f_constraint_current)
                      << std::endl;

            break;
          }
        }
        // std::cout << "AL iterations: " << iter - 1000 << std::endl;
        //(*a_interface)(i,j) = parabola; // final interface

        // #ifdef USE_MPI
        //       IRL::serializeAndPack(parabola, &interface_local);
        //       }
        //       count++;
        // #else
        // transform parabola back to original cell
        (*a_interface)(i, j) =
            IRL2D::ComputeTransformedParabola(UnitRectangle, parabola, false);

        //(*a_interface)(i, j) = parabola;
        // #endif
      }
    }
  }
  std::cout << "Cells not converged = " << count_noconv << std::endl;

  // #ifdef USE_MPI
  //   std::vector<int> proc_count(size);
  //   for (int r = 0; r < size; r++){
  //     proc_count[r] = size_parabola * (proc_offset[r + 1] - proc_offset[r]);
  //     proc_offset[r] = size_parabola * proc_offset[r];
  //   }

  //   interface_global.resize(size_parabola * nmixed_global);
  //   MPI_Allgatherv(interface_local.data(), size_parabola * nmixed_local,
  //   MPI_BYTE,
  //                  interface_global.data(), proc_count.data(),
  //                  proc_offset.data(), MPI_BYTE, MPI_COMM_WORLD);

  //   for (int i = mesh.imin(); i <= mesh.imax(); ++i){
  //     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j){
  //       const double liquid_volume_fraction = a_liquid_moments(i,j).m0() /
  //       cell_volume; if (liquid_volume_fraction >=
  //       IRL::global_constants::VF_LOW &&
  //           liquid_volume_fraction <= IRL::global_constants::VF_HIGH){
  //             IRL2D::Parabola parabola;
  //             IRL::unpackAndStore(&parabola, &interface_global);
  //             (*a_interface)(i,j) = parabola;
  //       }
  //     }
  //   }

  // #endif

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// optimization using nlopt
// ---------------------------------------------------------------
struct OptData {
  const IRL2D::BezierList m_cell;
  IRL2D::Parabola parabola_guess;
  IRL2D::Moments lms, gms;  // moments
  IRL2D::Vec lcs, gcs;      // target centroids
  IRL2D::Mat lm2s, gm2s;    // target second moments
  double lvfs;              // target liq vol frac
  IRL2D::Vec x0, x1;        // cell bounds
  double length_scale;
  IRL2D::Parabola parabola_opt;  // parabola obtained from optimization

  // parabola params
  double m_coeff;
  IRL2D::Vec m_datum;
  IRL2D::Mat m_frame;

  // constructor
  OptData(const IRL2D::BezierList& cell, const IRL2D::Parabola& interface_guess,
          const IRL2D::Moments& liq_moments, const IRL2D::Moments& gas_moments)
      : m_cell(cell),
        parabola_guess(interface_guess),
        lms(liq_moments),
        gms(gas_moments) {
    lvfs = liq_moments.m0() / IRL2D::ComputeArea(cell);
    lcs = liq_moments.m1() / liq_moments.m0();
    gcs = gas_moments.m1() / gas_moments.m0();
    RecenterMoments(&lms, lcs);
    RecenterMoments(&gms, gcs);
    lm2s = lms.m2();
    gm2s = gms.m2();
    x0 = cell[0].first;
    x1 = cell[2].first;
    length_scale = std::sqrt(IRL2D::ComputeArea(cell));
  }

  void setframe() {
    m_coeff = parabola_guess.coeff();
    m_frame = parabola_guess.frame();
    m_datum = parabola_guess.datum();
    const double maxkdx = 4.0;
    const double kdx = 2.0 * m_coeff * length_scale;
    if (std::abs(kdx) > maxkdx) {
      m_coeff = 0.0;
    }
  }

  std::vector<double> getMomentResiduals(const bool& scaled) {
    IRL2D::Moments lm_opt = IRL2D::ComputeMoments(m_cell, parabola_opt);
    IRL2D::Moments gm_opt = IRL2D::ComputeMoments(m_cell) - lm_opt;
    double lvf_opt = lm_opt.m0() / IRL2D::ComputeArea(m_cell);
    IRL2D::Vec lc_opt = lm_opt.m1() / lm_opt.m0();
    IRL2D::Vec gc_opt = gm_opt.m1() / gm_opt.m0();
    RecenterMoments(&lm_opt, lcs);
    RecenterMoments(&gm_opt, gcs);
    IRL2D::Mat lm2_opt = lm_opt.m2();
    IRL2D::Mat gm2_opt = gm_opt.m2();

    // calculating residuals
    std::vector<double> residual(10);

    if (scaled == true) {
      residual[0] = std::abs(lc_opt[0] - lcs[0]) / length_scale;
      residual[1] = std::abs(lc_opt[1] - lcs[1]) / length_scale;
      residual[2] = std::abs(gc_opt[0] - gcs[0]) / length_scale;
      residual[3] = std::abs(gc_opt[1] - gcs[1]) / length_scale;
      residual[4] = std::abs(lm2_opt[0][0] - lm2s[0][0]) /
                    (lm_opt.m0() * std::pow(length_scale, 2.0));
      residual[5] = std::abs(lm2_opt[1][1] - lm2s[1][1]) /
                    (lm_opt.m0() * std::pow(length_scale, 2.0));
      residual[6] = std::abs(lm2_opt[1][0] - lm2s[1][0]) /
                    (lm_opt.m0() * std::pow(length_scale, 2.0));
      residual[7] = std::abs(gm2_opt[0][0] - gm2s[0][0]) /
                    (gm_opt.m0() * std::pow(length_scale, 2.0));
      residual[8] = std::abs(gm2_opt[1][1] - gm2s[1][1]) /
                    (gm_opt.m0() * std::pow(length_scale, 2.0));
      residual[9] = std::abs(gm2_opt[1][0] - gm2s[1][0]) /
                    (gm_opt.m0() * std::pow(length_scale, 2.0));
    } else {
      residual[0] = std::abs(lc_opt[0] - lcs[0]);
      residual[1] = std::abs(lc_opt[1] - lcs[1]);
      residual[2] = std::abs(gc_opt[0] - gcs[0]);
      residual[3] = std::abs(gc_opt[1] - gcs[1]);
      residual[4] = std::abs(lm2_opt[0][0] - lm2s[0][0]);
      residual[5] = std::abs(lm2_opt[1][1] - lm2s[1][1]);
      residual[6] = std::abs(lm2_opt[1][0] - lm2s[1][0]);
      residual[7] = std::abs(gm2_opt[0][0] - gm2s[0][0]);
      residual[8] = std::abs(gm2_opt[1][1] - gm2s[1][1]);
      residual[9] = std::abs(gm2_opt[1][0] - gm2s[1][0]);
    }

    return residual;
  }
};

IRL2D::Parabola getParabola(const std::vector<double>& x, const OptData& optd) {
  auto rotation = IRL2D::ReferenceFrame(x[0]);
  auto new_frame = IRL2D::ReferenceFrame(rotation * optd.m_frame[0],
                                         rotation * optd.m_frame[1]);
  double new_coeff = optd.m_coeff + x[1] / optd.length_scale;
  IRL2D::Vec new_datum =
      optd.m_datum + IRL2D::Vec(x[2], x[3]) * optd.length_scale;
  return IRL2D::Parabola(new_datum, new_frame, new_coeff);
}

// objective function
double objective(const std::vector<double>& x, std::vector<double>& grad,
                 void* data) {
  auto* optd = static_cast<OptData*>(data);

  // compute quantities based on parabola guess
  IRL2D::Parabola p = getParabola(x, *optd);
  IRL2D::Moments lm = IRL2D::ComputeMoments(optd->m_cell, p);
  IRL2D::Moments gm = IRL2D::ComputeMoments(optd->m_cell) - lm;
  IRL2D::Vec lc = lm.m1() / lm.m0();
  IRL2D::Vec gc = gm.m1() / gm.m0();
  RecenterMoments(&lm, optd->lcs);
  RecenterMoments(&gm, optd->gcs);
  IRL2D::Mat lm2 = lm.m2();
  IRL2D::Mat gm2 = gm.m2();

  //  cost function
  IRL2D::Vec d_lc = (lc - optd->lcs) / optd->length_scale;
  IRL2D::Vec d_gc = (gc - optd->gcs) / optd->length_scale;
  const double eps = 1e-12;  // incase zeroth moment is almost zero
  IRL2D::Mat d_lm2 =
      (lm2 - optd->lm2s) *
      (1.0 / (lm.m0() * std::pow(optd->length_scale, 2.0) + eps));
  IRL2D::Mat d_gm2 =
      (gm2 - optd->gm2s) *
      (1.0 / (gm.m0() * std::pow(optd->length_scale, 2.0) + eps));

  double eq = std::pow(d_lc[0], 2.0) + std::pow(d_lc[1], 2.0) +
              std::pow(d_gc[0], 2.0) + std::pow(d_gc[1], 2.0) +
              std::pow(d_lm2[0][0], 2.0) + std::pow(d_lm2[1][0], 2.0) +
              std::pow(d_lm2[1][1], 2.0) + std::pow(d_gm2[0][0], 2.0) +
              std::pow(d_gm2[1][0], 2.0) + std::pow(d_gm2[1][1], 2.0);

  return eq;
}

// volume fraction constraint
double vf_constraint(const std::vector<double>& x, std::vector<double>& grad,
                     void* data) {
  auto* optd = static_cast<OptData*>(data);
  IRL2D::Parabola p = getParabola(x, *optd);
  double vfrac = IRL2D::ComputeVFrac(optd->m_cell, p);
  return (vfrac - optd->lvfs);
}

// confining datum within cell (inequalities)
double datum_xmin(const std::vector<double>& x, std::vector<double>& grad,
                  void* data) {
  auto* optd = static_cast<OptData*>(data);
  return (optd->m_datum[0] + x[2] * optd->length_scale) - optd->x0[0];  // left
}

double datum_xmax(const std::vector<double>& x, std::vector<double>& grad,
                  void* data) {
  auto* optd = static_cast<OptData*>(data);
  return optd->x1[0] - (optd->m_datum[0] + x[2] * optd->length_scale);  // right
}

double datum_ymin(const std::vector<double>& x, std::vector<double>& grad,
                  void* data) {
  auto* optd = static_cast<OptData*>(data);
  return (optd->m_datum[1] + x[3] * optd->length_scale) - optd->x0[1];  // top
}

double datum_ymax(const std::vector<double>& x, std::vector<double>& grad,
                  void* data) {
  auto* optd = static_cast<OptData*>(data);
  return optd->x1[1] -
         (optd->m_datum[1] + x[3] * optd->length_scale);  // bottom
}

// inline void enforce_bounds_order(double& low, double& up) {
//   if (low > up) std::swap(low, up);
// }

void NLOPT::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                              const Data<IRL2D::Moments>& a_gas_moments,
                              const double a_dt, const Data<double>& a_U,
                              const Data<double>& a_V,
                              Data<IRL2D::Parabola>* a_interface) {
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);

  // storing PLIC
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();

  // flag for mixed cells
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }

  // csv for storing moment residuals
  std::string dir = "/home/parinht2/Documents/testing code/moment_residuals";
  std::string filepath = dir + "/mom_residual.csv";
  std::ofstream csvfile(filepath);

  // header
  csvfile << "i,j,Vol Frac,Liq xc,Liq yc,Gas xc,Gas yc,"
          << "Liq M2[0][0],Liq M2[1][1],Liq M2[1][0],"
          << "Gas M2[0][0],Gas M2[1][1],Gas M2[1][0]\n";

  // parameters for particle method
  int N = 5;
  double Hp = 3.0;
  double h = mesh.dx();
  double eta = 0.5;

  // particle method for coefficient guess
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (band(i, j) == 1) {
        IRL2D::Parabola target_interface = plic(i, j);
        IRL2D::BezierList target_cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        IRL2D::Vec n_target = target_interface.frame()[1];
        std::vector<IRL2D::Parabola> all_interfaces;
        std::vector<IRL2D::BezierList> all_cells;
        for (int jj = -2; jj <= 2; jj++) {
          for (int ii = -2; ii <= 2; ii++) {
            if (band(ii + i, jj + j) == 1) {
              IRL2D::Vec n_neighbor = plic(ii + i, jj + j).frame()[1];
              if (n_neighbor * n_target <= -0.5) {
                continue;  // filtering based on normals
              }
              all_interfaces.push_back(plic(ii + i, jj + j));
              all_cells.push_back(IRL2D::RectangleFromBounds(
                  IRL2D::Vec(mesh.x(ii + i), mesh.y(jj + j)),
                  IRL2D::Vec(mesh.x(ii + i + 1), mesh.y(jj + j + 1))));
            }
          }
        }
        // center point of target plic (to use as datum)
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(target_cell, target_interface, true);
        target_interface.datum() =
            (clipped_plic[0].first + clipped_plic[1].first) / 2.0;
        target_interface.coeff() =
            0.5 * IRL2D::getCurvature(target_interface, target_cell,
                                      all_interfaces, all_cells, N, Hp, h, eta);

        const double INF = std::numeric_limits<double>::infinity();
        OptData optFunctor(target_cell, target_interface,
                           a_liquid_moments(i, j), a_gas_moments(i, j));
        optFunctor.setframe();

        // bounds for variables
        double rot_low = -M_PI, rot_up = M_PI;
        // double coeff_low = -INF, coeff_up = INF;
        double coeff_margin =
            std::max(0.1 * std::abs(target_interface.coeff()), 1e-6);
        double coeff_low = -coeff_margin * optFunctor.length_scale;
        double coeff_up = coeff_margin * optFunctor.length_scale;
        double dx_low = (optFunctor.x0[0] - optFunctor.m_datum[0]) /
                        optFunctor.length_scale;
        double dx_up = (optFunctor.x1[0] - optFunctor.m_datum[0]) /
                       optFunctor.length_scale;
        double dy_low = (optFunctor.x0[1] - optFunctor.m_datum[1]) /
                        optFunctor.length_scale;
        double dy_up = (optFunctor.x1[1] - optFunctor.m_datum[1]) /
                       optFunctor.length_scale;
        // enforce_bounds_order(dx_low, dx_up);
        // enforce_bounds_order(dy_low, dy_up);
        std::vector<double> lb = {rot_low, coeff_low, dx_low, dy_low};
        std::vector<double> ub = {rot_up, coeff_up, dx_up, dy_up};

        nlopt::opt opt(nlopt::LN_COBYLA, 4);  // 4 independent variables
        // nlopt::opt opt(nlopt::GN_ISRES, 4);
        opt.set_min_objective(objective, &optFunctor);
        opt.add_equality_constraint(vf_constraint, &optFunctor, 1e-14);
        // opt.add_inequality_constraint(datum_xmin, &optFunctor, 1e-8);
        // opt.add_inequality_constraint(datum_xmax, &optFunctor, 1e-8);
        // opt.add_inequality_constraint(datum_ymin, &optFunctor, 1e-8);
        // opt.add_inequality_constraint(datum_ymax, &optFunctor, 1e-8);
        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(ub);

        // varying tolerance types
        // opt.set_xtol_rel(1e-12);
        opt.set_xtol_abs(1e-12);
        opt.set_ftol_abs(1e-12);
        opt.set_maxeval(1000);

        std::vector<double> x = {0.0, 0.0, 0.0, 0.0};
        double obj_func_val;  // value of objective function for optimal x
        nlopt::result result = opt.optimize(x, obj_func_val);
        // std::cout << opt.get_numevals() << std::endl;

        (*a_interface)(i, j) = getParabola(x, optFunctor);
        (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(
            target_cell, (*a_interface)(i, j), optFunctor.lvfs);
        optFunctor.parabola_opt = (*a_interface)(i, j);

        // diagnostics
        // ------------------------------------------------------------------------
        bool outputScaledResiduals = false;
        std::vector<double> residuals =
            optFunctor.getMomentResiduals(outputScaledResiduals);
        // std::cout << "(" << i << "," << j << ") --> " <<
        // optFunctor.lvfs - IRL2D::ComputeVFrac(target_cell,
        // (*a_interface)(i,j)) << " --> " << obj_func_val << " --> "; std::cout
        // << "[ "; for (auto r : residuals){
        //   std::cout << r << " " ;
        // }
        // std::cout << "]" << std::endl;

        // writing diagnostics to csv file
        csvfile << i << "," << j << ","
                << optFunctor.lvfs -
                       IRL2D::ComputeVFrac(target_cell, (*a_interface)(i, j))
                << ",";
        for (auto r : residuals) {
          csvfile << std::scientific << std::setprecision(6) << r << ",";
        }
        csvfile << "\n";
      }
    }
  }
  csvfile.close();

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// least squares fitting of circles
// -----------------------------------------------------------------------

struct InterfaceData {
  bool mixed = false;
  int xIndex, yIndex;
  double vf;
  IRL2D::Vec a, b, center;  // start, end, and midpoint
  IRL2D::BezierList rectangle;
};

// Pratt fit
void Pratt::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                              const Data<IRL2D::Moments>& a_gas_moments,
                              const double a_dt, const Data<double>& a_U,
                              const Data<double>& a_V,
                              Data<IRL2D::Parabola>* a_interface) {
  // plic reconstruction
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();
  const double h = mesh.dx();

  // storing interface data
  Data<InterfaceData> plicData(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        plicData(i, j).mixed = true;
        plicData(i, j).vf = liquid_volume_fraction;
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        plicData(i, j).rectangle = cell;
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(cell, plic(i, j), true);
        plicData(i, j).a = clipped_plic[0].first;
        plicData(i, j).b = clipped_plic[1].first;
        plicData(i, j).center = (plicData(i, j).a + plicData(i, j).b) / 2.0;
        plicData(i, j).xIndex = i;
        plicData(i, j).yIndex = j;
      }
    }
  }

  // Pratt circle fit for curvature (parabola coefficient)
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (plicData(i, j).mixed == true) {
        IRL2D::Vec pref = plicData(i, j).center;
        IRL2D::Vec nref = plic(i, j).frame()[1];
        std::vector<double> vfw, dw, vfracs, nw;
        std::vector<IRL2D::Vec> ploc, nloc;
        std::vector<std::pair<IRL2D::Vec, IRL2D::Vec>> line_seg_endpoints;
        for (int ii = -2; ii <= 2; ii++) {
          for (int jj = -2; jj <= 2; jj++) {
            if (plicData(i + ii, j + jj).mixed == true) {
              line_seg_endpoints.push_back(
                  {plicData(i + ii, j + jj).a, plicData(i + ii, j + jj).b});
              vfracs.push_back(plicData(i + ii, j + jj).vf);
              ploc.push_back(plicData(i + ii, j + jj).center);
              nloc.push_back(plic(i + ii, j + jj).frame()[1]);
            }
          }
        }
        std::vector<IRL2D::Vec> pts = IRL2D::generatePoints(line_seg_endpoints);
        // computing weights
        for (int k = 0; k < line_seg_endpoints.size(); k++) {
          double vf_weight = IRL2D::getVfracWeight(vfracs[k]);
          double d_weight = IRL2D::getDistanceWeight(pref, ploc[k], h);
          double n_weight = IRL2D::getNormalWeight(nref, nloc[k]);
          // n_weight = 1.0;
          // double n_weight = IRL2D::getNormalGradWeight(nref, nloc[k], pref,
          // ploc[k]);
          for (int p = 0; p < (pts.size() / line_seg_endpoints.size()); p++) {
            vfw.push_back(vf_weight);
            dw.push_back(d_weight);
            nw.push_back(n_weight);
          }
        }
        // Pratt's parabola
        IRL2D::Parabola Prattparabola = IRL2D::getPrattParabola_localframe(
            pts, vfw, dw, nw, plic(i, j).frame(), plicData(i, j).center);
        // vf matching
        (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(
            plicData(i, j).rectangle, Prattparabola, plicData(i, j).vf);
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// Taubin fit
void Taubin::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  // plic reconstruction
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();
  const double h = mesh.dx();

  // storing interface data
  Data<InterfaceData> plicData(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        plicData(i, j).mixed = true;
        plicData(i, j).vf = liquid_volume_fraction;
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        plicData(i, j).rectangle = cell;
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(cell, plic(i, j), true);
        plicData(i, j).a = clipped_plic[0].first;
        plicData(i, j).b = clipped_plic[1].first;
        plicData(i, j).center = (plicData(i, j).a + plicData(i, j).b) / 2.0;
        plicData(i, j).xIndex = i;
        plicData(i, j).yIndex = j;
      }
    }
  }

  // Taubin circle fit for curvature (parabola coefficient)
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (plicData(i, j).mixed == true) {
        IRL2D::Vec pref = plicData(i, j).center;
        IRL2D::Vec nref = plic(i, j).frame()[1];
        std::vector<double> vfw, dw, vfracs, nw;
        std::vector<IRL2D::Vec> ploc, nloc;
        std::vector<std::pair<IRL2D::Vec, IRL2D::Vec>> line_seg_endpoints;
        for (int ii = -2; ii <= 2; ii++) {
          for (int jj = -2; jj <= 2; jj++) {
            if (plicData(i + ii, j + jj).mixed == true) {
              line_seg_endpoints.push_back(
                  {plicData(i + ii, j + jj).a, plicData(i + ii, j + jj).b});
              vfracs.push_back(plicData(i + ii, j + jj).vf);
              ploc.push_back(plicData(i + ii, j + jj).center);
              nloc.push_back(plic(i + ii, j + jj).frame()[1]);
            }
          }
        }
        std::vector<IRL2D::Vec> pts = IRL2D::generatePoints(line_seg_endpoints);
        // computing weights
        for (int k = 0; k < line_seg_endpoints.size(); k++) {
          double vf_weight = IRL2D::getVfracWeight(vfracs[k]);
          double d_weight = IRL2D::getDistanceWeight(pref, ploc[k], h);
          double n_weight = IRL2D::getNormalWeight(nref, nloc[k]);
          for (int p = 0; p < (pts.size() / line_seg_endpoints.size()); p++) {
            vfw.push_back(vf_weight);
            dw.push_back(d_weight);
            nw.push_back(n_weight);
          }
        }
        // Pratt's parabola
        IRL2D::Parabola Taubinparabola = IRL2D::getTaubinParabola_localframe(
            pts, vfw, dw, nw, plic(i, j).frame(), plicData(i, j).center);
        // vf matching
        (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(
            plicData(i, j).rectangle, Taubinparabola, plicData(i, j).vf);
        // (*a_interface)(i, j) = Taubinparabola;
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// LMA fit
struct LMAFunctor {
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // variables
  const std::vector<IRL2D::Vec>& points;
  const std::vector<double>&vfw, dw, nw;

  // constructor
  LMAFunctor(const std::vector<IRL2D::Vec>& pts,
             const std::vector<double>& vf_w, const std::vector<double>& d_w,
             const std::vector<double>& n_w)
      : points(pts), vfw(vf_w), dw(d_w), nw(n_w) {}

  int inputs() const { return 3; }  // A, D, theta
  int values() const { return static_cast<int>(points.size()); }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    double A = x(0), D = x(1), theta = x(2);
    double E = std::sqrt(1 + 4 * A * D);
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    for (int i = 0; i < values(); ++i) {
      double xi = points[i].x(), yi = points[i].y();
      double zi = xi * xi + yi * yi;
      double ui = xi * cos_theta + yi * sin_theta;
      double Pi = A * zi + E * ui + D;
      double Qi = std::sqrt(1.0 + 4.0 * A * Pi);
      fvec(i) = 2.0 * Pi / (1.0 + Qi) * std::sqrt(vfw[i] * dw[i] * nw[i]);
    }
    return 0;
  }

  // analytical Jacobian
  int df(const InputType& x, JacobianType& fjac) const {
    double A = x(0), D = x(1), theta = x(2);
    double E = std::sqrt(1.0 + 4.0 * A * D);
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    for (int i = 0; i < values(); ++i) {
      double xi = points[i].x(), yi = points[i].y();
      double zi = xi * xi + yi * yi;
      double ui = xi * cos_theta + yi * sin_theta;
      double Pi = A * zi + E * ui + D;
      double Qi = std::sqrt(1 + 4 * A * Pi);
      double di = 2.0 * Pi / (1.0 + Qi);
      double Ri = 2.0 * (1.0 - A * di / Qi) / (Qi + 1.0);

      double wi = std::sqrt(vfw[i] * dw[i] * nw[i]);  // weights

      fjac(i, 0) =
          wi * ((zi + 2.0 * D * ui / E) * Ri - (di * di / Qi));         // ∂d/∂A
      fjac(i, 1) = wi * ((2.0 * A * ui / E + 1.0) * Ri);                // ∂d/∂D
      fjac(i, 2) = wi * ((-xi * sin_theta + yi * cos_theta) * E * Ri);  // ∂d/∂θ
    }
    return 0;
  }
};

void LMA::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                            const Data<IRL2D::Moments>& a_gas_moments,
                            const double a_dt, const Data<double>& a_U,
                            const Data<double>& a_V,
                            Data<IRL2D::Parabola>* a_interface) {
  // plic reconstruction
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();
  const double h = mesh.dx();

  // storing interface data
  Data<InterfaceData> plicData(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        plicData(i, j).mixed = true;
        plicData(i, j).vf = liquid_volume_fraction;
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        plicData(i, j).rectangle = cell;
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(cell, plic(i, j), true);
        plicData(i, j).a = clipped_plic[0].first;
        plicData(i, j).b = clipped_plic[1].first;
        plicData(i, j).center = (plicData(i, j).a + plicData(i, j).b) / 2.0;
        plicData(i, j).xIndex = i;
        plicData(i, j).yIndex = j;
      }
    }
  }

  // LMA fit
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (plicData(i, j).mixed == true) {
        IRL2D::Vec pref = plicData(i, j).center;
        IRL2D::Vec nref = plic(i, j).frame()[1];
        std::vector<double> vfw, dw, vfracs, nw;
        std::vector<IRL2D::Vec> ploc, nloc;
        std::vector<std::pair<IRL2D::Vec, IRL2D::Vec>> line_seg_endpoints;
        for (int ii = -2; ii <= 2; ii++) {
          for (int jj = -2; jj <= 2; jj++) {
            if (plicData(i + ii, j + jj).mixed == true) {
              line_seg_endpoints.push_back(
                  {plicData(i + ii, j + jj).a, plicData(i + ii, j + jj).b});
              vfracs.push_back(plicData(i + ii, j + jj).vf);
              ploc.push_back(plicData(i + ii, j + jj).center);
              nloc.push_back(plic(i + ii, j + jj).frame()[1]);
            }
          }
        }
        std::vector<IRL2D::Vec> pts = IRL2D::generatePoints(line_seg_endpoints);
        // computing weights
        for (int k = 0; k < line_seg_endpoints.size(); k++) {
          double vf_weight = IRL2D::getVfracWeight(vfracs[k]);
          double d_weight = IRL2D::getDistanceWeight(pref, ploc[k], h);
          double n_weight = IRL2D::getNormalWeight(nref, nloc[k]);
          for (int p = 0; p < (pts.size() / line_seg_endpoints.size()); p++) {
            vfw.push_back(vf_weight);
            dw.push_back(d_weight);
            nw.push_back(n_weight);
          }
        }
        // initial guess for A, D and theta
        // std::vector<double> Pratt_params = IRL2D::getPrattParams(pts, vfw,
        // dw, nw);
        std::vector<double> Pratt_params =
            IRL2D::getTaubinParams(pts, vfw, dw, nw);
        Eigen::VectorXd x(3);
        x << Pratt_params[0], Pratt_params[1], Pratt_params[2];

        // Levenberg-Marquardt
        LMAFunctor functor(pts, vfw, dw, nw);
        Eigen::LevenbergMarquardt<LMAFunctor> lm(functor);
        lm.parameters.maxfev = 1000;
        lm.parameters.xtol = 1e-12;
        lm.minimize(x);

        // extracting cell center and radius
        double A = x(0), D = x(1), theta = x(2);
        double E = std::sqrt(1.0 + 4.0 * A * D);
        double B = E * std::cos(theta);
        double C = E * std::sin(theta);
        IRL2D::Vec circle_center = {-B / (2.0 * A), -C / (2.0 * A)};
        double radius = std::fabs(1.0 / (2.0 * A));

        // LMA parabola coefficient
        IRL2D::Parabola LMA_parabola;
        LMA_parabola.coeff() = 0.5 / radius;

        // reference frame
        bool flip_coeff = false;
        LMA_parabola.frame() =
            IRL2D::estimateFrame(circle_center, plicData(i, j).center, radius,
                                 plic(i, j).frame()[1], flip_coeff);
        if (flip_coeff == true) {
          LMA_parabola.coeff() = -LMA_parabola.coeff();
        }

        // datum
        IRL2D::Vec direction = plicData(i, j).center - circle_center;
        direction.normalize();
        LMA_parabola.datum() = circle_center + radius * direction;

        // vf matching
        (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(
            plicData(i, j).rectangle, LMA_parabola, plicData(i, j).vf);
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// variable support radius
void Pratt2::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  // plic reconstruction
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();
  const double h = mesh.dx();

  // storing interface data
  Data<InterfaceData> plicData(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        plicData(i, j).mixed = true;
        plicData(i, j).vf = liquid_volume_fraction;
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        plicData(i, j).rectangle = cell;
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(cell, plic(i, j), true);
        plicData(i, j).a = clipped_plic[0].first;
        plicData(i, j).b = clipped_plic[1].first;
        plicData(i, j).center = (plicData(i, j).a + plicData(i, j).b) / 2.0;
        plicData(i, j).xIndex = i;
        plicData(i, j).yIndex = j;
      }
    }
  }

  // Pratt circle fit for curvature (parabola coefficient)
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (plicData(i, j).mixed == true) {
        IRL2D::Vec pref = plicData(i, j).center;
        IRL2D::Vec nref = plic(i, j).frame()[1];
        std::vector<double> vfw, dw, vfracs, nw;
        std::vector<IRL2D::Vec> ploc, nloc;
        std::vector<std::pair<IRL2D::Vec, IRL2D::Vec>> line_seg_endpoints;
        for (int ii = -2; ii <= 2; ii++) {
          for (int jj = -2; jj <= 2; jj++) {
            if (plicData(i + ii, j + jj).mixed == true) {
              line_seg_endpoints.push_back(
                  {plicData(i + ii, j + jj).a, plicData(i + ii, j + jj).b});
              vfracs.push_back(plicData(i + ii, j + jj).vf);
              ploc.push_back(plicData(i + ii, j + jj).center);
              nloc.push_back(plic(i + ii, j + jj).frame()[1]);
            }
          }
        }
        std::vector<IRL2D::Vec> pts = IRL2D::generatePoints(line_seg_endpoints);

        // varying support radius
        std::vector<double> delta = {1.5, 2.5, 3.5};
        double min_centroid_error = std::numeric_limits<double>::max();
        IRL2D::Parabola best_parabola;
        for (int d = 0; d < delta.size(); d++) {
          // computing weights
          for (int k = 0; k < line_seg_endpoints.size(); k++) {
            double vf_weight = IRL2D::getVfracWeight(vfracs[k]);
            double d_weight = IRL2D::DistanceWeight(pref, ploc[k], h, delta[d]);
            double n_weight = IRL2D::getNormalWeight(nref, nloc[k]);
            // n_weight = 1.0;
            for (int p = 0; p < (pts.size() / line_seg_endpoints.size()); p++) {
              vfw.push_back(vf_weight);
              dw.push_back(d_weight);
              nw.push_back(n_weight);
            }
          }
          IRL2D::Parabola Prattparabola = IRL2D::getPrattParabola_localframe(
              pts, vfw, dw, nw, plic(i, j).frame(), plicData(i, j).center);
          Prattparabola = IRL2D::MatchToVolumeFraction(
              plicData(i, j).rectangle, Prattparabola, plicData(i, j).vf);

          IRL2D::Vec trueCentroid =
              a_liquid_moments(i, j).m1() / a_liquid_moments(i, j).m0();
          IRL2D::Moments PrattMoments =
              IRL2D::ComputeMoments(plicData(i, j).rectangle, Prattparabola);
          IRL2D::Vec computedCentroid = PrattMoments.m1() / PrattMoments.m0();
          double centroid_dist_error =
              (trueCentroid - computedCentroid).magnitude() / h;

          if (centroid_dist_error < min_centroid_error) {
            min_centroid_error = centroid_dist_error;
            best_parabola = Prattparabola;
          }
        }
        (*a_interface)(i, j) = best_parabola;
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// iterative circle fit
void Pratt3::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  // plic reconstruction
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);

  const BasicMesh& mesh = a_U.getMesh();
  const double h = mesh.dx();

  // params for iterative fit
  int fit_iter = 0;
  const int max_fit_iter = 1;

  while (fit_iter < max_fit_iter) {
    fit_iter++;
    Data<IRL2D::Parabola> interface = *a_interface;

    // interface data
    Data<InterfaceData> interfaceData(&mesh);
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        const double liquid_volume_fraction =
            (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          interfaceData(i, j).mixed = true;
          interfaceData(i, j).vf = liquid_volume_fraction;
          IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
              IRL2D::Vec(mesh.x(i), mesh.y(j)),
              IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
          interfaceData(i, j).rectangle = cell;
          IRL2D::BezierList clipped_plic =
              IRL2D::ParabolaClip(cell, interface(i, j), true);
          interfaceData(i, j).a = clipped_plic[0].first;
          interfaceData(i, j).b = clipped_plic[1].first;
          if (fit_iter == 1) {
            interfaceData(i, j).center =
                (interfaceData(i, j).a + interfaceData(i, j).b) / 2.0;
          } else {
            interfaceData(i, j).center = IRL2D::getParabolaCenter(
                {interfaceData(i, j).a, interfaceData(i, j).b},
                interface(i, j));
          }
          interfaceData(i, j).xIndex = i;
          interfaceData(i, j).yIndex = j;
        }
      }
    }

    // least squares circle fit
    for (int i = mesh.imin(); i <= mesh.imax(); i++) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
        if (interfaceData(i, j).mixed == true) {
          IRL2D::Vec pref = interfaceData(i, j).center;
          IRL2D::Vec nref = interface(i, j).frame()[1];
          std::vector<double> vfw, dw, vfracs, nw;
          std::vector<IRL2D::Vec> ploc, nloc;
          std::vector<std::pair<IRL2D::Vec, IRL2D::Vec>> line_seg_endpoints;
          std::vector<IRL2D::Parabola> interfaces;
          for (int ii = -2; ii <= 2; ii++) {
            for (int jj = -2; jj <= 2; jj++) {
              if (interfaceData(i + ii, j + jj).mixed == true) {
                line_seg_endpoints.push_back({interfaceData(i + ii, j + jj).a,
                                              interfaceData(i + ii, j + jj).b});
                vfracs.push_back(interfaceData(i + ii, j + jj).vf);
                ploc.push_back(interfaceData(i + ii, j + jj).center);
                nloc.push_back(interface(i + ii, j + jj).frame()[1]);
                interfaces.push_back(interface(i + ii, j + jj));
              }
            }
          }
          std::vector<IRL2D::Vec> pts;
          if (fit_iter == 1) {
            pts = IRL2D::generatePoints(line_seg_endpoints);
          } else {
            pts = IRL2D::generateParabolaPoints(line_seg_endpoints, interfaces);
          }

          // computing weights
          for (int k = 0; k < line_seg_endpoints.size(); k++) {
            double vf_weight = IRL2D::getVfracWeight(vfracs[k]);
            double d_weight = IRL2D::getDistanceWeight(pref, ploc[k], h);
            double n_weight = IRL2D::getNormalWeight(nref, nloc[k]);
            n_weight = 1.0;
            // double n_weight = IRL2D::getNormalGradWeight(nref, nloc[k], pref,
            // ploc[k]);
            for (int p = 0; p < (pts.size() / line_seg_endpoints.size()); p++) {
              vfw.push_back(vf_weight);
              dw.push_back(d_weight);
              nw.push_back(n_weight);
            }
          }
          // Pratt's parabola
          IRL2D::Parabola Prattparabola = IRL2D::getPrattParabola_localframe(
              pts, vfw, dw, nw, interface(i, j).frame(),
              interfaceData(i, j).center);
          // vf matching
          (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(
              interfaceData(i, j).rectangle, Prattparabola,
              interfaceData(i, j).vf);
        }
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// cubic approximation

struct interfaceInfo {
  bool mixed = false;
  double vf, nx, ny;
  IRL2D::Vec a, b, center;
  IRL2D::BezierList rectangle;
};

Eigen::MatrixXd pseudoInverse(const Eigen::MatrixXd& A, double tol = 1e-9) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const auto& S = svd.singularValues();
  Eigen::MatrixXd S_inv =
      Eigen::MatrixXd::Zero(svd.matrixV().cols(), svd.matrixU().cols());
  for (int i = 0; i < S.size(); ++i) {
    if (S(i) > tol) S_inv(i, i) = 1.0 / S(i);
  }
  return svd.matrixV() * S_inv * svd.matrixU().transpose();
}

double f(double x, double A, double B, double C, double D) {
  return A * x * x * x + B * x * x + C * x + D;
}
double fprime(double x, double A, double B, double C) {
  return 3.0 * A * x * x + 2.0 * B * x + C;
}
double fdoubleprime(double x, double A, double B) {
  return 6.0 * A * x + 2.0 * B;
}

struct ClosestPointFunctor {
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // variables
  double A, B, C, D;

  // constructor
  ClosestPointFunctor(const double& A_, const double& B_, const double& C_,
                      const double& D_)
      : A(A_), B(B_), C(C_), D(D_) {}

  int inputs() const { return 1; }
  int values() const { return 1; }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    double xp = x(0);
    double yp = A * xp * xp * xp + B * xp * xp + C * xp + D;
    ;
    fvec(0) = std::sqrt(xp * xp + yp * yp);
    return 0;
  }
};

void Cubic::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                              const Data<IRL2D::Moments>& a_gas_moments,
                              const double a_dt, const Data<double>& a_U,
                              const Data<double>& a_V,
                              Data<IRL2D::Parabola>* a_interface) {
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);
  Data<IRL2D::Parabola> plic = *a_interface;

  const BasicMesh& mesh = a_U.getMesh();

  // local <-> global
  auto globalToLocal = [&](const IRL2D::Vec& p,
                           const IRL2D::Parabola& interface) -> IRL2D::Vec {
    IRL2D::Vec t = interface.frame()[0], n = interface.frame()[1];
    IRL2D::Vec ploc = {(p.x() - interface.datum().x()) * t.x() +
                           (p.y() - interface.datum().y()) * t.y(),
                       (p.x() - interface.datum().x()) * n.x() +
                           (p.y() - interface.datum().y()) * n.y()};
    return ploc;
  };
  auto localToGlobal = [&](const IRL2D::Vec& ploc,
                           const IRL2D::Parabola& interface) -> IRL2D::Vec {
    IRL2D::Vec t = interface.frame()[0], n = interface.frame()[1];
    IRL2D::Vec p = {
        interface.datum().x() + t.x() * ploc.x() + n.x() * ploc.y(),
        interface.datum().y() + t.y() * ploc.x() + n.y() * ploc.y()};
    return p;
  };
  auto vectorGlobalToLocal =
      [&](const IRL2D::Vec& v, const IRL2D::Parabola& interface) -> IRL2D::Vec {
    IRL2D::Vec t = interface.frame()[0];
    IRL2D::Vec n = interface.frame()[1];
    IRL2D::Vec vloc = {v.x() * t.x() + v.y() * t.y(),
                       v.x() * n.x() + v.y() * n.y()};
    return vloc;
  };
  auto vectorLocalToGlobal =
      [&](const IRL2D::Vec& vloc,
          const IRL2D::Parabola& interface) -> IRL2D::Vec {
    IRL2D::Vec t = interface.frame()[0], n = interface.frame()[1];
    IRL2D::Vec v = {t.x() * vloc.x() + n.x() * vloc.y(),
                    t.y() * vloc.x() + n.y() * vloc.y()};
    return v;
  };

  // interface information
  Data<interfaceInfo> plicData(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        plicData(i, j).mixed = true;
        plicData(i, j).vf = liquid_volume_fraction;
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        plicData(i, j).rectangle = cell;
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(cell, plic(i, j), true);
        plicData(i, j).a = clipped_plic[0].first;
        plicData(i, j).b = clipped_plic[1].first;
        plicData(i, j).center = (plicData(i, j).a + plicData(i, j).b) / 2.0;
        plicData(i, j).nx = plic(i, j).frame()[1][0];
        plicData(i, j).ny = plic(i, j).frame()[1][1];
        plic(i, j).datum() =
            plicData(i, j).center;  // for coordinate transformation
      }
    }
  }

  const int num_points = 10;  // per plic
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (plicData(i, j).mixed == true) {
        IRL2D::Vec pref = plicData(i, j).center;
        IRL2D::Vec pref_local = globalToLocal(pref, plic(i, j));
        IRL2D::Vec nref_local =
            vectorGlobalToLocal(plic(i, j).frame()[1], plic(i, j));
        std::vector<double> nxs, nys, weights;
        std::vector<IRL2D::Vec> points;
        for (int ii = -2; ii <= 2; ii++) {
          for (int jj = -2; jj <= 2; jj++) {
            if (plicData(i + ii, j + jj).mixed == true) {
              IRL2D::Vec ploc_local =
                  globalToLocal(plicData(i + ii, j + jj).center, plic(i, j));
              IRL2D::Vec n_local = vectorGlobalToLocal(
                  plic(i + ii, j + jj).frame()[1], plic(i, j));
              // std::cout << n_local << std::endl;
              IRL2D::Vec a_local =
                  globalToLocal(plicData(i + ii, j + jj).a, plic(i, j));
              IRL2D::Vec b_local =
                  globalToLocal(plicData(i + ii, j + jj).b, plic(i, j));
              std::vector<IRL2D::Vec> pts =
                  IRL2D::getPoints({a_local, b_local}, num_points);
              double vf_weight =
                  IRL2D::getVfracWeight(plicData(i + ii, j + jj).vf);
              double d_weight =
                  IRL2D::getDistanceWeight(pref_local, ploc_local, mesh.dx());
              double n_weight = IRL2D::getNormalWeight(nref_local, n_local);
              // vf_weight = 1.0; d_weight =1.0;
              for (int k = 0; k < pts.size(); k++) {
                points.push_back(pts[k]);
                weights.push_back(vf_weight * d_weight * n_weight);
                nxs.push_back(n_local[0]);
                nys.push_back(n_local[1]);
              }
            }
          }
        }

        // printing points
        // IRL2D::Vec p_global;
        // if (i == 23 && j == 51){
        // std::cout << "x = [ ";
        // for (const auto& p : points){
        //   p_global = localToGlobal(p, plic(i,j));
        //   std::cout << p_global.x() << " ";
        // }
        // std::cout << "];" << std::endl;
        // std::cout << "y = [ ";
        // for (const auto& p : points){
        //   p_global = localToGlobal(p, plic(i,j));
        //   std::cout << p_global.y() << " ";
        // }
        // std::cout << "];" << std::endl;
        // std::cout << plic(i,j).datum() << std::endl;
        // std::cout << plic(i,j).frame() << std::endl;
        // std::cout << "x_loc = [ ";
        // for (const auto& p : points){
        //   std::cout << p.x() << " ";
        // }
        // std::cout << "];" << std::endl;
        // std::cout << "y_loc = [ ";
        // for (const auto& p : points){
        //   std::cout << p.y() << " ";
        // }
        // std::cout << "];" << std::endl;
        // }

        // least squares fit to find coeffs
        const int N = points.size();
        Eigen::MatrixXd W(2 * N, 2 * N);
        W.setZero();
        Eigen::MatrixXd U(2 * N, 4);
        Eigen::MatrixXd z(2 * N, 1);
        Eigen::MatrixXd A(4, 4);
        Eigen::MatrixXd b(4, 1);
        for (int k = 0; k < N; k++) {
          W(2 * k, 2 * k) = std::sqrt(weights[k]);
          W(2 * k + 1, 2 * k + 1) = std::sqrt(weights[k]);
          double xi = points[k].x();
          U.row(2 * k) << xi * xi * xi, xi * xi, xi, 1.0;
          // U.row(2*k + 1) << 3.0 * xi * xi, 2.0 * xi, 1.0, 0.0;
          U.row(2 * k + 1) << 3.0 * xi * xi * nys[k], 2.0 * xi * nys[k],
              1.0 * nys[k], 0.0;
          z(2 * k, 0) = points[k].y();
          // z(2*k + 1, 0) = -nxs[k]/nys[k];
          z(2 * k + 1, 0) = -nxs[k];
        }
        A = U.transpose() * W.transpose() * W * U;
        b = U.transpose() * W.transpose() * W * z;
        // Eigen::Vector4d coeffs = A.ldlt().solve(b);
        Eigen::Vector4d coeffs = A.colPivHouseholderQr().solve(b);

        // fit without normals
        // const int N = points.size();
        // Eigen::MatrixXd W(N, N);  W.setZero();
        // Eigen::MatrixXd U(N, 4);
        // Eigen::MatrixXd z(N, 1);
        // Eigen::MatrixXd A(4, 4);
        // Eigen::MatrixXd b(4, 1);
        // for (int k = 0; k < N; k++){
        //   W(k, k) = std::sqrt(weights[k]);
        //   double xi = points[k].x();
        //   U.row(k) << xi * xi * xi, xi * xi, xi, 1.0;
        //   z(k, 0) = points[k].y();
        // }
        // A = U.transpose() * W.transpose() * W * U;
        // b = U.transpose() * W.transpose() * W * z;
        // Eigen::Vector4d coeffs = A.colPivHouseholderQr().solve(b);

        // if (i == 23 && j == 51){
        //   std::cout << "A = " << coeffs(0) << ";" << std::endl;
        //   std::cout << "B = " << coeffs(1) << ";" << std::endl;
        //   std::cout << "C = " << coeffs(2) << ";" << std::endl;
        //   std::cout << "D = " << coeffs(3) << ";" << std::endl;
        // }

        // bool a_solution_exists = (A*coeffs).isApprox(b, 1e-6);
        // std::cout << a_solution_exists << std::endl;
        // std::cout << "coefficients = " << coeffs << std::endl;

        // pseduo inverse for finding coefficients
        // auto coeffs = pseudoInverse(W*U) * W * z;

        // finding point closest to curve in local frame
        ClosestPointFunctor functor(coeffs(0), coeffs(1), coeffs(2), coeffs(3));
        Eigen::NumericalDiff<ClosestPointFunctor> numDiff(functor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<ClosestPointFunctor>,
                                  double>
            lm(numDiff);
        Eigen::VectorXd x(1);
        x(0) = 0.0;
        lm.parameters.maxfev = 1000;
        lm.parameters.xtol = 1e-12;
        lm.minimize(x);
        double x_star = x(0);
        double y_star = f(x_star, coeffs(0), coeffs(1), coeffs(2), coeffs(3));
        IRL2D::Parabola cubic_interface;
        cubic_interface.datum() = localToGlobal({x_star, y_star}, plic(i, j));
        // if (i == 23 && j == 51){
        //   std::cout << "datum = [" << cubic_interface.datum()[0] << " "
        //             << cubic_interface.datum()[1] << "];" << std::endl;
        // }

        // curvature
        double fp = fprime(x_star, coeffs(0), coeffs(1), coeffs(2));
        double fpp = fdoubleprime(x_star, coeffs(0), coeffs(1));
        // double curvature = std::abs(fpp) / std::pow(1.0 + fp * fp, 1.5);
        double curvature = -fpp / std::pow(1.0 + fp * fp, 1.5);
        cubic_interface.coeff() = 0.5 * curvature;
        // std::cout << "curvature = " << curvature << std::endl;

        // reference frame
        IRL2D::Vec normal_local = {-fp, 1.0};
        normal_local.normalize();
        IRL2D::Vec normal = vectorLocalToGlobal(normal_local, plic(i, j));
        if ((normal * plic(i, j).frame()[1]) < 0) {
          normal *= -1.0;
          cubic_interface.coeff() *= -1.0;
        }
        IRL2D::Vec tangent = {normal.y(), -normal.x()};
        cubic_interface.frame() = {tangent, normal};

        // cubic_interface.frame() = plic(i,j).frame();

        // vf matching
        (*a_interface)(i, j) = IRL2D::MatchToVolumeFraction(
            plicData(i, j).rectangle, cubic_interface, plicData(i, j).vf);
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// partition of unity
// ----------------------------------------------------------------------------------------

void PU::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                           const Data<IRL2D::Moments>& a_gas_moments,
                           const double a_dt, const Data<double>& a_U,
                           const Data<double>& a_V,
                           Data<IRL2D::Parabola>* a_interface) {
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);
  Data<IRL2D::Parabola> plic = *a_interface;
  const BasicMesh& mesh = a_U.getMesh();
  std::vector<IRL2D::Vec> centroids;
  std::vector<IRL2D::Vec> normals;
  Data<bool> mixed(&mesh);
  Data<IRL2D::Vec> plic_center(&mesh);
  Data<std::pair<IRL2D::Vec, IRL2D::Vec>> end_points(&mesh);

  // writing plic interface to csv
  // std::string dir = "/home/parinht2/Documents/testing
  // code/partition_of_unity"; std::string filepath = dir + "/interface.csv";
  // std::string filepath = dir + "/implicit_frame.csv";
  // std::ofstream csvfile(filepath);
  // csvfile << "ax,ay,bx,by,tx,ty,nx,ny\n";
  // csvfile << "x,y,tx,ty,nx,ny\n";

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      mixed(i, j) = false;
      const double lvf = (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (lvf >= IRL::global_constants::VF_LOW &&
          lvf <= IRL::global_constants::VF_HIGH) {
        mixed(i, j) = true;
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(cell, plic(i, j), true);
        IRL2D::Vec center =
            (clipped_plic[0].first + clipped_plic[1].first) / 2.0;
        centroids.push_back(center);
        plic_center(i, j) = center;
        end_points(i, j) = {clipped_plic[0].first, clipped_plic[1].first};
        normals.push_back(plic(i, j).frame()[1]);
        // csvfile << clipped_plic[0].first[0] << "," <<
        // clipped_plic[0].first[1] << ","
        //         << clipped_plic[1].first[0] << "," <<
        //         clipped_plic[1].first[1] << ","
        //         << plic(i,j).frame()[0][0] << "," << plic(i,j).frame()[0][1]
        //         << ","
        //         << plic(i,j).frame()[1][0] << "," << plic(i,j).frame()[1][1]
        //         << "\n";
      }
    }
  }

  // projecting all centroids on implicit surface
  // std::vector<IRL2D::Vec> projected_points(centroids.size());
  // for (int i = 0; i < projected_points.size(); i++){
  //   IRL2D::Vec x = centroids[i];
  //   projected_points[i] = IRL2D::projectToImplicitSurface(x, centroids,
  //   normals, mesh.dx());
  //   // std::cout << projected_points[i] << std::endl;
  // }
  // for (int i = mesh.imin(); i <= mesh.imax(); i++){
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); j++){
  //     if (mixed(i,j) == true){
  //       IRL2D::Parabola PU_parabola =
  //       IRL2D::getPU_interface(plic_center(i,j), centroids, normals,
  //       mesh.dx());
  //       // csvfile << PU_parabola.datum()[0] << "," << PU_parabola.datum()[1]
  //       << ","
  //       //         << PU_parabola.frame()[0][0] << "," <<
  //       PU_parabola.frame()[0][1] << ","
  //       //         << PU_parabola.frame()[1][0] << "," <<
  //       PU_parabola.frame()[1][1] << "\n";
  //       // std::cout << PU_parabola.frame() << std::endl;
  //     }
  //   }
  // }

  // std::cout << "x = [ ";
  // for (int i = 0; i < projected_points.size(); i++){
  //   std::cout << projected_points[i][0] << " ";
  // }
  // std::cout << "];" << std::endl;

  // std::cout << "y = [ ";
  // for (int i = 0; i < projected_points.size(); i++){
  //   std::cout << projected_points[i][1] << " ";
  // }
  // std::cout << "];" << std::endl;

  // parabolic interface using PU
  const double kernel_size = 2.5 * mesh.dx();
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (mixed(i, j) == true) {
        bool usePlane = false;
        IRL2D::Parabola PU_parabola = IRL2D::getPU_interface(
            plic_center(i, j), centroids, normals, kernel_size, usePlane);

        if (usePlane) {
          // std::cout << "(i,j) = " << i << " , " << j << std::endl;
          // // outputting plic data
          // std::string dir =
          //     "/home/parinht2/Desktop/partition_of_unity/2d_interface";
          // std::string filepath = dir + "/nan_data.csv";
          // std::ofstream csvfile(filepath);
          // csvfile << "ax,ay,bx,by,tx,ty,nx,ny\n";
          // for (int ii = mesh.imin(); ii <= mesh.imax(); ii++) {
          //   for (int jj = mesh.jmin(); jj <= mesh.jmax(); jj++) {
          //     if (mixed(ii, jj) == true) {
          //       csvfile << end_points(ii, jj).first[0] << ","
          //               << end_points(ii, jj).first[1] << ","
          //               << end_points(ii, jj).second[0] << ","
          //               << end_points(ii, jj).second[1] << ","
          //               << plic(ii, jj).frame()[0][0] << ","
          //               << plic(ii, jj).frame()[0][1] << ","
          //               << plic(ii, jj).frame()[1][0] << ","
          //               << plic(ii, jj).frame()[1][1] << "\n";
          //     }
          //   }
          // }
          continue;  // use LVIRA
        }

        // return LVIRA with PU
        // std::vector<IRL2D::Vec> cen = {plic_center(i,j)};
        // std::vector<IRL2D::Vec> nor = {plic(i,j).frame()[1]};
        // IRL2D::Parabola PU_parabola =
        // IRL2D::getPU_interface(plic_center(i,j), cen, nor, kernel_size);

        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        double vfrac = (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
        (*a_interface)(i, j) =
            IRL2D::MatchToVolumeFraction(cell, PU_parabola, vfrac);
        // (*a_interface)(i, j) = PU_parabola;
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// partition of unity with plic ------------------------------------------------

void PUPLIC::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                           a_interface);
  Data<IRL2D::Parabola> plic = *a_interface;  // previous plic
  const BasicMesh& mesh = a_U.getMesh();
  std::vector<IRL2D::Vec> centroids;
  std::vector<IRL2D::Vec> normals;
  Data<bool> mixed(&mesh);
  Data<IRL2D::Vec> plic_center(&mesh);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      mixed(i, j) = false;
      const double lvf = (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (lvf >= IRL::global_constants::VF_LOW &&
          lvf <= IRL::global_constants::VF_HIGH) {
        mixed(i, j) = true;
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        IRL2D::BezierList clipped_plic =
            IRL2D::ParabolaClip(cell, plic(i, j), true);
        IRL2D::Vec center =
            (clipped_plic[0].first + clipped_plic[1].first) / 2.0;
        centroids.push_back(center);
        plic_center(i, j) = center;
        normals.push_back(plic(i, j).frame()[1]);
      }
    }
  }

  // plic approximation using PU
  const double kernel_size = 2.5 * mesh.dx();
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (mixed(i, j) == true) {
        bool usePlane = false;
        IRL2D::Parabola PU_interface = plic(i, j);
        IRL2D::Vec x_proj = IRL2D::projectToImplicitSurface(
            plic_center(i, j), centroids, normals, kernel_size, usePlane);
        IRL2D::ImplicitSurface IS(centroids, normals, kernel_size);
        double Fx = IS.Fx(x_proj);
        double Fy = IS.Fy(x_proj);
        IRL2D::Vec normal = {Fx, Fy};
        normal.normalize();
        IRL2D::Vec tangent = {normal.y(), -normal.x()};
        PU_interface.frame() = {tangent, normal};

        if (usePlane) {
          continue;  // use LVIRA
        }

        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(mesh.x(i), mesh.y(j)),
            IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
        double vfrac = (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
        (*a_interface)(i, j) =
            IRL2D::MatchToVolumeFraction(cell, PU_interface, vfrac);
      }
    }
  }

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

// ---------------------------------------------------------------------------------------------------------

void correctInterfaceBorders(Data<IRL2D::Parabola>* a_interface) {
  const BasicMesh& mesh = (*a_interface).getMesh();
  // Fix distances in reconstruction for periodic boundary

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      IRL2D::Vec& datum = (*a_interface)(i, j).datum();
      datum[0] -= mesh.lx();
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      IRL2D::Vec& datum = (*a_interface)(i, j).datum();
      datum[0] += mesh.lx();
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      IRL2D::Vec& datum = (*a_interface)(i, j).datum();
      datum[1] -= mesh.ly();
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      IRL2D::Vec& datum = (*a_interface)(i, j).datum();
      datum[1] += mesh.ly();
    }
  }
}
