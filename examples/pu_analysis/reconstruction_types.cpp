// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include "examples/pu_analysis/reconstruction_types.h"

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
#include "examples/pu_analysis/basic_mesh.h"
#include "examples/pu_analysis/data.h"
#include "examples/pu_analysis/vof_advection.h"

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
  } else if (a_reconstruction_method == "LVIRAQ") {
    LVIRAQ::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "Taubin") {
    Taubin::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "PU") {
    PU::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                          a_interface);
  } else if (a_reconstruction_method == "PUIterate") {
    PUIterate::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U,
                                 a_V, a_interface);
  } else if (a_reconstruction_method == "PUPLIC") {
    PUPLIC::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: ELVIRA, LVIRA, Jibben, Particle, "
                 "LVIRAQ, Taubin, PU, PUIterate. \n";
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

// least squares fitting of circles
// -----------------------------------------------------------------------

struct InterfaceData {
  bool mixed = false;
  int xIndex, yIndex;
  double vf;
  IRL2D::Vec a, b, center;  // start, end, and midpoint
  IRL2D::BezierList rectangle;
};

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
            // vfw.push_back(vf_weight);
            vfw.push_back(1.0);
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
      }
    }
  }

  // parabolic interface using PU
  const double kernel_size = 2.5 * mesh.dx();
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (mixed(i, j) == true) {
        bool usePlane = false;
        IRL2D::Parabola PU_parabola = IRL2D::getPU_interface(
            plic_center(i, j), centroids, normals, kernel_size, usePlane);

        if (usePlane) {
          continue;  // use LVIRA
        }

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

// iterating partition of unity -
// ----------------------------------------------------------------------------------------
void PUIterate::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                                  const Data<IRL2D::Moments>& a_gas_moments,
                                  const double a_dt, const Data<double>& a_U,
                                  const Data<double>& a_V,
                                  Data<IRL2D::Parabola>* a_interface) {
  Data<IRL2D::Parabola> pu_interface = *a_interface;
  const BasicMesh& mesh = a_U.getMesh();
  std::vector<IRL2D::Vec> centroids, datums, normals, tangents;
  std::vector<double> a_list;
  Data<bool> mixed(&mesh);
  Data<IRL2D::Vec> interface_centroid(&mesh);

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
        IRL2D::BezierList clipped_interface =
            IRL2D::ParabolaClip(cell, pu_interface(i, j), true);
        IRL2D::Vec start = clipped_interface[0].first;
        IRL2D::Vec end = clipped_interface[1].first;
        IRL2D::Vec centroid = IRL2D::centroidParabolaSegmentAnalytic(
            pu_interface(i, j).datum(), pu_interface(i, j).frame()[0],
            pu_interface(i, j).frame()[1], pu_interface(i, j).coeff(), start,
            end);
        centroids.push_back(centroid);
        interface_centroid(i, j) = centroid;

        datums.push_back(pu_interface(i, j).datum());
        normals.push_back(pu_interface(i, j).frame()[1]);
        tangents.push_back(pu_interface(i, j).frame()[0]);
        a_list.push_back(pu_interface(i, j).coeff());
      }
    }
  }

  // parabolic interface using PU
  const double kernel_size = 2.5 * mesh.dx();
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      if (mixed(i, j) == true) {
        bool usePlane = false;
        IRL2D::Parabola PU_parabola = IRL2D::getPU_parabola_interface(
            interface_centroid(i, j), centroids, datums, tangents, normals,
            a_list, kernel_size, usePlane);

        if (usePlane) {
          continue;  // use previous interface
        }

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

// iterating partition of unity using PLIC
// ----------------------------------------------------------------------------------------
void PUPLIC::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                               const Data<IRL2D::Moments>& a_gas_moments,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V,
                               Data<IRL2D::Parabola>* a_interface) {
  // LVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
  //                          a_interface);
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
