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
  if (a_reconstruction_method == "ELVIRA") {
    ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "LVIRAQ") {
    LVIRAQ::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: ELVIRA, LVIRAQ. \n";
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
        const auto planar_separator = reconstructionWithELVIRA2D(neighborhood);
        const auto normal = planar_separator[0].normal();
        const auto frame = IRL2D::ReferenceFrame{{normal[1], -normal[0]},
                                                 {normal[0], normal[1]}};
        const auto datum = IRL2D::Vec{mesh.xm(i), mesh.ym(j)};
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
        IRL2D::ReferenceFrame({rotation * m_frame[0], rotation * m_frame[1]});
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
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);

  // Now fit parabola LVIRA-style
  const BasicMesh& mesh = a_U.getMesh();

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double vfrac = a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (vfrac < IRL::global_constants::VF_LOW) {
        (*a_interface)(i, j).markAsAlwaysBelow();
      } else if (vfrac > IRL::global_constants::VF_HIGH) {
        (*a_interface)(i, j).markAsAlwaysAbove();
      } else {
        // Fill stencil of moments
        std::array<std::array<double, 3>, 3> vfracs;
        std::array<std::array<IRL2D::BezierList, 3>, 3> cells;
        for (int jj = 0; jj < 3; ++jj) {
          for (int ii = 0; ii < 3; ++ii) {
            vfracs[ii][jj] = a_liquid_moments(i + ii - 1, j + jj - 1).m0() /
                             mesh.cell_volume();
            const auto x0 = IRL2D::Vec(mesh.x(i + ii - 1), mesh.y(j + jj - 1));
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
        Eigen::VectorXd x(2);
        x.setZero();
        Eigen::LevenbergMarquardtSpace::Status status =
            LVIRAQ_LM.minimizeInit(x);
        do {
          status = LVIRAQ_LM.minimizeOneStep(x);
        } while (status == Eigen::LevenbergMarquardtSpace::Running);
        (*a_interface)(i, j) = myLVIRAQFunctor.getparabola(x);
      }
    }
  }
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

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
