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
  } else if (a_reconstruction_method == "LVIRAQ") {
    LVIRAQ::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "MOF1") {
    MOF1::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "MOF2") {
    MOF2::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                              a_interface);
  } else if (a_reconstruction_method == "MOF2AL"){
    MOF2AL::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
      a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: ELVIRA, LVIRA, LVIRAQ, MOF1, MOF2, MOF2AL. \n";
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

// Functor for MOF
struct MOF1Functor{
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum{
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

  MOF1Functor(const IRL2D::BezierList& cell, const IRL2D::Vec& liq_centroid_star,
             const IRL2D::Vec& gas_centroid_star, const double liq_f_star)
    : m_cell(cell), m_liq_centroid_star(liq_centroid_star),
      m_gas_centroid_star(gas_centroid_star), m_liq_f_star(liq_f_star) {
      m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cell));
  }
  
  void setframe(const IRL2D::Parabola& guess_plane){
    m_coeff = 0.0;
    m_frame = guess_plane.frame();
    m_datum = IRL2D::ComputeMoments(m_cell).m1() / IRL2D::ComputeArea(m_cell);
    //m_datum = guess_plane.datum();
  }

  const IRL2D::Parabola getplane(const Eigen::VectorXd& x) const {
    const auto rotation = IRL2D::ReferenceFrame(x(0));
    const auto new_frame = IRL2D::ReferenceFrame(rotation * m_frame[0], rotation * m_frame[1]);
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
                            Data<IRL2D::Parabola>* a_interface){
  
  // initial guess
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);
  
  const BasicMesh& mesh = a_U.getMesh();

  #ifdef USE_MPI
    const double cell_volume = mesh.cell_volume();
    int nmixed_global = 0;
    for (int i = mesh.imin(); i <= mesh.imax(); ++i){
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j){
        const double liquid_volume_fraction = a_liquid_moments(i,j).m0() / cell_volume;
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH){
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
    IRL::serializeAndPack(dummy_par , &dummy_buffer);
    const int size_parabola = dummy_buffer.size();

    int nmixed_local = std::max(nmixed_global / size , 1);
    std::vector<int> proc_offset(size + 1);
    proc_offset[0] = 0;
    for (int r = 0; r < size; r++){
      proc_offset[r + 1] = proc_offset[r] + nmixed_local;
    }
    proc_offset[size] = nmixed_global;
    for (int r = 1; r < size + 1; r++){
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
      const double liquid_volume_fraction = a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {

  #ifdef USE_MPI
        if (count >= proc_offset[rank] && proc_offset[rank + 1]) {
  #endif
        
        IRL2D::Vec x0 = IRL2D::Vec(mesh.x(i) , mesh.y(j));
        IRL2D::Vec x1 = IRL2D::Vec(mesh.x(i+1) , mesh.y(j+1));
        IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(x0, x1);
        
        // matching volume fraction and centroid
        IRL2D::Vec liq_centroid_star = a_liquid_moments(i,j).m1() / a_liquid_moments(i,j).m0();
        IRL2D::Vec gas_centroid_star = a_gas_moments(i,j).m1() / a_gas_moments(i,j).m0();
        MOF1Functor myMOFFunctor(rectangle, liq_centroid_star, gas_centroid_star, liquid_volume_fraction);
        myMOFFunctor.setframe((*a_interface)(i, j));
        Eigen::NumericalDiff<MOF1Functor> numericalDiffMyFunctor(myMOFFunctor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF1Functor>, double> lm(numericalDiffMyFunctor);
        Eigen::VectorXd x(1);
        x.setZero();
        Eigen::LevenbergMarquardtSpace::Status status =
              lm.minimizeInit(x);
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
    for (int r = 0; r < size; r++){
      proc_count[r] = size_parabola * (proc_offset[r + 1] - proc_offset[r]);
      proc_offset[r] = size_parabola * proc_offset[r];
    }

    interface_global.resize(size_parabola * nmixed_global);
    MPI_Allgatherv(interface_local.data(), size_parabola * nmixed_local, MPI_BYTE,
                   interface_global.data(), proc_count.data(), proc_offset.data(),
                   MPI_BYTE, MPI_COMM_WORLD);

    for (int i = mesh.imin(); i <= mesh.imax(); ++i){
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j){
        const double liquid_volume_fraction = a_liquid_moments(i,j).m0() / cell_volume;
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH){
              IRL2D::Parabola parabola;
              IRL::unpackAndStore(&parabola, &interface_global);
              (*a_interface)(i,j) = parabola;
        }
      }
    }
                            
  #endif

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);

}

// Functor for MOF2
struct MOF2Functor{
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum{
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // variables
  const int m_inputs, m_values;
  const IRL2D::BezierList& m_cell;
  IRL2D::Moments& m_liq_moments;
  IRL2D::Moments& m_gas_moments;
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
              IRL2D::Moments& liq_moments, IRL2D::Moments& gas_moments)
    : m_inputs(inputs),
      m_values(values),
      m_cell(cell),
      m_liq_moments(liq_moments),
      m_gas_moments(gas_moments),
      m_liq_f_star(liq_moments.m0() / IRL2D::ComputeArea(cell)),
      m_liq_centroid_star(liq_moments.m1() / liq_moments.m0()),
      m_gas_centroid_star(gas_moments.m1() / gas_moments.m0()) {  
      //m_datum = ( (1-m_liq_f_star) * m_liq_centroid_star + m_liq_f_star * m_gas_centroid_star );
      m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cell));
      RecenterMoments(&m_liq_moments, m_liq_centroid_star);
      RecenterMoments(&m_gas_moments, m_gas_centroid_star);
      m_liq_M2_star = m_liq_moments.m2();
      m_gas_M2_star = m_gas_moments.m2();
  }
  
  void setframe(const IRL2D::Parabola& guess_parabola){
    m_coeff = guess_parabola.coeff();
    m_frame = guess_parabola.frame();
    m_datum = guess_parabola.datum();
    //m_frame = IRL2D::Mat(IRL2D::Vec(1.0,0.0), IRL2D::Vec(0.0,1.0));
    //m_datum = IRL2D::ComputeMoments(m_cell).m1() / IRL2D::ComputeArea(m_cell);

    // check for curvature
    const double maxkdx = 4.0;
    const double kdx = 2.0 * m_coeff * m_length_scale;
    if (std::abs(kdx) > maxkdx){
      m_coeff = 0.0; // plane
    }
  }

  const IRL2D::Parabola getparabola(const Eigen::VectorXd& x) const {
    const auto rotation = IRL2D::ReferenceFrame(x(0));
    const auto new_frame = IRL2D::ReferenceFrame(rotation * m_frame[0], rotation * m_frame[1]);
    const double new_coeff = m_coeff + x(1) / m_length_scale;
    const auto parabola = IRL2D::Parabola(m_datum, new_frame, new_coeff);
    return IRL2D::MatchToVolumeFraction(m_cell, parabola, m_liq_f_star);
  }

  mutable int iteration = 0; // counting iterations for convergence
  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto parabola = this->getparabola(x);
    IRL2D::Moments liq_mom = IRL2D::ComputeMoments(m_cell, parabola);
    IRL2D::Moments gas_mom = IRL2D::ComputeMoments(m_cell) - liq_mom;
    IRL2D::Vec liq_centroid_h = liq_mom.m1() / liq_mom.m0();
    IRL2D::Vec gas_centroid_h = gas_mom.m1() / gas_mom.m0();
    RecenterMoments(&liq_mom, m_liq_centroid_star);
    RecenterMoments(&gas_mom, m_gas_centroid_star);
    IRL2D::Mat liq_M2_h = liq_mom.m2();
    IRL2D::Mat gas_M2_h = gas_mom.m2();
    fvec(0) = (m_liq_centroid_star[0] - liq_centroid_h[0]) / m_length_scale;
    fvec(1) = (m_liq_centroid_star[1] - liq_centroid_h[1]) / m_length_scale;
    fvec(2) = (m_gas_centroid_star[0] - gas_centroid_h[0]) / m_length_scale;
    fvec(3) = (m_gas_centroid_star[1] - gas_centroid_h[1]) / m_length_scale;
    fvec(4) = (m_liq_M2_star[0][0] - liq_M2_h[0][0]) / (liq_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(5) = (m_liq_M2_star[1][0] - liq_M2_h[1][0]) / (liq_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(6) = (m_liq_M2_star[1][1] - liq_M2_h[1][1]) / (liq_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(7) = (m_gas_M2_star[0][0] - gas_M2_h[0][0]) / (gas_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(8) = (m_gas_M2_star[1][0] - gas_M2_h[1][0]) / (gas_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(9) = (m_gas_M2_star[1][1] - gas_M2_h[1][1]) / (gas_mom.m0()*std::pow(m_length_scale, 2.0));

    //scaling based on Shashkov's paper
    // fvec(0) = std::pow((liq_centroid_h[0] - m_liq_centroid_star[0]), 2.0) / m_liq_moments.m0();
    // fvec(1) = std::pow((liq_centroid_h[1] - m_liq_centroid_star[1]), 2.0) / m_liq_moments.m0();
    // fvec(2) = std::pow((gas_centroid_h[0] - m_gas_centroid_star[0]), 2.0) / m_gas_moments.m0();
    // fvec(3) = std::pow((gas_centroid_h[1] - m_gas_centroid_star[1]), 2.0) / m_gas_moments.m0();
    // fvec(4) = std::pow((liq_M2_h[0][0] - m_liq_M2_star[0][0]), 2.0) / (std::pow(m_liq_M2_star[0][0],2.0)) ; 
    // fvec(5) = std::pow((liq_M2_h[1][1] - m_liq_M2_star[1][1]), 2.0) / (std::pow(m_liq_M2_star[1][1],2.0));
    // fvec(6) = std::pow((liq_M2_h[1][0] - m_liq_M2_star[1][0]), 2.0) / (std::pow(m_liq_M2_star[0][0],2.0) + std::pow(m_liq_M2_star[1][1],2.0));
    // fvec(7) = std::pow((gas_M2_h[0][0] - m_gas_M2_star[0][0]), 2.0) / (std::pow(m_gas_M2_star[0][0],2.0));
    // fvec(8) = std::pow((gas_M2_h[1][1] - m_gas_M2_star[1][1]), 2.0) / (std::pow(m_gas_M2_star[1][1],2.0));
    // fvec(9) = std::pow((gas_M2_h[1][0] - m_gas_M2_star[1][0]), 2.0) / (std::pow(m_gas_M2_star[0][0],2.0) + std::pow(m_gas_M2_star[1][1],2.0));

    // Penalty to prevent kappa * dx > 4
    const double mu = 50.0;
    const double kdx = 2.0 * m_coeff * m_length_scale + x(1);
    const double maxkdx = 4.0;
    fvec(10) = mu * std::max(0.0, std::abs(kdx) - maxkdx);
    
    std::cout << "Iteration: " << iteration << ", ||fvec|| = " << fvec.norm() <<std::endl;

    if (iteration == 0){
      std::cout << fvec.transpose() << std::endl;
      std::cout << fvec.norm() << std::endl;
    }

    iteration++;
    
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
                            Data<IRL2D::Parabola>* a_interface){
  
  // initial guess
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);
  
  const BasicMesh& mesh = a_U.getMesh();

  // #ifdef USE_MPI
  //   const double cell_volume = mesh.cell_volume();
  //   int nmixed_global = 0;
  //   for (int i = mesh.imin(); i <= mesh.imax(); ++i){
  //     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j){
  //       const double liquid_volume_fraction = a_liquid_moments(i,j).m0() / cell_volume;
  //       if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
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

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction = a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH && i == 23 && j == 51) {

  // #ifdef USE_MPI
  //       if (count >= proc_offset[rank] && proc_offset[rank + 1]) {
  // #endif
        
        IRL2D::Vec x0 = IRL2D::Vec(mesh.x(i) , mesh.y(j));
        IRL2D::Vec x1 = IRL2D::Vec(mesh.x(i+1) , mesh.y(j+1));
        IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(x0, x1);
        
        auto liq_moments = a_liquid_moments(i,j);
        auto gas_moments = a_gas_moments(i,j);
        MOF2Functor myMOFFunctor(4, 16, rectangle, liq_moments, gas_moments);
        myMOFFunctor.setframe((*a_interface)(i, j));
        Eigen::NumericalDiff<MOF2Functor> numericalDiffMyFunctor(myMOFFunctor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF2Functor>, double> lm(numericalDiffMyFunctor);
        Eigen::VectorXd x(4);
        x.setZero();
        lm.parameters.ftol = 1e-14;
        lm.parameters.xtol = 1e-14;
        //  lm.parameters.maxfev = 500;
        Eigen::LevenbergMarquardtSpace::Status status =
              lm.minimizeInit(x);
        do {
          status = lm.minimizeOneStep(x);
        } while (status == Eigen::LevenbergMarquardtSpace::Running);
        const auto parabola = IRL2D::MatchToVolumeFractionBisection(
            rectangle, myMOFFunctor.getparabola(x), liquid_volume_fraction, 500);
  // #ifdef USE_MPI
  //       IRL::serializeAndPack(parabola, &interface_local);
  //       }
  //       count++;
  // #else       
        (*a_interface)(i, j) = parabola;
  //#endif       
      }
    }
  }

  // #ifdef USE_MPI
  //   std::vector<int> proc_count(size);
  //   for (int r = 0; r < size; r++){
  //     proc_count[r] = size_parabola * (proc_offset[r + 1] - proc_offset[r]);
  //     proc_offset[r] = size_parabola * proc_offset[r];
  //   }

  //   interface_global.resize(size_parabola * nmixed_global);
  //   MPI_Allgatherv(interface_local.data(), size_parabola * nmixed_local, MPI_BYTE,
  //                  interface_global.data(), proc_count.data(), proc_offset.data(),
  //                  MPI_BYTE, MPI_COMM_WORLD);

  //   for (int i = mesh.imin(); i <= mesh.imax(); ++i){
  //     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j){
  //       const double liquid_volume_fraction = a_liquid_moments(i,j).m0() / cell_volume;
  //       if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
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

// Functor for MOF2
struct MOF2AugmentedLagrangianFunctor{
  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  enum{
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
  IRL2D::Vec m_cell_centroid;
  double m_coeff;
  double m_length_scale;
  std::vector<double> m_lambda;
  std::vector<double> m_mu;

  // constructor
  MOF2AugmentedLagrangianFunctor(int inputs, int values, const IRL2D::BezierList& cell,
                                 const IRL2D::Moments& liq_moments, const IRL2D::Moments& gas_moments,
                                const std::vector<double>& lambda, const std::vector<double>& mu)
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
      //m_datum = ( (1-m_liq_f_star) * m_liq_centroid_star + m_liq_f_star * m_gas_centroid_star );
      m_length_scale = std::sqrt(IRL2D::ComputeArea(m_cell));
      RecenterMoments(&m_liq_moments, m_liq_centroid_star);
      RecenterMoments(&m_gas_moments, m_gas_centroid_star);
      m_liq_M2_star = m_liq_moments.m2();
      m_gas_M2_star = m_gas_moments.m2();
      m_cell_centroid = IRL2D::ComputeMoments(m_cell).m1() / IRL2D::ComputeArea(m_cell);
  }
  
  void setframe(const IRL2D::Parabola& guess_parabola){
    m_coeff = guess_parabola.coeff();
    m_frame = guess_parabola.frame();
    m_datum = guess_parabola.datum();
    //m_frame = IRL2D::Mat(IRL2D::Vec(1.0,0.0), IRL2D::Vec(0.0,1.0));
    //m_datum = IRL2D::ComputeMoments(m_cell).m1() / IRL2D::ComputeArea(m_cell);

    // check for curvature
    const double maxkdx = 4.0;
    const double kdx = 2.0 * m_coeff * m_length_scale;
    if (std::abs(kdx) > maxkdx){
      m_coeff = 0.0; // plane
    }
  }

  // x(0): theta    x(1): alpha   x(2): datum[0]    x(3): datum[1]
  const IRL2D::Parabola getparabola(const Eigen::VectorXd& x) const {
    const auto rotation = IRL2D::ReferenceFrame(x(0));
    const auto new_frame = IRL2D::ReferenceFrame(rotation * m_frame[0], rotation * m_frame[1]);
    const double new_coeff = m_coeff + x(1) / m_length_scale;
    const auto new_datum = IRL2D::Vec(m_datum[0] + x(2)*m_length_scale , m_datum[1] + x(3)*m_length_scale);
    const auto parabola = IRL2D::Parabola(new_datum, new_frame, new_coeff);
    return parabola;
  }

  double getVolfrac (const Eigen::VectorXd& x) const{
    const auto parabola = this->getparabola(x);
    const auto moments = IRL2D::ComputeMoments(m_cell, parabola);
    return (moments.m0() / IRL2D::ComputeArea(m_cell));
  }

  // std::vector<double> getIneqFunction(const Eigen::VectorXd& x) const{
  //   std::vector<double> IneqFunctions;
  //   const double dx = m_length_scale;
  //   const double x0 = (m_cell[0].first)[0], x1 = (m_cell[1].first)[0], 
  //                y0 = (m_cell[0].first)[1], y1 = (m_cell[2].first)[1];
  //   std::vector<double> bounds = { (x0 - (m_datum[0] + x(2))) , ((m_datum[0] + x(2)) - x1) , 
  //                                  (y0 - (m_datum[1] + x(3))) , ((m_datum[1] + x(3)) - y1) };
  //   for (int bound = 0; bound < bounds.size(); bound++){
  //     IneqFunctions[bound] = std::max(0.0 , bounds[bound]);
  //   }             
  //   return IneqFunctions;
  // }

  mutable int iteration = 0; // counting iterations for convergence
  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto parabola = this->getparabola(x);
    IRL2D::Moments liq_mom = IRL2D::ComputeMoments(m_cell, parabola);
    IRL2D::Moments gas_mom = IRL2D::ComputeMoments(m_cell) - liq_mom;
    IRL2D::Vec liq_centroid_h = liq_mom.m1() / liq_mom.m0();
    IRL2D::Vec gas_centroid_h = gas_mom.m1() / gas_mom.m0();
    RecenterMoments(&liq_mom, m_liq_centroid_star);
    RecenterMoments(&gas_mom, m_gas_centroid_star);
    IRL2D::Mat liq_M2_h = liq_mom.m2();
    IRL2D::Mat gas_M2_h = gas_mom.m2();

    // centroids
    fvec(0) = (m_liq_centroid_star[0] - liq_centroid_h[0]) / m_length_scale;
    fvec(1) = (m_liq_centroid_star[1] - liq_centroid_h[1]) / m_length_scale;
    fvec(2) = (m_gas_centroid_star[0] - gas_centroid_h[0]) / m_length_scale;
    fvec(3) = (m_gas_centroid_star[1] - gas_centroid_h[1]) / m_length_scale;

    // second moments
    fvec(4) = (m_liq_M2_star[0][0] - liq_M2_h[0][0]) / (liq_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(5) = (m_liq_M2_star[1][0] - liq_M2_h[1][0]) / (liq_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(6) = (m_liq_M2_star[1][1] - liq_M2_h[1][1]) / (liq_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(7) = (m_gas_M2_star[0][0] - gas_M2_h[0][0]) / (gas_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(8) = (m_gas_M2_star[1][0] - gas_M2_h[1][0]) / (gas_mom.m0()*std::pow(m_length_scale, 2.0));
    fvec(9) = (m_gas_M2_star[1][1] - gas_M2_h[1][1]) / (gas_mom.m0()*std::pow(m_length_scale, 2.0));

    // volume fraction constraint
    double liq_f_h = liq_mom.m0() / IRL2D::ComputeArea(m_cell);
    double f_constraint = liq_f_h - m_liq_f_star;
    fvec(10) = std::sqrt(m_mu[0]) * f_constraint + 0.5 / (std::sqrt(m_mu[0])) * m_lambda[0];

    // std::cout << iteration++ << ". residual = [" << fvec(0) << " " << fvec(1) << " " << fvec(2) << " " << fvec(3) << " " << fvec(4) << " "
    //           << fvec(5) << " " << fvec(6) << " " << fvec(7) << " " << fvec(8) << " " << fvec(9) << " " << fvec(10) << "] " 
    //           << "mu = " << m_mu[0] << " lambda = " << m_lambda[0] << " norm = " << fvec.norm() << std::endl;
    // confining the datum within the cell
    // const double dx = m_length_scale;
    // const double x0 = (m_cell[0].first)[0], x1 = (m_cell[1].first)[0], 
    //              y0 = (m_cell[0].first)[1], y1 = (m_cell[2].first)[1];
    //std::vector<double> bounds = { -(dx/2.0 + x(2)) , (x(2) - dx/2.0) , -(dx/2.0 + x(3)) , (x(3) - dx/2.0) };
    // std::vector<double> bounds = { (x0 - (m_datum[0] + x(2))) , ((m_datum[0] + x(2)) - x1) , 
    //                                (y0 - (m_datum[1] + x(3))) , ((m_datum[1] + x(3)) - y1) };
    // for (int bound = 0; bound < bounds.size(); bound++){
    //   fvec(10 + (bound + 1)) = std::sqrt(m_mu[bound + 1]) * std::max(0.0 , bounds[bound]) + 
    //                            0.5 / (std::sqrt(m_mu[bound + 1])) * m_lambda[bound + 1];
    // }

    // std::vector<double> IneqFunctions = this->getIneqFunction(x);
    // for (int bound = 0; bound < IneqFunctions.size(); bound++){
    //   fvec(10 + (bound + 1)) = std::sqrt(m_mu[bound + 1]) * IneqFunctions[bound] + 
    //                            0.5 / (std::sqrt(m_mu[bound + 1])) * m_lambda[bound + 1];
    // }

    // // curvature inequality
    // const double kdx = 2.0 * m_coeff * m_length_scale + x(1);
    // const double maxkdx = 4.0;
    // fvec(15) = std::sqrt(m_mu[5]) * std::max(0.0, std::abs(kdx) - maxkdx); + 
    //            0.5 / (std::sqrt(m_mu[5])) * m_lambda[5];
  
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    this->errorvec(x, fvec);
    return 0;
  }

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

}; 


void MOF2AL::getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                            const Data<IRL2D::Moments>& a_gas_moments,
                            const double a_dt, const Data<double>& a_U,
                            const Data<double>& a_V,
                            Data<IRL2D::Parabola>* a_interface){
  
  // initial guess
  ELVIRA::getReconstruction(a_liquid_moments, a_gas_moments, a_dt, a_U, a_V,
                            a_interface);
  
  const BasicMesh& mesh = a_U.getMesh();

  #ifdef USE_MPI
    const double cell_volume = mesh.cell_volume();
    int nmixed_global = 0;
    for (int i = mesh.imin(); i <= mesh.imax(); ++i){
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j){
        const double liquid_volume_fraction = a_liquid_moments(i,j).m0() / cell_volume;
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH){
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
    IRL::serializeAndPack(dummy_par , &dummy_buffer);
    const int size_parabola = dummy_buffer.size();

    int nmixed_local = std::max(nmixed_global / size , 1);
    std::vector<int> proc_offset(size + 1);
    proc_offset[0] = 0;
    for (int r = 0; r < size; r++){
      proc_offset[r + 1] = proc_offset[r] + nmixed_local;
    }
    proc_offset[size] = nmixed_global;
    for (int r = 1; r < size + 1; r++){
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
      const double liquid_volume_fraction = a_liquid_moments(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH){ // && i == 23 && j == 51) {

  #ifdef USE_MPI
        if (count >= proc_offset[rank] && proc_offset[rank + 1]) {
  #endif
        
        // current cell
        IRL2D::Vec x0 = IRL2D::Vec(mesh.x(i) , mesh.y(j));
        IRL2D::Vec x1 = IRL2D::Vec(mesh.x(i+1) , mesh.y(j+1));
        IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(x0, x1);
        
        // moments
        const auto liq_moments = a_liquid_moments(i,j);
        const auto gas_moments = a_gas_moments(i,j);
        
        // parameters for augmented lagrangian
        std::vector<double> lambda(1, 0.0); // CHANGE 1 -> 6 with inequalities
        std::vector<double> mu(1, 1.0);
        double mu_max = 100.0;
        double tol = 1e-12;
        double f_constraint_current = 1;
        double f_constraint_prev;
        int iter = 0;
        int max_iter = 1000;

        // initial guess for first iteration
        auto guess_interface = (*a_interface)(i,j);
        IRL2D::Parabola parabola;
        
        // for LM solver
        int num_inputs = 4;
        int num_eq = 11;

        Eigen::VectorXd x(num_inputs);
        x.setZero();

        while (std::abs(f_constraint_current) > tol){

          iter++;

          //std::cout << "-----------------Augmented Lagrangian Iteration: " << iter << "-------------------" << std::endl;

          // minimization variables
          // Eigen::VectorXd x(num_inputs);
          // x.setZero();

          // setting up LM solver
          MOF2AugmentedLagrangianFunctor myMOFFunctor(num_inputs, num_eq, rectangle, 
                                                      liq_moments, gas_moments,
                                                      lambda, mu);
          myMOFFunctor.setframe(guess_interface); // guess for interface
          Eigen::NumericalDiff<MOF2AugmentedLagrangianFunctor> numericalDiffMyFunctor(myMOFFunctor);
          Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOF2AugmentedLagrangianFunctor>, double> lm(numericalDiffMyFunctor);

          //LM solver parameters
          lm.parameters.ftol = 1e-12;
          lm.parameters.xtol = 1e-12;
          
          auto x_prev = x; // storing prev iter for AL multiplier updates
          lm.minimize(x); // minimization

          // compute volume fraction with new and prev x 
          f_constraint_current = myMOFFunctor.getVolfrac(x) - liquid_volume_fraction;
          f_constraint_prev = myMOFFunctor.getVolfrac(x_prev) - liquid_volume_fraction;
          
          // update Lagrangian multiplier for constraint
          lambda[0] += 2.0 * mu[0] * f_constraint_current;

          // update penalty parameter for constraint
          if (std::abs(f_constraint_current) < 0.25 * std::abs(f_constraint_prev)){
            //mu[0] = mu[0];
            mu[0] = std::min(mu[0], mu_max);
          } else {
            //mu[0] *= 2.0;
            mu[0] = std::min(2*mu[0], mu_max);
          }

          // updating parameters for inequalities
          // std::vector<double> IneqFunctions = myMOFFunctor.getIneqFunction(x);
          // for (int ineq = 0 ; ineq < IneqFunctions.size(); ineq++){
          //   // Lagrangian multiplier update
          //   lambda[ineq + 1] += 2.0 * mu[ineq + 1] * IneqFunctions[ineq];

          //   // penalty term update
          // }

          // updating interface guess for next AL iteration
          parabola = myMOFFunctor.getparabola(x); // new parabola
          //guess_interface = parabola;
          
          // break if it reaches max iter
          if (iter == max_iter){
            break;
          }

        }
        //std::cout << "AL iterations: " << iter << std::endl;
        //(*a_interface)(i,j) = parabola; // final interface

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
    for (int r = 0; r < size; r++){
      proc_count[r] = size_parabola * (proc_offset[r + 1] - proc_offset[r]);
      proc_offset[r] = size_parabola * proc_offset[r];
    }

    interface_global.resize(size_parabola * nmixed_global);
    MPI_Allgatherv(interface_local.data(), size_parabola * nmixed_local, MPI_BYTE,
                   interface_global.data(), proc_count.data(), proc_offset.data(),
                   MPI_BYTE, MPI_COMM_WORLD);

    for (int i = mesh.imin(); i <= mesh.imax(); ++i){
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j){
        const double liquid_volume_fraction = a_liquid_moments(i,j).m0() / cell_volume;
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH){
              IRL2D::Parabola parabola;
              IRL::unpackAndStore(&parabola, &interface_global);
              (*a_interface)(i,j) = parabola;
        }
      }
    }
                            
  #endif

  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);

}

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
