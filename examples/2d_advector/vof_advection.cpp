// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <mpi.h>
#include <iostream>

#include "examples/2d_advector/edge_transport.h"
#include "examples/2d_advector/functions.h"
#include "examples/2d_advector/irl2d.h"

#include "examples/2d_advector/reconstruction_types.h"
#include "examples/2d_advector/rotation_2d.h"
#include "examples/2d_advector/vof_advection.h"
#include "examples/2d_advector/vtk.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection_amr.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"

void resetMoments(const Data<IRL::LocalizedParaboloidLink<double>>&
                      a_link_localized_paraboloid,
                  Data<IRL::GeneralMoments3D<2>>* a_liquid_moments,
                  Data<IRL::GeneralMoments3D<2>>* a_gas_moments) {
  const BasicMesh& mesh = a_link_localized_paraboloid.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        auto moments_exact =
            IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
                cell,
                a_link_localized_paraboloid(i, j, k).getNextReconstruction());
        auto moments_quad = IRL::getVolumeMoments<
            IRL::SeparatedMoments<IRL::GeneralMoments3D<2>>>(
            cell, a_link_localized_paraboloid(i, j, k).getNextReconstruction());
        for (int l = 1; l < 4; ++l) {
          (*a_liquid_moments)(i, j, k)[l] = moments_exact[0].centroid()[l - 1];
          (*a_gas_moments)(i, j, k)[l] = moments_exact[1].centroid()[l - 1];
        }
        for (int l = 4; l < 10; ++l) {
          (*a_liquid_moments)(i, j, k)[l] = moments_quad[0][l];
          (*a_gas_moments)(i, j, k)[l] = moments_quad[1][l];
        }
      }
    }
  }
}

void resetCentroids(const Data<IRL::LocalizedParaboloidLink<double>>&
                        a_link_localized_paraboloid,
                    Data<IRL::Pt>* a_liquid_centroid,
                    Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_link_localized_paraboloid.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        auto moments = IRL::getNormalizedVolumeMoments<
            IRL::SeparatedMoments<IRL::VolumeMoments>>(
            cell, a_link_localized_paraboloid(i, j, k).getNextReconstruction());
        (*a_liquid_centroid)(i, j, k) = moments[0].centroid();
        (*a_gas_centroid)(i, j, k) = moments[1].centroid();
      }
    }
  }
}

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid) {
  IRL::LocalizedParaboloidLink<double>* neighbor_ptr;
  // Provide mesh connectivity information.
  IRL::UnsignedIndex_t unique_id = 0;

  for (int i = a_mesh.imino(); i <= a_mesh.imaxo(); ++i) {
    for (int j = a_mesh.jmino(); j <= a_mesh.jmaxo(); ++j) {
      for (int k = a_mesh.kmino(); k <= a_mesh.kmaxo(); ++k) {
        (*a_link_localized_paraboloid)(i, j, k).setId(unique_id);
        neighbor_ptr = i - 1 < a_mesh.imino()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i - 1, j, k);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            0, neighbor_ptr);
        neighbor_ptr = i + 1 > a_mesh.imaxo()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i + 1, j, k);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            1, neighbor_ptr);
        neighbor_ptr = j - 1 < a_mesh.jmino()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i, j - 1, k);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            2, neighbor_ptr);
        neighbor_ptr = j + 1 > a_mesh.jmaxo()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i, j + 1, k);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            3, neighbor_ptr);
        neighbor_ptr = k - 1 < a_mesh.kmino()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i, j, k - 1);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            4, neighbor_ptr);
        neighbor_ptr = k + 1 > a_mesh.kmaxo()
                           ? nullptr
                           : &(*a_link_localized_paraboloid)(i, j, k + 1);
        (*a_link_localized_paraboloid)(i, j, k).setEdgeConnectivity(
            5, neighbor_ptr);
        ++unique_id;
      }
    }
  }
}

std::array<int, 3> getIndexFromTag(const BasicMesh& a_mesh,
                                   const IRL::UnsignedIndex_t a_tag) {
  std::array<int, 3> indices;
  auto int_tag = static_cast<int64_t>(a_tag);
  indices[0] = static_cast<int>(int_tag / (a_mesh.getNzo() * a_mesh.getNyo()));
  indices[1] = static_cast<int>(
      (int_tag - indices[0] * (a_mesh.getNzo() * a_mesh.getNyo())) /
      a_mesh.getNzo());
  indices[2] =
      static_cast<int>(int_tag - a_mesh.getNzo() * indices[1] -
                       (indices[0] * (a_mesh.getNzo() * a_mesh.getNyo())));

  return {indices[0] + a_mesh.imino(), indices[1] + a_mesh.jmino(),
          indices[2] + a_mesh.kmino()};
}

void advectVOF(
    const std::string& a_simulation_type, const std::string& a_advection_method,
    const std::string& a_reconstruction_method, const double a_dt,
    const double a_time, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<IRL::GeneralMoments3D<2>>* a_liquid_moments,
    Data<IRL::GeneralMoments3D<2>>* a_gas_moments,
    Data<IRL::Paraboloid>* a_interface) {
  if (a_advection_method == "SemiLagrangianCorrected") {
    SemiLagrangianCorrected::advectVOF(
        a_simulation_type, a_reconstruction_method, a_dt, a_time, a_U, a_V, a_W,
        a_link_localized_paraboloid, a_liquid_moments, a_gas_moments);
  } else {
    std::cout << "Unknown advection method of : " << a_advection_method << '\n';
    std::cout << "Value entries are: SemiLagrangianCorrected. \n";
    std::exit(-1);
  }
}

void SemiLagrangianCorrected::advectVOF(
    const std::string& a_simulation_type,
    const std::string& a_reconstruction_method, const double a_dt,
    const double a_time, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<IRL::GeneralMoments3D<2>>* a_liquid_moments,
    Data<IRL::GeneralMoments3D<2>>* a_gas_moments) {
  const BasicMesh& mesh = a_liquid_moments->getMesh();

  const std::array<double, 3> (*getExactVelocity)(const IRL::Pt& a_location,
                                                  const double a_time);
  const std::array<std::array<double, 3>, 3> (*getExactGradient)(
      const IRL::Pt& a_location, const double a_time);

  const IRL2D::Vec (*getExactVelocity2D)(const double t, const IRL2D::Vec& P);
  const IRL2D::Mat (*getExactGradient2D)(const double t, const IRL2D::Vec& P);

  if (a_simulation_type == "Rotation2D") {
    getExactVelocity = Rotation2D::getExactVelocity;
    getExactGradient = Rotation2D::getExactVelocityGradient;
    getExactVelocity2D = Rotation2D::getExactVelocity2D;
    getExactGradient2D = Rotation2D::getExactVelocityGradient2D;
  }

  resetMoments(*a_link_localized_paraboloid, a_liquid_moments, a_gas_moments);

  Data<double> U_face(&mesh);
  Data<double> V_face(&mesh);
  Data<double> W_face(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax() + 1; ++k) {
        U_face(i, j, k) = 0.5 * (a_U(i, j, k) + a_U(i - 1, j, k));
        V_face(i, j, k) = 0.5 * (a_V(i, j, k) + a_V(i, j - 1, k));
        W_face(i, j, k) = 0.5 * (a_W(i, j, k) + a_W(i, j, k - 1));
      }
    }
  }

  Data<IRL::Normal> gradU(&mesh);
  Data<IRL::Normal> gradV(&mesh);
  Data<IRL::Normal> gradW(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        gradU(i, j, k)[0] = (-a_U(i + 2, j, k) + 8.0 * a_U(i + 1, j, k) -
                             8.0 * a_U(i - 1, j, k) + a_U(i - 2, j, k)) /
                            (12.0 * mesh.dx());
        gradU(i, j, k)[1] = (-a_U(i, j + 2, k) + 8.0 * a_U(i, j + 1, k) -
                             8.0 * a_U(i, j - 1, k) + a_U(i, j - 2, k)) /
                            (12.0 * mesh.dy());
        gradU(i, j, k)[2] = (-a_U(i, j, k + 2) + 8.0 * a_U(i, j, k + 1) -
                             8.0 * a_U(i, j, k - 1) + a_U(i, j, k - 2)) /
                            (12.0 * mesh.dz());
        gradV(i, j, k)[0] = (-a_V(i + 2, j, k) + 8.0 * a_V(i + 1, j, k) -
                             8.0 * a_V(i - 1, j, k) + a_V(i - 2, j, k)) /
                            (12.0 * mesh.dx());
        gradV(i, j, k)[1] = (-a_V(i, j + 2, k) + 8.0 * a_V(i, j + 1, k) -
                             8.0 * a_V(i, j - 1, k) + a_V(i, j - 2, k)) /
                            (12.0 * mesh.dy());
        gradV(i, j, k)[2] = (-a_V(i, j, k + 2) + 8.0 * a_V(i, j, k + 1) -
                             8.0 * a_V(i, j, k - 1) + a_V(i, j, k - 2)) /
                            (12.0 * mesh.dz());
        gradW(i, j, k)[0] = (-a_W(i + 2, j, k) + 8.0 * a_W(i + 1, j, k) -
                             8.0 * a_W(i - 1, j, k) + a_W(i - 2, j, k)) /
                            (12.0 * mesh.dx());
        gradW(i, j, k)[1] = (-a_W(i, j + 2, k) + 8.0 * a_W(i, j + 1, k) -
                             8.0 * a_W(i, j - 1, k) + a_W(i, j - 2, k)) /
                            (12.0 * mesh.dy());
        gradW(i, j, k)[2] = (-a_W(i, j, k + 2) + 8.0 * a_W(i, j, k + 1) -
                             8.0 * a_W(i, j, k - 1) + a_W(i, j, k - 2)) /
                            (12.0 * mesh.dz());
      }
    }
  }
  gradU.updateBorder();
  gradV.updateBorder();
  gradW.updateBorder();

  using MomentType = IRL::GeneralMoments3D<2>;

  // Allocate storage for face fluxes
  Data<IRL::SeparatedMoments<MomentType>> face_flux[3] = {
      Data<IRL::SeparatedMoments<MomentType>>(&mesh),
      Data<IRL::SeparatedMoments<MomentType>>(&mesh),
      Data<IRL::SeparatedMoments<MomentType>>(&mesh)};

  // Calculate CFL
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        CFL = std::fmax(CFL,
                        std::fmax(a_U(i, j, k) * a_dt / mesh.dx(),
                                  std::fmax(a_V(i, j, k) * a_dt / mesh.dy(),
                                            a_W(i, j, k) * a_dt / mesh.dz())));
      }
    }
  }

  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        band(i, j, k) = 0;
        const double liquid_volume_fraction =
            (*a_liquid_moments)(i, j, k)[0] / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          band(i, j, k) = 1;
        }
      }
    }
  }
  band.updateBorder();

  const int nlayers = static_cast<int>(std::ceil(CFL));
  for (int n = 0; n < nlayers; ++n) {
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
          if (band(i, j, k) == 0) {
            for (int ii = -1; ii <= 1; ++ii) {
              for (int jj = -1; jj <= 1; ++jj) {
                for (int kk = -1; kk <= 1; ++kk) {
                  if (band(i + ii, j + jj, k + kk) == n + 1) {
                    band(i, j, k) = n + 2;
                  }
                }
              }
            }
          }
        }
      }
    }
    band.updateBorder();
  }

  int nmixed_global = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax() + 1; ++k) {
        if (band(i, j, k) > 0 || band(i - 1, j, k) > 0 ||
            band(i, j - 1, k) > 0 || band(i, j, k - 1) > 0) {
          nmixed_global++;
        }
      }
    }
  }

  const int size_moments = 3 * 2 * 10;
  std::vector<double> flux_global(size_moments * nmixed_global);

  int count = 0, count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax() + 1; ++k) {
        for (int dim = 0; dim < 3; ++dim) {
          // Store face flux
          (face_flux[dim])(i, j, k) =
              IRL::SeparatedMoments<MomentType>::fromScalarConstant(0.0);
        }
        if (1) {
          auto cell = IRL2D::RectangleListFromBounds(mesh.x(i), mesh.x(i + 1),
                                                     mesh.y(j), mesh.y(j + 1));
          IRL2D::Print(cell);

          auto datum = IRL2D::Vec(0.0, 0.0);
          auto frame = IRL2D::ReferenceFrame(-M_PI / 3.0);
          auto parabola = IRL2D::Parabola(datum, frame, -1.0);
          std::cout << "Parabola  : " << parabola << std::endl;

          // cell[0].second.y() += 0.1;
          // cell[1].second.y() -= 2.0;
          // cell[2].second.y() -= 0.1;
          // cell[3].second.y() += 2.0;
          IRL2D::ToVTK(cell, "cell");

          auto clipped_cell = IRL2D::ParabolaClip(cell, parabola);
          IRL2D::ToVTK(clipped_cell, "clipped_cell");

          clipped_cell = IRL2D::ClipByRectangleAndParabola(
              cell, IRL2D::Vec(-0.4, -0.4), IRL2D::Vec(0.4, 0.4), parabola);

          IRL2D::ToVTK(clipped_cell, "clipped_by_rectangle_and_parabola");
          std::cout << "Clipped cell has " << clipped_cell.size() << " curves"
                    << " and volume = " << IRL2D::CellVolume(clipped_cell)
                    << std::endl;

          cell = IRL2D::RectangleListFromBounds(-0.4, 0.4, -0.4, 0.4);
          const bool correct_fluxes = true;

          std::vector<IRL2D::BezierList> flux_cell;
          double area = IntegrateFlux(cell[0].first, cell[1].first, a_dt,
                                      a_time, getExactVelocity2D);
          flux_cell.push_back(IRL2D::CreateFluxCell(
              cell[0].first, cell[1].first, -a_dt, a_time, getExactVelocity2D,
              getExactGradient2D, correct_fluxes, area));
          area = IntegrateFlux(cell[1].first, cell[2].first, a_dt, a_time,
                               getExactVelocity2D);
          flux_cell.push_back(IRL2D::CreateFluxCell(
              cell[1].first, cell[2].first, -a_dt, a_time, getExactVelocity2D,
              getExactGradient2D, correct_fluxes, area));
          area = IntegrateFlux(cell[2].first, cell[3].first, a_dt, a_time,
                               getExactVelocity2D);
          flux_cell.push_back(IRL2D::CreateFluxCell(
              cell[2].first, cell[3].first, -a_dt, a_time, getExactVelocity2D,
              getExactGradient2D, correct_fluxes, area));
          area = IntegrateFlux(cell[3].first, cell[0].first, a_dt, a_time,
                               getExactVelocity2D);
          flux_cell.push_back(IRL2D::CreateFluxCell(
              cell[3].first, cell[0].first, -a_dt, a_time, getExactVelocity2D,
              getExactGradient2D, correct_fluxes, area));
          IRL2D::ToVTK(flux_cell, "flux_cell");

          std::cout << "Divergence = " << std::scientific
                    << std::setprecision(6)
                    << (IRL2D::CellVolume(flux_cell[0]) +
                        IRL2D::CellVolume(flux_cell[2]) +
                        IRL2D::CellVolume(flux_cell[1]) +
                        IRL2D::CellVolume(flux_cell[3]))
                    << std::endl;

          exit(0);
        }
        if (band(i, j, k) > 0 || band(i - 1, j, k) > 0 ||
            band(i, j - 1, k) > 0 || band(i, j, k - 1) > 0) {
          auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
              IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));

          //////////////////////////////////////////////////////
          //////////////////////////////////////////////////////
          //////////////////////////////////////////////////////

          // ASHISH::Point p00 = {mesh.x(i), mesh.y(j)},
          //               p10 = {mesh.x(i + 1), mesh.y(j)},
          //               p11 = {mesh.x(i + 1), mesh.y(j + 1)},
          //               p01 = {mesh.x(i), mesh.y(j + 1)};
          // ASHISH::PointAndControl pc0 = {
          //     {p00.x, p00.y}, {.5 * (p00.x + p10.x), .5 * (p00.y + p10.y)}};
          // ASHISH::PointAndControl pc1 = {
          //     {p10.x, p10.y}, {.5 * (p10.x + p11.x), .5 * (p10.y + p11.y)}};
          // ASHISH::PointAndControl pc2 = {
          //     {p11.x, p11.y}, {.5 * (p11.x + p01.x), .5 * (p11.y + p01.y)}};
          // ASHISH::PointAndControl pc3 = {
          //     {p01.x, p01.y}, {.5 * (p01.x + p00.x), .5 * (p01.y + p00.y)}};

          // ASHISH::BezierList list_cell{pc0, pc1, pc2, pc3};

          // ASHISH::Parabola = {0., 0., 0., .5 * M_PI};

          // ASHISH::Point_N p00 = {{mesh.x(i), mesh.y(j)}, {0., 0.}};
          // ASHISH::Point_N p10 = {{mesh.x(i + 1), mesh.y(j)}, {0., 0.}};
          // ASHISH::Point_N p11 = {{mesh.x(i + 1), mesh.y(j + 1)}, {0., 0.}};
          // ASHISH::Point_N p01 = {{mesh.x(i), mesh.y(j + 1)}, {0., 0.}};

          // std::array<std::vector<ASHISH ::Point_N>, 4> edges;
          // ASHISH::EdgeTransport::transport_edge_with_tangents(
          //     p00, p10, a_dt, a_time, getExactVelocity2D, getExactGradient2D,
          //     0, edges[0]);
          // ASHISH::EdgeTransport::transport_edge_with_tangents(
          //     p10, p11, a_dt, a_time, getExactVelocity2D, getExactGradient2D,
          //     0, edges[1]);
          // ASHISH::EdgeTransport::transport_edge_with_tangents(
          //     p11, p01, a_dt, a_time, getExactVelocity2D, getExactGradient2D,
          //     0, edges[2]);
          // ASHISH::EdgeTransport::transport_edge_with_tangents(
          //     p01, p00, a_dt, a_time, getExactVelocity2D, getExactGradient2D,
          //     0, edges[3]);

          // for (auto edge : edges) {
          //   for (int n = 0; n < edge.size() - 1; n++) {
          //     auto P0 = edge[n];
          //     auto P1 = edge[n + 1];
          //     P1.flip_dir();
          //     std::cout << "Point   " << n << " = ( " << edge[n].P.x << ", "
          //               << edge[n].P.y << " )" << std::endl;
          //     auto inter = ASHISH::EdgeTransport::rays_intersect(P0, P1);
          //     if (inter.first) {
          //       double t0 = inter.second;
          //       ASHISH ::Point ctrl = {P0.P.x + P0.N.u * t0,
          //                             P0.P.y + P0.N.v * t0};
          //       std::cout << "Control   = ( " << ctrl.x << ", " << ctrl.y
          //                 << " )" << std::endl;
          //     } else {
          //       std::cout << "Could not find control point..." << std::endl;
          //     }
          //     P1.flip_dir();
          //   }
          //   std::cout << "Point   " << edge.size() - 1 << " = ( "
          //             << edge[edge.size() - 1].P.x << ", "
          //             << edge[edge.size() - 1].P.y << " )" << std::endl;
          // }

          //////////////////////////////////////////////////////
          //////////////////////////////////////////////////////
          //////////////////////////////////////////////////////

          // Get the back projected Dodecahedron.
          IRL::Dodecahedron transported_cell;
          for (IRL::UnsignedIndex_t n = 0; n < 8; ++n) {
            transported_cell[n] =
                back_project_vertex(cell[n], a_simulation_type, -a_dt, a_time);
            // back_project_vertex(cell[n], -a_dt, a_U, a_V, a_W);
          }
          // Initialize face flux hexahedra
          IRL::CappedDodecahedron face_cell[3];
          IRL::Pt face_center_pt =
              0.25 * (cell[7] + cell[4] + cell[5] + cell[6]);
          IRL::Pt correction_pt = back_project_vertex(
              face_center_pt, a_simulation_type, -a_dt, a_time);
          // back_project_vertex(face_center_pt, -a_dt, a_U, a_V, a_W);
          face_cell[0] = IRL::CappedDodecahedron(
              {cell[7], cell[4], cell[5], cell[6], transported_cell[7],
               transported_cell[4], transported_cell[5], transported_cell[6],
               correction_pt});
          face_center_pt = 0.25 * (cell[0] + cell[4] + cell[7] + cell[3]);
          correction_pt = back_project_vertex(face_center_pt, a_simulation_type,
                                              -a_dt, a_time);
          // back_project_vertex(face_center_pt, -a_dt, a_U, a_V, a_W);
          face_cell[1] = IRL::CappedDodecahedron(
              {cell[0], cell[4], cell[7], cell[3], transported_cell[0],
               transported_cell[4], transported_cell[7], transported_cell[3],
               correction_pt});
          face_center_pt = 0.25 * (cell[5] + cell[4] + cell[0] + cell[1]);
          correction_pt = back_project_vertex(face_center_pt, a_simulation_type,
                                              -a_dt, a_time);
          // back_project_vertex(face_center_pt, -a_dt, a_U, a_V, a_W);
          face_cell[2] = IRL::CappedDodecahedron(
              {cell[5], cell[4], cell[0], cell[1], transported_cell[5],
               transported_cell[4], transported_cell[0], transported_cell[1],
               correction_pt});
          // Tack on the corrective 9th vertex for each face hexahedra
          face_cell[0].adjustCapToMatchVolume(a_dt * U_face(i, j, k) *
                                              mesh.dy() * mesh.dz());
          face_cell[1].adjustCapToMatchVolume(a_dt * V_face(i, j, k) *
                                              mesh.dx() * mesh.dz());
          face_cell[2].adjustCapToMatchVolume(a_dt * W_face(i, j, k) *
                                              mesh.dx() * mesh.dy());

          for (int dim = 0; dim < 3; ++dim) {
            // Store face flux
            (face_flux[dim])(i, j, k) =
                IRL::getVolumeMoments<IRL::SeparatedMoments<MomentType>>(
                    face_cell[dim], (*a_link_localized_paraboloid)(i, j, k));
            auto flux_m0_m1 = IRL::getVolumeMoments<
                IRL::SeparatedMoments<IRL::VolumeMoments>>(
                face_cell[dim], (*a_link_localized_paraboloid)(i, j, k));
            (face_flux[dim])(i, j, k)[0][0] = flux_m0_m1[0].volume();
            (face_flux[dim])(i, j, k)[0][1] = flux_m0_m1[0].centroid()[0];
            (face_flux[dim])(i, j, k)[0][2] = flux_m0_m1[0].centroid()[1];
            (face_flux[dim])(i, j, k)[0][3] = flux_m0_m1[0].centroid()[2];
            (face_flux[dim])(i, j, k)[1][0] = flux_m0_m1[1].volume();
            (face_flux[dim])(i, j, k)[1][1] = flux_m0_m1[1].centroid()[0];
            (face_flux[dim])(i, j, k)[1][2] = flux_m0_m1[1].centroid()[1];
            (face_flux[dim])(i, j, k)[1][3] = flux_m0_m1[1].centroid()[2];
            for (int m = 0; m < 10; m++) {
              flux_global[count_local++] = (face_flux[dim])(i, j, k)[0][m];
              flux_global[count_local++] = (face_flux[dim])(i, j, k)[1][m];
            }
          }
          count++;
        }
      }
    }
  }

  // face_flux[0].updateBorder();
  // face_flux[1].updateBorder();
  // face_flux[2].updateBorder();

  count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax() + 1; ++k) {
        if (band(i, j, k) > 0 || band(i - 1, j, k) > 0 ||
            band(i, j - 1, k) > 0 || band(i, j, k - 1) > 0) {
          for (int dim = 0; dim < 3; ++dim) {
            for (int m = 0; m < 10; m++) {
              (face_flux[dim])(i, j, k)[0][m] = flux_global[count_local++];
              (face_flux[dim])(i, j, k)[1][m] = flux_global[count_local++];
            }
          }
        }
      }
    }
  }

  // Now calculate VOF from the face fluxes.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto cell_moments =
            IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
                cell, IRL::Paraboloid::createAlwaysAbove());
        const double cell_volume = mesh.cell_volume();
        const double liquid_volume_fraction =
            (*a_liquid_moments)(i, j, k)[0] / cell_volume;
        if (band(i, j, k) > 0) {
          // Compute moments in transported cell image
          (*a_liquid_moments)(i, j, k) =
              (*a_liquid_moments)(i, j, k) + (face_flux[0])(i, j, k)[0] -
              (face_flux[0])(i + 1, j, k)[0] + (face_flux[1])(i, j, k)[0] -
              (face_flux[1])(i, j + 1, k)[0] + (face_flux[2])(i, j, k)[0] -
              (face_flux[2])(i, j, k + 1)[0];
          (*a_gas_moments)(i, j, k) =
              (*a_gas_moments)(i, j, k) + (face_flux[0])(i, j, k)[1] -
              (face_flux[0])(i + 1, j, k)[1] + (face_flux[1])(i, j, k)[1] -
              (face_flux[1])(i, j + 1, k)[1] + (face_flux[2])(i, j, k)[1] -
              (face_flux[2])(i, j, k + 1)[1];

          const std::array<Data<IRL::GeneralMoments3D<2>>*, 2> moment_list(
              {a_liquid_moments, a_gas_moments});

          for (int m = 0; m < 2; m++) {
            const auto moments = moment_list[m];
            // Extract moment in tensor notation
            //  1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, x^2 y, ...
            const double M0 = (*moments)(i, j, k)[0];
            const Eigen::Matrix<double, 3, 1> M1{(*moments)(i, j, k)[1],
                                                 (*moments)(i, j, k)[2],
                                                 (*moments)(i, j, k)[3]};
            const Eigen::Matrix<double, 3, 3> M2{
                {(*moments)(i, j, k)[4], (*moments)(i, j, k)[5],
                 (*moments)(i, j, k)[6]},
                {(*moments)(i, j, k)[5], (*moments)(i, j, k)[7],
                 (*moments)(i, j, k)[8]},
                {(*moments)(i, j, k)[6], (*moments)(i, j, k)[8],
                 (*moments)(i, j, k)[9]}};

            // Correct moment 0
            double M0_final = M0;
            if (M0 < 0.0) {
              M0_final = 0.0;
            } else if (M0 > cell_volume) {
              M0_final = cell_volume;
            }
            (*moments)(i, j, k)[0] = M0_final;

            auto x0_k1 = IRL::Pt(M1(0), M1(1), M1(2)) / IRL::safelyEpsilon(M0);
            const auto cell_centroid =
                IRL::Pt(mesh.xm(i), mesh.ym(j), mesh.zm(k));
            const auto img_centroid = back_project_vertex(
                cell_centroid, a_simulation_type, -a_dt, a_time);
            const double dist = IRL::magnitude(img_centroid - x0_k1);
            // if (dist > 2.0 * std::sqrt(3.0) * mesh.dx()) x0_k1 =
            // img_centroid;
            const auto u0_k1 = IRL::Normal::fromRawDoublePointer(
                getExactVelocity(x0_k1, a_time).data());
            const auto gradu0_k1 = getExactGradient(x0_k1, a_time);
            const auto x0_k2 = IRL::Pt(x0_k1 + a_dt * u0_k1 * 0.5);
            const auto u0_k2 = IRL::Normal::fromRawDoublePointer(
                getExactVelocity(x0_k2, a_time + 0.5 * a_dt).data());
            const auto gradu0_k2 = getExactGradient(x0_k2, a_time + 0.5 * a_dt);
            const auto x0_k3 = IRL::Pt(x0_k1 + a_dt * u0_k2 * 0.5);
            const auto u0_k3 = IRL::Normal::fromRawDoublePointer(
                getExactVelocity(x0_k3, a_time + 0.5 * a_dt).data());
            const auto gradu0_k3 = getExactGradient(x0_k3, a_time + 0.5 * a_dt);
            const auto x0_k4 = IRL::Pt(x0_k1 + a_dt * u0_k3);
            const auto u0_k4 = IRL::Normal::fromRawDoublePointer(
                getExactVelocity(x0_k4, a_time + a_dt).data());
            const auto gradu0_k4 = getExactGradient(x0_k4, a_time + a_dt);

            Eigen::Matrix<double, 3, 1> X0_k1, X0_k2, X0_k3, X0_k4;
            Eigen::Matrix<double, 3, 1> U0_k1, U0_k2, U0_k3, U0_k4;
            Eigen::Matrix<double, 3, 3> gradU0_k1, gradU0_k2, gradU0_k3,
                gradU0_k4;

            for (int ii = 0; ii < 3; ii++) {
              X0_k1(ii) = x0_k1[ii];
              X0_k2(ii) = x0_k2[ii];
              X0_k3(ii) = x0_k3[ii];
              X0_k4(ii) = x0_k4[ii];
              U0_k1(ii) = u0_k1[ii];
              U0_k2(ii) = u0_k2[ii];
              U0_k3(ii) = u0_k3[ii];
              U0_k4(ii) = u0_k4[ii];
              for (int jj = 0; jj < 3; jj++) {
                gradU0_k1(ii, jj) = gradu0_k1[ii][jj];
                gradU0_k2(ii, jj) = gradu0_k2[ii][jj];
                gradU0_k3(ii, jj) = gradu0_k3[ii][jj];
                gradU0_k4(ii, jj) = gradu0_k4[ii][jj];
              }
            }

            const auto M1_final =
                M1 + a_dt * M0_final *
                         (U0_k1 + 2.0 * U0_k2 + 2.0 * U0_k3 + U0_k4) / 6.0;

            const auto M2t0_k1 = M0_final * X0_k1 * U0_k1.transpose();
            const auto M2t0_k2 = M0_final * X0_k2 * U0_k2.transpose();
            const auto M2t0_k3 = M0_final * X0_k3 * U0_k3.transpose();
            const auto M2t0_k4 = M0_final * X0_k4 * U0_k4.transpose();
            const auto M2t1_k1 =
                -M0_final * (X0_k1 * X0_k1.transpose()) * gradU0_k1.transpose();
            const auto M2t1_k2 =
                -M0_final * (X0_k2 * X0_k2.transpose()) * gradU0_k2.transpose();
            const auto M2t1_k3 =
                -M0_final * (X0_k3 * X0_k3.transpose()) * gradU0_k3.transpose();
            const auto M2t1_k4 =
                -M0_final * (X0_k4 * X0_k4.transpose()) * gradU0_k4.transpose();
            const auto M2t2_k1 = M2 * gradU0_k1.transpose();
            const auto M2t2_k2 = (M2 + a_dt *
                                           (M2t0_k1 + M2t0_k1.transpose() +
                                            M2t1_k1 + M2t1_k1.transpose() +
                                            M2t2_k1 + M2t2_k1.transpose()) /
                                           2.0) *
                                 gradU0_k2.transpose();
            const auto M2t2_k3 = (M2 + a_dt *
                                           (M2t0_k2 + M2t0_k2.transpose() +
                                            M2t1_k2 + M2t1_k2.transpose() +
                                            M2t2_k2 + M2t2_k2.transpose()) /
                                           2.0) *
                                 gradU0_k3.transpose();
            const auto M2t2_k4 = (M2 + a_dt * (M2t0_k3 + M2t0_k3.transpose() +
                                               M2t1_k3 + M2t1_k3.transpose() +
                                               M2t2_k3 + M2t2_k3.transpose())) *
                                 gradU0_k4.transpose();
            const auto M2_rk4 =
                a_dt *
                ((M2t0_k1 + 2.0 * M2t0_k2 + 2.0 * M2t0_k3 + M2t0_k4) +
                 (M2t1_k1 + 2.0 * M2t1_k2 + 2.0 * M2t1_k3 + M2t1_k4) +
                 (M2t2_k1 + 2.0 * M2t2_k2 + 2.0 * M2t2_k3 + M2t2_k4)) /
                6.0;
            const auto M2_final = M2 + M2_rk4 + M2_rk4.transpose();

            (*moments)(i, j, k)[1] = M1_final(0);
            (*moments)(i, j, k)[2] = M1_final(1);
            (*moments)(i, j, k)[3] = M1_final(2);
            (*moments)(i, j, k)[4] = M2_final(0, 0);
            (*moments)(i, j, k)[5] = M2_final(0, 1);
            (*moments)(i, j, k)[6] = M2_final(0, 2);
            (*moments)(i, j, k)[7] = M2_final(1, 1);
            (*moments)(i, j, k)[8] = M2_final(1, 2);
            (*moments)(i, j, k)[9] = M2_final(2, 2);
          }
        } else if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
                   liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
          const double clipped_vfrac =
              std::max(0.0, std::min(1.0, liquid_volume_fraction));
          (*a_liquid_moments)(i, j, k) = clipped_vfrac * cell_moments;
          (*a_gas_moments)(i, j, k) = (1.0 - clipped_vfrac) * cell_moments;
        }
      }
    }
  }
  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();
  correctMomentLocations(a_liquid_moments);
  correctMomentLocations(a_gas_moments);
}

void correctMomentLocations(Data<IRL::GeneralMoments3D<2>>* a_liquid_moments) {
  const BasicMesh& mesh = (*a_liquid_moments).getMesh();
  // Fix distance to recreate volume fraction

  // Extract moment in tensor notation
  Data<Eigen::Matrix<double, 3, 1>> Shift(&mesh);
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        Shift(i, j, k).setZero();
      }
    }
  }

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        Shift(i, j, k)(0) -= mesh.lx();
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        Shift(i, j, k)(0) += mesh.lx();
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        Shift(i, j, k)(1) -= mesh.ly();
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        Shift(i, j, k)(1) += mesh.ly();
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        Shift(i, j, k)(2) -= mesh.lz();
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        Shift(i, j, k)(2) += mesh.lz();
      }
    }
  }

  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        const double M0 = (*a_liquid_moments)(i, j, k)[0];
        const Eigen::Matrix<double, 3, 1> M1{(*a_liquid_moments)(i, j, k)[1],
                                             (*a_liquid_moments)(i, j, k)[2],
                                             (*a_liquid_moments)(i, j, k)[3]};
        const Eigen::Matrix<double, 3, 3> M2{
            {(*a_liquid_moments)(i, j, k)[4], (*a_liquid_moments)(i, j, k)[5],
             (*a_liquid_moments)(i, j, k)[6]},
            {(*a_liquid_moments)(i, j, k)[5], (*a_liquid_moments)(i, j, k)[7],
             (*a_liquid_moments)(i, j, k)[8]},
            {(*a_liquid_moments)(i, j, k)[6], (*a_liquid_moments)(i, j, k)[8],
             (*a_liquid_moments)(i, j, k)[9]}};

        auto M1_final = M1 + M0 * Shift(i, j, k);
        auto M2_final = M2 + M1 * Shift(i, j, k).transpose() +
                        Shift(i, j, k) * M1.transpose() +
                        M0 * Shift(i, j, k) * Shift(i, j, k).transpose();

        (*a_liquid_moments)(i, j, k)[1] = M1_final(0);
        (*a_liquid_moments)(i, j, k)[2] = M1_final(1);
        (*a_liquid_moments)(i, j, k)[3] = M1_final(2);
        (*a_liquid_moments)(i, j, k)[4] = M2_final(0, 0);
        (*a_liquid_moments)(i, j, k)[5] = M2_final(0, 1);
        (*a_liquid_moments)(i, j, k)[6] = M2_final(0, 2);
        (*a_liquid_moments)(i, j, k)[7] = M2_final(1, 1);
        (*a_liquid_moments)(i, j, k)[8] = M2_final(1, 2);
        (*a_liquid_moments)(i, j, k)[9] = M2_final(2, 2);
      }
    }
  }
}

void correctCentroidLocation(Data<IRL::Pt>* a_liquid_centroid,
                             Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = (*a_liquid_centroid).getMesh();
  // Fix distance to recreate volume fraction

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[0] - mesh.lx();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[0] - mesh.lx();
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[0] + mesh.lx();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[0] + mesh.lx();
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[1] - mesh.ly();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[1] - mesh.ly();
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[1] + mesh.ly();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[1] + mesh.ly();
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[2] - mesh.lz();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[2] - mesh.lz();
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        (*a_liquid_centroid)(i, j, k) =
            (*a_liquid_centroid)(i, j, k)[2] + mesh.lz();
        (*a_gas_centroid)(i, j, k) = (*a_gas_centroid)(i, j, k)[2] + mesh.lz();
      }
    }
  }
}
