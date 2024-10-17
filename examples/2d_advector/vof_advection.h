// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_

#include <string>

#include "examples/2d_advector/rotation_2d.h"
#include "irl/geometry/general/pt.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

#include "examples/2d_advector/data.h"

void resetCentroids(const Data<IRL::LocalizedParaboloidLink<double>>&
                        a_link_localized_paraboloid,
                    Data<IRL::Pt>* a_liquid_centroid,
                    Data<IRL::Pt>* a_gas_centroid);

void resetMoments(const Data<IRL::LocalizedParaboloidLink<double>>&
                      a_link_localized_paraboloid,
                  Data<IRL::GeneralMoments3D<2>>* a_liquid_moments,
                  Data<IRL::GeneralMoments3D<2>>* a_gas_moments);

std::array<int, 3> getIndexFromTag(const BasicMesh& a_mesh,
                                   const IRL::UnsignedIndex_t a_tag);

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid);

void advectVOF(
    const std::string& a_simulation_type, const std::string& a_advection_method,
    const std::string& a_reconstruction_method, const double a_dt,
    const double a_time, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<IRL::GeneralMoments3D<2>>* a_liquid_moments,
    Data<IRL::GeneralMoments3D<2>>* a_gas_moments,
    Data<IRL::Paraboloid>* a_interface);

struct SemiLagrangianCorrected {
  static void advectVOF(
      const std::string& a_simulation_type,
      const std::string& a_reconstruction_method, const double a_dt,
      const double a_time, const Data<double>& a_U, const Data<double>& a_V,
      const Data<double>& a_W,
      Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
      Data<IRL::GeneralMoments3D<2>>* a_liquid_moments,
      Data<IRL::GeneralMoments3D<2>>* a_gas_moments);
};

inline IRL::Vec3<double> getVelocity(const IRL::Pt& a_location,
                                     const Data<double>& a_U,
                                     const Data<double>& a_V,
                                     const Data<double>& a_W);

inline IRL::Pt back_project_vertex(const IRL::Pt& a_initial_pt,
                                   const double a_dt, const Data<double>& a_U,
                                   const Data<double>& a_V,
                                   const Data<double>& a_W);

inline IRL::Pt back_project_vertex(const IRL::Pt& a_initial_pt,
                                   const std::string& a_simulation_type,
                                   const double a_dt, const double a_time);

void correctCentroidLocation(Data<IRL::Pt>* a_liquid_centroid,
                             Data<IRL::Pt>* a_gas_centroid);

void correctMomentLocations(Data<IRL::GeneralMoments3D<2>>* a_liquid_moments);

// ************************************************
//     Inlined functions below this
// ************************************************
// Spatial RK4 projection of a point in a velocity field..
inline IRL::Pt back_project_vertex(const IRL::Pt& a_initial_pt,
                                   const double a_dt, const Data<double>& a_U,
                                   const Data<double>& a_V,
                                   const Data<double>& a_W) {
  auto v1 = getVelocity(a_initial_pt, a_U, a_V, a_W);
  // return a_initial_pt + IRL::Pt::fromVec3(a_dt * v1);
  auto v2 = getVelocity(a_initial_pt + IRL::Pt::fromVec3(0.5 * a_dt * v1), a_U,
                        a_V, a_W);
  auto v3 = getVelocity(a_initial_pt + IRL::Pt::fromVec3(0.5 * a_dt * v2), a_U,
                        a_V, a_W);
  auto v4 =
      getVelocity(a_initial_pt + IRL::Pt::fromVec3(a_dt * v3), a_U, a_V, a_W);
  return a_initial_pt +
         IRL::Pt::fromVec3(a_dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0);
}

inline IRL::Pt back_project_vertex(const IRL::Pt& a_initial_pt,
                                   const std::string& a_simulation_type,
                                   const double a_dt, const double a_time) {
  const std::array<double, 3> (*getExactVelocity)(const IRL::Pt& a_location,
                                                  const double a_time);
  if (a_simulation_type == "Rotation2D") {
    getExactVelocity = Rotation2D::getExactVelocity;
  }

  auto v1 = IRL::Vec3<double>::fromRawDoublePointer(
      getExactVelocity(a_initial_pt, a_time).data());
  auto v2 = IRL::Vec3<double>::fromRawDoublePointer(
      getExactVelocity(a_initial_pt + IRL::Pt::fromVec3(0.5 * a_dt * v1),
                       a_time + 0.5 * a_dt)
          .data());
  auto v3 = IRL::Vec3<double>::fromRawDoublePointer(
      getExactVelocity(a_initial_pt + IRL::Pt::fromVec3(0.5 * a_dt * v2),
                       a_time + 0.5 * a_dt)
          .data());
  auto v4 = IRL::Vec3<double>::fromRawDoublePointer(
      getExactVelocity(a_initial_pt + IRL::Pt::fromVec3(a_dt * v3),
                       a_time + a_dt)
          .data());
  return a_initial_pt +
         IRL::Pt::fromVec3(a_dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0);
}

inline IRL::Vec3<double> getVelocity(const IRL::Pt& a_location,
                                     const Data<double>& a_U,
                                     const Data<double>& a_V,
                                     const Data<double>& a_W) {
  return IRL::Vec3<double>(a_U.interpolate(a_location),
                           a_V.interpolate(a_location),
                           a_W.interpolate(a_location));
}

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_
