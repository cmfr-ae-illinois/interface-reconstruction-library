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

#include "examples/paraboloid_advector/deformation_3d.h"
#include "examples/paraboloid_advector/rotation_3d.h"
#include "examples/paraboloid_advector/stagnation_3d.h"
#include "examples/paraboloid_advector/translation_3d.h"
#include "irl/geometry/general/pt.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

#include "examples/paraboloid_advector/data.h"

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

// struct Split {
//   static void advectVOF(
//       const std::string& a_reconstruction_method, const double a_dt,
//       const Data<double>& a_U, const Data<double>& a_V, const Data<double>&
//       a_W, Data<IRL::LocalizedParaboloidLink<double>>*
//       a_link_localized_paraboloid, Data<IRL::GeneralMoments3D<2>>*
//       a_liquid_moments, Data<IRL::GeneralMoments3D<2>>* a_gas_moments,
//       Data<IRL::Paraboloid>* a_interface);
// };

// struct FullLagrangian {
//   static void advectVOF(
//       const std::string& a_reconstruction_method, const double a_dt,
//       const Data<double>& a_U, const Data<double>& a_V, const Data<double>&
//       a_W, Data<IRL::LocalizedParaboloidLink<double>>*
//       a_link_localized_paraboloid, Data<IRL::GeneralMoments3D<2>>*
//       a_liquid_moments, Data<IRL::GeneralMoments3D<2>>* a_gas_moments);
// };

// struct SemiLagrangian {
//   static void advectVOF(
//       const std::string& a_reconstruction_method, const double a_dt,
//       const Data<double>& a_U, const Data<double>& a_V, const Data<double>&
//       a_W, Data<IRL::LocalizedParaboloidLink<double>>*
//       a_link_localized_paraboloid, Data<IRL::GeneralMoments3D<2>>*
//       a_liquid_moments, Data<IRL::GeneralMoments3D<2>>* a_gas_moments);
// };

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

inline IRL::Pt back_project_vertex(const IRL::Pt& a_loc,
                                   const std::string& a_simulation_type,
                                   const double a_dt, const double a_time) {
  const Eigen::Vector3d (*getExactVelocity)(const Eigen::Vector3d&,
                                            const double);
  if (a_simulation_type == "Deformation3D")
    getExactVelocity = Deformation3D::getExactVelocity;
  else if (a_simulation_type == "Translation3D")
    getExactVelocity = Translation3D::getExactVelocity;
  else if (a_simulation_type == "Rotation3D")
    getExactVelocity = Rotation3D::getExactVelocity;
  else if (a_simulation_type == "Stagnation3D")
    getExactVelocity = Stagnation3D::getExactVelocity;
  const Eigen::Vector3d loc{a_loc[0], a_loc[1], a_loc[2]};
  const auto v1 = getExactVelocity(loc, a_time);
  const auto v2 = getExactVelocity(loc + 0.5 * a_dt * v1, a_time + 0.5 * a_dt);
  const auto v3 = getExactVelocity(loc + 0.5 * a_dt * v2, a_time + 0.5 * a_dt);
  const auto v4 = getExactVelocity(loc + a_dt * v3, a_time + a_dt);
  const auto f = a_dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0;
  return IRL::Pt(a_loc[0] + f[0], a_loc[1] + f[1], a_loc[2] + f[2]);
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
