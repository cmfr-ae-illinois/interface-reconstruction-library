// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_SURFACE_ADVECTOR_VOF_ADVECTION_H_
#define EXAMPLES_SURFACE_ADVECTOR_VOF_ADVECTION_H_

#include <string>

#include "examples/surface_advector/data.h"
#include "examples/surface_advector/reconstruction_types.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron.h"
#include "irl/geometry/polyhedrons/polyhedron_24.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"

void resetMoments(
    const Data<IRL::LocalizedSeparatorVariantLink>& a_link_localized_interface,
    Data<IRL::VolumeMoments>* a_liq_moments,
    Data<IRL::VolumeMoments>* a_gas_moments);

std::array<int, 3> getIndexFromTag(const BasicMesh& a_mesh,
                                   const IRL::UnsignedIndex_t a_tag);

void advectVOF(
    const std::string& a_advection_method,
    const std::string& a_reconstruction_method, const double a_dt,
    const double a_time, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::SeparatorVariant>* a_interface,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface,
    Data<IRL::VolumeMoments>* a_liq_moments,
    Data<IRL::VolumeMoments>* a_gas_moments, Data<double>* a_surfactant_mass);

struct FullLagrangianCorrected {
  static void advectVOF(
      const std::string& a_reconstruction_method, const double a_dt,
      const double a_time, const Data<double>& a_U, const Data<double>& a_V,
      const Data<double>& a_W, Data<IRL::SeparatorVariant>* a_interface,
      Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface,
      Data<IRL::VolumeMoments>* a_liq_moments,
      Data<IRL::VolumeMoments>* a_gas_moments, Data<double>* a_surfactant_mass);
};

inline IRL::Vec3<double> getVelocity(const IRL::Pt& a_location,
                                     const Data<double>& a_U,
                                     const Data<double>& a_V,
                                     const Data<double>& a_W);

inline IRL::Pt project_vertex(const IRL::Pt& a_initial_pt, const double a_dt,
                              const Data<double>& a_U, const Data<double>& a_V,
                              const Data<double>& a_W);

void correctMomentLocations(Data<IRL::VolumeMoments>* a_liq_moments,
                            Data<IRL::VolumeMoments>* a_gas_moments);

// ************************************************
//     Inlined functions below this
// ************************************************
// Interpolate velocity
inline IRL::Vec3<double> getVelocity(const IRL::Pt& a_location,
                                     const Data<double>& a_U,
                                     const Data<double>& a_V,
                                     const Data<double>& a_W) {
  return IRL::Vec3<double>(a_U.interpolate(a_location),
                           a_V.interpolate(a_location),
                           a_W.interpolate(a_location));
}
// Spatial RK4 projection of a point in a velocity field..
inline IRL::Pt project_vertex(const IRL::Pt& a_initial_pt, const double a_dt,
                              const Data<double>& a_U, const Data<double>& a_V,
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

// Lookup tables for construction of flux-corrected Poly24
static const std::array<std::array<int, 4>, 6> flux_id_table = {{{4, 5, 6, 7},
                                                                 {0, 1, 2, 3},
                                                                 {1, 5, 4, 0},
                                                                 {2, 6, 7, 3},
                                                                 {6, 5, 1, 2},
                                                                 {7, 4, 0, 3}}};
static const std::array<int, 6> face_center_id_table = {{13, 8, 9, 11, 10, 12}};

#endif  // EXAMPLES_SURFACE_ADVECTOR_VOF_ADVECTION_H_
