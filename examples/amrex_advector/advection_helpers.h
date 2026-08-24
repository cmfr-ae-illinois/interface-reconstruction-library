// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_H_
#define EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_H_

#include "irl/amrex/sepunion_multifab.h"

#include "examples/amrex_advector/case_deformation_3d.h"
#include "examples/amrex_advector/case_rotation_3d.h"
#include "examples/amrex_advector/case_translation_3d.h"

using namespace amrex;

enum class VelocityFieldType {
  Interpolated,
  Rotation,
  Translation,
  Deformation
};

// general velocity accessor (vx, vy, vz is only used if analytical velocity is
// not used)
IRL::Vec3<double> GetVelocity(const IRL::Pt& pt, const double time,
                              const VelocityFieldType velocity_field_type,
                              Array4<Real const> const& vx,
                              Array4<Real const> const& vy,
                              Array4<Real const> const& vz, const Box& bx,
                              const Geometry& a_geom);

inline IRL::Vec3<double> GetInterpolatedVelocity(const IRL::Pt& pt,
                                                 Array4<Real const> const& vx,
                                                 Array4<Real const> const& vy,
                                                 Array4<Real const> const& vz,
                                                 const Box& bx,
                                                 const Geometry& a_geom);

Eigen::Vector3d GetVelocity(const Eigen::Vector3d& pt, const double time,
                            const VelocityFieldType velocity_field_type,
                            Array4<Real const> const& vx,
                            Array4<Real const> const& vy,
                            Array4<Real const> const& vz, const Box& bx,
                            const Geometry& a_geom);

Eigen::Matrix3d GetVelocityGradient(const Eigen::Vector3d& pt,
                                    const double time,
                                    const VelocityFieldType velocity_field_type,
                                    Array4<Real const> const& vx,
                                    Array4<Real const> const& vy,
                                    Array4<Real const> const& vz, const Box& bx,
                                    const Geometry& a_geom);

inline IRL::Pt ProjectVertex(const IRL::Pt& pt, const double dt,
                             const double time,
                             const VelocityFieldType velocity_field_type,
                             Array4<Real const> const& vx,
                             Array4<Real const> const& vy,
                             Array4<Real const> const& vz, const Box& bx,
                             const Geometry& a_geom);

#include "examples/amrex_advector/advection_helpers.tpp"

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_H_
