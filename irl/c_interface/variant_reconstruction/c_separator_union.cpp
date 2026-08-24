// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/variant_reconstruction/c_separator_union.h"

#include <iostream>

#include "irl/interface_reconstruction_methods/volume_fraction_matching.h"
#include "irl/parameters/constants.h"

extern "C" {

void c_SeparatorUnion_setToOnePlane_raw(IRL::SeparatorUnion& a_self,
                                        const double* a_normal,
                                        const double* a_distance) {
  assert(a_normal != nullptr);
  assert(a_distance != nullptr);
  a_self.setToPlane(
      IRL::Plane(IRL::Normal::fromRawDoublePointer(a_normal), *a_distance));
}

void c_SeparatorUnion_setToFull_raw(IRL::SeparatorUnion& a_self) {
  a_self.setToFull();
}

void c_SeparatorUnion_setToEmpty_raw(IRL::SeparatorUnion& a_self) {
  a_self.setToEmpty();
}

void c_SeparatorUnion_getPlane_raw(const IRL::SeparatorUnion& a_self,
                                   const int* a_index, double* a_plane_listed) {
  assert(a_index != nullptr);
  assert(*a_index >= 0);
  assert(*a_index <= 1);
  assert(a_plane_listed != nullptr);
  const IRL::Plane& plane =
      a_self.getPlane(static_cast<IRL::UnsignedIndex_t>(*a_index));
  a_plane_listed[0] = plane.normal()[0];
  a_plane_listed[1] = plane.normal()[1];
  a_plane_listed[2] = plane.normal()[2];
  a_plane_listed[3] = plane.distance();
}

void c_SeparatorUnion_reflect_raw(IRL::SeparatorUnion& a_self,
                                  const IRL::SeparatorUnion& a_reflected_sep,
                                  const int* a_dir, const double* a_loc) {
  assert(a_dir != nullptr);
  assert(*a_dir >= 0);
  assert(*a_dir <= 2);
  assert(a_loc != nullptr);
  a_self.reflect(a_reflected_sep, *a_dir, *a_loc);
}

bool c_SeparatorUnion_isOnePlane_raw(const IRL::SeparatorUnion& a_self) {
  return a_self.type() == IRL::SeparatorUnion::SeparatorType::OnePlane;
}

bool c_SeparatorUnion_isParaboloid_raw(const IRL::SeparatorUnion& a_self) {
  return a_self.type() == IRL::SeparatorUnion::SeparatorType::Paraboloid;
}

bool c_SeparatorUnion_isTypeDefined_raw(const IRL::SeparatorUnion& a_self) {
  return (a_self.type() == IRL::SeparatorUnion::SeparatorType::OnePlane) ||
         (a_self.type() == IRL::SeparatorUnion::SeparatorType::TwoPlanes) ||
         (a_self.type() == IRL::SeparatorUnion::SeparatorType::Paraboloid) ||
         (a_self.type() == IRL::SeparatorUnion::SeparatorType::Cylinder);
}

bool c_SeparatorUnion_isFull_raw(const IRL::SeparatorUnion& a_self) {
  return a_self.isFull();
}

bool c_SeparatorUnion_isEmpty_raw(const IRL::SeparatorUnion& a_self) {
  return a_self.isEmpty();
}

double c_SeparatorUnion_getMeanCurvature_raw(IRL::SeparatorUnion& a_self) {
  if (a_self.type() == IRL::SeparatorUnion::SeparatorType::Paraboloid) {
    const auto& aligned_paraboloid =
        a_self.getParaboloid().getAlignedParaboloid();
    return 2.0 * aligned_paraboloid.a() + 2.0 * aligned_paraboloid.b();
  } else if (a_self.type() == IRL::SeparatorUnion::SeparatorType::Cylinder) {
    const auto& aligned_cylinder = a_self.getCylinder().getAlignedCylinder();
    return aligned_cylinder.b() / aligned_cylinder.r();
  }
  return 0.0;
}

}  // end extern C
