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

void c_SeparatorUnion_new(c_SeparatorUnion* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::SeparatorUnion;
}

void c_SeparatorUnion_delete(c_SeparatorUnion* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_SeparatorUnion_setNumberOfPlanes(c_SeparatorUnion* a_self,
                                        const int* a_number_to_set) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_number_to_set >= 0);
  assert(*a_number_to_set <= 2);
  if (*a_number_to_set == 2) {
    a_self->obj_ptr->type() = IRL::SeparatorUnion::SeparatorType::TwoPlanes;
  } else {
    a_self->obj_ptr->type() = IRL::SeparatorUnion::SeparatorType::OnePlane;
  }
}

void c_SeparatorUnion_setPlane(c_SeparatorUnion* a_self,
                               const int* a_plane_index_to_set,
                               const double* a_normal,
                               const double* a_distance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_plane_index_to_set >= 0);
  assert(*a_plane_index_to_set <= 1);
  a_self->obj_ptr->setPlane(
      *a_plane_index_to_set,
      IRL::Plane(IRL::Normal::fromRawDoublePointer(a_normal), *a_distance));
}

void c_SeparatorUnion_copy(c_SeparatorUnion* a_self,
                           const c_SeparatorUnion* a_other_planar_separator) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_other_planar_separator != nullptr);
  assert(a_other_planar_separator->obj_ptr != nullptr);
  (*a_self->obj_ptr) = (*a_other_planar_separator->obj_ptr);
}

int c_SeparatorUnion_getNumberOfPlanes(const c_SeparatorUnion* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  if (a_self->obj_ptr->type() == IRL::SeparatorUnion::SeparatorType::OnePlane) {
    return 1;
  } else if (a_self->obj_ptr->type() ==
             IRL::SeparatorUnion::SeparatorType::TwoPlanes) {
    return 2;
  } else {
    return 0;
  }
}

void c_SeparatorUnion_getPlane(c_SeparatorUnion* a_self, const int* a_index,
                               double* a_plane_listed) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(*a_index <= 1);
  const IRL::Plane& plane =
      a_self->obj_ptr->getPlane(static_cast<IRL::UnsignedIndex_t>(*a_index));
  a_plane_listed[0] = plane.normal()[0];
  a_plane_listed[1] = plane.normal()[1];
  a_plane_listed[2] = plane.normal()[2];
  a_plane_listed[3] = plane.distance();
}

bool c_SeparatorUnion_isFlipped(const c_SeparatorUnion* a_self) {
  assert(a_self != nullptr);
  return false;
}

void c_SeparatorUnion_printToScreen(const c_SeparatorUnion* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << *(a_self->obj_ptr);
}

void c_SeparatorUnion_shift(c_SeparatorUnion* a_self, const double* a_shift) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_shift != nullptr);
  a_self->obj_ptr->shift(IRL::Pt::fromRawDoublePointer(a_shift));
}

}  // end extern C
