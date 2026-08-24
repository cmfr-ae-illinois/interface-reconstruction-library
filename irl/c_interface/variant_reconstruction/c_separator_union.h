// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_SEPARATOR_UNION_H_
#define IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_SEPARATOR_UNION_H_

#include "irl/variant_reconstruction/separator_union.h"

extern "C" {

void c_SeparatorUnion_setToOnePlane_raw(IRL::SeparatorUnion& a_self,
                                        const double* a_normal,
                                        const double* a_distance);

void c_SeparatorUnion_setToFull_raw(IRL::SeparatorUnion& a_self);

void c_SeparatorUnion_setToEmpty_raw(IRL::SeparatorUnion& a_self);

void c_SeparatorUnion_getPlane_raw(const IRL::SeparatorUnion& a_self,
                                   const int* a_index, double* a_plane_listed);

void c_SeparatorUnion_reflect_raw(IRL::SeparatorUnion& a_self,
                                  const IRL::SeparatorUnion& a_reflected_sep,
                                  const int* a_dir, const double* a_loc);

bool c_SeparatorUnion_isOnePlane_raw(const IRL::SeparatorUnion& a_self);

bool c_SeparatorUnion_isParaboloid_raw(const IRL::SeparatorUnion& a_self);

bool c_SeparatorUnion_isTypeDefined_raw(const IRL::SeparatorUnion& a_self);

bool c_SeparatorUnion_isFull_raw(const IRL::SeparatorUnion& a_self);

bool c_SeparatorUnion_isEmpty_raw(const IRL::SeparatorUnion& a_self);

double c_SeparatorUnion_getMeanCurvature_raw(IRL::SeparatorUnion& a_self);

}  // end extern C

#endif  // IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_SEPARATOR_UNION_H_
