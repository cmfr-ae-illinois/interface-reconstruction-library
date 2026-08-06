// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/variant_reconstruction/c_separator_variant.h"
#include <iostream>
#include "irl/c_interface/paraboloid_reconstruction/c_paraboloid.h"

#include "irl/interface_reconstruction_methods/volume_fraction_matching.h"
#include "irl/parameters/constants.h"

extern "C" {

void c_SeparatorVariant_new(c_SeparatorVariant* a_self) {
  assert(a_self != nullptr);
  // assert(a_self->obj_ptr == nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::SeparatorVariant;
}

void c_SeparatorVariant_newFromObjectAllocationServer(
    c_SeparatorVariant* a_self,
    c_ObjServer_SeparatorVariant* a_object_allocation_server) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
}

void c_SeparatorVariant_delete(c_SeparatorVariant* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_SeparatorVariant_setNumberOfPlanes(c_SeparatorVariant* a_self,
                                          const int* a_number_to_set) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_number_to_set >= 0);
  a_self->obj_ptr->setToPlanarSeparator();
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_self->obj_ptr)) {
    separator->setNumberOfPlanes(
        static_cast<IRL::UnsignedIndex_t>(*a_number_to_set));
  }
}

void c_SeparatorVariant_setPlane(c_SeparatorVariant* a_self,
                                 const int* a_plane_index_to_set,
                                 const double* a_normal,
                                 const double* a_distance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_plane_index_to_set >= 0);
  a_self->obj_ptr->setToPlanarSeparator();
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_self->obj_ptr)) {
    (*separator)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index_to_set)] =
        IRL::Plane(IRL::Normal::fromRawDoublePointer(a_normal), *a_distance);
  }
}

void c_SeparatorVariant_copy(
    c_SeparatorVariant* a_self,
    const c_SeparatorVariant* a_other_planar_separator) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_other_planar_separator != nullptr);
  assert(a_other_planar_separator->obj_ptr != nullptr);
  (*a_self->obj_ptr) = (*a_other_planar_separator->obj_ptr);
}

int c_SeparatorVariant_getNumberOfPlanes(const c_SeparatorVariant* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_self->obj_ptr)) {
    return static_cast<int>(separator->getNumberOfPlanes());
  } else {
    return 0;
  }
}

void c_SeparatorVariant_getPlane(c_SeparatorVariant* a_self, const int* a_index,
                                 double* a_plane_listed) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_self->obj_ptr)) {
    assert(static_cast<IRL::UnsignedIndex_t>(*a_index) <
           separator->getNumberOfPlanes());
    a_plane_listed[0] =
        (*separator)[static_cast<IRL::UnsignedIndex_t>(*a_index)].normal()[0];
    a_plane_listed[1] =
        (*separator)[static_cast<IRL::UnsignedIndex_t>(*a_index)].normal()[1];
    a_plane_listed[2] =
        (*separator)[static_cast<IRL::UnsignedIndex_t>(*a_index)].normal()[2];
    a_plane_listed[3] =
        (*separator)[static_cast<IRL::UnsignedIndex_t>(*a_index)].distance();
  }
}

void c_SeparatorVariant_setPrincipalCurvatures(c_SeparatorVariant* a_self,
                                               double* a_curvatures) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_curvatures != nullptr);
  a_self->obj_ptr->setPrincipalCurvatures(a_curvatures[0], a_curvatures[1]);
}

void c_SeparatorVariant_getPrincipalCurvatures(c_SeparatorVariant* a_self,
                                               double* a_curvatures) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_curvatures != nullptr);
  const auto curv_pair = a_self->obj_ptr->getPrincipalCurvatures();
  a_curvatures[0] = curv_pair.first;
  a_curvatures[1] = curv_pair.second;
}

void c_SeparatorVariant_getParaboloid(c_SeparatorVariant* a_self,
                                      double* a_paraboloid_listed) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  if (IRL::Paraboloid* separator =
          std::get_if<IRL::Paraboloid>(a_self->obj_ptr)) {
    // Datum
    a_paraboloid_listed[0] = (*separator).getDatum()[0];
    a_paraboloid_listed[1] = (*separator).getDatum()[1];
    a_paraboloid_listed[2] = (*separator).getDatum()[2];
    // std::cout << "Datum Z\n" << (*separator).getDatum()[2] << "\n";
    // Reference Frame
    a_paraboloid_listed[3] = (*separator).getReferenceFrame()[0][0];
    a_paraboloid_listed[4] = (*separator).getReferenceFrame()[0][1];
    a_paraboloid_listed[5] = (*separator).getReferenceFrame()[0][2];
    a_paraboloid_listed[6] = (*separator).getReferenceFrame()[1][0];
    a_paraboloid_listed[7] = (*separator).getReferenceFrame()[1][1];
    a_paraboloid_listed[8] = (*separator).getReferenceFrame()[1][2];
    a_paraboloid_listed[9] = (*separator).getReferenceFrame()[2][0];
    a_paraboloid_listed[10] = (*separator).getReferenceFrame()[2][1];
    a_paraboloid_listed[11] = (*separator).getReferenceFrame()[2][2];
    // Aligned Paraboloid
    a_paraboloid_listed[12] = (*separator).getAlignedParaboloid().a();
    a_paraboloid_listed[13] = (*separator).getAlignedParaboloid().b();
  }
}

void c_SeparatorVariant_getParaboloidObject(c_SeparatorVariant* a_self,
                                            c_Paraboloid* a_para) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  if (IRL::Paraboloid* separator =
          std::get_if<IRL::Paraboloid>(a_self->obj_ptr)) {
    a_para->obj_ptr = separator;
  }
}

void c_SeparatorVariant_setParaboloid(
    c_SeparatorVariant* a_self, const double* a_datum, const double* a_normal1,
    const double* a_normal2, const double* a_normal3, const double* a_coeff_a,
    const double* a_coeff_b) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->setToParaboloid();
  IRL::Pt datum = IRL::Pt::fromRawDoublePointer(a_datum);

  IRL::Normal n1 = IRL::Normal::fromRawDoublePointer(a_normal1);
  IRL::Normal n2 = IRL::Normal::fromRawDoublePointer(a_normal2);
  IRL::Normal n3 = IRL::Normal::fromRawDoublePointer(a_normal3);
  IRL::ReferenceFrame RF = IRL::ReferenceFrame(n1, n2, n3);

  if (IRL::Paraboloid* separator =
          std::get_if<IRL::Paraboloid>(a_self->obj_ptr)) {
    (*separator) = IRL::Paraboloid(datum, RF, *a_coeff_a, *a_coeff_b);
  }
}

bool c_SeparatorVariant_isFlipped(const c_SeparatorVariant* a_self) {
  assert(a_self != nullptr);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_self->obj_ptr)) {
    return separator->isFlipped();
  } else if (IRL::Paraboloid* paraboloid =
                 std::get_if<IRL::Paraboloid>(a_self->obj_ptr)) {
    return paraboloid->isFlipped();
  } else if (IRL::Cylinder* cylinder =
                 std::get_if<IRL::Cylinder>(a_self->obj_ptr)) {
    return cylinder->isFlipped();
  } else {
    throw std::runtime_error("Variant type unknown");
  }
}

bool c_SeparatorVariant_isPlane(const c_SeparatorVariant* a_self) {
  assert(a_self != nullptr);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_self->obj_ptr)) {
    return true;
  } else {
    return false;
  }
}

bool c_SeparatorVariant_isParaboloid(const c_SeparatorVariant* a_self) {
  assert(a_self != nullptr);
  if (IRL::Paraboloid* paraboloid =
          std::get_if<IRL::Paraboloid>(a_self->obj_ptr)) {
    return true;
  } else {
    return false;
  }
}

void c_SeparatorVariant_printToScreen(const c_SeparatorVariant* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_self->obj_ptr)) {
    std::cout << (*separator);
  } else if (IRL::Paraboloid* paraboloid =
                 std::get_if<IRL::Paraboloid>(a_self->obj_ptr)) {
    std::cout << (*paraboloid);
  } else if (IRL::Cylinder* cylinder =
                 std::get_if<IRL::Cylinder>(a_self->obj_ptr)) {
    std::cout << (*cylinder);
  } else {
    throw std::runtime_error("Variant type unknown");
  }
}

void c_SeparatorVariant_shift(c_SeparatorVariant* a_self,
                              const double* a_shift) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_shift != nullptr);
  const IRL::Pt shift = IRL::Pt::fromRawDoublePointer(a_shift);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_self->obj_ptr)) {
    for (auto& plane : *separator) {
      plane.distance() += plane.normal() * shift;
    }
  } else if (IRL::Paraboloid* paraboloid =
                 std::get_if<IRL::Paraboloid>(a_self->obj_ptr)) {
    const IRL::Pt& datum = paraboloid->getDatum();
    paraboloid->setDatum(datum + shift);
  } else if (IRL::Cylinder* cylinder =
                 std::get_if<IRL::Cylinder>(a_self->obj_ptr)) {
    const IRL::Pt& datum = cylinder->getDatum();
    cylinder->setDatum(datum + shift);
  } else {
    throw std::runtime_error("Variant type unknown");
  }
}

}  // end extern C
