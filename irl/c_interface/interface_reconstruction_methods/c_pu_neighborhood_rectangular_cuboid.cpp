#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"

#include <cassert>

extern "C" {

void c_PUNeigh_RectCub_new(c_PUNeigh_RectCub* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PUNeighborhood<IRL::RectangularCuboid>;
}

void c_PUNeigh_RectCub_delete(c_PUNeigh_RectCub* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_PUNeigh_RectCub_setSize(c_PUNeigh_RectCub* a_self, const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->resize(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_PUNeigh_RectCub_reserve(c_PUNeigh_RectCub* a_self, const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->reserve(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_PUNeigh_RectCub_setMember(c_PUNeigh_RectCub* a_self, const int* a_index,
                                 const double* __restrict__ a_centroid,
                                 const double* a_weight,
                                 const c_SeparatorVariant* a_separator,
                                 const double* a_scalar) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_centroid != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
  IRL::Pt centroid = IRL::Pt::fromRawDoublePointer(a_centroid);
  a_self->obj_ptr->setMember(static_cast<IRL::UnsignedIndex_t>(*a_index),
                             &centroid, a_separator->obj_ptr, *a_weight,
                             *a_scalar);
}

void c_PUNeigh_RectCub_addMember(c_PUNeigh_RectCub* a_self,
                                 const double* __restrict__ a_centroid,
                                 const double* a_weight,
                                 const c_SeparatorVariant* a_separator,
                                 const double* a_scalar) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_centroid != nullptr);

  IRL::Pt centroid = IRL::Pt::fromRawDoublePointer(a_centroid);
  a_self->obj_ptr->addMember(&centroid, a_separator->obj_ptr, *a_weight,
                             *a_scalar);
}

void c_PUNeigh_RectCub_emptyNeighborhood(c_PUNeigh_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->emptyNeighborhood();
}

void c_PUNeigh_RectCub_setCenterOfStencil(c_PUNeigh_RectCub* a_self,
                                          const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  a_self->obj_ptr->setCenterOfStencil(
      static_cast<IRL::UnsignedIndex_t>(*a_index));
}
}