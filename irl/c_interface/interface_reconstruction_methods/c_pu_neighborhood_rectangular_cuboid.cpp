#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"

#include <cassert>

extern "C" {

void c_PUSTNeigh_RectCub_new(c_PUSTNeigh_RectCub* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PUSTNeighborhood<IRL::RectangularCuboid>;
}

void c_PUSTNeigh_RectCub_delete(c_PUSTNeigh_RectCub* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_PUSTNeigh_RectCub_setSize(c_PUSTNeigh_RectCub* a_self,
                                 const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->!= nullptr);
  a_self->obj_ptr->resize(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_PUSTNeigh_RectCub_reserve(c_PUSTNeigh_RectCub* a_self,
                                 const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->!= nullptr);
  a_self->obj_ptr->reserve(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_PUSTNeigh_RectCub_setMember(c_PUSTNeigh_RectCub* a_self,
                                   const int* a_index,
                                   const double* __restrict__ a_centroid,
                                   const c_SeparatorVariant* a_separator) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_centroid != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
  IRL::Pt centroid = IRL::Pt::fromRawDoublePointer(a_centroid);
  a_self->obj_ptr->setMember(static_cast<IRL::UnsignedIndex_t>(*a_index),
                             &centroid, a_separator->obj_ptr);
}

void c_PUSTNeigh_RectCub_addMember(c_PUSTNeigh_RectCub* a_self,
                                   const double* __restrict__ a_centroid,
                                   const c_SeparatorVariant* a_separator) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_centroid != nullptr);

  IRL::Pt centroid = IRL::Pt::fromRawDoublePointer(a_centroid);
  a_self->obj_ptr->addMember(&centroid, a_separator->obj_ptr);
}

void c_PUSTNeigh_RectCub_emptyNeighborhood(c_PUSTNeigh_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->emptyNeighborhood();
}
}