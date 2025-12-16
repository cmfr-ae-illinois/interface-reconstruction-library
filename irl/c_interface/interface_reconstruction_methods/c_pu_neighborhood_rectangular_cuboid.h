
#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_

#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/c_interface/variant_reconstruction/c_separator_variant.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/variant_reconstruction/separator_variant.h"

extern "C" {

struct c_PUNeigh_RectCub {
  IRL::PUNeighborhood<IRL::RectangularCuboid>* obj_ptr = nullptr;
};

void c_PUNeigh_RectCub_new(c_PUNeigh_RectCub* a_self);

void c_PUNeigh_RectCub_delete(c_PUNeigh_RectCub* a_self);

void c_PUNeigh_RectCub_setSize(c_PUNeigh_RectCub* a_self,
                                 const int* a_size);

void c_PUNeigh_RectCub_reserve(c_PUNeigh_RectCub* a_self,
                                 const int* a_size);

void c_PUNeigh_RectCub_setMember(c_PUNeigh_RectCub* a_self,
                                   const int* a_index,
                                   const double* __restrict__ a_centroid,
                                   const double* a_weight,
                                   const c_SeparatorVariant* a_separator);

void c_PUNeigh_RectCub_addMember(c_PUNeigh_RectCub* a_self,
                                   const double* __restrict__ a_centroid,
                                   const double* a_weight,
                                   const c_SeparatorVariant* a_separator);

void c_PUNeigh_RectCub_emptyNeighborhood(c_PUNeigh_RectCub* a_self);
}

#endif