
#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_

#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/c_interface/variant_reconstruction/c_separator_variant.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/variant_reconstruction/separator_variant.h"

extern "C" {

struct c_PUSTNeigh_RectCub {
  IRL::PUSTNeighborhood<IRL::RectangularCuboid>* obj_ptr = nullptr;
};

void c_PUSTNeigh_RectCub_new(c_PUSTNeigh_RectCub* a_self);

void c_PUSTNeigh_RectCub_delete(c_PUSTNeigh_RectCub* a_self);

void c_PUSTNeigh_RectCub_setSize(c_PUSTNeigh_RectCub* a_self,
                                 const int* a_size);

void c_PUSTNeigh_RectCub_setMember(c_PUSTNeigh_RectCub* a_self,
                                   const int* a_index,
                                   const double* __restrict__ a_centroid,
                                   const c_SeparatorVariant* a_separator);

void c_PUSTNeigh_RectCub_addMember(c_PUSTNeigh_RectCub* a_self,
                                   const double* __restrict__ a_centroid,
                                   const c_SeparatorVariant* a_separator);

void c_PUSTNeigh_RectCub_emptyNeighborhood(c_PUSTNeigh_RectCub* a_self);
}

#endif