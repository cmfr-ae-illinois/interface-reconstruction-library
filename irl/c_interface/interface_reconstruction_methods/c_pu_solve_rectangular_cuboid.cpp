#include "irl/c_interface/interface_reconstruction_methods/c_pu_solve_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"
#include "irl/geometry/general/normal.h"
#include "irl/parameters/defined_types.h"

#include <cassert>

extern "C" {

void c_PUST_RectCub_new(c_PUST_RectCub* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PUST<IRL::RectangularCuboid>;
}

void c_PUST_RectCub_delete(c_PUST_RectCub* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_PUST_RectCub_solveEdge(c_PUST_RectCub* a_self, double STCoeff,
                              double* P0, double* P1, double* a_force) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);
  IRL::Pt P1temp = IRL::Pt::fromRawDoublePointer(P1);

  IRL::Normal force = a_self->obj_ptr->solveEdge(STCoeff, P0temp, P1temp);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    a_force[n] = force[n];
  }
}
}