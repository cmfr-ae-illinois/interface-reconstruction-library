#include "irl/c_interface/interface_reconstruction_methods/c_pust_rectangular_cuboid.h"
#include <cassert>
#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"
#include "irl/geometry/general/normal.h"
#include "irl/parameters/defined_types.h"

extern "C" {

void c_PUST_RectCub_new(c_PUST_RectCub* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PUST<IRL::RectangularCuboid>;
}

void c_PUST_RectCub_delete(c_PUST_RectCub* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_PUST_RectCub_setNeighborhood(c_PUST_RectCub* a_self,
                                    c_PUNeigh_RectCub* a_neighborhood) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->setNeighborhood(*a_neighborhood->obj_ptr);
}

void c_PUST_RectCub_setKernelSize(c_PUST_RectCub* a_self, double* kernel_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->setKernelSize(*kernel_size);
}

void c_PUST_RectCub_solveFace(c_PUST_RectCub* a_self, double* STCoeff,
                              double* P0, double* P1, double* P2, double* P3,
                              double* delta, double* Pressure,
                              double* Marangoni, double* a_force) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);
  IRL::Pt P1temp = IRL::Pt::fromRawDoublePointer(P1);
  IRL::Pt P2temp = IRL::Pt::fromRawDoublePointer(P2);
  IRL::Pt P3temp = IRL::Pt::fromRawDoublePointer(P3);

  IRL::Normal MarangoniTemp = IRL::Normal::fromRawDoublePointer(Marangoni);

  IRL::Normal force =
      a_self->obj_ptr->solveFace(*STCoeff, P0temp, P1temp, P2temp, P3temp,
                                 *delta, *Pressure, MarangoniTemp);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(a_force + n) = force[n];
  }
}

void c_PUST_RectCub_solveEdge(c_PUST_RectCub* a_self, double* STCoeff,
                              double* P0, double* P1, double* delta,
                              double* Pressure, double* Marangoni,
                              double* a_force) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);
  IRL::Pt P1temp = IRL::Pt::fromRawDoublePointer(P1);

  IRL::Normal MarangoniTemp = IRL::Normal::fromRawDoublePointer(Marangoni);

  IRL::Normal force = a_self->obj_ptr->solveEdge(
      *STCoeff, P0temp, P1temp, *delta, *Pressure, MarangoniTemp);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(a_force + n) = force[n];
  }
}

void c_PUST_RectCub_getValue(c_PUST_RectCub* a_self, double* x, double* y,
                             double* z, double* delta, double* value) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  IRL::Pt P = {*x, *y, *z};
  *value = a_self->obj_ptr->getPU(P);
}

void c_PUST_RectCub_getWeight(c_PUST_RectCub* a_self, double* x, double* y,
                              double* z, double* delta, double* weight) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  IRL::Pt P = {*x, *y, *z};
  *weight = a_self->obj_ptr->getTotalWeight(P);
}

void c_PUST_RectCub_projectToPU(c_PUST_RectCub* a_self, double* P0, double* dx,
                                bool* success, double* Pout) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);

  IRL::Pt Pouttemp = a_self->obj_ptr->projectOntoPU(P0temp, *dx, *success);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(Pout + n) = Pouttemp[n];
  }
}

void c_PUST_RectCub_getMeanCurvature(c_PUST_RectCub* a_self, double* x,
                                     double* y, double* z, double* delta,
                                     double* value) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt Ptemp = {*(x), *(y), *(z)};
  double temp = a_self->obj_ptr->getMeanCurvature(Ptemp);
  *value = temp;
}

// Get Normal for Both
void c_PUST_RectCub_getNormal(c_PUST_RectCub* a_self, double* x, double* y,
                              double* z, double* delta, double* normal) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  IRL::Pt Ptemp = {*(x), *(y), *(z)};
  IRL::Normal N = a_self->obj_ptr->getNormal(Ptemp);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(normal + n) = N[n];
  }
}

// Debug
void c_PUST_RectCub_printSolver(c_PUST_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->printSolver();
}
}