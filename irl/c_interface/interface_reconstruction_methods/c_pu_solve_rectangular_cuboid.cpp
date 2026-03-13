#include "irl/c_interface/interface_reconstruction_methods/c_pu_solve_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"
#include "irl/geometry/general/normal.h"
#include "irl/parameters/defined_types.h"

#include <cassert>

extern "C" {

void c_PU_RectCub_new(c_PU_RectCub* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PU<IRL::RectangularCuboid>;
}

void c_PU_RectCub_delete(c_PU_RectCub* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_PU_RectCub_setNeighborhood(c_PU_RectCub* a_self,
                                  c_PUNeigh_RectCub* a_neighborhood) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->setNeighborhood(*a_neighborhood->obj_ptr);
}

void c_PU_RectCub_setThreshold(c_PU_RectCub* a_self, double* a_threshold) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->setThreshold(*a_threshold);
}

void c_PU_RectCub_solveFace(c_PU_RectCub* a_self, double* STCoeff, double* P0,
                            double* P1, double* P2, double* P3, double* delta,
                            double* Pressure, double* Marangoni,
                            double* a_force) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);
  IRL::Pt P1temp = IRL::Pt::fromRawDoublePointer(P1);
  IRL::Pt P1temp = IRL::Pt::fromRawDoublePointer(P2);
  IRL::Pt P1temp = IRL::Pt::fromRawDoublePointer(P3);

  IRL::Normal MarangoniTemp = IRL::Normal::fromRawDoublePointer(Marangoni);

  IRL::Normal force = a_self->obj_ptr->solveEdge(
      *STCoeff, P0temp, P1temp, *delta, *Pressure, MarangoniTemp);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(a_force + n) = force[n];
  }
}

void c_PU_RectCub_solveFace(c_PU_RectCub* a_self, double* STCoeff, double* P0,
                            double* P1, double* delta, double* Pressure,
                            double* Marangoni, double* a_force) {
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

void c_PU_RectCub_getValue(c_PU_RectCub* a_self, double* x, double* y,
                           double* z, double* delta, double* value) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  *value = a_self->obj_ptr->getValue(*x, *y, *z, *delta);
}

void c_PU_RectCub_getTangent(c_PU_RectCub* a_self, double* x, double* y,
                             double* z, double* delta, double* tangent) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Normal T = a_self->obj_ptr->getTangent(*x, *y, *z, *delta);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(tangent + n) = T[n];
  }
}

void c_PU_RectCub_getWeight(c_PU_RectCub* a_self, double* x, double* y,
                            double* z, double* delta, double* weight) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  *weight = a_self->obj_ptr->getWeight(*x, *y, *z, *delta);
}

// Cylinder Versions
void c_PU_RectCub_solveEdgeCylinder(c_PU_RectCub* a_self, double* STCoeff,
                                    double* P0, double* P1, double* radius,
                                    double* center, double* delta,
                                    double* a_force) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);
  IRL::Pt P1temp = IRL::Pt::fromRawDoublePointer(P1);

  IRL::Pt CenterTemp = IRL::Pt::fromRawDoublePointer(center);
  IRL::Normal force = a_self->obj_ptr->solveEdgeCylinder(
      *STCoeff, P0temp, P1temp, *radius, CenterTemp, *delta);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(a_force + n) = force[n];
  }
}

void c_PU_RectCub_getValueCylinder(c_PU_RectCub* a_self, double* x, double* y,
                                   double* z, double* radius, double* center,
                                   double* value) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt C = IRL::Pt::fromRawDoublePointer(center);

  *value = a_self->obj_ptr->getValueCylinder(*x, *y, *z, *radius, C);
}

void c_PU_RectCub_getTangentCylinder(c_PU_RectCub* a_self, double* x, double* y,
                                     double* z, double* radius, double* center,
                                     double* tangent) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt C = {*(center), *(center + 1), *(center + 2)};
  IRL::Normal T = a_self->obj_ptr->getTangentCylinder(*x, *y, *z, *radius, C);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(tangent + n) = T[n];
  }
}

void c_PU_RectCub_getWeightCylinder(c_PU_RectCub* a_self, double* x, double* y,
                                    double* z, double* radius, double* center,
                                    double* weight) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt C = {*(center), *(center + 1), *(center + 2)};
  *weight = a_self->obj_ptr->getWeightCylinder(*x, *y, *z, *radius, C);
}
// Debug
void c_PU_RectCub_printSolver(c_PU_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->printSolver();
}
}