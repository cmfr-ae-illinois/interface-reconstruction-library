#include "irl/c_interface/interface_reconstruction_methods/c_pu_solve_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"
#include "irl/geometry/general/normal.h"
#include "irl/parameters/defined_types.h"

#include <cassert>

extern "C" {

void c_PUST_OLD_RectCub_new(c_PUST_OLD_RectCub* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PUST_OLD<IRL::RectangularCuboid>;
}

void c_PUST_OLD_RectCub_delete(c_PUST_OLD_RectCub* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_PUST_OLD_RectCub_setNeighborhood(c_PUST_OLD_RectCub* a_self,
                                        c_PUNeigh_RectCub* a_neighborhood) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->setNeighborhood(*a_neighborhood->obj_ptr);
}

void c_PUST_OLD_RectCub_setThreshold(c_PUST_OLD_RectCub* a_self,
                                     double* a_threshold) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->setThreshold(*a_threshold);
}

void c_PUST_OLD_RectCub_solveFace(c_PUST_OLD_RectCub* a_self, double* STCoeff,
                                  double* P0, double* P1, double* P2,
                                  double* P3, double* delta, double* Pressure,
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

void c_PUST_OLD_RectCub_solveFaceEllipsoid(
    c_PUST_OLD_RectCub* a_self, double* STCoeff, double* P0, double* P1,
    double* P2, double* P3, double* column1, double* column2, double* column3,
    double* center, double* Pressure, double* Marangoni, double* a_force) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);
  IRL::Pt P1temp = IRL::Pt::fromRawDoublePointer(P1);
  IRL::Pt P2temp = IRL::Pt::fromRawDoublePointer(P2);
  IRL::Pt P3temp = IRL::Pt::fromRawDoublePointer(P3);
  IRL::Normal column1Temp = IRL::Normal::fromRawDoublePointer(column1);
  IRL::Normal column2Temp = IRL::Normal::fromRawDoublePointer(column2);
  IRL::Normal column3Temp = IRL::Normal::fromRawDoublePointer(column3);
  IRL::Pt centerTemp = IRL::Pt::fromRawDoublePointer(center);
  IRL::Normal MarangoniTemp = IRL::Normal::fromRawDoublePointer(Marangoni);

  IRL::Normal force = a_self->obj_ptr->solveFaceEllipsoid(
      *STCoeff, P0temp, P1temp, P2temp, P3temp, column1Temp, column2Temp,
      column3Temp, centerTemp, *Pressure, MarangoniTemp);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(a_force + n) = force[n];
  }
}

void c_PUST_OLD_RectCub_solveEdge(c_PUST_OLD_RectCub* a_self, double* STCoeff,
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

void c_PUST_OLD_RectCub_getValue(c_PUST_OLD_RectCub* a_self, double* x,
                                 double* y, double* z, double* delta,
                                 double* value) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  *value = a_self->obj_ptr->getValue(*x, *y, *z, *delta);
}

void c_PUST_OLD_RectCub_getTangent(c_PUST_OLD_RectCub* a_self, double* x,
                                   double* y, double* z, double* delta,
                                   double* tangent) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Normal T = a_self->obj_ptr->getTangent(*x, *y, *z, *delta);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(tangent + n) = T[n];
  }
}

void c_PUST_OLD_RectCub_getWeight(c_PUST_OLD_RectCub* a_self, double* x,
                                  double* y, double* z, double* delta,
                                  double* weight) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  *weight = a_self->obj_ptr->getWeight(*x, *y, *z, *delta);
}

// Cylinder Versions
void c_PUST_OLD_RectCub_solveEdgeCylinder(c_PUST_OLD_RectCub* a_self,
                                          double* STCoeff, double* P0,
                                          double* P1, double* radius,
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

void c_PUST_OLD_RectCub_getValueCylinder(c_PUST_OLD_RectCub* a_self, double* x,
                                         double* y, double* z, double* radius,
                                         double* center, double* value) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt C = IRL::Pt::fromRawDoublePointer(center);

  *value = a_self->obj_ptr->getValueCylinder(*x, *y, *z, *radius, C);
}

void c_PUST_OLD_RectCub_getTangentCylinder(c_PUST_OLD_RectCub* a_self,
                                           double* x, double* y, double* z,
                                           double* radius, double* center,
                                           double* tangent) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt C = {*(center), *(center + 1), *(center + 2)};
  IRL::Normal T = a_self->obj_ptr->getTangentCylinder(*x, *y, *z, *radius, C);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(tangent + n) = T[n];
  }
}

void c_PUST_OLD_RectCub_getWeightCylinder(c_PUST_OLD_RectCub* a_self, double* x,
                                          double* y, double* z, double* radius,
                                          double* center, double* weight) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt C = {*(center), *(center + 1), *(center + 2)};
  *weight = a_self->obj_ptr->getWeightCylinder(*x, *y, *z, *radius, C);
}

void c_PUST_OLD_RectCub_projectToPU(c_PUST_OLD_RectCub* a_self, double* P0,
                                    double* delta, double* Pout) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);
  IRL::Pt Pouttemp = a_self->obj_ptr->projectOntoPU(P0temp, *delta);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(Pout + n) = Pouttemp[n];
  }
}

void c_PUST_OLD_RectCub_projectToEllipsoid(c_PUST_OLD_RectCub* a_self,
                                           double* P0, double* column1,
                                           double* column2, double* column3,
                                           double* center, double* Pout) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt P0temp = IRL::Pt::fromRawDoublePointer(P0);
  IRL::Normal col1 = IRL::Normal::fromRawDoublePointer(column1);
  IRL::Normal col2 = IRL::Normal::fromRawDoublePointer(column2);
  IRL::Normal col3 = IRL::Normal::fromRawDoublePointer(column3);
  IRL::Pt centertemp = IRL::Pt::fromRawDoublePointer(center);

  IRL::Pt Pouttemp = a_self->obj_ptr->projectOntoEllipsoid(P0temp, col1, col2,
                                                           col3, centertemp);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(Pout + n) = Pouttemp[n];
  }
}

void c_PUST_OLD_RectCub_getCurvature(c_PUST_OLD_RectCub* a_self, double* x,
                                     double* y, double* z, double* delta,
                                     double* value) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt Ptemp = {*(x), *(y), *(z)};
  double temp = a_self->obj_ptr->getCurvature(Ptemp, *delta);
  *value = temp;
}

void c_PUST_OLD_RectCub_getMeanCurvatureEllipsoid(
    c_PUST_OLD_RectCub* a_self, double* x, double* y, double* z,
    double* column1, double* column2, double* column3, double* center,
    double* value) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Pt Ptemp = {*(x), *(y), *(z)};
  IRL::Normal col1 = IRL::Normal::fromRawDoublePointer(column1);
  IRL::Normal col2 = IRL::Normal::fromRawDoublePointer(column2);
  IRL::Normal col3 = IRL::Normal::fromRawDoublePointer(column3);
  IRL::Pt centertemp = IRL::Pt::fromRawDoublePointer(center);
  double temp = a_self->obj_ptr->getMeanCurvatureEllipsoid(Ptemp, col1, col2,
                                                           col3, centertemp);
  *value = temp;
}

// Get Normal for Both
void c_PUST_OLD_RectCub_getNormal(c_PUST_OLD_RectCub* a_self, double* x,
                                  double* y, double* z, double* delta,
                                  double* normal) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Normal N = a_self->obj_ptr->getNormal(*x, *y, *z, *delta);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(normal + n) = N[n];
  }
}

void c_PUST_OLD_RectCub_getNormalEllipsoid(c_PUST_OLD_RectCub* a_self,
                                           double* x, double* y, double* z,
                                           double* column1, double* column2,
                                           double* column3, double* center,
                                           double* normal) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  IRL::Normal col1 = IRL::Normal::fromRawDoublePointer(column1);
  IRL::Normal col2 = IRL::Normal::fromRawDoublePointer(column2);
  IRL::Normal col3 = IRL::Normal::fromRawDoublePointer(column3);
  IRL::Pt centertemp = IRL::Pt::fromRawDoublePointer(center);

  IRL::Normal N = a_self->obj_ptr->getNormalEllipsoid(*x, *y, *z, col1, col2,
                                                      col3, centertemp);

  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    *(normal + n) = N[n];
  }
}
// Debug
void c_PUST_OLD_RectCub_printSolver(c_PUST_OLD_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  a_self->obj_ptr->printSolver();
}
}