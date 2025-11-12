
#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_SOLVE_RECTANGULAR_CUBOID_SEPARATOR_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_SOLVE_RECTANGULAR_CUBOID_SEPARATOR_H_

#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_solve_rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/interface_reconstruction_methods/pu_solve.h"
#include "irl/variant_reconstruction/separator_variant.h"
extern "C" {

struct c_PUST_RectCub {
  IRL::PUST<IRL::RectangularCuboid>* obj_ptr = nullptr;
};

void c_PUST_RectCub_new(c_PUST_RectCub* a_self);

void c_PUST_RectCub_delete(c_PUST_RectCub* a_self);

void c_PUST_RectCub_setNeighborhood(c_PUST_RectCub* a_self,
                                    c_PUSTNeigh_RectCub* a_neighborhood);

void c_PUST_RectCub_setTreshold(c_PUST_RectCub* a_self, double* a_tresh);

void c_PUST_RectCub_solveEdge(c_PUST_RectCub* a_self, double* STCoeff,
                              double* P0, double* P1, double* delta,
                              double* Pressure, double* Marangoni,
                              double* a_force);

void c_PUST_RectCub_getValue(c_PUST_RectCub* a_self, double* x, double* y,
                             double* z, double* delta, double* value);

void c_PUST_RectCub_getTangent(c_PUST_RectCub* a_self, double* x, double* y,
                               double* z, double* delta, double* tangent);

void c_PUST_RectCub_getWeight(c_PUST_RectCub* a_self, double* x, double* y,
                              double* z, double* delta, double* value);

// Cylinder Versions
void c_PUST_RectCub_solveEdgeCylinder(c_PUST_RectCub* a_self, double* STCoeff,
                                      double* P0, double* P1, double* radius,
                                      double* center, double* delta,
                                      double* a_force);

void c_PUST_RectCub_getValueCylinder(c_PUST_RectCub* a_self, double* x,
                                     double* y, double* z, double* radius,
                                     double* center, double* value);

void c_PUST_RectCub_getTangentCylinder(c_PUST_RectCub* a_self, double* x,
                                       double* y, double* z, double* radius,
                                       double* center, double* tangent);

void c_PUST_RectCub_getWeightCylinder(c_PUST_RectCub* a_self, double* x,
                                      double* y, double* z, double* radius,
                                      double* center, double* value);

// Debug
void c_PUST_RectCub_printSolver(c_PUST_RectCub* a_self);
}

#endif