
#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_SOLVE_RECTANGULAR_CUBOID_SEPARATOR_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_SOLVE_RECTANGULAR_CUBOID_SEPARATOR_H_

#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_solve_rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/interface_reconstruction_methods/pu_solve.h"
#include "irl/variant_reconstruction/separator_variant.h"
extern "C" {

struct c_PU_RectCub {
  IRL::PU<IRL::RectangularCuboid>* obj_ptr = nullptr;
};

void c_PU_RectCub_new(c_PU_RectCub* a_self);

void c_PU_RectCub_delete(c_PU_RectCub* a_self);

void c_PU_RectCub_setNeighborhood(c_PU_RectCub* a_self,
                                  c_PUNeigh_RectCub* a_neighborhood);

void c_PU_RectCub_setTreshold(c_PU_RectCub* a_self, double* a_tresh);

void c_PU_RectCub_solveEdge(c_PU_RectCub* a_self, double* STCoeff, double* P0,
                            double* P1, double* delta, double* Pressure,
                            double* Marangoni, double* a_force);

void c_PU_RectCub_solveFace(c_PU_RectCub* a_self, double* STCoeff, double* P0,
                            double* P1, double* P2, double* P3, double* delta,
                            double* Pressure, double* Marangoni,
                            double* a_force);

void c_PU_RectCub_getValue(c_PU_RectCub* a_self, double* x, double* y,
                           double* z, double* delta, double* value);

void c_PU_RectCub_getTangent(c_PU_RectCub* a_self, double* x, double* y,
                             double* z, double* delta, double* tangent);

void c_PU_RectCub_getWeight(c_PU_RectCub* a_self, double* x, double* y,
                            double* z, double* delta, double* value);

void c_PU_RectCub_getCurvature(c_PU_RectCub* a_self, double* x, double* y,
                               double* z, double* delta, double* value);

// Cylinder Versions
void c_PU_RectCub_solveEdgeCylinder(c_PU_RectCub* a_self, double* STCoeff,
                                    double* P0, double* P1, double* radius,
                                    double* center, double* delta,
                                    double* a_force);

void c_PU_RectCub_getValueCylinder(c_PU_RectCub* a_self, double* x, double* y,
                                   double* z, double* radius, double* center,
                                   double* value);

void c_PU_RectCub_getTangentCylinder(c_PU_RectCub* a_self, double* x, double* y,
                                     double* z, double* radius, double* center,
                                     double* tangent);

void c_PU_RectCub_getWeightCylinder(c_PU_RectCub* a_self, double* x, double* y,
                                    double* z, double* radius, double* center,
                                    double* value);
void c_PU_RectCub_projectToPU(c_PU_RectCub* a_self, double* P0, double* delta,
                              double* Pout);
void c_PU_RectCub_projectToEllipsoid(c_PU_RectCub* a_self, double* P0,
                                     double* column1, double* column2,
                                     double* column3, double* center,
                                     double* Pout);
void c_PU_RectCub_getMeanCurvatureEllipsoid(c_PU_RectCub* a_self, double* x,
                                            double* y, double* z,
                                            double* column1, double* column2,
                                            double* column3, double* center,
                                            double* value);
// Get Normal for Both
void c_PU_RectCub_getNormal(c_PU_RectCub* a_self, double* x, double* y,
                            double* z, double* delta, double* normal);
void c_PU_RectCub_getNormalEllipsoid(c_PU_RectCub* a_self, double* x, double* y,
                                     double* z, double* column1,
                                     double* column2, double* column3,
                                     double* center, double* normal);

// Updated Solve Face
void c_PU_RectCub_solveFaceEllipsoid(c_PU_RectCub* a_self, double* STCoeff,
                                     double* P0, double* P1, double* P2,
                                     double* P3, double* column1,
                                     double* column2, double* column3,
                                     double* center, double* Pressure,
                                     double* Marangoni, double* a_force);
// Debug
void c_PU_RectCub_printSolver(c_PU_RectCub* a_self);
}

#endif