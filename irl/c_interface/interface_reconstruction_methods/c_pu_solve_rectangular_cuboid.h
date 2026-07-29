
#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_SOLVE_RECTANGULAR_CUBOID_SEPARATOR_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_SOLVE_RECTANGULAR_CUBOID_SEPARATOR_H_

#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood_rectangular_cuboid.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_solve_rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/interface_reconstruction_methods/pu_solve.h"
#include "irl/variant_reconstruction/separator_variant.h"
extern "C" {

struct c_PUST_OLD_RectCub {
  IRL::PUST_OLD<IRL::RectangularCuboid>* obj_ptr = nullptr;
};

void c_PUST_OLD_RectCub_new(c_PUST_OLD_RectCub* a_self);

void c_PUST_OLD_RectCub_delete(c_PUST_OLD_RectCub* a_self);

void c_PUST_OLD_RectCub_setNeighborhood(c_PUST_OLD_RectCub* a_self,
                                        c_PUNeigh_RectCub* a_neighborhood);

void c_PUST_OLD_RectCub_setTreshold(c_PUST_OLD_RectCub* a_self,
                                    double* a_tresh);

void c_PUST_OLD_RectCub_solveEdge(c_PUST_OLD_RectCub* a_self, double* STCoeff,
                                  double* P0, double* P1, double* delta,
                                  double* Pressure, double* Marangoni,
                                  double* a_force);

void c_PUST_OLD_RectCub_solveFace(c_PUST_OLD_RectCub* a_self, double* STCoeff,
                                  double* P0, double* P1, double* P2,
                                  double* P3, double* delta, double* Pressure,
                                  double* Marangoni, double* a_force);

void c_PUST_OLD_RectCub_getValue(c_PUST_OLD_RectCub* a_self, double* x,
                                 double* y, double* z, double* delta,
                                 double* value);

void c_PUST_OLD_RectCub_getTangent(c_PUST_OLD_RectCub* a_self, double* x,
                                   double* y, double* z, double* delta,
                                   double* tangent);

void c_PUST_OLD_RectCub_getWeight(c_PUST_OLD_RectCub* a_self, double* x,
                                  double* y, double* z, double* delta,
                                  double* value);

void c_PUST_OLD_RectCub_getCurvature(c_PUST_OLD_RectCub* a_self, double* x,
                                     double* y, double* z, double* delta,
                                     double* value);

// Cylinder Versions
void c_PUST_OLD_RectCub_solveEdgeCylinder(c_PUST_OLD_RectCub* a_self,
                                          double* STCoeff, double* P0,
                                          double* P1, double* radius,
                                          double* center, double* delta,
                                          double* a_force);

void c_PUST_OLD_RectCub_getValueCylinder(c_PUST_OLD_RectCub* a_self, double* x,
                                         double* y, double* z, double* radius,
                                         double* center, double* value);

void c_PUST_OLD_RectCub_getTangentCylinder(c_PUST_OLD_RectCub* a_self,
                                           double* x, double* y, double* z,
                                           double* radius, double* center,
                                           double* tangent);

void c_PUST_OLD_RectCub_getWeightCylinder(c_PUST_OLD_RectCub* a_self, double* x,
                                          double* y, double* z, double* radius,
                                          double* center, double* value);
void c_PUST_OLD_RectCub_projectToPU(c_PUST_OLD_RectCub* a_self, double* P0,
                                    double* delta, double* Pout);
void c_PUST_OLD_RectCub_projectToEllipsoid(c_PUST_OLD_RectCub* a_self,
                                           double* P0, double* column1,
                                           double* column2, double* column3,
                                           double* center, double* Pout);
void c_PUST_OLD_RectCub_getMeanCurvatureEllipsoid(
    c_PUST_OLD_RectCub* a_self, double* x, double* y, double* z,
    double* column1, double* column2, double* column3, double* center,
    double* value);
// Get Normal for Both
void c_PUST_OLD_RectCub_getNormal(c_PUST_OLD_RectCub* a_self, double* x,
                                  double* y, double* z, double* delta,
                                  double* normal);
void c_PUST_OLD_RectCub_getNormalEllipsoid(c_PUST_OLD_RectCub* a_self,
                                           double* x, double* y, double* z,
                                           double* column1, double* column2,
                                           double* column3, double* center,
                                           double* normal);

// Updated Solve Face
void c_PUST_OLD_RectCub_solveFaceEllipsoid(
    c_PUST_OLD_RectCub* a_self, double* STCoeff, double* P0, double* P1,
    double* P2, double* P3, double* column1, double* column2, double* column3,
    double* center, double* Pressure, double* Marangoni, double* a_force);
// Debug
void c_PUST_OLD_RectCub_printSolver(c_PUST_OLD_RectCub* a_self);
}

#endif