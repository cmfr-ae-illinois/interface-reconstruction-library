
#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_SOLVE_RECTANGULAR_CUBOID_SEPARATOR_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_SOLVE_RECTANGULAR_CUBOID_SEPARATOR_H_

#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/interface_reconstruction_methods/pu_solve.h"
#include "irl/variant_reconstruction/separator_variant.h"

extern "C" {

    struct c_PUST_RectCub {
        IRL::PUST<IRL::RectangularCuboid>* obj_ptr = nullptr;
    }

    void c_PUST_RectCub_new(c_PUST_RectCub* a_self);

    void c_PUST_RectCub_delete(c_PUST_RectCub* a_self);

    void c_PUST_RectCub_solveEdge(c_PUST_RectCub* a_self, double STCoeff, double* P0, double* P1, double* a_force);
}

#endif