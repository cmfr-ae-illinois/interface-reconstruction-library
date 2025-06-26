#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_H_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_H_


#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "examples/PUSurfaceTension/pu_neighborhood.h"

#include "irl/moments/cell_collection.h"
#include "irl/moments/cell_grouped_moments.h"

#include "irl/variant_reconstruction/separator_variant.h"
#include "irl/geometry/general/normal.h"
namespace IRL {
    /// This file contains all the functions
    /// needed to calculate the surface tension
    /// Forces.

    // Here we are going to make a static class for the Wendland function
    // This is because instead of having to make different objects for different 
    // detla, x_eval values, we will just use the Wendland class functions.
    // This means we will need to pass in delta,x_eval as parameters, but that should not be a problem maybe?
    class Wendland {
        // Note that for these functions, xi is the center, delta is the radius,
        // and x_eval is the location we are evaluating our function. 
        public:
            // Wendland Function Eval
            static double phi(Normal xi, double delta, Normal x_eval);

            // Wendland x derivative eval
            static double dphidx(Normal xi, double delta, Normal x_eval);

            // Wendland y derivative eval
            static double dphidy(Normal xi, double delta, Normal x_eval);

            // Wendland xx 2nd derivative eval
            static double ddphidxx(Normal xi, double delta, Normal x_eval);

            // Wendland yy 2nd derivative eval
            static double ddphidyy(Normal xi, double delta, Normal x_eval);

            // Wendland xy 2nd derivative eval
            static double ddphidxy(Normal xi, double delta, Normal x_eval);

            // Disallow Instance Creation
            Wendland() = delete;
    };


    class ImplicitSurface {
        private:
            const std::vector<Normal> centroids, normals;
            const double kernel_size;
        public:
            // Constructor
            ImplicitSurface(const std::vector<Normal>& centroids_, const std::vector<Normal>& normals_, const double& kernel_size_);
            // Function Eval
            double F(Normal x); // Change to Points
            // x derivative eval
            double Fx(Normal x);
            // y derivative eval
            double Fy(Normal x);
            // Hession Eval
            std::vector<double> HessianTerms(Normal x);
            // Given an initial guess x0, finding the nearest point for when F(x)=0
            Normal projectToImplicitSurface(const Normal& x0, bool& usePlane); // Can add in options for max_iter and tol later
            // Find intersection between the implicit curve and a provided line. 
            std::vector<Normal> intersectEdge(const Normal& x0, const Normal& x1, const int& Npartitions);
    };

    template<class CellType>
    class PUST {
        private:
            const PUSTNeighborhood<CellType> stencil;
        public:
            // Constructor
            PUST(const PUSTNeighborhood<CellType> stencil_);
            // Solve Method - Returns the surface tension vector in center cell
            Normal solve(void);

    };
} // End Namespace IRL

#include "examples/PUSurfaceTension/pu_solve.tpp"

#endif