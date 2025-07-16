#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_H_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_H_


#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <tuple>

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
            static double phi(Pt xi, double delta, Pt x_eval);

            // Wendland x derivative eval
            static double dphidx(Pt xi, double delta, Pt x_eval);

            // Wendland y derivative eval
            static double dphidy(Pt xi, double delta, Pt x_eval);

            // Wendland xx 2nd derivative eval
            static double ddphidxx(Pt xi, double delta, Pt x_eval);

            // Wendland yy 2nd derivative eval
            static double ddphidyy(Pt xi, double delta, Pt x_eval);

            // Wendland xy 2nd derivative eval
            static double ddphidxy(Pt xi, double delta, Pt x_eval);

            // Disallow Instance Creation
            Wendland() = delete;
    };

    // template<class SeparatorType>
    class ImplicitSurface {
        private:
            const std::vector<Pt> centroids;
            const std::vector<SeparatorVariant> separators;
            const double kernel_size;
        public:
            // Constructor
            ImplicitSurface(const std::vector<Pt>& centroids_, const std::vector<SeparatorVariant>& separators, const double& kernel_size_);
            // Function Eval
            double F(Pt& x); // Change to Points
            // x derivative eval
            double Fx(Pt& x);
            // y derivative eval
            double Fy(Pt& x);
            // Get Value,Gradient,Hessian
            std::tuple<double,Eigen::Vector3d,Eigen::Matrix3d>
                getValueAndGradAndHessian(Pt& x);
            // Hession Eval
            std::vector<double> HessianTerms(Pt& x);
            // Given an initial guess x0, finding the nearest point for when F(x)=0
            Normal projectToImplicitSurface(const Pt& x0, bool& usePlane); // Can add in options for max_iter and tol later
            // Find intersection between the implicit curve and a provided line. 
            std::vector<Pt> intersectEdge(const Pt& x0, const Pt& x1, const int& Npartitions);
    };

    template<class CellType>
    class PUST {
        private:
            const PUSTNeighborhood<CellType> stencil_m;

        public:
            // Constructor
            PUST(const PUSTNeighborhood<CellType> stencil_);
            // Takes Neighborhood and Returns the Implicit Surface
            ImplicitSurface neighborhoodToImplicitSurface(double delta);
            // Solve Method - Returns the surface tension vector in center cell
            Normal solve(double STCoeff);
            
    };


} // End Namespace IRL

#include "irl/conservative_surface_tension/pu_solve.tpp"

#endif