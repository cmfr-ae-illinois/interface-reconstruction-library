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

    // template<class SeparatorType>
    class PUImplicitSurface {
        private:
            const std::vector<Pt> centroids;
            const std::vector<SeparatorVariant> separators;
            const double kernel_size;
        public:
            // Constructor
            PUImplicitSurface(const std::vector<Pt>& centroids_, const std::vector<SeparatorVariant>& separators, const double& kernel_size_);
            // Function Eval
            double F(Pt& x); // Change to Points
            // Get Value,Gradient,Hessian
            std::tuple<double,Eigen::Vector3d,Eigen::Matrix3d>
                getValueAndGradAndHessian(Pt& x);
            // Hession Eval
            std::vector<double> HessianTerms(Pt& x);
            // Given an initial guess x0, finding the nearest point for when F(x)=0
            Normal projectToImplicitSurface(const Pt& x0, bool& usePlane); // Can add in options for max_iter and tol later
            // Find intersection between the implicit curve and a provided line. 
            std::vector<Pt> intersectEdge(const Pt& x0, const Pt& x1, const int& Npartitions);

            // Find the Tangent and Curvature at the point
            Normal getTangent(Pt& x);
            double getCurvature(Pt& x);
    };

    template<class CellType>
    class PUST {
        private:
            const PUSTNeighborhood<CellType> stencil_m;

        public:
            // Constructor
            PUST(const PUSTNeighborhood<CellType> stencil_);
            // Takes Neighborhood and Returns the Implicit Surface
            PUImplicitSurface neighborhoodToImplicitSurface(double delta);
            // Solve Method - Returns the surface tension vector in center cell
            std::vector<double> solve(double STCoeff,int direction);
            
    };


} // End Namespace IRL

#include "irl/interface_reconstruction_methods/pu_solve.tpp"

#endif