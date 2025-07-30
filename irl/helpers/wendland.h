#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_H_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_H_

namespace IRL {
     // Here we are going to make a static class for the Wendland function
    // This is because instead of having to make different objects for different 
    // detla, x_eval values, we will just use the Wendland class functions.
    // This means we will need to pass in delta,x_eval as parameters, but that should not be a problem maybe?
    class Wendland {
        // Note that for these functions, xi is the center, delta is the radius,
        // and x_eval is the location we are evaluating our function. 
        public:
            // New
            // Compute the Radius
            static double computeR(Pt xi, Pt x_eval);
            // Compute Zeroth Derivative (Function)
            static double eval(double r, double delta);
            // Compute First Derivative
            static double firstDer(double r, double delta);
            // Compute Second Derivative
            static double secondDer(double r, double delta);

            // Compute the Wendland Function
            static std::tuple<double,Eigen::Vector3d,Eigen::Matrix3d>
                evaluateValGradHessian(Pt xi, double delta, Pt x_eval);

            // Disallow Instance Creation
            Wendland() = delete;
    };
} // End Namespace




#include "irl/helpers/wendland.cpp"

#endif