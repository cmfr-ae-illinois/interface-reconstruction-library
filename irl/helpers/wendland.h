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
} // End Namespace




#include "irl/helpers/wendland.cpp"

#endif