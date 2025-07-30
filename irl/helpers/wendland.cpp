#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_

#include <vector>
#include <limits>
#include <tuple>

#include "examples/PUSurfaceTension/pu_neighborhood.h"
#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/moments/general_moments.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

namespace IRL {
    // ============== Wendland Class Functions
    double Wendland::phi(Pt xi, double delta, Pt x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]));
        if (r >= delta) return 0.0;
        double s = 1.0 - r / delta;
        return std::pow(s, 4) * (4.0 * r / delta + 1.0);
    }

    double Wendland::dphidx(Pt xi, double delta, Pt x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) );
        if (r >= delta) return 0.0;
        return 1.0 / std::pow(delta, 5.0) * 20.0 * (x_eval[0] - xi[0]) * std::pow(r - delta, 3.0);
    }

    double Wendland::dphidy(Pt xi,double delta,Pt x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) );
        if (r >= delta) return 0.0;
        return 1.0 / std::pow(delta, 5.0) * 20.0 * (x_eval[1] - xi[1]) * std::pow(r - delta, 3.0);
  
    }

    double Wendland::ddphidxx(Pt xi,double delta,Pt x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) );
        if (r >= delta ) return 0.0;
        if (r < 1e-12) r = 1.0; 
        double num = 20.0 * std::pow(r - delta, 2.0) * ( (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) * (r - delta) 
                                                        + x_eval[0] * x_eval[0] * (4.0 * r - delta) 
                                                        + xi[0] * xi[0] * (4.0 * r - delta)
                                                        + 2.0 * x_eval[0] * xi[0] * (-4.0 * r + delta) );
        double denom = std::pow(r, 2.0) * std::pow(delta, 5.0);
        return num/denom;
    }

    double Wendland::ddphidyy(Pt xi,double delta,Pt x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) );
        if (r >= delta ) return 0.0;
        if (r < 1e-12) r = 1.0;
        double num = 20.0 * std::pow(r - delta, 2.0) * ( x_eval[0] * x_eval[0] * (r - delta) 
                                                    + xi[0] * xi[0] * (r - delta) 
                                                    + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) * (4.0 * r - delta)
                                                    + 2.0 * x_eval[0] * xi[0] * (-r + delta) );
        double denom = std::pow(r, 2.0) * std::pow(delta, 5.0);
        return num/denom;
    }

    double Wendland::ddphidxy(Pt xi,double delta,Pt x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) );
        if (r >= delta ) return 0.0;
        if (r < 1e-12) r = 1.0;
        double num = 60.0 * (x_eval[0] - xi[0]) * (x_eval[1] - xi[1]) * std::pow(r - delta, 2.0);
        double denom = r * std::pow(delta, 5.0);
        return num/denom;
    }
}

#endif

