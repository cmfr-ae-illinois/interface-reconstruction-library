#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_

#include <vector>
#include <limits>
#include <tuple>

#include "examples/PUSurfaceTension/pu_neighborhood.h"
#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/moments/general_moments.h"
#include "irl/variant_reconstruction/separator_variant.h"
#include "irl/helpers/wendland.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>


namespace IRL {
    

    // =================== Implicit Surface Class Functions
    // template <class SeparatorType>
    PUImplicitSurface::PUImplicitSurface(const std::vector<Pt>& centroids_, 
                                    const std::vector<SeparatorVariant>& separators_,
                                    const double& kernel_size_)
                                    : centroids(centroids_), separators(separators_), kernel_size(kernel_size_) {}
    
    
    // template <class SeparatorType>
    double PUImplicitSurface::F(Pt& x) {
        double num = 0.0;
        double denom = 0.0;
        for(int i=0; i < centroids.size(); i++) {
            double wi = Wendland::phi(centroids[i],kernel_size,x);
            auto variantRet = separators[i].getSignedDistanceAndGradAndHessianSep(x,centroids[i]);
            double Fi = std::get<0>(variantRet);
            
            num += wi * Fi;
            denom += wi;
        }
        denom = safelyEpsilon(denom);
        return num/ denom;
    }
    
    std::tuple<double,Eigen::Vector3d,Eigen::Matrix3d> PUImplicitSurface::getValueAndGradAndHessian(Pt& x) {
        // =====Implicit Function Value Variables
        double F_phi_sum = 0.0; // Sum of Fi*phi_i
        double phi_sum   = 0.0; // sum of phi_i

        // ===== Gradient Value Variables
        // x Gradient
        double dphidx_F_sum = 0.0; // Sum of dphi_dx * F_i
        double phi_dFdx_sum = 0.0; // Sum of phi*dFi_dx
        double dphidx_sum   = 0.0; // Sum of dphi_dx
        // y Gradient
        double dphidy_F_sum = 0.0;
        double phi_dFdy_sum = 0.0;
        double dphidy_sum   = 0.0;

        // ===== Hessian Value Variables
        // x-y mixed derivative
        double dphi_dxdy_F_sum = 0.0;
        double dphidx_dFdy_sum = 0.0;
        double dphidy_dFdx_sum = 0.0;
        double phi_dF_dxdy_sum = 0.0;
        double dphi_dxdy_sum = 0.0;
    
        // x-x second derivative
        double dphi_dxdx_F_sum = 0.0;
        double dphidx_dFdx_sum = 0.0; // NOTICE: THIS WILL BE DOUBLED SINCE IT REPEATS
        double phi_dF_dxdx_sum = 0.0;
        double dphi_dxdx_sum = 0.0;
        // y-y second derivative
        double dphi_dydy_F_sum = 0.0;
        double dphidy_dFdy_sum = 0.0; // NOTICE: THIS WILL BE DOUBLED SINCE IT REPEATS IN EQUATION0
        double phi_dF_dydy_sum = 0.0;
        double dphi_dydy_sum = 0.0;
        
        // Calculate
        for(int i = 0; i < centroids.size(); ++i) {
            // Distance Weights
            double phi = Wendland::phi(centroids[i],kernel_size,x);
            double dphidx = Wendland::dphidx(centroids[i],kernel_size,x);
            double dphidy = Wendland::dphidy(centroids[i],kernel_size,x);
            double ddphidxx = Wendland::ddphidxx(centroids[i],kernel_size,x);
            double ddphidyy = Wendland::ddphidyy(centroids[i],kernel_size,x);
            double ddphidxy = Wendland::ddphidxy(centroids[i],kernel_size,x);

            // Function Values
            auto valGradHess = separators[i].getSignedDistanceAndGradAndHessianSep(x, centroids[i]);
            double Fi = std::get<0>(valGradHess);
            Eigen::Vector3d gradFi = std::get<1>(valGradHess);
            Eigen::Matrix3d hessFi = std::get<2>(valGradHess);

            double Fix = gradFi(0);
            double Fiy = gradFi(1);

            double Fixx = hessFi(0,0);
            double Fixy = hessFi(0,1);
            double Fiyy = hessFi(1,1);

            // Add To Values
            F_phi_sum += Fi*phi;
            phi_sum += phi;

            dphidx_F_sum += dphidx*Fi;
            phi_dFdx_sum += phi*Fix;
            dphidx_sum += dphidx;

            dphidy_F_sum += dphidy*Fi;
            phi_dFdy_sum += phi*Fiy;
            dphidy_sum += dphidy;

            dphi_dxdy_F_sum += ddphidxy*Fi;
            dphidx_dFdy_sum += dphidx*Fiy;
            dphidy_dFdx_sum += dphidy*Fix;
            phi_dF_dxdy_sum += phi*Fixy;
            dphi_dxdy_sum += ddphidxy;

            dphi_dxdx_F_sum += ddphidxx*Fi;
            dphidx_dFdx_sum += dphidx*Fix;
            phi_dF_dxdx_sum += phi*Fixx;
            dphi_dxdx_sum += ddphidxx;

            dphi_dydy_F_sum += ddphidyy*Fi;
            dphidy_dFdy_sum += dphidy*Fiy;
            phi_dF_dydy_sum += phi*Fiyy;
            dphi_dydy_sum += ddphidyy;
        }
        // Value
        const double F = F_phi_sum / phi_sum;
        // Gradient 
        const double Fx = ((dphidx_F_sum+phi_dFdx_sum)*phi_sum - F_phi_sum * dphidx_sum) / (phi_sum * phi_sum);
        const double Fy = ((dphidy_F_sum+phi_dFdy_sum)*phi_sum - F_phi_sum * dphidy_sum) / (phi_sum * phi_sum);
        Eigen::Vector3d gradF = Eigen::Vector3d(Fx,Fy,0);
        // Hessian
        // xx derivative
        const double Term1InsideLine1xx = (dphi_dxdx_F_sum + 2.0*dphidx_dFdx_sum + phi_dF_dxdx_sum) * phi_sum;
        const double Term1InsideLine2xx = (dphidx_F_sum+phi_dFdx_sum)*dphidx_sum - (dphidx_F_sum+phi_dFdx_sum)*dphidx_sum; // Goes to zero for xx,yy
        const double Term1InsideLine3xx = -F_phi_sum*dphi_dxdx_sum;
        const double term1xx = (Term1InsideLine1xx + Term1InsideLine2xx + Term1InsideLine3xx)/(phi_sum*phi_sum);

        const double term2xx = Fx * (2.0*dphidx_sum) / phi_sum;
        const double Fxx = term1xx - term2xx;

        // yy Deriative
        const double Term1InsideLine1yy = (dphi_dydy_F_sum + 2.0*dphidy_dFdy_sum + phi_dF_dydy_sum) * phi_sum;
        const double Term1InsideLine2yy = (dphidy_F_sum+phi_dFdy_sum)*dphidy_sum - (dphidy_F_sum+phi_dFdy_sum)*dphidy_sum; // Goes to zero for xx,yy
        const double Term1InsideLine3yy = -F_phi_sum*dphi_dydy_sum;
        const double term1yy = (Term1InsideLine1yy + Term1InsideLine2yy + Term1InsideLine3yy)/(phi_sum*phi_sum);

        const double term2yy = Fy * (2.0*dphidy_sum) / phi_sum;
        const double Fyy = term1yy - term2yy;

        // Mixed Deriative
        const double Term1InsideLine1xy = (dphi_dxdy_F_sum + dphidx_dFdy_sum + dphidy_dFdx_sum+phi_dF_dxdy_sum) * phi_sum;
        const double Term1InsideLine2xy = (dphidx_F_sum+phi_dFdx_sum)*dphidy_sum - (dphidy_F_sum+phi_dFdy_sum) * dphidx_sum;
        const double Term1InsideLine3xy = - F_phi_sum*dphi_dxdy_sum;
        const double term1xy = (Term1InsideLine1xy + Term1InsideLine2xy + Term1InsideLine3xy)/(phi_sum*phi_sum);

        const double term2xy = Fx*(2.0*dphidy_sum)/phi_sum;
        const double term2xy_other = Fy*(2.0*dphidx_sum)/phi_sum;

        const double Fxy = term1xy - term2xy;
        Eigen::Matrix3d hessF = Eigen::Matrix3d::Zero();
        hessF(0,0) = Fxx;
        hessF(0,1) = Fxy;
        hessF(1,0) = Fxy;
        hessF(1,1) = Fyy;
        return std::make_tuple(F,gradF,hessF);
    }
    




    // template <class SeparatorType>
    std::vector<double> PUImplicitSurface::HessianTerms(Pt& x) {
        double N   = 0.0;
        double Nx  = 0.0;
        double Ny  = 0.0;
        double Nxx = 0.0;
        double Nyy = 0.0;
        double Nxy = 0.0;
        double D   = 0.0;
        double Dx  = 0.0;
        double Dy  = 0.0;
        double Dxx = 0.0;
        double Dyy = 0.0;
        double Dxy = 0.0;

        for (int i = 0; i < centroids.size(); i++){
            // phi derivatives
            double phi_i      = Wendland::phi(centroids[i],kernel_size,x);
            double dphidx_i   = Wendland::dphidx(centroids[i],kernel_size,x);
            double dphidy_i   = Wendland::dphidy(centroids[i],kernel_size,x);
            double ddphidxx_i   = Wendland::ddphidxx(centroids[i],kernel_size,x);
            double ddphidyy_i   = Wendland::ddphidyy(centroids[i],kernel_size,x);
            double ddphidxy_i   = Wendland::ddphidxy(centroids[i],kernel_size,x);
            double F_i = 0.0;
            double dFdx_i = 0.0;
            double dFdy_i = 0.0;
            double ddFdxx_i = 0.0;
            double ddFdyy_i = 0.0;
            double ddFdxy_i = 0.0;
            // F derivatives
            if(const IRL::PlanarSeparator* sepPtr =
                    std::get_if<IRL::PlanarSeparator>(&(separators[i]))) { // If Plane
                // std::cout << "Plane Detected\n";
                Normal n = (*sepPtr)[0].normal();
                F_i = n[0] * (x[0] - centroids[i][0]) +
                        n[1] * (x[1] - centroids[i][1]); 
                dFdx_i = n[0];
                dFdy_i = n[1];
                // ddFdxx_i = 0.0;
                // ddFdyy_i = 0.0;
                // ddFdxy_i = 0.0;
            }
            // Normal n = separators[i][0].normal();
            // double F_i = n[0] * (x[0] - centroids[i][0]) +
            //             n[1] * (x[1] - centroids[i][1]); 
            // double dFdx_i = n[0];
            // double dFdy_i = n[1];
            // double ddFdxx_i = 0.0;
            // double ddFdyy_i = 0.0;
            // double ddFdxy_i = 0.0;

            // numerator terms
            N   += phi_i * F_i;
            Nx  += dphidx_i * F_i + phi_i * dFdx_i;
            Ny  += dphidy_i * F_i + phi_i * dFdy_i;
            Nxx += F_i * ddphidxx_i + phi_i * ddFdxx_i + 2.0 * dphidx_i * dFdx_i;
            Nyy += F_i * ddphidyy_i + phi_i * ddFdyy_i + 2.0 * dphidy_i * dFdy_i;
            Nxy += F_i * ddphidxy_i + phi_i * ddFdxy_i + dphidx_i * dFdy_i + dphidy_i * dFdx_i;

            // denominator terms
            D   += phi_i;
            Dx  += dphidx_i;
            Dy  += dphidy_i;
            Dxx += ddphidxx_i;
            Dyy += ddphidyy_i;
            Dxy += ddphidxy_i;
        }
        double Fxx = (D * (Nxx * D - N * Dxx) - 2.0 * Dx * (Nx * D - N * Dx)) / (D * D * D);
        double Fyy = (D * (Nyy * D - N * Dyy) - 2.0 * Dy * (Ny * D - N * Dy)) / (D * D * D);
        double Fxy = (D * (Nxy * D + Nx * Dy - Ny * Dx - N * Dxy) - 2.0 * Dy * (Nx * D - N * Dx)) / (D * D * D);
        
        return {Fxx, Fyy, Fxy};
    }
    // template <class SeparatorType>
    Normal PUImplicitSurface::projectToImplicitSurface(const Pt& x0, bool& usePlane) {
        int max_iter = 200;
        double tol = 1e-12;
        Pt x = x0;

        // std::ostringstream diagnostics;
        // diagnostics << "x0 = " << x0 << std::endl;
        bool nan_detected = false;

        for(int i =0; i < max_iter; i++) {
            double F = this->F(x);
            double Fx = this->Fx(x);
            double Fy = this->Fy(x);
            Pt gradF = Pt(Fx,Fy,0);
            double denom = 1/safelyEpsilon(Fx*Fx+Fy*Fy);
            Pt change = denom*F*gradF;
            x = x - change;

            if((std::isnan(x[0])) || (std::isnan(x[1]) || (std::isnan(x[2])))  && !nan_detected) {
                nan_detected = true;
                usePlane = true;
                std::cout << "========== NaN detected at iteration " << i+1 << " ==========\n";
                // std::cout << diagnostics.str(); 
            }
            double magChange = std::sqrt(change[0]*change[0] + change[1]*change[1] + change[2]*change[2]);
            if(std::fabs(magChange) < tol && !nan_detected) break;
        }
        return x;
    }

    // template <class SeparatorType>
    std::vector<Pt> PUImplicitSurface::intersectEdge(const Pt& x0, const Pt& x1, const int& Npartitions) {
        // Split the domain into segments
        std::vector<Pt> sampleLocations = {};
        // At these locations, calculate the function value
        std::vector<double> values = {};
        // Also get the sign of these values
        std::vector<double> signs = {};
        for(int i = 0; i < Npartitions+1;i++) {
            Pt temp = (1-static_cast<double>(i)/static_cast<double>(Npartitions)) * x0 + (static_cast<double>(i)/static_cast<double>(Npartitions)) * x1;
            sampleLocations.push_back(temp);
            
            double val = this->F(temp);
            values.push_back(val);

            double sgn = (0.0 < val) - (val < 0.0);
            signs.push_back(sgn);
        }

        // Loop over all the partitions. If the signs are different, do a bisection method to find the root
        std::vector<Pt> intersections = {};
        // std::cout << intersections.size() << std::endl;
        Pt upperX;
        Pt lowerX;
        Pt midX;
        double upperVal;
        double lowerVal;
        double midVal;

        double tol = 1e-12;
        double max_iters = 200;

        for(int i = 0; i < Npartitions; i++) {
            if(signs[i] == 0 || signs[i+1] ==0 ) { // At least one root
                if(signs[i] == 0) {
                    std::cout << "Zero Intersection Found" << std::endl;
                    intersections.push_back(sampleLocations[i]);
                }
                if(signs[i+1] == 0 && i + 1 == Npartitions) { 
                    std::cout << "Zero Intersection Found" << std::endl;
                    intersections.push_back(sampleLocations[i+1]);
                }
            } else if(signs[i] != signs[i+1]) { // Different signs, root somewhere
                // Decide which Side is upper and lower
                if(signs[i] == 1) { // Left  is upper value
                    upperX = sampleLocations[i];
                    upperVal = values[i];

                    lowerX = sampleLocations[i+1];
                    lowerVal = values[i+1];
                } else {
                    upperX = sampleLocations[i+1];
                    upperVal = values[i+1];

                    lowerX = sampleLocations[i];
                    lowerVal = values[i];
                }
                // Apply Bisection Method
                for(int j =0; j < max_iters; j++) { // Do until you reach max iters
                    midX = 0.5*(upperX + lowerX);
                    std::cout << "Current Point = "  << midX[0] << ","<< midX[1] << ","<< midX[2] << std::endl;
                    midVal = this-> F(midX);
                    std::cout << "Current Val = "  << midVal<<std::endl;
                    if(midVal < tol) {
                        break;
                    } else if(midVal > 0.0) {
                        upperX = midX;
                    } else {
                        lowerX = midX;
                    }
                } // End Bisection
                // Add Intersection
                intersections.push_back(midX);
                
            }
        } // End loop over partitions
        return intersections;
    }

    Normal PUImplicitSurface::getTangent(Pt& x) {
        auto holdsGrad = this->getValueAndGradAndHessian(x);
        auto gradF = std::get<1>(holdsGrad);
        double Fx = gradF(0);
        double Fy = gradF(1);

        return Normal(-Fy,Fx,0.0);
    }

    double PUImplicitSurface::getCurvature(Pt& x) {
        auto holdsGradAndHessian = this->getValueAndGradAndHessian(x);
        auto gradF = std::get<1>(holdsGradAndHessian);
        auto hessF = std::get<2>(holdsGradAndHessian);

        double Fxx = hessF(0,0);
        double Fyy = hessF(1,1);
        double Fxy = hessF(0,1);
        double Fx = gradF(0);
        double Fy = gradF(1);

        double numer = Fxx*Fy*Fy - 2*Fxy*Fx*Fy + Fx*Fx*Fyy;
        double magGradF = std::sqrt(Fx*Fx+Fy*Fy);
        double denom = magGradF*magGradF*magGradF;

        double kz = -numer/denom;

        return kz;
    }   

    // ============== Solver Methods
    template <class CellType>
    PUST<CellType>::PUST(const PUSTNeighborhood<CellType> stencil_) : stencil_m(stencil_){
        this->surface_m = this ->neighborhoodToImplicitSurface(5.0);
    }

    // void getIntersectionPts(const IRL::Polygon& a_polygon,
    //                     const IRL::Plane& a_cutting_plane,
    //                     IRL::StackVector<IRL::Pt, 2>* a_intersection_pts) {
    //     a_intersection_pts->resize(0);
    //     double distance = a_cutting_plane.signedDistanceToPoint(a_polygon[0]);
    //     for (int n = 0; n < a_polygon.getNumberOfVertices() - 1; ++n) {
    //         double next_distance =
    //             a_cutting_plane.signedDistanceToPoint(a_polygon[n + 1]);
    //         if (distance * next_distance < 0.0) {
    //         a_intersection_pts->push_back(IRL::Pt::fromEdgeIntersection(
    //             a_polygon[n], distance, a_polygon[n + 1], next_distance));
    //         if (a_intersection_pts->size() == 2) {
    //             break;
    //         }
    //         }
    //         distance = next_distance;
    //     }s
    // }
    
    template <class CellType>
    PUImplicitSurface PUST<CellType>::neighborhoodToImplicitSurface(double delta) {
        const auto centroids = stencil_m.getCentroids();
        const auto separators = stencil_m.getSeparators();

        return PUImplicitSurface(centroids,separators,delta);
    }

    template <class CellType>
    std::vector<double> PUST<CellType>::solve(double STCoeff,int direction) {
        if(constexpr(std::is_same_v<CellType,RectangularCuboid>)) {
            // Direction should be 0 for x, 1 for y
            
            // First, Make the Implicit Edge
            PUImplicitSurface s = this->neighborhoodToImplicitSurface(5.0);

            // Below Here is the Intersection
            Pt P0,P1;
            std::vector<Pt> inters;

            CellType c = stencil_m.getCenterCell(); // Should be an x cell or a y cell
            Pt BL = c.getLowerLimits();
            Pt TR = c.getUpperLimits();
            std::cout << "Lower Point: " << BL[0] << "," << BL[1] << "," << BL[2] << "\n";
            std::cout << "Upper Point: " << TR[0] << "," << TR[1] << "," << TR[2] << "\n";
            Pt BR = Pt(TR[0],BL[1],TR[2]);
            Pt TL = Pt(BL[0],TR[1],BL[2]);
            std::vector<Pt> poly = {BR,TR,TL,BL,BR};

            Normal tangent;
            std::vector<double> surfaceTensionComponents = {0,0,0,0}; 
            // sigma_xx_i, sigma_xy_j+1/2, sigma_xx_i-1, sigma_xy_j-1/2 for x cell
            // sigma_yx,i+1/2, sigma_yy_j, sigma_yx_i-1/2, sigma_yy_j-1 for y cell
            // Assume proper cell is set as center cell.

            double D;
            Pt dP;
            for(int i = 0; i < poly.size()-1; i++) {
                P0 = poly[i];
                P1 = poly[i+1];
                inters = s.intersectEdge(P0,P1,10);
                dP = P1-P0;
                D = std::sqrt(dP[0]*dP[0] + dP[1]*dP[1]);

                if(inters.size() > 0) {
                    for(int j = 0; j < inters.size(); j++) {
                        tangent = s.getTangent(inters[j]);
                        
                        surfaceTensionComponents[i] += STCoeff*tangent[direction]/D;
                    }
                }
            }

            return surfaceTensionComponents;
        } else {
            throw std::invalid_argument("Solve Failure: CellType not RectangularCuboid");
        }
    }
                        
} // End Namespace IRL


#endif

