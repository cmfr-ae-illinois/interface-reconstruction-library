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
            double wi = std::get<0>(Wendland::evaluateValGradHessian(centroid[i],kernel_size,x));
            auto variantRet = PUImplicitSurface::getSignedDistanceAndGradAndHessianSep(x,centroids[i],&separators[i]);
            double Fi = std::get<0>(variantRet);
            
            num += wi * Fi;
            denom += wi;
        }
        denom = safelyEpsilon(denom);
        return num/ denom;
    }

    std::tuple<double,Eigen::Vector3d,Eigen::Matrix3d>
      PUImplicitSurface::getSignedDistanceAndGradAndHessianSep(const Pt& a_pt, const Pt& a_centroid,const SeparatorVariant* a_sepPtr) const {
        const Pt x = a_pt - a_centroid;
        double F;
        Eigen::Vector3d gradF;
        Eigen::Matrix3d hessF;
        if(const auto sepPtr = std::get_if<PlanarSeparator>(a_sepPtr)) {
            // std::cout << "Variant Plane Detected\n";
            if(sepPtr->getNumberOfPlanes() > 0){
            const Plane plane = (*sepPtr)[0];
            const Normal n = plane.normal();
            F = n * x;
            gradF = Eigen::Vector3d(n[0],n[1],n[2]);
            hessF = Eigen::Matrix3d::Zero();
            }
        } else if (const auto sepPtr = std::get_if<Paraboloid>(a_sepPtr)) {
            // std::cout << "Variant Paraboloid Detected\n";
            const ReferenceFrame frame = sepPtr->getReferenceFrame();
            const double a = sepPtr->getAlignedParaboloid().a();
            const double b = sepPtr->getAlignedParaboloid().b();
            // Move to local frame
            const Pt tmp = a_pt - sepPtr->getDatum();
            Pt xloc;
            for(int d = 0; d < 3; ++d) {
                xloc[d] = frame[d]*tmp;
            }
            const Eigen::Vector3d e0(frame[0][0], frame[0][1], frame[0][2]);
            const Eigen::Vector3d e1(frame[1][0], frame[1][1], frame[1][2]);
            const Eigen::Vector3d e2(frame[2][0], frame[2][1], frame[2][2]);

            // Taubin Distance Norm,grad,hessian
            const double dist_norm =
                    1. / std::sqrt(1. + 4. * a * a * xloc[0] * xloc[0] +
                                    4. * b * b * xloc[1] * xloc[1]);
            const Eigen::Vector3d grad_dist_norm =
                -4. * dist_norm * dist_norm * dist_norm *
                (a * a * e0 * xloc[0] + b * b * e1 * xloc[1]);
            const Eigen::Matrix3d hess_dist_norm =
                4. * dist_norm * dist_norm * dist_norm * dist_norm * dist_norm *
                (a * a *
                    (8. * a * a * xloc[0] * xloc[0] -
                    4. * b * b * xloc[1] * xloc[1] - 1.) *
                    e0 * e0.transpose() +
                b * b *
                    (8. * b * b * xloc[1] * xloc[1] -
                    4. * a * a * xloc[0] * xloc[0] - 1.) *
                    e1 * e1.transpose() +
                12. * a * a * b * b * xloc[0] * xloc[1] *
                    (e1 * e0.transpose() + e0 * e1.transpose()));
            
            // Compute Algebraic Distance, grad, hessian
            const double F_alg =
                    xloc[2] + a * xloc[0] * xloc[0] + b * xloc[1] * xloc[1];
            const Eigen::Vector3d grad_F_alg =
                e2 + 2. * (a * e0 * xloc[0] + b * e1 * xloc[1]);
            const Eigen::Matrix3d hess_F_alg =
                2. * (a * e0 * e0.transpose() + b * e1 * e1.transpose());
            
            // Compute Signed Distance,grad,hessian
            F = F_alg * dist_norm;
            gradF =
                F_alg * grad_dist_norm + grad_F_alg * dist_norm;
            hessF =
                F_alg * hess_dist_norm +
                grad_F_alg * grad_dist_norm.transpose() +
                grad_dist_norm * grad_F_alg.transpose() +
                dist_norm * hess_F_alg;
        } else {
            throw std::runtime_error("No signed distance for Variant Type");
        }
        return std::make_tuple(F,gradF,hessF);
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
            auto retWend = Wendland::evaluateValGradHessian(centroids[i],kernel_size,x);
            double phi = std::get<0>(retWend);
            Eigen::Vector3d gradPhi = std::get<1>(retWend);
            Eigen::Matrix3d hessPhi = std::get<2>(retWend);
            double dphidx = gradPhi(0);
            double dphidy = gradPhi(1);
            double ddphidxx = hessPhi(0,0);
            double ddphidyy = hessPhi(1,1);
            double ddphidxy = hessPhi(1,0);

            // Function Values
            auto valGradHess = PUImplicitSurface::getSignedDistanceAndGradAndHessianSep(x, centroids[i],&separators[i]);
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
                    midVal = this->F(midX);
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

