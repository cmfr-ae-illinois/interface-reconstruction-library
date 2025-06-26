#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_

#include <vector>

namespace IRL {
    // ============== Wendland Class Functions
    double Wendland::phi(Normal xi, double delta, Normal x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]));
        if (r >= delta) return 0.0;
        double s = 1.0 - r / delta;
        return std::pow(s, 4) * (4.0 * r / delta + 1.0);
    }

    double Wendland::dphidx(Normal xi, double delta, Normal x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) );
        if (r >= delta) return 0.0;
        return 1.0 / std::pow(delta, 5.0) * 20.0 * (x_eval[0] - xi[0]) * std::pow(r - delta, 3.0);
    }

    double Wendland::dphidy(Normal xi,double delta,Normal x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) );
        if (r >= delta) return 0.0;
        return 1.0 / std::pow(delta, 5.0) * 20.0 * (x_eval[1] - xi[1]) * std::pow(r - delta, 3.0);
  
    }

    double Wendland::ddphidxx(Normal xi,double delta,Normal x_eval) {
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

    double Wendland::ddphidyy(Normal xi,double delta,Normal x_eval) {
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

    double Wendland::ddphidxy(Normal xi,double delta,Normal x_eval) {
        double r = std::sqrt( (x_eval[0] - xi[0]) * (x_eval[0] - xi[0]) 
                        + (x_eval[1] - xi[1]) * (x_eval[1] - xi[1]) );
        if (r >= delta ) return 0.0;
        if (r < 1e-12) r = 1.0;
        double num = 60.0 * (x_eval[0] - xi[0]) * (x_eval[1] - xi[1]) * std::pow(r - delta, 2.0);
        double denom = r * std::pow(delta, 5.0);
        return num/denom;
    }

    // =================== Implicit Surface Class Functions
    ImplicitSurface::ImplicitSurface(const std::vector<Normal>& centroids_, 
                                    const std::vector<Normal>& normals_,
                                    const double& kernel_size_)
                                    : centroids(centroids_), normals(normals_), kernel_size(kernel_size_) {}
    
    double ImplicitSurface::F(Normal x) {
        double num = 0.0;
        double denom = 0.0;
        for(int i=0; i < centroids.size(); i++) {
            double wi = Wendland::phi(centroids[i],kernel_size,x);
            double Fi = normals[i][0] * (x[0] - centroids[i][0]) + 
                        normals[i][1] * (x[1] - centroids[i][1]);

            num += wi * Fi;
            denom += wi;
        }

        if(denom < 1e-12) {
            std::cout << "Sum of weights is too small, denominator = " << denom << std::endl;
        }
        std::cout << num / denom << "return value \n";
        return num/ denom;
    }

    double ImplicitSurface::Fx(Normal x) {
        double sum_phi      = 0.0;
        double sum_phi_F    = 0.0;
        double sum_dphidx   = 0.0;
        double sum_dphidx_F = 0.0;
        double sum_phi_dFdx = 0.0;

        for(int i = 0; i < centroids.size(); i++) {
            // phi derivatives
            double phi_i      = Wendland::phi(centroids[i],kernel_size,x);
            double dphidx_i   = Wendland::dphidx(centroids[i],kernel_size,x);
            // Fi derivatives
            Normal n = normals[i];
            double F_i = n[0] * (x[0] - centroids[i][0]) +
                        n[1] * (x[1] - centroids[i][1]);
            double dFdx_i = n[0];
            // terms for Fx
            sum_phi       += phi_i;
            sum_phi_F     += phi_i * F_i;
            sum_dphidx    += dphidx_i;
            sum_dphidx_F  += dphidx_i * F_i;
            sum_phi_dFdx  += phi_i * dFdx_i;
        }

        return ((sum_dphidx_F + sum_phi_dFdx) * sum_phi - sum_phi_F * sum_dphidx) / (sum_phi * sum_phi);
    }

    double ImplicitSurface::Fy(Normal x) {
        double sum_phi      = 0.0;
        double sum_phi_F    = 0.0;
        double sum_dphidy   = 0.0;
        double sum_dphidy_F = 0.0;
        double sum_phi_dFdy = 0.0;

        for (int i = 0; i < centroids.size(); i++){
            // phi derivatives
            double phi_i      = Wendland::phi(centroids[i],kernel_size,x);
            double dphidy_i   = Wendland::dphidy(centroids[i],kernel_size,x);
            // Fi derivatives
            Normal n = normals[i];
            double F_i = n[0] * (x[0] - centroids[i][0]) +
                        n[1] * (x[1] - centroids[i][1]);
            double dFdy_i = n[1];
            // terms for Fy
            sum_phi       += phi_i;
            sum_phi_F     += phi_i * F_i;
            sum_dphidy    += dphidy_i;
            sum_dphidy_F  += dphidy_i * F_i;
            sum_phi_dFdy  += phi_i * dFdy_i;
        }
        return ((sum_dphidy_F + sum_phi_dFdy) * sum_phi - sum_phi_F * sum_dphidy) / (sum_phi * sum_phi);
    }

    std::vector<double> ImplicitSurface::HessianTerms(Normal x) {
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

            // F derivatives
            Normal n = normals[i];
            double F_i = n[0] * (x[0] - centroids[i][0]) +
                        n[1] * (x[1] - centroids[i][1]); 
            double dFdx_i = n[0];
            double dFdy_i = n[1];
            double ddFdxx_i = 0.0;
            double ddFdyy_i = 0.0;
            double ddFdxy_i = 0.0;

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

    Normal ImplicitSurface::projectToImplicitSurface(const Normal& x0, bool& usePlane) {
        int max_iter = 200;
        double tol = 1e-12;
        Normal x = x0;

        // std::ostringstream diagnostics;
        // diagnostics << "x0 = " << x0 << std::endl;
        bool nan_detected = false;

        for(int i =0; i < max_iter; i++) {
            double F = this->F(x);
            double Fx = this->Fx(x);
            double Fy = this->Fy(x);
            Normal gradF = Normal(Fx,Fy,0);

            Normal change = F*gradF/(Fx*Fx + Fy*Fy);
            x = x - change;

            if(std::isnan(x.calculateMagnitude()) && !nan_detected) {
                nan_detected = true;
                usePlane = true;
                std::cout << "========== NaN detected at iteration " << i+1 << " ==========\n";
                // std::cout << diagnostics.str(); 
            }

            if(std::fabs(change.calculateMagnitude()) < tol && !nan_detected) break;
        }
        return x;
    }

    std::vector<Normal> ImplicitSurface::intersectEdge(const Normal& x0, const Normal& x1, const int& Npartitions) {
        // Split the domain into segments
        std::vector<Normal> sampleLocations = {};
        // At these locations, calculate the function value
        std::vector<double> values = {};
        // Also get the sign of these values
        std::vector<double> signs = {};
        for(int i = 0; i < Npartitions+1;i++) {
            Normal temp = (1-double(i)/Npartitions) * x0 + (double(i)/Npartitions) * x1;
            sampleLocations.push_back(temp);
            
            double val = this->F(temp);
            values.push_back(val);

            double sgn = (0.0 < val) - (val < 0.0);
            signs.push_back(sgn);
        }

        // Loop over all the partitions. If the signs are different, do a bisection method to find the root
        std::vector<Normal> intersections = {};
        std::cout << intersections.size() << std::endl;
        Normal upperX;
        Normal lowerX;
        Normal midX;
        double upperVal;
        double lowerVal;
        double midVal;

        double tol = 1e-12;
        double max_iters = 200;

        for(int i = 0; i < Npartitions; i++) {
            if(signs[i] == 0 || signs[i+1] ==0 ) { // At least one root
                if(signs[i] == 0) intersections.push_back(sampleLocations[i]);
                if(signs[i+1] == 0) intersections.push_back(sampleLocations[i+1]);
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

    // ============== Solver Methods

} // End Namespace IRL


#endif

