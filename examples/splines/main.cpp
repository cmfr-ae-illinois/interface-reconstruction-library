#include <iostream>
#include "Spline.h"

int main() {
    std :: cout << "Hello World!\n";

    std::vector<std::vector<double>> V1 = {{1,2},{3,4},{4,5},{6,7},{8,9}};
    std::vector<std::vector<double>> V2 = {{2,3},{4,5}};
    std::vector<double> KV = {0,0,0,1,2,3,4,4,5,5,5};
    int p = 2;
    std::vector<double> W = {5,6,7,8};
    
    Spline SP1 = Spline(V1,KV,W);
    Spline SP2 = Spline(V2,W,KV);
    
    
    for(int k = 0; k < 8;k++){
        std :: cout << "Spline number ";
        std :: cout << k;
        std :: cout << " Result\n";
        std::vector<std::vector<double>> coeffs = SP1.BasisCoefficients(k);
        std::vector<std::vector<double>> bounds = SP1.BasisCoefficientBounds(k);
        for(int i = 0;i < coeffs.size(); i++) {
            for(int j = 0; j <coeffs[0].size();j++) {
                std::cout << coeffs[i][j];
                std::cout <<",";
            }
            std::cout << "      ";
            for(int j = 0; j <bounds[0].size();j++) {
                std::cout << bounds[i][j];
                std::cout <<",";
            }
            std::cout << "\n";
        }
    }






    return 0;
}