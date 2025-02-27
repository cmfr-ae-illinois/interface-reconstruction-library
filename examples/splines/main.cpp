#include <iostream>
#include "Spline.h"

int main() {
    std :: cout << "Hello World!\n";

    std::vector<std::vector<double>> V1 = {{1,1},{1,2},{1,3},{1,4},{1,5},{1,6},{1,7},{1,8},{1,9}};
    std::vector<std::vector<double>> V2 = {{2,3},{4,5}};
    std::vector<double> KV = {0,0,0,1,2,3,4,4,5,5,5};
    int p = 2;
    std::vector<double> W = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    
    Spline SP1 = Spline(V1,KV,W);
    Spline SP2 = Spline(V2,W,KV);
    std::vector<std::vector<double>> Dcoeffs = SP1.CurveCoefficients();
    Dcoeffs = SP1.getYCoeffs();
    std::vector<std::vector<double>> spans = SP1.getSpans();
    std::cout << spans.size() << "\n";
    for(int i = 0; i<Dcoeffs.size();i++) {
        for(int j = 0; j<Dcoeffs[0].size();j++) {
            std :: cout << Dcoeffs[i][j] <<",";
        }
        std::cout << "\n";
    }
    
    std::cout << SP1.getArcLength();
    return 0;
}