#include <iostream>
#include "Spline.h"

int main() {
    std :: cout << "Hello World!\n";

    std::vector<std::vector<double>> V1 = {{1,1},{1,2},{1,3},{1,4},{1,5},{1,6},{1,7},{1,8},{1,9}};
    std::vector<std::vector<double>> V2 = {{2,3},{4,5}};
    std::vector<double> KV = {0,0,0,1,2,3,4,4,5,5,5};
    int p = 2;
    std::vector<double> W = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    
    for(int i = 0; i < KV.size(); i++) {
        std::cout << KV[i] << ",";
    }
    std::cout << "\n";
    
    auto spline = Spline(V1, KV, W);

    std::cout << "Arc length = " << spline.getArcLength() << std::endl;

    spline.saveToVTK("spline_test");
    
    return 0;
}