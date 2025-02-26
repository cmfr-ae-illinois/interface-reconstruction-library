#include <iostream>
#include "Spline.h"

int main() {
    std :: cout << "Hello World!\n";

    std::vector<std::vector<double>> V1 = {{1,2},{3,4},{4,5},{6,7},{8,9}};
    std::vector<std::vector<double>> V2 = {{2,3},{4,5}};
    std::vector<double> KV = {0,0,0,0.25,0.5,0.75,1,1,1};
    std::vector<double> W = {5,6,7,8};
    
    Spline SP1 = Spline(V1,KV,W);
    Spline SP2 = Spline(V2,W,KV);

    double res = SP1.BBasisFunction(2,0,0.125);
    std :: cout << res;
    res = SP1.BBasisFunction(2,0,0.4);
    std :: cout << res;
    
    return 0;
}