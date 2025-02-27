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
    
    
    std::vector<double> breakpoints = SP1.getBreakpoints();
    std::vector<std::vector<double>> spans = SP1.getSpans();
    for(int i = 0;i<breakpoints.size();i++) {
        std::cout << breakpoints[i];
        std::cout << ",";
    }
    std :: cout << "\n";
    for(int i = 0;i<spans.size();i++) {
        for(int j = 0; j< spans[0].size();j++) {
            std::cout << spans[i][j] << ",";
        }
        std::cout << "\n";
    }






    return 0;
}