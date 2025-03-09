#include <iostream>
#include "Spline.h"
#include <cmath>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

int main() {
    std :: cout << "Hello World!\n";
    
    std::vector<std::vector<double>> Q = {{1,0},{0,1},{-1,0},{0,-1},{1,0}};
    std::vector<std::vector<double>> T = {{0,1},{-1,0},{0,-1},{1,0},{0,1}};
    Spline test = Spline::LocalRQuadInterp(Q,T);
    
    std::vector<std::vector<double>> CP = test.getControlPoints();
    std::vector<double> W = test.getWeights();
    std::vector<double> U = test.getKnotVector();

    std::cout<<"\n\n\n";
    // Print Result
    std::cout<<"Knot Vector: \n";
    for(int i = 0; i < U.size();i++) {
        std::cout<<U[i] << ",";
    }
    std::cout <<"\n";
    std::cout<<"Weight Vector: \n";
    for(int i = 0; i < W.size();i++) {
        std::cout<<W[i] << ",";
    }
    std::cout <<"\n";
    std::cout<<"Control Points: \n";
    for(int i = 0; i < CP.size();i++) {
        for(int j = 0; j < CP[0].size();j++) {
            std::cout << CP[i][j] << ",";
        }
        std::cout <<"\n";
    }
    std::vector<double> uset = {0};
    int N = 10;
    for(int i = 1; i <= N; i++) {
        uset.insert(uset.end(),i/N);
    }
    std::vector<std::vector<double>> curvePoints = test.makeRationalQuadCurve(uset);
    std::vector<double> distances (curvePoints.size(),0);
    std::cout<<"Distance from Origin";
    for(int i = 0;i < curvePoints.size(); i++) {
        for(int j = 0;j < curvePoints[0].size();j++) {
            distances[i] += pow(curvePoints[i][j],2);
        }
        distances[i] = sqrt(distances[i]);
        std::cout << distances[i] << ",";
    }


    return 0;
}