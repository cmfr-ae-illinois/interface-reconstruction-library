#include "Spline.h"

// Will Delete, for testing
float Spline::add(float x, float y) {
    float add;
    add = x-y;
    return add;
}

// Constructors
Spline::Spline(std::vector<std::vector<double>> CP,std::vector<double> KV,std::vector<double> W) {
    Spline::ControlPoints = CP;
    Spline::KnotVector = KV;
    Spline::Weights = W;
}


// Static Methods ***********************


// Dynamic Methods **************************
double Spline::BBasisFunction(int i, int p,double u) {
    double numer1 = u - KnotVector[i];
    double numer2 = KnotVector[i+p+1]-u;
    double denom1 = KnotVector[i+p] - KnotVector[i];
    double denom2 = KnotVector[i+p+1]-KnotVector[i+1];
    double P;
    double Q1;
    double Q2;

    if(p == 0){
        std :: cout << "Base Case";
        std :: cout << "\n";
        if(KnotVector[i] <= u && u <= KnotVector[i+1]) {
            P = 1;
        } else {
            P = 0;
        }
    } else {
        std :: cout << "P curr = ";
        std :: cout << p;
        std :: cout << "\n";
        if(denom1 == 0) {
            Q1 = 0;
        } else {
            Q1 = numer1/denom1;
        }

        if(denom2 == 0) {
            Q2 = 0;
        } else {
            Q2 = numer2/denom2;
        }

        double P  = Q1*BBasisFunction(i,p-1,u)+Q2*BBasisFunction(i+1,p-1,u);
    }
    return P;
}









// Getters
std::vector<std::vector<double>> Spline::getControlPoints() {
    return ControlPoints;
}
std::vector<double> Spline::getKnotVector() {
    return KnotVector;
}
std::vector<double> Spline::getWeights() {
    return Weights;
}
// Setters
void Spline::setControlPoints(std::vector<std::vector<double>> input) {
    ControlPoints = input;
}
void Spline::setKnotVector(std::vector<double> input) {
    KnotVector = input;
}
void Spline::setWeights(std::vector<double> input){
    Weights = input;
}