#include "Spline.h"
#include <cmath>
#include <vector>

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
        if(KnotVector[i] <= u && u <= KnotVector[i+1]) {
            P = 1;
        } else {
            P = 0;
        }
    } else {
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
        P  = Q1*BBasisFunction(i,p-1,u)+Q2*BBasisFunction(i+1,p-1,u);
    } 
    return P;
}

// Assumes Quadratic NURBs
std::vector<std::vector<double>> Spline::BasisCoefficients(int i) { 
    std::vector<std::vector<double>> coeffs = {{0,0,0},{0,0,0},{0,0,0}};
    // Make sure i is within possible bounds
    if(i+4 > KnotVector.size()){
        // Error Case, once I learn to throw errors
        coeffs = {{-1}};
        return coeffs;
    }
    // Denominators for all terms
    const double denom1 = (KnotVector[i+2]-KnotVector[i])*(KnotVector[i+1]-KnotVector[i]);
    const double denom2 = (KnotVector[i+2]-KnotVector[i])*(KnotVector[i+2]-KnotVector[i+1]);
    const double denom3 = (KnotVector[i+3]-KnotVector[i+1])*(KnotVector[i+2]-KnotVector[i+1]);
    const double denom4 = (KnotVector[i+3]-KnotVector[i+1])*(KnotVector[i+3]-KnotVector[i+2]);
    // Quadratic Term
    const double numer1Q = BBasisFunction(i,0,KnotVector[i]);
    const double numer2Q = -BBasisFunction(i+1,0,KnotVector[i+1]);
    const double numer3Q = -BBasisFunction(i+1,0,KnotVector[i+1]);
    const double numer4Q = BBasisFunction(i+2,0,KnotVector[i+2]);
    // Linear Terms
    const double numer1L = -2*KnotVector[i]*BBasisFunction(i,0,KnotVector[i]);
    const double numer2L = (KnotVector[i+2]+KnotVector[i])*BBasisFunction(i+1,0,KnotVector[i+1]);
    const double numer3L = (KnotVector[i+3]+KnotVector[i+1])*BBasisFunction(i+1,0,KnotVector[i+1]);
    const double numer4L = -2*KnotVector[i+3]*BBasisFunction(i+2,0,KnotVector[i+2]);
    // Constant Term
    const double numer1C = pow(KnotVector[i],2)*BBasisFunction(i,0,KnotVector[i]);
    const double numer2C = -KnotVector[i]*KnotVector[i+2]*BBasisFunction(i+1,0,KnotVector[i+1]);
    const double numer3C = -KnotVector[i+1]*KnotVector[i+3]*BBasisFunction(i+1,0,KnotVector[i+1]);
    const double numer4C = pow(KnotVector[i+3],2)*BBasisFunction(i+2,0,KnotVector[i+2]);
    // Set up Coefficients and BOunds
    // First Span
    if(denom1 == 0){
        coeffs[0][0] = 0;
        coeffs[0][1] = 0;
        coeffs[0][2] = 0;
    } else {
        coeffs[0][0] = numer1Q/denom1;
        coeffs[0][1] = numer1L/denom1;
        coeffs[0][2] = numer1C/denom1;
    }
    // Second Span
    if(denom2 == 0) {
        coeffs[1][0] = 0;
        coeffs[1][1] = 0;
        coeffs[1][2] = 0;
    } else {
        coeffs[1][0] = numer2Q/denom2 + numer3Q/denom3;
        coeffs[1][1] = numer2L/denom2 + numer3L/denom3;
        coeffs[1][2] = numer2C/denom2 + numer3C/denom3;
    }
    // Third Span
    if(denom4 == 0) {
        coeffs[2][0] = 0;
        coeffs[2][1] = 0;
        coeffs[2][2] = 0;
    } else {
        coeffs[2][0] = numer4Q/denom4;
        coeffs[2][1] = numer4L/denom4;
        coeffs[2][2] = numer4C/denom4;
    }
    return coeffs;
}

std::vector<std::vector<double>> Spline::BasisCoefficientBounds(int i) {
    std::vector<std::vector<double>> bounds = {{0,0},{0,0},{0,0}};
    // First Span
    bounds[0][0] = KnotVector[i];
    bounds[0][1] = KnotVector[i+1];
    // Second Span
    bounds[1][0] = KnotVector[i+1];
    bounds[1][1] = KnotVector[i+2];
    // Third Span 
    bounds[2][0] = KnotVector[i+2];
    bounds[2][1] = KnotVector[i+3];

    return bounds;
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