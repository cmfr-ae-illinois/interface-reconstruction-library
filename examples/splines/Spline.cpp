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

    this->makeBreakpoints();
    this->makeSpans();
    this->CurveCoefficients();
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

// Currently for 2D
std::vector<std::vector<double>> Spline::CurveCoefficients() {
    int spansSize = spans.size();

    numerCoeffsX = {{0,0,0}};
    numerCoeffsY = {{0,0,0}};
    denomCoeffs = {{0,0,0}};
    // Loop over basis functions:
    for(int i = 0;i < KnotVector.size()-3;i++){
        std::vector<std::vector<double>> coeffs = this->BasisCoefficients(i);
        std::vector<std::vector<double>> bounds = this->BasisCoefficientBounds(i);

        // Loop Over spans, add coefficients to appropriate span 
        for(int j = 0;j < spans.size();j++) {
            if(i == 0 && j > 0){
                numerCoeffsX.insert(numerCoeffsX.end(),{0,0,0});
                numerCoeffsY.insert(numerCoeffsY.end(),{0,0,0});
                denomCoeffs.insert(denomCoeffs.end(),{0,0,0});
            }
            // First Span
            if(bounds[0][0] == spans[j][0] && bounds[0][1] == spans[j][1]) { 
                for(int k = 0;k < numerCoeffsX[j].size();k++) {
                    denomCoeffs[j][k] += coeffs[0][k]*Weights[i];
                    numerCoeffsX[j][k] += coeffs[0][k]*Weights[i]*ControlPoints[i][0];
                    numerCoeffsY[j][k] += coeffs[0][k]*Weights[i]*ControlPoints[i][1];
                }
            }
            // Second Span
            if(bounds[1][0] == spans[j][0] && bounds[1][1] == spans[j][1]) { 
                for(int k = 0;k < numerCoeffsX[j].size();k++) {
                    denomCoeffs[j][k] += coeffs[1][k]*Weights[i];
                    numerCoeffsX[j][k] += coeffs[1][k]*Weights[i]*ControlPoints[i][0];
                    numerCoeffsY[j][k] += coeffs[1][k]*Weights[i]*ControlPoints[i][1];
                }
            }
            // Third Span
            if(bounds[2][0] == spans[j][0] && bounds[2][1] == spans[j][1]) { 
                for(int k = 0;k < numerCoeffsX[j].size();k++) {
                    denomCoeffs[j][k] += coeffs[2][k]*Weights[i];
                    numerCoeffsX[j][k] += coeffs[2][k]*Weights[i]*ControlPoints[i][0];
                    numerCoeffsY[j][k] += coeffs[2][k]*Weights[i]*ControlPoints[i][1];
                }
            }

        }
    }

    return denomCoeffs;
}

std::vector<double> Spline::makeBreakpoints() {
    breakpoints = {KnotVector[0]};
    int lenBreak = 1;
    for(int i = 1; i < KnotVector.size();i++) {
        if(breakpoints[lenBreak-1] != KnotVector[i]) {
            breakpoints.insert(breakpoints.end(),KnotVector[i]);
            lenBreak++;
        }
    }
    return breakpoints;
}

std::vector<std::vector<double>> Spline::makeSpans() {
    spans = {{breakpoints[0],breakpoints[1]}};
    for(int i = 1; i< breakpoints.size()-1;i++) {
        spans.insert(spans.end(),{breakpoints[i],breakpoints[i+1]});
    }
    return spans;
}

double Spline::getArcLength() { // Testing Needed
    double nudge = 1e-12;
    // Three Point Quadrature
    // std::vector<double> GaussPoints = {-sqrt(3/5),0,sqrt(3/5)};
    // std::vector<double> GaussWeights = {5/9,8/9,5/9};
    // Five Point Quadrature
    std::vector<double> GaussPoints = {-sqrt(5+2*sqrt(10/7))/3,-sqrt(5-2*sqrt(10/7))/3,0,sqrt(5-2*sqrt(10/7))/3,sqrt(5+2*sqrt(10/7))/3};
    std::vector<double> GaussWeights = {(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900};
    double AL = 0;
    double u1;
    double u2;
    for(int i = 0;i < breakpoints.size()-1;i++) {
        // Determine integrating bounds for segment
        if(i == 0) {
            u1 = breakpoints[i];
            u2 = breakpoints[i+1]-nudge;
        } else if(i == breakpoints.size()-1) {
            u1 = breakpoints[i]+nudge;
            u2 = breakpoints[i+1];
        } else {
            u1 = breakpoints[i]+nudge;
            u2 = breakpoints[i+1]-nudge;
        }

        // Change Endpoints
        double m = (u2-u1)/2;
        double bm = (u1+u2)/2;
        for(int j = 0; j < GaussPoints.size(); j++) {
            double ueff = m*GaussPoints[j] + bm;

            // Coefficients
            double a = numerCoeffsX[i][0];
            double b = numerCoeffsX[i][1];
            double c = numerCoeffsX[i][2];

            double d = numerCoeffsY[i][0];
            double e = numerCoeffsY[i][1];
            double f = numerCoeffsY[i][2];

            double alpha = denomCoeffs[i][0];
            double beta = denomCoeffs[i][1];
            double gamma = denomCoeffs[i][2];

            // Derivative Coefficients
            double ax = a*beta-b*alpha;
            double bx = 2*(a*gamma-alpha*c);
            double cx = b*gamma-beta*c;

            double ay = d*beta-e*alpha;
            double by = 2*(d*gamma-alpha*f);
            double cy = e*gamma-beta*f;

            // Normal Vector Magnitude
            double nx = (ax*pow(ueff,2)+bx*ueff+cx)/pow(alpha*pow(ueff,2)+beta*ueff+gamma,2);
            double ny = (ay*pow(ueff,2)+by*ueff+cy)/pow(alpha*pow(ueff,2)+beta*ueff+gamma,2);

            double n = sqrt(pow(nx,2)+pow(ny,2));

            AL += m*GaussWeights[j]*n;
        }
    }
    return AL;
}

double Spline::getCurvature(double u) {
    
    // Find Span we are in
    int spanIndex = -1;
    if(u > 1) {// Too Large Error
        spanIndex = -1; 
    }else if(u < 0) {// Negative Error
        spanIndex = -2;
    } else if(u == 1) { // Edge Case
        spanIndex = breakpoints.size()-2;
    } else {
        for(int i = 0; i <breakpoints.size()-1;i++) {
            if(u >= breakpoints[i] && u < breakpoints[j+1]) {
                spanIndex = i;
                break;
            }
        }
    }

    // 0th Derivative Coefficients
    double a = numerCoeffsX[spanIndex][0];
    double b = numerCoeffsX[spanIndex][1];
    double c = numerCoeffsX[spanIndex][2];

    double d = numerCoeffsY[spanIndex][0];
    double e = numerCoeffsY[spanIndex][1];
    double f = numerCoeffsY[spanIndex][2];

    double alpha = denomCoeffs[spanIndex][0];
    double beta = denomCoeffs[spanIndex][1];
    double gamma = denomCoeffs[spanIndex][2];

    // 1st Derivative Coefficients
    double ax = a*beta-b*alpha;
    double bx = 2*(a*gamma-alpha*c);
    double cx = b*gamma-beta*c;

    double ay = d*beta-e*alpha;
    double by = 2*(d*gamma-alpha*f);
    double cy = e*gamma-beta*f;

    // First Derivatives 
    double xp = (ax*pow(u,2)+bx*u+cx)/pow(alpha*pow(u,2)+beta*u+gamma,2);
    double yp = (ay*pow(u,2)+by*u+cy)/pow(alpha*pow(u,2)+beta*u+gamma,2);

    // Second Derivatives

    // X
    double denom = pow(alpha*pow(u,2)+beta*u+gamma,4);
    double term1 = (2*ax*u+bx)*pow(alpha*pow(u,2)+beta*u+gamma,2);
    double term2 = 2*(alpha*pow(u,2)+beta*u+gamma)*(2*alpha*u+beta)*(ax*pow(u,2)+bx*u+cx);
    double xpp = (term1-term2)/denom;

    // Y
    term1 = (2*ay*u+by)*pow(alpha*pow(u,2)+beta*u+gamma,2);
    term2 = 2*(alpha*pow(u,2)+beta*u+gamma)*(2*alpha*u+beta)*(ay*pow(u,2)+by*u+cy);
    double ypp = (term1-term2)/denom;   

    // Curvature
    term1 = xp*ypp-yp*xpp;
    term2 = pow(pow(xp,2)+pow(yp,2),3/2);

    double k = term1/term2;
    return k;
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
std::vector<double> Spline::getBreakpoints() {
    return breakpoints;
}
std::vector<std::vector<double>> Spline::getSpans() {
    return spans;
}
std::vector<std::vector<double>> Spline::getXCoeffs() {
    return numerCoeffsX;
}
std::vector<std::vector<double>> Spline::getYCoeffs() {
    return numerCoeffsY;
}
std::vector<std::vector<double>> Spline::getDCoeffs() {
    return denomCoeffs;
}

// Setters
void Spline::setControlPoints(std::vector<std::vector<double>> input) {
    ControlPoints = input;
}
void Spline::setKnotVector(std::vector<double> input) {
    KnotVector = input;
    this->makeBreakpoints();
    this->makeSpans();
}
void Spline::setWeights(std::vector<double> input){
    Weights = input;
    
}