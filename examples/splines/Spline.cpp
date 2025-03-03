#include "Spline.h"
#include <cmath>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

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
std::vector<double> Spline::BesselTangentUVec(std::vector<std::vector<double>> Q) {// Testing Needed
    int n = Q.size()-1;
    double d = 0;
    std::vector<double> dset = {0};
    // Make Distance Vector
    for(int i =1; i <Q.size();i++){
        double mag = 0;
        for(int j = 0;j < Q[0].size();j++) {
            mag += pow(Q[i][j]-Q[i-1][j],2);
        } 
        mag = sqrt(mag);
        if(i == 1) {
            dset[0] == mag;
        } else {
            dset.insert(dset.end(),mag);
        }
        d += mag;
    }
    // Make U from Distance Vector
    std::vector<double> ubar = {0};
    for(int i = 1; i < n+1; i++) {
        ubar.insert(ubar.end(),ubar[i-1]+dset[i-1]/d);
    }

    return ubar;
}

std::vector<std::vector<double>> Spline::BesselTangents(std::vector<std::vector<double>> Q) { // Testing Needed
    int n = Q.size()-1;

    std::vector<double> ubar = BesselTangentUVec(Q);

    std::vector<double> dUbar = {0};
    
    std::vector<std::vector<double>> q = {{0,0}};
    for(int i = 0;i<n;i++) {
        if(i == 0) {
            dUbar[i] = ubar[i+1]-ubar[i];
            for(int j = 0; j < Q[0].size(); j++) {
                q[i][j] = Q[i+1][j]-Q[i][j];
            }
        } else {
            dUbar.insert(dUbar.end(),ubar[i+1]-ubar[i]);
            std::vector<double> temp = {0,1};
            for(int j = 0; j < temp.size(); j++) {
                temp[j] = Q[i+1][j]-Q[i][j];
            }
            q.insert(q.end(),temp);
        }
    }
    for(int i = 0; i < q.size(); i++) {
        for(int j = 0; j < q[0].size(); j++) {
            q[i][j] = q[i][j]/dUbar[i];
        }
    }

    // Calculate Alphas
    std::vector<double> alpha = {0};
    for(int i = 1; i<n;i++) { 
        double numer = dUbar[i-1];
        double denom = dUbar[i-1]-dUbar[i];

        alpha.insert(alpha.end(),numer/denom);
    }

    // Calculate D
    std::vector<std::vector<double>> D = {{0,0}};
    std::vector<double> temp = {0,0};
    for(int i = 1; i < n; i++) {
        temp = {0,0};
        for(int j = 0; j <temp.size();j++) {
            temp[j] = (1-alpha[i])*q[i-1][j]+alpha[i]*q[i][j];
        }
        
        D.insert(D.end(),temp);
    }
    // D[0]
    temp = {0,0};
    for(int j = 0; j <temp.size();j++) {
        temp[j] = 2*q[0][j] - D[1][j];
    }
    D[0] = temp;
    // D[end]
    temp = {0,0};
    for(int j = 0; j <temp.size();j++) {
        temp[j] = 2*q[q.size()-1][j] - D[D.size()-1][j];
    }
    D.insert(D.end(),temp);

    //  Loop through and normalize
    for(int i = 0; i < D.size(); i++) {
        double mag = 0;
        for(int j = 0; j<D[0].size();j++) {
            mag += pow(D[i][j],2);
        }
        mag = sqrt(mag);

        for(int j = 0; j<D[0].size();j++) {
            D[i][j] /= mag; 
        }
    }

    // Return
    return D;
}

//Assumes 2D
std::vector<std::vector<double>> Spline::solvePointTangentIntersection(std::vector<double> Q1, 
                                                          std::vector<double> Q2,
                                                          std::vector<double> T1,
                                                          std::vector<double> T2) { // Testing Needed
        double detA = T2[0]*T1[1]-T2[1]*T1[0];
        double dQx = Q1[0]-Q2[0];
        double dQy = Q1[1]-Q2[1];
        
        double g1 = (T2[1]*dQx-T2[0]*dQy)/detA;
        double g2 = (T1[1]*dQx-T2[0]*dQy)/detA;
        std::vector<double> R = {0,0};
        std::vector<double> g = {g1,g2};
        for(int i = 0;i <R.size();i++){
            R[i] = Q1[i] + g1*T1[i];
        }
        return {R,g};
    }

//Assumes 2D
double Spline::makeWeight(std::vector<double> Qkm, std::vector<double> Rk, std::vector<double> Qk) { // Testing Needed
    std::vector<double> dQmR = {0,0};
    std::vector<double> dQR = {0,0};
    // Set up Dot Product
    double dQmRnorm = 0;
    double dQRnorm = 0;
    double dotProd = 0;
    for(int i = 0;i < dQR.size(); i++) {
        dQmR[i] = Rk[i]-Qkm[i];
        dQmRnorm += pow(dQmR[i],2);

        dQR[i] = Qk[i]-Rk[i];
        dQRnorm += pow(dQR[i],2);

        dotProd += dQmR[i]*dQR[i];
    }
    dotProd /= dQmRnorm*dQRnorm;
    double tolerance = 1e-6;
    double weight = -1;
    if(abs(dotProd-1) < tolerance) { // Collinear
        weight = 1;
    } else {
        if(abs(dQmRnorm-dQRnorm) < tolerance) { // Isosceles
            std::vector<double> M = {(Qkm[0]+Qk[0])/2,(Qkm[1]+Qk[1])/2};
            double eVal = 0;
            for(int i = 0;i<M.size();i++) {
                M[i] = Qk[i] - M[i];
                eVal += pow(M[i],2);
            }
            weight = eVal/dQmRnorm;
        } else { // Other Case
            std::vector<double> M = {(Qkm[0]+Qk[0])/2,(Qkm[1]+Qk[1])/2};
            // S1
            std::vector<double> MR = {0,0};
            double MRnorm = 0;
            std::vector<double> MQm = {0,0};
            double MQmnorm = 0;
            std::vector<double> QRm = {0,0};
            double QRmnorm = 0;
            // S2
            std::vector<double> MQ = {0,0};
            double MQnorm = 0;
            std::vector<double> QR = {0,0};
            double QRnorm = 0;
            for(int i = 0; i < MR.size(); i++) {
                // S1
                MR[i] = Rk[i] - M[i];
                MRnorm += pow(MR[i],2);
                MQm[i] = M[i] - Qkm[i];
                MQmnorm += pow(MQm[i],2);
                QRm[i] = Rk[i] - Qkm[i];
                QRmnorm += pow(QRm[i],2);
                // S2
                MQ[i] = M[i] - Qk[i];
                MQnorm += pow(MQ[i],2);
                QR[i] = Rk[i] - Qk[i];
                QRnorm += pow(QR[i],2);
            }
            // normalize
            MRnorm = sqrt(MRnorm);
            MQmnorm = sqrt(MQmnorm);
            QRmnorm = sqrt(QRmnorm);
            MQnorm = sqrt(MQnorm);
            QRnorm = sqrt(QRnorm);

            for(int i = 0; i<MR.size();i++){
                // S1
                MR[i] /= MRnorm;
                MQm[i] /= MQmnorm;
                QRm[i] /= QRmnorm;
                // S2
                MQ[i] /= MQnorm;
                QR[i] /= QRnorm;   
            }
            // Bisection Vectors
            std::vector<double> BiVec1 = {0,0};
            double BiVec1norm = 0;
            std::vector<double> BiVec2 = {0,0};
            double BiVec2norm = 0;
            for(int i = 0; i < MR.size(); i++) {
                // S1
                BiVec1[i] = QRm[i] + MQm[i];
                BiVec1norm += BiVec1[i];
                // S2
                BiVec2[i] = MQ[i]+QR[i];
                BiVec2norm += BiVec2[i];
            }
            // Normalize
            BiVec1norm = sqrt(BiVec1norm);
            BiVec2norm = sqrt(BiVec2norm);
            for(int i = 0; i < MR.size(); i++) {
                // S1
                BiVec1[i] /= BiVec1norm;
                // S2
                BiVec2[i] /= BiVec2norm;
            }
            // Solve for S1,S2
            std::vector<double> S1 = Spline::solvePointTangentIntersection(Qkm,BiVec1,M,MR)[0];
            std::vector<double> S2 = Spline::solvePointTangentIntersection(Qk,BiVec2,M,MR)[0];
            // Get S
            std::vector<double> S = {(S1[0]+S2[0])/2,(S1[1]+S2[1])/2};
            // Make Weight
            double numer = 0;
            double denom = 0;
            for(int i = 0; i < S.size(); i++) {
                numer += pow(M[i]-S[i],2);
                denom += pow(Rk[i]-S[i],2);
            }
            weight = sqrt(numer/denom);
        }
    }
    return weight;
}

Spline Spline::LocalRQuadInterp(std::vector<std::vector<double>> Q,std::vector<std::vector<double>> T) { // Testing and Debugging Needed
    int n = Q.size()-1;
    std::vector<std::vector<double>> CPoints = {Q[1]};
    std::vector<std::vector<double>> gamma = {{0,0}};
    std::vector<double> weight = {1};
    
    for(int i =0;i<=n;i++) { // Loop over consecutive points.

        // Involved Tangents and Points
        std::vector<double> Tk = T[i+1];
        double normTk = 0;

        std::vector<double> Tkm = T[i];
        double normTkm = 0;

        std::vector<double> Qk = Q[i+1];
        std::vector<double> Qkm = Q[i];

        // Normalize Tangents
        for(int j = 0; j < Tk.size();j++) {
            normTk += pow(Tk[j],2);
            normTkm += pow(Tkm[j],2);
        }
        normTk = sqrt(normTk);
        normTkm = sqrt(normTkm);
        for(int j = 0; j < Tk.size();j++) {
            Tk[j] /= normTk;
            Tkm[j] /= normTkm;
        }

        // Dot Product of Tangents
        double DotProductValue = 0;
        for(int j = 0; j < Tk.size();j++) {
            DotProductValue += Tk[j]*Tkm[j];
        }

        // Check for Parallel
        if(abs(abs(DotProductValue)-1) <= 1e-6) { // Parallel Case
            std::vector<double> dQVec = {Qk[0] - Qkm[0],Qk[1] - Qkm[1]};
            double dQnorm = sqrt(pow(dQVec[0],2) + pow(dQVec[1],2));
            dQVec[0] /= dQnorm;
            dQVec[1] /= dQnorm;

            // Dot Products
            double Dk = dQVec[0]*Tk[0]+dQVec[1]*Tk[1];
            if(abs(abs(Dk)-1) <= 1e-6) { // Straight Line Case
                std::vector<double> Rk = {(Qkm[0]+Qk[0])/2,(Qkm[1]+Qk[1])/2};
                CPoints.insert(CPoints.end(),Rk);
                CPoints.insert(CPoints.end(),Qk);

                double w1 = makeWeight(Qkm,Rk,Qk);
                weight.insert(weight.end(),w1);
                weight.insert(weight.end(),1);
            } else { // Parallel Tangents, not Parallel to Chord
                double gammak = 0.5*dQnorm;
                double gammakp = gammak;
                std::vector<double> Rkp = {Qkm[0]+gammak*Tkm[0],Qkm[1]+gammak*Tkm[1]};
                std::vector<double> Rkpp = {Qk[0]-gammakp*Tk[0],Qk[1]-gammakp*Tk[1]};
                std::vector<double> Qkp = {(gammak*Rkpp[0]+gammakp*Rkp[0])/(gammak+gammakp),
                                           (gammak*Rkpp[1]+gammakp*Rkp[1])/(gammak+gammakp)};
                
                double w1 = makeWeight(Qkm,Rkp,Qkp);
                double w2 = makeWeight(Qkp,Rkpp,Qk);
                weight.insert(weight.end(),{w1,1,w2,1});
                CPoints.insert(CPoints.end(),{Rkp,Qkp,Rkpp,Qk});
            }
        } else { // Not Parllel
            std::vector<std::vector<double>> sol = solvePointTangentIntersection(Qkm,Qk,Tkm,Tk);
            std::vector<double> Rk1 = sol[0];
            std::vector<double> g = sol[1];

            if(g[0] <= 1e-12 || g[1] >= -1e12) { // Solution Exists, but not in bounds (Inflection Point or 180 turn)
                std::vector<double> dQVec = {Qk[0] - Qkm[0],Qk[1] - Qkm[1]};
                double dQnorm = sqrt(pow(dQVec[0],2) + pow(dQVec[1],2));
                dQVec[0] /= dQnorm;
                dQVec[1] /= dQnorm;

                // Dot Product
                double Dk = dQVec[0]*Tk[0]+dQVec[1]*Tk[1];
                double Dkm = dQVec[0]*Tkm[0]+dQVec[1]*Tkm[1];

                //Angles
                double thetak = abs(acos(Dk));
                double thetakm = abs(acos(Dkm));

                double a = 2/3; // Tunable Parameter

                double numer = dQnorm;
                double denom1 = 2*(1+a*cos(thetak)+(1-a)*cos(thetakm));
                double denom2 = 2*(1+a*cos(thetakm)+(1-a)*cos(thetak));

                double gammak = numer/denom1;
                double gammakp = numer/denom2;

                std::vector<double> Rkp = {Qkm[0]+gammak*Tkm[0],Qkm[1]+gammak*Tkm[1]};
                std::vector<double> Rkpp = {Qk[0]-gammakp*Tk[0],Qk[1]-gammakp*Tk[1]};
                std::vector<double> Qkp = {(gammak*Rkpp[0]+gammakp*Rkp[0])/(gammak+gammakp),
                                           (gammak*Rkpp[1]+gammakp*Rkp[1])/(gammak+gammakp)};
                double w1 = makeWeight(Qkm,Rkp,Qkp);
                double w2 = makeWeight(Qkp,Rkpp,Qk);
                weight.insert(weight.end(),{w1,1,w2,1});

                CPoints.insert(CPoints.end(),{Rkp,Qkp,Rkpp,Qk});
            } else { // Normal Case
                double w1 = makeWeight(Qkm,Rk1,Qk);

                weight.insert(weight.end(),{w1,1});
                CPoints.insert(CPoints.end(),{Rk1,Qk});
            }
        }
    }
    // At this point, weight and CPoints are made. The last thing we need is the knot vector
    // Split the Q and R sets
    std::vector<std::vector<double>> Qset = {CPoints[0]};
    std::vector<std::vector<double>> Rset = {CPoints[1]};
    for(int i = 2; i < CPoints.size();i+=2){ 
        Qset.insert(Qset.end(),CPoints[i]); // Get every other control point
        Rset.insert(Rset.end(),CPoints[i+1]); // Get every other other control point.
    }
    n = Qset.size()-1;
    // Make ubar Helper
    std::vector<double> ubar= {0,1};
    for(int i = 2;i < n+1; i++) {
        double coeff = ubar[i-1]-ubar[i-2];

        double numer = sqrt(pow(Rset[i-1][0]-Qset[i-1][0],2)+pow(Rset[i-1][1]-Qset[i-1][1],2));
        double denom = sqrt(pow(Qset[i-1][0]-Rset[i-2][0],2)+pow(Qset[i-1][1]-Rset[i-2][1],2));

        ubar.insert(ubar.end(),ubar[i-1]+coeff*numer/denom);
    }
    double un = ubar[ubar.size()-1];
    // Make Knot Vector
    std::vector<double> U = {0,0,0};
    for(int i = 1;i < ubar.size();i++) {
        U.insert(U.end(),{ubar[i]/un,ubar[i]/un});
    }
    U.insert(U.end(),1);
    // Finally Make the Spline Object
    Spline ret = Spline(CPoints,U,weight);
    return ret;
}

// Dynamic Methods **************************
// Assuming 2D
int Spline::findSpan(double u) {
    int spanIndex = -1;
    if(u > 1) {// Too Large Error
        spanIndex = -1; 
    }else if(u < 0) {// Negative Error
        spanIndex = -2;
    } else if(u == 1) { // Edge Case
        spanIndex = breakpoints.size()-2;
    } else {
        for(int i = 0; i <breakpoints.size()-1;i++) {
            if(u >= breakpoints[i] && u < breakpoints[i+1]) {
                spanIndex = i;
                break;
            }
        }
    }
    return spanIndex;
}

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

double Spline::getSurfaceEnergy() { // Testing Needed
    double nudge = 1e-12;
    // Three Point Quadrature
    // std::vector<double> GaussPoints = {-sqrt(3/5),0,sqrt(3/5)};
    // std::vector<double> GaussWeights = {5/9,8/9,5/9};
    // Five Point Quadrature
    std::vector<double> GaussPoints = {-sqrt(5+2*sqrt(10/7))/3,-sqrt(5-2*sqrt(10/7))/3,0,sqrt(5-2*sqrt(10/7))/3,sqrt(5+2*sqrt(10/7))/3};
    std::vector<double> GaussWeights = {(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900};
    double Ek = 0;
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
            
            // Curvature
            double k = this->getCurvature(ueff);
            Ek += m*GaussWeights[j]*n*abs(k); // Unsigned Surface Energy is the real one
        }
    }
    return Ek;
}

double Spline::getCurvature(double u) { // Testing Needed
    
    // Find Span we are in
    int spanIndex = findSpan(u);

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

std::vector<std::vector<double>> Spline::makeRationalQuadCurve(std::vector<double> uset) { // Testing Needed
    std::vector<std::vector<double>> curve = {{0,0}};
    for(int i = 0; i < uset.size();i++) { // Loop over all u values
        std::vector<double> numer = {0,0};
        double denom = 0;
        std::vector<double> val = {0,0};
        for(int j = 0; j < ControlPoints.size();j++){ // Loop over control points
            for(int k = 0; k < ControlPoints[0].size();k++) { // Loop over contorl point indices
                numer[k] += BBasisFunction(j,2,uset[i])*Weights[j]*ControlPoints[j][k];
                denom += BBasisFunction(j,2,uset[i])*Weights[j];
            }
        }

        // Now that we have numer and denom, calculate value
        for(int j = 0; j < numer.size();j++) {
            val[j] = numer[j]/denom;
        }

        // Now insert val
        if(i == 0) {
            curve[i] = val;
        } else {
            curve.insert(curve.end(),val);
        }
    }
    return curve;
}

double Spline::integratedSpline(double u) {// Testing Needed
    // Find Span we are in
    int spanIndex = findSpan(u);

    // Get Appropriate Coefficients
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

    // Numerator Coefficients
    double term1 = a*ay;
    double term2 = a*by+b*ay;
    double term3 = a*cy+b*by+c*ay;
    double term4 = b*cy+c*by;
    double term5 = c*cy;

    // Rename to match Mathematica
    double a2 = term1;
    double b2 = term2;
    double c2 = term3;
    double d2 = term4;
    double e2 = term5;
    double f2 = alpha;
    double g = beta;
    double h = gamma;

    double N1B = c2*pow(f2,2)*pow(g,3) - b2*f2*pow(g,4) + a2*pow(g,5) +
                2*c2*pow(f2,3)*g*h + 5*b2*pow(f2,2)*pow(g,2)*h - 8*a2*f2*pow(g,3)*h
                - 16*b2*pow(f2,3)*pow(h,2) + 22*a2*pow(f2,2)*g*pow(h,2) +
                2*c2*pow(f2,3)*pow(g,2)*u - 2*a2*f2*pow(g,4)*u + 4*c2*pow(f2,4)*h*u
                - 6*b2*pow(f2,3)*g*h*u + 16*a2*pow(f2,2)*pow(g,2)*h*u -
                20*a2*pow(f2,3)*pow(h,2)*u + 6*e2*pow(f2,4)*(g + 2*f2*u) -
                3*d2*pow(f2,3)*g*(g + 2*f2*u);

    double D1B = (pow(f2,3)*pow(pow(g,2) - 4*f2*h,2)*(h + u*(g + f2*u)));

    double N2B = (c2*pow(f2,2)*g*h - b2*f2*pow(g,2)*h + a2*pow(g,3)*h +
                2*b2*pow(f2,2)*pow(h,2) - 3*a2*f2*g*pow(h,2) +
                c2*pow(f2,2)*pow(g,2)*u - b2*f2*pow(g,3)*u + a2*pow(g,4)*u -
                2*c2*pow(f2,3)*h*u + 3*b2*pow(f2,2)*g*h*u - 4*a2*f2*pow(g,2)*h*u +
                2*a2*pow(f2,2)*pow(h,2)*u + e2*pow(f2,3)*(g + 2*f2*u) -
                d2*pow(f2,3)*(2*h + g*u));
    
    double D2B = (pow(f2,3)*(-pow(g,2) + 4*f2*h)*pow(h + u*(g + f2*u),2));

    double N3B = (4*(6*e2*pow(f2,2) - 3*d2*f2*g + c2*pow(g,2) + 2*c2*f2*h - 3*b2*g*h +
                6*a2*pow(h,2))*atan((g + 2*f2*u)/sqrt(-pow(g,2) + 4*f2*h)));

    double D3B = pow(-pow(g,2) + 4*f2*h,2.5);

    double value = 0.5*(N1B/D1B + N2B/D2B + N3B/D3B);

    return value;
}

double Spline::getArea() { // Testing Needed
    double Area = 0;
    double nudge = 1e-8;
    // Loop over each breakpoint, find area of that section, then add them together.
    double uL;
    double uR;
    for(int i =0; i <breakpoints.size()-1; i++) {
        uL = breakpoints[i]; // Left Breakpoint
        uR = breakpoints[i+1]; // Right Breakpoint

        uL += nudge;
        uR -= nudge;

        Area += integratedSpline(uR) - integratedSpline(uL);
    }
    return Area;
}

std::vector<double> Spline::lineCurveIntersection(std::vector<double> P1, std::vector<double> P2) { // Testing Needed
    std::vector<double> tangent = {P2[0]-P1[0],P2[1]-P1[1]};
    std::vector<double> normal = {-tangent[1],tangent[0]};
    int numSpans = breakpoints.size()-1;
    std::vector<double> uIntersections = {0};
    std::vector<double> sIntersections = {0};
    double a = normal[0];
    double b = normal[1];
    double c = -(normal[0]*P1[0]+normal[1]*P1[1]);
    int count = 0;
    for(int i = 0; i < numSpans; i ++) { // Loop Over Spans
        // Quadratic Coefficients
        double aPoly = a*numerCoeffsX[i][0] + b*numerCoeffsY[i][0]+c*denomCoeffs[i][0];
        double bPoly = a*numerCoeffsX[i][1] + b*numerCoeffsY[i][1]+c*denomCoeffs[i][1];
        double cPoly = a*numerCoeffsX[i][2] + b*numerCoeffsY[i][2]+c*denomCoeffs[i][2];
        // Solve Quadratic
        double discrim = pow(bPoly,2)-4*aPoly*cPoly;
        // Filter Values for real, in span on curve,and in spans of line
        if(discrim >= 0) { // Real Values Only
            double u1 = (-bPoly + sqrt(discrim))/(2*aPoly);
            double u2 = (-bPoly - sqrt(discrim))/(2*aPoly);
            std::vector<double> Inters = {u1,u2};
            std::vector<std::vector<double>> pointInter = this->makeRationalQuadCurve(Inters);
            // Check u1
            for(int j = 0; j < Inters.size();j++) {
                u1 = Inters[j];
                if(u1 >= spans[i][0] && u1 <= spans[i][1]) { // Within Span of curve
                    std::vector<double> Pcurr = pointInter[j];
                    std::vector<double> dP1 = {Pcurr[0]-P1[0],Pcurr[1]-P1[1]};
                    std::vector<double> dP2 = {Pcurr[0]-P2[0],Pcurr[1]-P2[1]};
                    std::vector<double> dP = {P1[0]-P2[0],P1[1]-P2[1]};
                    double lineCheck;
                    if(abs(dP[0]) > abs(dP[1])) { // Checks if point between endpoints of line
                        lineCheck = abs(dP1[0])+abs(dP2[0])-abs(dP[0]);
                    } else {
                        lineCheck = abs(dP1[1])+abs(dP2[1])-abs(dP[1]);
                    }
                    // If point on line, add
                    if(lineCheck <= 1e-8) {
                        if(count == 0){ // If first one, replace first value
                            uIntersections[0] = u1;
                            sIntersections[0] = std::min(abs(dP1[0]/dP[0]),abs(dP1[1]/dP[1]));
                            count++;
                        } else { // else, add to end
                            uIntersections.insert(uIntersections.end(),u1);
                            sIntersections.insert(sIntersections.end(),std::min(abs(dP1[0]/dP[0]),abs(dP1[1]/dP[1])));
                        }
                    }
                }
            }
        }
        return uIntersections;
    }
    // At this point we have finished looping over the spans. We now reorder the values to be along the line
    // This section of code was taken from https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
    std::vector<size_t> idx(sIntersections.size());
    std:: iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&sIntersections](size_t i1, size_t i2) {return sIntersections[i1] < sIntersections[i2];});
    // Now back to me, rearrange uInters to be along the line
    std::vector<double> ret(uIntersections.size());
    for(int i = 0; i< idx.size(); i++) {
        ret[i] = uIntersections[idx[i]];
    }   
    return ret;
}

std::vector<std::vector<double>> Spline::getParameterLoop(std::vector<std::vector<double>> square) { // Testing Needed
    std::vector<double> P1 = square[0];
    std::vector<double> P2;
    int I = 0;
    for(int i = 0; i < ControlPoints.size();i++) { // Loop over control points
        if(ControlPoints[i][0] > ControlPoints[I][0]) { // Find Furthest Right Point
            I = i;
        }
    }
    P2 = ControlPoints[I];
    P2[0] += 1; // Move a little to the right
    std::vector<double> intersections = lineCurveIntersection(P1,P2);
    
    std::vector<double> parameter;
    std::vector<double> indicator; // 1 = corner (outside), 2 = corner (inside), 3 = Curve hit (entry), 4 = curve hit (Exit)
    // If number of intersections even, outside. If odd, inside.
    if(intersections.size() % 2 == 0) {
        parameter = {0};
        indicator = {1};
    } else {
        parameter = {0};
        indicator = {2};
    }

    for(int i = 0; i < square.size()-1;i++) { // Square includes start twice, so we have to exclude the last one
        P1 = square[i];
        P2 = square[i+1];

        intersections = lineCurveIntersection(P1,P2);

        // Assign Intersections
        for(int j = 0;j < intersections.size();j++) {
            parameter.insert(parameter.end(),intersections[j]);
            double lastInd = indicator[indicator.size()-1];
            if(lastInd == 1) {
                indicator.insert(indicator.end(),3);
            } else if(lastInd == 2) {
                indicator.insert(indicator.end(),4);
            } else if(lastInd == 3) {
                indicator.insert(indicator.end(),4);
            } else if(lastInd == 4) {
                indicator.insert(indicator.end(),3);
            }
        }

        // Assign Endpoint
        parameter.insert(parameter.end(),(i+1)%4);
        int lastInd = indicator[indicator.size()-1];
        if(lastInd == 1) {
            indicator.insert(indicator.end(),1);
        } else if(lastInd == 2) {
            indicator.insert(indicator.end(),2);
        } else if(lastInd == 3) {
            indicator.insert(indicator.end(),2);
        } else if(lastInd == 4) {
            indicator.insert(indicator.end(),1);
        }
    }
    // Pop 1s
    std::vector<double> tempParameters = {0};
    std::vector<double> tempIndicators = {0};
    int count = 0;
    for(int i = 0; i < indicator.size();i++) {
        if(indicator[i] != 1) {
            if(count == 0) {
                tempParameters[0] = parameter[i];
                tempIndicators[0] = indicator[i];
                count++;
            } else {
                tempParameters.insert(tempParameters.end(),parameter[i]);
                tempIndicators.insert(tempIndicators.end(),indicator[i]);
            }
        }
    }
    if(indicator[0] == 1 && count != 0) {
        tempParameters.insert(tempParameters.end(),tempParameters[0]);
        tempParameters.insert(tempParameters.end(),tempParameters[0]);
    }

    return {parameter,indicator};
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