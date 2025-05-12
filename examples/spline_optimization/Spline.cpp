#include "./Spline.h"
#include <cmath>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <complex>
#include <math.h>
#include <array>
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
// TEST ALL STATIC METHODS ***************************************************************************************
// Packing for Optimization 
std::vector<double> Spline::pack(std::vector<std::vector<double>> Q,std::vector<std::vector<double>> T) {
    int Nq = Q.size();
    int Nt = T.size();

    if(Nq == Nt) {// Needed
        std::vector<double> ret(2*(Nq-1)+(Nt-1));
        for(int i = 0;i < Nq-1; i++) {
            ret[3*i] = Q[i][0];
            ret[3*i+1] = Q[i][1];
            ret[3*i+2] = double(std::atan2(T[i][1],T[i][0])); // Convert into angle, from -pi to pi.
        }
        // std::cout << "Size ret = " << ret.size() << "\n [";
        // for(int i = 0; i < ret.size(); i++) {
        //     std::cout << ret[i] << ",";
        // }
        // std::cout << "]\n";
        return ret;
    } else {
        return {-1};
    }
    return {0};
}
// Unpacking from Optimization to Interpolation
std::vector<std::vector<std::vector<double>>> Spline::unpack(std::vector<double> V) {
    int N = V.size();
    int Nq = int(N/3)+1;
    // std::cout << "Nq = " << Nq << "\n";
    if(N%4 == 0) { // Required
        std::vector<std::vector<double>> Q(Nq,std::vector<double> (2));
        std::vector<std::vector<double>> T(Nq,std::vector<double> (2));
        for(int i = 0; i < Q.size(); i++) {
            Q[i][0] = V[3*i];
            Q[i][1] = V[3*i+1];
            T[i][0] = std::cos(V[3*i+2]);
            T[i][1] = std::sin(V[3*i+2]);
        }
        // Make sure we are periodic
        Q[Nq-1][0] = Q[0][0];
        Q[Nq-1][1] = Q[0][1];
        T[Nq-1][0] = T[0][0];
        T[Nq-1][1] = T[0][1];
        return {Q,T};
    } else {
        return {{{-1}}};
    }
    return {{{0}}};
    
}

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
                                                          std::vector<double> T2) {
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
double Spline::makeWeight(std::vector<double> Qkm, std::vector<double> Rk, std::vector<double> Qk) { 
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
    dQmRnorm = sqrt(dQmRnorm);
    dQRnorm = sqrt(dQRnorm);

    dotProd /= dQmRnorm*dQRnorm;
    double tolerance = 1e-6;
    double weight = -1;
    if(fabs(dotProd-1) < tolerance) { // Collinear
        // std::cout << "Collinear\n";
        weight = 0;
    } else {
        if(fabs(dQmRnorm-dQRnorm) < tolerance) { // Isosceles
            // std::cout << "Isosceles\n";
            std::vector<double> M = {(Qkm[0]+Qk[0])/2,(Qkm[1]+Qk[1])/2};
            double eVal = 0;
            for(int i = 0;i<M.size();i++) {
                M[i] = Qk[i] - M[i];
                eVal += pow(M[i],2);
            }
            eVal = sqrt(eVal);
            weight = eVal/dQmRnorm;
        } else { // Other Case
            // std::cout << "Other\n";
            std::vector<double> M = {(Qkm[0]+Qk[0])/2,(Qkm[1]+Qk[1])/2};
            // std::cout << "Qk = " << Qk[0] << "," << Qk[1] << "\n";
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
                BiVec1norm += pow(BiVec1[i],2);
                // S2
                BiVec2[i] = MQ[i]+QR[i];
                BiVec2norm += pow(BiVec2[i],2);
            }
            // Normalize
            BiVec1norm = sqrt(BiVec1norm);
            BiVec2norm = sqrt(BiVec2norm);
            // std::cout <<"BiVec1 = " <<BiVec1[0] << "," << BiVec1[1] << "\n";
            // std::cout << "BiVec1norm = " << BiVec1norm << "\n";
            // std::cout << "BiVec2norm = " << BiVec2norm << "\n";
            for(int i = 0; i < BiVec1.size(); i++) {
                // S1
                BiVec1[i] /= BiVec1norm;
                // S2
                BiVec2[i] /= BiVec2norm;
            }
            // Solve for S1,S2
            std::vector<double> S1 = Spline::solvePointTangentIntersection(Qkm,M,BiVec1,MR)[0];
            std::vector<double> S2 = Spline::solvePointTangentIntersection(Qk,M,BiVec2,MR)[0];
            // Get S
            std::vector<double> S = {(S1[0]+S2[0])/2,(S1[1]+S2[1])/2};
            // std::cout <<"S = " <<S[0] << "," << S[1] << "\n";
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

Spline Spline::LocalRQuadInterp(std::vector<std::vector<double>> Q,std::vector<std::vector<double>> T) { // Should be working (Tested)
    int n = Q.size()-1;
    std::vector<std::vector<double>> CPoints = {Q[0]};
    std::vector<std::vector<double>> gamma = {{0,0}};
    std::vector<double> weight = {1};
    
    for(int i =0;i<n;i++) { // Loop over consecutive points.
        // std::cout << "i = " << i << "\n";

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
        // std::cout << "===================== i = " << i << " ================================\n";
        // std::cout << "Current Qk:\n" << Qk[0] << "\n" << Qk[1] << "\n"; 
        // std::cout << "Current Qkm:\n" << Qkm[0] << "\n" << Qkm[1] << "\n"; 
        // std::cout << "Current Tk:\n" << Tk[0] << "\n" << Tk[1] << "\n"; 
        // std::cout << "Current Tkm:\n" << Tkm[0] << "\n" << Tkm[1] << "\n"; 
        if(fabs(fabs(DotProductValue)-1) <= 1e-6) { // Parallel Case
            // std::cout << "Parallel Case\n";
            std::vector<double> dQVec = {Qk[0] - Qkm[0],Qk[1] - Qkm[1]};
            double dQnorm = sqrt(pow(dQVec[0],2) + pow(dQVec[1],2));
            dQVec[0] /= dQnorm;
            dQVec[1] /= dQnorm;

            // Dot Products
            double Dk = dQVec[0]*Tk[0]+dQVec[1]*Tk[1];
            if(fabs(fabs(Dk)-1) <= 1e-6) { // Straight Line Case
                std::vector<double> Rk = {(Qkm[0]+Qk[0])/2,(Qkm[1]+Qk[1])/2};
                // std::cout << "i = " << i << "\n";
                // std::cout <<"Rk = \n" << Rk[0]<< "\n" << Rk[1] << "\n";
                CPoints.insert(CPoints.end(),Rk);
                CPoints.insert(CPoints.end(),Qk);

                double w1 = makeWeight(Qkm,Rk,Qk);
                // std::cout <<"w1 = \n" << w1<< "\n";
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
            // std::cout << "Not Parallel\n";
            std::vector<std::vector<double>> sol = solvePointTangentIntersection(Qkm,Qk,Tkm,Tk);
            std::vector<double> Rk1 = sol[0];
            std::vector<double> g = sol[1];

            if(g[0] <= 1e-12 || g[1] >= -1e12) { // Solution Exists, but not in bounds (Inflection Point or 180 turn)
                // std::cout << "Inflection or 180\n";
                // std::cout << "Qk = " << Qk[0] << "," << Qk[1] << "\n";
                // std::cout << "Qkm = " << Qkm[0] << "," << Qkm[1] << "\n";
                std::vector<double> dQVec = {Qk[0] - Qkm[0],Qk[1] - Qkm[1]};
                double dQnorm = sqrt(pow(dQVec[0],2) + pow(dQVec[1],2));
                // std::cout << "dQnorm = " << dQnorm << "\n"; 
                dQVec[0] /= dQnorm;
                dQVec[1] /= dQnorm;

                // Dot Product
                double Dk = dQVec[0]*Tk[0]+dQVec[1]*Tk[1];
                double Dkm = dQVec[0]*Tkm[0]+dQVec[1]*Tkm[1];
                // std::cout << "Dk = " << Dk << "\n";
                // std::cout << "Dkm = " << Dkm<< "\n";
                //Angles
                double thetak = fabs(acos(Dk));
                double thetakm = fabs(acos(Dkm));
                // std::cout << "thetak = " << thetak<< "\n";
                // std::cout << "thetakm = " << thetakm<< "\n";

                double a = 2.0/3.0; // Tunable Parameter

                double numer = dQnorm;
                double denom1 = 2*(1+a*cos(thetak)+(1-a)*cos(thetakm));
                double denom2 = 2*(1+a*cos(thetakm)+(1-a)*cos(thetak));
                // std::cout << "denom1 = " << (1-a)*cos(thetakm)<< "\n";
                // std::cout << "denom2 = " << denom2<< "\n";
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
                // std::cout << "Normal\n";
                double w1 = makeWeight(Qkm,Rk1,Qk);

                weight.insert(weight.end(),{w1,1});
                CPoints.insert(CPoints.end(),{Rk1,Qk});
            }
        }
    }
    // std::cout << "Loop Done\n";
    // At this point, weight and CPoints are made. The last thing we need is the knot vector
    // Split the Q and R sets
    std::vector<std::vector<double>> Qset = {CPoints[0]};
    std::vector<std::vector<double>> Rset = {CPoints[1]};
    // std::cout << "CHECK QSET and RSET Creation **********************************\n";
    for(int i = 2; i < CPoints.size()-1;i+=2){ 
        Qset.insert(Qset.end(),CPoints[i]); // Get every other control point
        Rset.insert(Rset.end(),CPoints[i+1]); // Get every other other control point.
        // std::cout << "Qset added: " << CPoints[i][0] <<","<<CPoints[i][1] <<"\n";
        // std::cout << "Rset added: " << CPoints[i+1][0] <<","<<CPoints[i+1][1] <<"\n";
    }
    n = Qset.size()-1;
    // Make ubar Helper
    std::vector<double> ubar= {0.0,1.0};
    for(int i = 2;i <= n+1; i++) {
        double coeff = ubar[i-1]-ubar[i-2];

        double numer = sqrt(pow(Rset[i-1][0]-Qset[i-1][0],2)+pow(Rset[i-1][1]-Qset[i-1][1],2));
        double denom = sqrt(pow(Qset[i-1][0]-Rset[i-2][0],2)+pow(Qset[i-1][1]-Rset[i-2][1],2));

        ubar.insert(ubar.end(),ubar[i-1]+coeff*numer/denom);
    }
    double un = ubar[ubar.size()-1];
    // Make Knot Vector
    std::vector<double> U = {0,0,0};
    for(int i = 1;i < ubar.size();i++) {
        U.insert(U.end(),{double(i+1)/ubar.size(),double(i+1)/ubar.size()});
    }
    U.insert(U.end(),1);
    
    // std::cout << " Knot Vector Inside: ";
    // for(int i = 0; i < ubar.size();i++) {
    //     std::cout << ubar[i] << ",";
    // }
    Spline ret = Spline(CPoints,U,weight);
    return ret;
}

// Fabien Area Finding
std::array<double, 2> Spline::coeffsAreaExact(const double w) {
    const auto L = 1.0 / (w * w - 1.0);
    const auto S = (w < 1.0) ? sqrt(1.0 - w * w) : sqrt(w * w - 1.0);
    const auto T = (w < 1.0) ? atan((1.0 - w) / S) / S : atanh((w - 1.0) / S) / S;
    return {L * (0.5 - w * T), L * (0.5 * w * w - w * T)};
}

std::array<double, 2> Spline::coeffsAreaSeries(const double w) {
    std::array<double, 2> K;
    K.fill(double(0.0));
    double x = double(1.0);
    int i = 0;
    while (i <= 40) {
        for (int j = 0; j < 2; ++j) {
        double add_to_coeff = Spline::ASeries[i][j] * x;
        K[j] += add_to_coeff;
        }
        x *= w - 1.0;
        i++;
    }
    return K;
}

// Dynamic Methods **************************
// Assuming 2D
int Spline::findSpan(double u) {
    int spanIndex = -1;
    // std::cout << fabs(u-1) << "\n";
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
            P = 1.0;
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
            // std::cout << i << "," << j <<"\n";
            if(i == 0 && j > 0){
                numerCoeffsX.insert(numerCoeffsX.end(),{0,0,0});
                numerCoeffsY.insert(numerCoeffsY.end(),{0,0,0});
                denomCoeffs.insert(denomCoeffs.end(),{0,0,0});
            }
            // First Span
            if(bounds[0][0] == spans[j][0] && bounds[0][1] == spans[j][1]) { 
                // std::cout << "First Span\n";
                for(int k = 0;k < numerCoeffsX[j].size();k++) {
                    denomCoeffs[j][k] += coeffs[0][k]*Weights[i];
                    numerCoeffsX[j][k] += coeffs[0][k]*Weights[i]*ControlPoints[i][0];
                    numerCoeffsY[j][k] += coeffs[0][k]*Weights[i]*ControlPoints[i][1];
                }
            }
            // Second Span
            if(bounds[1][0] == spans[j][0] && bounds[1][1] == spans[j][1]) {
                // std::cout << "Second Span\n";
                for(int k = 0;k < numerCoeffsX[j].size();k++) {
                    denomCoeffs[j][k] += coeffs[1][k]*Weights[i];
                    numerCoeffsX[j][k] += coeffs[1][k]*Weights[i]*ControlPoints[i][0];
                    numerCoeffsY[j][k] += coeffs[1][k]*Weights[i]*ControlPoints[i][1];
                }
            }
            // Third Span
            if(bounds[2][0] == spans[j][0] && bounds[2][1] == spans[j][1]) { 
                // std::cout << "Third Span\n";
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
// TEST BELOW HERE ***********************************************************************************
double Spline::getArcLength() { 
    double nudge = 1e-12;
    // Three Point Quadrature
    // std::vector<double> GaussPoints = {-sqrt(3/5),0,sqrt(3/5)};
    // std::vector<double> GaussWeights = {5/9,8/9,5/9};
    // Five Point Quadrature
    std::vector<double> GaussPoints = {-sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0,
                                       -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0,
                                       0,
                                       sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0,
                                       sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0};
    std::vector<double> GaussWeights = {(322.0-13.0*sqrt(70.0))/900.0,
                                        (322.0+13.0*sqrt(70.0))/900.0,
                                        128.0/225.0,
                                        (322.0+13.0*sqrt(70.0))/900.0,
                                        (322.0-13.0*sqrt(70.0))/900.0};
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

        // std::cout << "u1: " << u1 << "\n u2: "<< u2 << std::endl;
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
            // std::cout << "Gauss: " << GaussPoints[j] <<"\n";
        }
        
    }
    
    return AL;
}

double Spline::getSurfaceEnergy() {
    double nudge = 1e-12;
    // Three Point Quadrature
    // std::vector<double> GaussPoints = {-sqrt(3/5),0,sqrt(3/5)};
    // std::vector<double> GaussWeights = {5/9,8/9,5/9};
    // Five Point Quadrature
    std::vector<double> GaussPoints = {-sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0,
        -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0,
        0,
        sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0,
        sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0};
    std::vector<double> GaussWeights = {(322.0-13.0*sqrt(70.0))/900.0,
         (322.0+13.0*sqrt(70.0))/900.0,
         128.0/225.0,
         (322.0+13.0*sqrt(70.0))/900.0,
         (322.0-13.0*sqrt(70.0))/900.0};
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
            Ek += m*GaussWeights[j]*n*std::abs(k); // Unsigned Surface Energy is the real one
        }
    }
    return Ek;
}

double Spline::getCurvature(double u) {
    
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
    term2 = pow(pow(xp,2)+pow(yp,2),3.0/2.0);

    double k = term1/term2;
    return k;
}

std::vector<std::vector<double>> Spline::makeRationalQuadCurve(std::vector<double> uset) {
    std::vector<std::vector<double>> curve = {{0,0}};
    for(int i = 0; i < uset.size();i++) { // Loop over all u values
        std::vector<double> numer = {0,0};
        double denom = 0;
        std::vector<double> val = {0,0};
        for(int j = 0; j < ControlPoints.size();j++){ // Loop over control points
            for(int k = 0; k < ControlPoints[0].size();k++) { // Loop over contorl point indices
                numer[k] += BBasisFunction(j,2,uset[i])*Weights[j]*ControlPoints[j][k]*2;
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

double Spline::integratedSpline(double u) {
    // Find Span we are in
    int spanIndex = findSpan(u);
    // std::cout << "Found Span\n" << spanIndex << "\n";
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

    // std::cout <<"\nParameters -" << g << "," << f2 << "," << u << "," << h << "\n";
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

    double determinant = -pow(g,2) + 4*f2*h;
    // std::cout << "Determinant = " << determinant << "\n";
    double N3B;
    double D3B;
    if(determinant >= 1e-12) { // Strictly Positive 
        // std::cout << "positive\n";
        N3B = (4*(6*e2*pow(f2,2) - 3*d2*f2*g + c2*pow(g,2) + 2*c2*f2*h - 3*b2*g*h +
                    6*a2*pow(h,2))*atan((g + 2*f2*u)/sqrt(determinant)));

        D3B = pow(determinant,2.5);
    } else if(determinant <= -1e-12) { // Strictly Negative
        // std::cout << "negative\n";
        N3B = (4*(6*e2*pow(f2,2) - 3*d2*f2*g + c2*pow(g,2) + 2*c2*f2*h - 3*b2*g*h +
                    6*a2*pow(h,2))*-1*atanh((g + 2*f2*u)/sqrt(-determinant)));

        D3B = pow(-determinant,2.5);
    } else { // Near Zero
        N3B = 0;
        D3B = 1;
    }

    if(fabs(D2B) < 1e-6) { // Near 0
        // D2B = 1;
        // N2B = 0;
        std::cout << "D2 Near 0\n";
    }

    if(fabs(D1B) < 1e-6) { // Near
        // D1B = 1;
        // N1B = 0;
        std::cout << "D1 Near 0\n";
    }
    // std::cout << "End \n";

    double value = 0.5*(N1B/D1B + N2B/D2B + N3B/D3B);

    return value;
}

double Spline::getArea() {
    // std::cout << "GETTING AREA\n";
    double Area = 0;
    double nudge = 1e-8;
    for (int i = 0; i < ControlPoints.size()-2; i += 2) {
        const auto x0 = ControlPoints[i][0], y0 = ControlPoints[i][1];
        const auto x1 = ControlPoints[i + 1][0], y1 = ControlPoints[i+1][1];
        const auto x2 = ControlPoints[(i + 2) % ControlPoints.size()][0],
                   y2 = ControlPoints[(i + 2) % ControlPoints.size()][1];
        const auto w = Weights[i + 1];
        if (w < 0.35 || w > 1.7) {
          auto K = Spline::coeffsAreaExact(w);
          Area -= (x0 * y2 - x2 * y0) * K[0] +
                (x1 * y0 + x2 * y1 - x0 * y1 - x1 * y2) * K[1];
        } else {
          auto K = Spline::coeffsAreaSeries(w);
          Area -= (x0 * y2 - x2 * y0) * K[0] +
                (x1 * y0 + x2 * y1 - x0 * y1 - x1 * y2) * K[1];
        }
    }
    // double Acheck = Area;
    // Area = 0;
    // nudge = 1e-8;
    // // Loop over each breakpoint, find area of that section, then add them together.
    // double uL;
    // double uR;
    // // std::cout << breakpoints.size() << "\n";
    // for(int i =0; i <breakpoints.size()-1; i++) {
    //     // std:: cout << "====================== i = " << i << " ========================\n";
    //     uL = breakpoints[i]; // Left Breakpoint
    //     uR = breakpoints[i+1]; // Right Breakpoint

    //     uL += nudge;
    //     uR -= nudge;

    //     // std :: cout << "uL = " << uL << "\n";
    //     // std :: cout << "uLVal = " << integratedSpline(uL) << "\n";
    //     // std :: cout << "uR = " << uR << "\n";
    //     // std :: cout << "uRVal = " << integratedSpline(uR) << "\n";

    //     Area += integratedSpline(uR) - integratedSpline(uL);
    //     // std::cout <<"Added to Area\n";
    //     // std:: cout << "Area = " << Area << "\n";
    // }
    // std:: cout << "Area = " << Area << "\n";
    // std:: cout << "ACheck = " << Acheck << "\n";
    // std::cout << "Difference = " << (Area-Acheck) <<"\n";
    return Area;
}

std::vector<double> Spline::lineCurveIntersection(std::vector<double> P1, std::vector<double> P2) { // Should be working (Tested)
    std::vector<double> tangent = {P2[0]-P1[0],P2[1]-P1[1]};
    std::vector<double> normal = {-tangent[1],tangent[0]};
    int numSpans = breakpoints.size()-1;
    std::vector<double> uIntersections = {0};
    std::vector<double> sIntersections = {0};
    double a = normal[0];
    double b = normal[1];
    double c = -(normal[0]*P1[0]+normal[1]*P1[1]);
    int count = 0;
    int numInt = 0;
    for(int i = 0; i < numSpans; i ++) { // Loop Over Spans
        // Quadratic Coefficients
        double aPoly = a*numerCoeffsX[i][0] + b*numerCoeffsY[i][0]+c*denomCoeffs[i][0];
        double bPoly = a*numerCoeffsX[i][1] + b*numerCoeffsY[i][1]+c*denomCoeffs[i][1];
        double cPoly = a*numerCoeffsX[i][2] + b*numerCoeffsY[i][2]+c*denomCoeffs[i][2];
        // Solve Quadratic
        double discrim = pow(bPoly,2)-4*aPoly*cPoly;
        
        // std::cout << "Discrim = " << discrim << "\n";
        // std::cout << "aPoly = " << aPoly << "\n";
        // std::cout << "bPoly = " << bPoly << "\n";
        // std::cout << "cPoly = " << cPoly << "\n";
        // Filter Values for real, in span on curve,and in spans of line
        double u1;
        double u2;
        if(discrim >= 0) { // Real Values Only
            std::vector<double> Inters;
            if(aPoly == 0) {
                u1 = -cPoly/bPoly;
                Inters = {u1};
            } else {
                u1 = (-bPoly + sqrt(discrim))/(2*aPoly);
                u2 = (-bPoly - sqrt(discrim))/(2*aPoly);
                Inters = {u1,u2};
            }

            std::vector<std::vector<double>> pointInter = this->makeRationalQuadCurve(Inters);
            // Check u1
            for(int j = 0; j < Inters.size();j++) {
                u1 = Inters[j];
                // std::cout << "span2= " << spans[i][1] << "\n"; 
                if(u1 >= (spans[i][0]-1e-12) && (u1 < (spans[i][1]-1e-12))) { // Within Span of curve, not including the endpoint of span.
                    std::vector<double> Pcurr = pointInter[j];
                    Pcurr = this->makeRationalQuadCurve({u1})[0];
                    std::vector<double> dP1 = {Pcurr[0]-P1[0],Pcurr[1]-P1[1]};
                    std::vector<double> dP2 = {Pcurr[0]-P2[0],Pcurr[1]-P2[1]};
                    std::vector<double> dP = {P1[0]-P2[0],P1[1]-P2[1]};
                    double lineCheck;
                    if(fabs(dP[0]) > fabs(dP[1])) { // Checks if point between endpoints of line
                        lineCheck = fabs(dP1[0])+fabs(dP2[0])-fabs(dP[0]);
                    } else {
                        lineCheck = fabs(dP1[1])+fabs(dP2[1])-fabs(dP[1]);
                    }
                    // If point on line, add
                    if(lineCheck <= 1e-8) {
                        if(count == 0){ // If first one, replace first value
                            uIntersections[0] = u1;
                            sIntersections[0] = std::min(fabs(dP1[0]/dP[0]),fabs(dP1[1]/dP[1]));
                            count++;
                        } else { // else, add to end
                            uIntersections.insert(uIntersections.end(),u1);
                            sIntersections.insert(sIntersections.end(),std::min(fabs(dP1[0]/dP[0]),fabs(dP1[1]/dP[1])));
                        }
                        numInt++;
                    }
                }
            }
        }
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

    // Finally, we will loop over u and only take unique values, since a line can only intersect a curve once at the same location
    std::vector<double> ret2 = {ret[0]};
    for(int i = 1; i < ret.size();i++) {
        int present = 0;
        for(int j = 0;j < ret2.size(); j++) {
            if(ret2[j] == ret[i]) {
                present = 1;
            }
        }

        if(present == 0) {
            ret2.insert(ret2.end(),ret[i]);
        }
    }
    if(numInt > 0){
        return ret2;
    } else {
        return {};
    }
    
}

std::vector<std::vector<double>> Spline::getParameterLoop(std::vector<std::vector<double>> square) { // Should be working (Tested)
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

        parameter.insert(parameter.end(),(i+1)%4);
        // std::cout << "\n Current INdicator \n";
        // for(int i = 0; i < indicator.size() ;i++ ){
        //     std :: cout << indicator[i] << ",";
        // }
        // std::cout << "\n Current Indicator \n";
        int lastInd = indicator[indicator.size()-1];
        // std::cout << "Last Ind = " << lastInd << "\n";
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
    if(indicator[0] == 1.0 && count != 0) {
        tempParameters.insert(tempParameters.end(),tempParameters[0]);
        tempIndicators.insert(tempIndicators.end(),tempIndicators[0]);
    }

    return {tempParameters,tempIndicators};
}

std::vector<double> Spline::getTangent(double u) {
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

    return {xp,yp};
}

std::vector<std::vector<double>> Spline::clippedBezier(double u1, double u2) {
    std::vector<std::vector<double>> Points;
    std::vector<std::vector<double>> tans;
    std::vector<double> curves; // Stores Curvatures

    Points = {this->makeRationalQuadCurve(
        {u1})[0]};              // Add Starting Point to Points (Exit Point)
    tans = {this->getTangent(u1)};  // Add starting point Tangent (Exit Point)
    curves = {this->getCurvature(u1)};
    // Find Spans, in order
    int spanIndex1 = this->findSpan(u1);
    int spanIndex2 = this->findSpan(u2);
    std::cout << "spanIndex1 = " << spanIndex1 << "\n";
    std::cout << "spanIndex2 = " << spanIndex2 << "\n";
    // if u1,u2 in the same span, we can go direct.
    // If u1,u2 are in different spans, we have to segment through each span
    if (spanIndex1 == spanIndex2) {  // Same Span 
        std::vector<double> Pend =
        this->makeRationalQuadCurve({u2})[0];
        std::vector<double> Tend = this->getTangent(u2);

        std::vector<double> Pstart = this->makeRationalQuadCurve({u1})[0];
        std::vector<double> Tstart = this->getTangent(u1);

        // Calculate Intersection
        std::vector<std::vector<double>> solution =
            Spline::solvePointTangentIntersection(Pstart, Pend, Tstart, Tend);
        std::vector<double> inter = solution[0];

        // Now, add intersection and end point to array
        Points.insert(Points.end(), inter);
        Points.insert(Points.end(), Pend);

        tans.insert(tans.end(),Tend);  // Put end tangent to back of this to move forward.
        curves.insert(curves.end(),0);
        curves.insert(curves.end(),this->getCurvature(u2));
    } else {  // Different Spans

      int numSpans = spans.size();
      int tempSpanIndex2 = spanIndex2;
      if (spanIndex2 < spanIndex1) {
        tempSpanIndex2 = spanIndex2 + numSpans;
      }

      // Get Set of Breakpoint Indices
      std::vector<int> breakIndexSet = {spanIndex1 + 1};
      for (int j = spanIndex1 + 2; j <= tempSpanIndex2; j++) {
        breakIndexSet.insert(breakIndexSet.end(), j);
      }

      // Mod them into range
      for (int j = 0; j < breakIndexSet.size(); j++) {
        breakIndexSet[j] = breakIndexSet[j] % numSpans;
      }
      // Calculate Break Values
      std::vector<double> breaks = {u1};
      for (int j = 0; j < breakIndexSet.size(); j++) {
        breaks.insert(breaks.end(), breakpoints[breakIndexSet[j]]);
      }

      breaks.insert(breaks.end(),
                    u2);  // This contains the parameter value of the start,
                          // all the breakpoints between
      // segments, then the end.
      std::cout << "breaks size = " << breaks.size()  << "\n";
      for (int j = 0; j < breaks.size() - 1; j++) {
        // For each of these points, add the intersection and end point, along
        // with the previous tangent
        std::vector<double> Pend =
        this->makeRationalQuadCurve({breaks[j + 1]})[0];
        std::vector<double> Tend = this->getTangent(breaks[j + 1]);

        std::vector<double> Pstart = this->makeRationalQuadCurve({breaks[j]})[0];
        std::vector<double> Tstart = this->getTangent(breaks[j]);
        std::cout << "j = " << j << "==================\n";
        std::cout << "Pstart = " << Pstart[0] << "," << Pstart[1] << "\n";
        std::cout << "Pend = " << Pend[0] << "," << Pend[1] << "\n";
        std::cout << "Tstart = " << Tstart[0] << "," << Tstart[1] << "\n";
        std::cout << "Tend = " << Tend[0] << "," << Tend[1] << "\n";

        // Calculate Intersection
        std::cout << "Find Intersection\n";
        std::vector<std::vector<double>> solution =
            Spline::solvePointTangentIntersection(Pstart, Pend, Tstart, Tend);
        std::cout << "End Intersection\n";
        std::vector<double> inter = solution[0];
        std::cout << "End inter\n";
        // Now, add intersection and end point to array
        Points.insert(Points.end(), inter);
        Points.insert(Points.end(), Pend);
        curves.insert(curves.end(),0); // Intersection Curvature (unknown, unused)
        curves.insert(curves.end(),this->getCurvature(breaks[j + 1])); // End Point Curvature (Will be used)
        tans.insert(tans.end(),Tend);  // Put end tangent to back of this to move forward.

        }
        std::cout << "End loop\n";
    }
    // From here, we will use this information to calculate the  weights.
    std::cout <<"======== POINTS RESULT ========\n";
    for(int i= 0; i < Points.size();i++) {
        std::cout << "(" << Points[i][0] << "," << Points[i][1] << ")\n";
    }
    std::cout <<"Start Weights\n";
    std::vector<double> weights = {1};  // First weight is always one (every other one will be too);
    for (int i = 0; i < Points.size() - 1;i += 2) {  
        // We go every 2 since the points go in triples (0,1,2 then
        // 2,3,4 then 4,5,6 etc.
        // Here is an example of how to use the Points array, and what the ordering
        // of it is
        std::cout <<"P0 Start\n";
        std::vector<double> P0 = Points[i];      // Start
        std::cout <<"P1 Start\n";
        std::vector<double> P1 = Points[i + 1];  // Intersection
        std::cout <<"P2 Start\n";
        std::vector<double> P2 = Points[i + 2];  // End
        std::cout <<"P2 End\n";
        // Since I don't fully understand the curvature method in paper yet, I am
        // just quickly doing this to get weights Since I know they are well
        // behaved. I will understand the method in the paper soon and then that can
        // be used In general case.

        // First, calculate midpoint of P0,P2
        // std::cout << "Get Mid \n";
        // std::vector<double> mid = {(P0[0] + P2[0]) / 2, (P0[1] + P2[1]) / 2};
        // std::cout << "Got Mid \n";
        // // Next, calculate intersection with spline (Shoulder Point)
        // std::cout << "uHit Getting\n";
        // double uHit = this->lineCurveIntersection(mid, P1)[0];
        // std::cout << "uHit = " << uHit << "\n";
        // std::vector<double> S = this->makeRationalQuadCurve({uHit})[0];
        // // Calculate Midpoint to Shoulder Point Distance
        // double MS = sqrt(pow(S[0] - mid[0], 2) + pow(S[1] - mid[1], 2));
        // // Calculate Shoulder Point to intersection distance
        // double SP1 = sqrt(pow(S[0] - P1[0], 2) + pow(S[1] - P1[1], 2));
        // // Calculate Weight (eq 7.32)
        // double w = MS / SP1;
        
        // std::cout << "weight = " << w << "\n";

        // ==========Calculate Weight Using Curvature
        // Calculate Triangle Area
        double a = sqrt(pow(P0[0] - P1[0], 2) + pow(P0[1] - P1[1], 2)); // P0 -> P1
        double b = sqrt(pow(P2[0] - P1[0], 2) + pow(P2[1] - P1[1], 2)); // P1 -> P2
        double c = sqrt(pow(P0[0] - P2[0], 2) + pow(P0[1] - P2[1], 2)); // P2 -> P0
        double s = (a+b+c)/2;
        double Atri = sqrt(s*(s-a)*(s-b)*(s-c));
        // l is the same as a in my case.
        double w2;
        double k0 = curves[i];
        std::cout << "k0 = " << k0 << "\n";
        std::cout << "Atri = " << Atri << "\n";
        // Catch Straight Line Case

        // ================================ WEIGHT PROBLEMS HERE ===============================
        if(fabs(k0) >= 1e-6){
            w2 = sqrt(Atri/(fabs(k0)*a*a*a)); // can go to NaN if k0 or a is 0.
        } else {
            w2 = 0;
        }
        std::cout << "weight = " << w2 << "\n";
        // Add Weight
        weights.insert(weights.end(), w2);
        weights.insert(weights.end(), 1);  // Always 1 at the start and end
        //
        // std::cout << "weight should = " << s.makeWeight(P0,P1,P2) << "\n";
    }
    std::cout <<"end Weights\n";
    
    // Get into form of [Px,Py,w] array
    std::vector<std::vector<double>> ret = {{Points[0][0],Points[0][1],weights[0]}};
    for(int i = 1; i < weights.size();i++) {
        ret.insert(ret.end(),{Points[i][0],Points[i][1],weights[i]});
    }
    std::cout << "============ CLIPPED SECTION INSIDE========== \n";
    // for(int i = 0; i < ret.size(); i++) {
    //     std::cout << "(" << ret[i][0] <<"," << ret[i][1] << ")" << "," << ret[i][2] << "\n";
    // }
    return ret;
}

double Spline::integrateSplineSquare(std::vector<std::vector<double>> square) { // Should be working (Tested)
    std::cout << "Begin Parameter Getting \n";
    std::vector<std::vector<double>> loop = this->getParameterLoop(square);
    std::cout << "End Parameter Getting \n";
    std::vector<double> parameter = loop[0];
    std::vector<double> indicator = loop[1];
    // std::cout << "Begin Parameter \n";
    // for(int i = 0; i < parameter.size(); i++) {
    //     std::cout << parameter[i] << ",";
    // }
    // std::cout << "\nEnd Parameter \n";

    // std::cout << "Begin Indicator \n";
    // for(int i = 0; i < indicator.size(); i++) {
    //     std::cout << indicator[i] << ",";
    // }
    // std::cout << "\nEnd Indicator \n";
    double Area = 0;
    // std::cout << "Indictator/Parameter size = " << indicator.size() << "\n";
    for(int i = 0;i<indicator.size()-1;i++) { // Loop over pairs of indicators/parameters
        // std::cout << "i = " << i << "\n";
        double ind1 = indicator[i];
        double ind2 = indicator[i+1];
        // Most indicator pairs are either 0 or impossible. As such, we only look at relevant cases.
        // std::cout << "Start Branches - ind1 = " << ind1 << " ind2 = " << ind2 << "\n";
        if(ind1 == 2 && ind2 == 2) { // 2 To 2 - Two Inside Corners - Integrate Square *********************************************
            std::cout << "Case 22\n";
            std::vector<double> P1 = square[parameter[i]];
            std::vector<double> P2 = square[parameter[i+1]];

            double integral = 0.5*(P2[0]*P1[1]-P1[0]*P2[1]);
            Area -= integral;
        } else if(ind1 == 2 && ind2 == 4) {// 2 To 4 - Inside Corner to Exit - Integrate Square ************************************
            std::cout << "Case 24\n";
            std::vector<double> P1 = square[parameter[i]];
            std::vector<double> P2 = makeRationalQuadCurve({parameter[i+1]})[0];

            double integral = 0.5*(P2[0]*P1[1]-P1[0]*P2[1]);
            Area -= integral;
        } else if(ind1 == 3 && ind2 == 2) {// 3 to 2 - Entry to Inside Corner - Integrate Square ***********************************
            std::cout << "Case 32\n";
            std::vector<double> P1 = makeRationalQuadCurve({parameter[i]})[0];
            std::vector<double> P2 = square[parameter[i+1]];

            double integral = 0.5*(P2[0]*P1[1]-P1[0]*P2[1]);
            Area -= integral;
        } else if(ind1 == 3 && ind2 == 4) {// 3 to 4 - Entry to Exit - Integrate Square ********************************************
            std::cout << "Case 34\n";
            std::vector<double> P1 = makeRationalQuadCurve({parameter[i]})[0];
            std::vector<double> P2 = makeRationalQuadCurve({parameter[i+1]})[0];

            double integral = 0.5*(P2[0]*P1[1]-P1[0]*P2[1]);
            Area -= integral;
        } else if(ind1 == 4 && ind2 == 3) {// 4 to 3 - Exit to Entry - Integrate Curve
            std::cout << "Case 43\n";
            // Get Parameter Values
            double u1 = parameter[i];
            double u2 = parameter[i+1];
            std::cout << "Start Clipping\n";
            std::vector<std::vector<double>> clippedSection = this->clippedBezier(u1,u2);
            double integral = 0.0;
            std::cout << "============ CLIPPED SECTION ========== \n";
            for(int i = 0; i < clippedSection.size(); i++) {
                std::cout << "(" << clippedSection[i][0] <<"," << clippedSection[i][1] << ")" << "," << clippedSection[i][2] << "\n";
            }
            // std::cout << "\n ========= INTEGRAL CALCULATION ========\n";
            for(int i = 0; i < clippedSection.size()-2; i += 2) {
                // std::cout << "Get Points \n";
                // std::cout << "Get Point 1\n";
                const auto x0 = clippedSection[i][0], y0 = clippedSection[i][1];
                // std::cout << "Get Point 2\n";
                const auto x1 = clippedSection[i+1][0], y1 = clippedSection[i+1][1];
                const auto w = clippedSection[i + 1][2];
                // std::cout << "Get Point 3\n";
                const auto x2 = clippedSection[i + 2][0],
                        y2 = clippedSection[i + 2][1];
                

                std::cout << "P0 = " << x0 << "," << y0 << ", w = " << clippedSection[i][2] << "\n";
                std::cout << "P1 = " << x1 << "," << y1 << ", w = " << clippedSection[i+1][2]<< "\n";
                std::cout << "P2 = " << x2 << "," << y2 << ", w = " << clippedSection[i+2][2]<< "\n";
                double dA;
                // Uses weights to calculate area here. ==========================================================================
                if (w < 0.35 || w > 1.7) {
                    auto K = Spline::coeffsAreaExact(w);
                    dA = (x0 * y2 - x2 * y0) * K[0] +
                        (x1 * y0 + x2 * y1 - x0 * y1 - x1 * y2) * K[1];
                } else {
                    auto K = Spline::coeffsAreaSeries(w);
                    dA = (x0 * y2 - x2 * y0) * K[0] +
                        (x1 * y0 + x2 * y1 - x0 * y1 - x1 * y2) * K[1];
                }
                integral -= dA;
                // std::cout << "dA = " << -dA << "\n";
            }
            // std::cout << "Case 5 \n";
            // std::cout << "Component Integral = " << integral << "\n";
            // std::cout << "\n4-3 Area Contribution: " << integral << "\n";
            Area += integral;
        }   
    }
    return Area;
}

void Spline::saveToVTK(const std::string& filename, const int nsamples){
    std::vector<double> uset(nsamples, 0.);
    for (int i = 0; i < nsamples; i++){
        uset[i] = KnotVector[0] + (KnotVector[KnotVector.size()-1] - KnotVector[0]) * static_cast<double>(i) / static_cast<double>(nsamples - 1);
    }
    const auto curve = this->makeRationalQuadCurve(uset);

    std::ofstream file;
    file.open(filename + std::string(".vtu"));
    file << "<VTKFile type=\"UnstructuredGrid\">\n<UnstructuredGrid>\n";
    file << "<Piece NumberOfPoints=\"" << nsamples << "\" NumberOfCells=\"" << 1 << "\">\n";
    file << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";
    for (int i = 0; i < nsamples; i++){
        file << std::scientific << std::setprecision(15) << curve[i][0] << " " << curve[i][1] << " 0. \n";
    }
    file << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    int count = 0;
    for (int i = 0; i < nsamples; i++){
        file << count++ << " ";
    }
    file << "\n</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    file << nsamples << " ";
    file << "\n</DataArray>\n";
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    file << "7 ";
    file << "\n</DataArray>\n";
    file << "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    file.close();
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

// Bebugging
void Spline::printControlPoints() {
    std::cout << "\nPrinting Control Points \n";
    for(int i = 0; i < ControlPoints.size(); i++ ){
        for(int j = 0; j <ControlPoints[0].size();j++) {
            std::cout << ControlPoints[i][j] << ",";
        }
        std::cout << "\n";
    }
    std::cout << "End Control Points \n";
}
void Spline::printKnotVector() {
    std::cout << "\nPrinting Knot Vector\n";
    for(int i = 0; i < KnotVector.size(); i++ ){
        std::cout << KnotVector[i] << ",";
    }
    std::cout << "\nEnd  Knot Vector \n";
}
void Spline::printWeights(){
    std::cout << "\nPrinting Weights\n";
    for(int i = 0; i < Weights.size(); i++ ){
        std::cout << Weights[i] << ",";
    }
    std::cout << "\nEnd Weights \n";
}
void Spline::printXCoeffs() {
    std::cout << "\nPrinting X Coeffs \n";
    for(int i = 0; i < numerCoeffsX.size(); i++ ){
        for(int j = 0; j <numerCoeffsX[0].size();j++) {
            std::cout << numerCoeffsX[i][j] << ",";
        }
        std::cout << "\n";
    }
    std::cout << "End X Coeffs \n";
}
void Spline::printYCoeffs() {
    std::cout << "\nPrinting Y Coeffs \n";
    for(int i = 0; i < numerCoeffsY.size(); i++ ){
        for(int j = 0; j <numerCoeffsY[0].size();j++) {
            std::cout << numerCoeffsY[i][j] << ",";
        }
        std::cout << "\n";
    }
    std::cout << "End Y Coeffs\n";
}
void Spline::printDCoeffs() {
    std::cout << "\nPrinting Denom Coeffs \n";
    for(int i = 0; i < denomCoeffs.size(); i++ ){
        for(int j = 0; j <denomCoeffs[0].size();j++) {
            std::cout << denomCoeffs[i][j] << ",";
        }
        std::cout << "\n";
    }
    std::cout << "End Denom Coeffs\n";
}

