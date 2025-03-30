// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com> ?????????????????
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/splines/Spline.h"
#include <iostream>
#include <vector>
#include "gtest/gtest.h" // I don't know why this is happening


// =================================== TESTS FOR GETTERS =======================================
// std::vector<std::vector<double>> Spline::getControlPoints()
TEST(GetterTests,ControlPoints) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};

    Spline spline = Spline(CP,KV,W);
    splineCP = spline.getControlPoints();
    ASSERT_EQ(CP.size(),splineCP.size()) << "Different Number of Control Points";
    ASSERT_EQ(CP[0].size(),splineCP[0].size())<< "Different Size of Control Points";

    for(int i = 0; i < CP.size(); i++) {
        for(int j = 0; j <CP[0].size(); j++) {
            EXPECT_EQ(splineCP[i][j],CP[i][j]) << "Control Point " << i << " differ in index " << j;
        }
    }


}
// std::vector<double> Spline::getKnotVector()
TEST(GetterTests,KnotVector) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};

    Spline spline = Spline(CP,KV,W);
    splineKV = spline.getKnotVector();
    ASSERT_EQ(KV.size(),splineKV.size()) << "Different Number of Knots";

    for(int i = 0; i < KV.size(); i++) {
        EXPECT_EQ(splineKV[i],KV[i]) << "Knot Vectors Differ at index " << i;
    }
}
// std::vector<double> Spline::getWeights() 
TEST(GetterTests,Weights) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};

    Spline spline = Spline(CP,KV,W);
    splineW = spline.getKnotVector();
    ASSERT_EQ(W.size(),splineW.size()) << "Different Number of Weights";

    for(int i = 0; i < W.size(); i++) {
        EXPECT_EQ(splineW[i],W[i]) << "Weights Differ at index " << i;
    }
}
// std::vector<double> Spline::getBreakpoints()
TEST(GetterTests,Breakpoints) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};
    std::vector<double> BP = {0,0.25,0.5,1};

    Spline spline = Spline(CP,KV,W);
    BPs = spline.getBreakpoints();

    ASSERT_EQ(BP.size(),BPS.size()) << "Different Breakpoints Size";

    for(int i = 0; i < BPs.size(); i++) {
        EXPECT_EQ(BPs[i],BP[i]) << "Breakpoints Differ at Index" << i;
    }
}
// std::vector<std::vector<double>> Spline::getSpans() 
TEST(GetterTests,Spans) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};
    std::vector<std::vector<double>> S = {{0,0.25},{0.25,0.5},{0.5,1}};

    Spline spline = Spline(CP,KV,W);
    Spans = spline.getSpans();

    ASSERT_EQ(S.size(),Spans.size()) << "Different Number of Control Points";
    ASSERT_EQ(S[0].size(),Spans[0].size())<< "Different Size of Control Points";

    for(int i = 0; i < S.size(); i++) {
        for(int j = 0; j <S[0].size(); j++) {
            EXPECT_EQ(Spans[i][j],S[i][j]) << "Spans " << i << " differ in index " << j;
        }
    }
}
// std::vector<std::vector<double>> Spline::getXCoeffs()
TEST(GetterTests,XCoeffs) {
    ASSERT_EQ(1,1)<< "Different Size of Control Points";

}
// std::vector<std::vector<double>> Spline::getYCoeffs()
TEST(GetterTests,YCoeffs) {
    ASSERT_EQ(1,1)<< "Different Size of Control Points";
}
// std::vector<std::vector<double>> Spline::getDCoeffs()
TEST(GetterTests,DCoeffs) {
    ASSERT_EQ(1,1)<< "Different Size of Control Points";
}
// =================================== TESTS FOR KNOT VECTOR OPERATIONS =======================================
// int Spline::findSpan(double u) 
TEST(KnotVectorOperations,SpanFinding) {

}
// std::vector<std::vector<double>> Spline::BasisCoefficients(int i)
TEST(KnotVectorOperations,BasisCoeffs) {
    
}
// std::vector<std::vector<double>> Spline::BasisCoefficientBounds(int i)
TEST(KnotVectorOperations,BasisBounds) {
    
}
// std::vector<double> Spline::makeBreakpoints() REDUNDANT TO GETTER
// TEST(KnotVectorOperations,BreakpointMaking) {
    
// }
// // std::vector<std::vector<double>> Spline::makeSpans() REDUNDANT TO GETTER
// TEST(KnotVectorOperations,SpanMaking) {

// }


// =================================== TESTS FOR SPLINE CONSTRUCTION =======================================
// std::vector<std::vector<double>> Spline::solvePointTangentIntersection(std::vector<double> Q1, 
//                                                           std::vector<double> Q2,
//                                                           std::vector<double> T1,
//                                                           std::vector<double> T2)

TEST(InterpolationMethods,PTInter) {
    // Case 1: Test Two Regular Sloped Lines
    std::vector<double> Q1 = {1.0,2.0};
    std::vector<double> Q2 = {5.0,1.0};
    std::vector<double> T1 = {1.0,1.0};
    std::vector<double> T2 = {-1.0,1.0};
    std::vector<double> ExactIntersection = {2.5,3.5}

    std::vector<std::vector<double>> Rg1 = Spline::solvePointTangentIntersection(Q1,Q2,T1,T2);
    std::vector<double> R1 = Rg1[0];

    for(int i = 0; i < ExactIntersection.size();i++) {
        EXPECT_EQ(ExactIntersection[i],R1[i]) << "Case 1 Intersection differs at index " << i;
    }

    // Case 2: One Vertical Line
    T1 = {0.0,1.0};
    ExactIntersection = {1.0,5.0};
    Rg1 = Spline::solvePointTangentIntersection(Q1,Q2,T1,T2);
    R1 = Rg1[0];

    for(int i = 0; i < ExactIntersection.size();i++) {
        EXPECT_EQ(ExactIntersection[i],R1[i]) << "Case 2 Intersection differs at index " << i;
    }

    // Case 3: One Horizontal Line
    T1 = {1.0,0};
    ExactIntersection = {5.0,1.0};
    Rg1 = Spline::solvePointTangentIntersection(Q1,Q2,T1,T2);
    R1 = Rg1[0];

    for(int i = 0; i < ExactIntersection.size();i++) {
        EXPECT_EQ(ExactIntersection[i],R1[i]) << "Case 2 Intersection differs at index " << i;
    }

    // Case 4: No Intersections
    // Not in Method yet

    // Case 5: Infinite Intersections
    // Not in Method yet
}

// double Spline::makeWeight(std::vector<double> Qkm, std::vector<double> Rk, std::vector<double> Qk)
TEST(InterpolationMethods,WeightCalculator) {
    // Collinear Case
    std::vector<double> Qkm = {1,1};
    std::vector<double> Rk = {4,1};
    std::vector<double> Qk = {7,1};
    double wExact = 0;
    double w = Spline::makeWeight(Qkm,Rk,Qk); 

    EXPECT_EQ(w,wExact) << "Collinear Weight Fail";

    // Isosceles Case
    Rk = {4,5};
    wExact = 0.6;
    w = Spline::makeWeight(Qkm,Rk,Qk); 
    EXPECT_EQ(w,wExact) << "Isosceles Weight Fail";

    // General Case
    Rk = {3,5};
    wExact = 0.597492359431871;
    w = Spline::makeWeight(Qkm,Rk,Qk);
    EXPECT_EQ(w,wExact) << "General Weight Fail";
}


// Spline Spline::LocalRQuadInterp(std::vector<std::vector<double>> Q,std::vector<std::vector<double>> T)
TEST(InterpolationMethods,Interpolator) {
    // Compare Control Points, Knot Vector, and Weights for Circle and Blob Cases
}

// =================================== TESTS FOR SPLINE PROPERTIES =======================================
// std::vector<std::vector<double>> Spline::CurveCoefficients()
TEST(SplineProperties,CurveCoefficients) {
    
}

// double Spline::getArcLength()
TEST(SplineProperties,ArcLength) {
    
}

// double Spline::getCurvature(double u)
TEST(SplineProperties,Curvature) {
    // Figure out Later
}

// double Spline::getSurfaceEnergy()
TEST(SplineProperties,SurfaceEnergy) {
    
}
// double Spline::getArea()
TEST(SplineProperties,Area) {

}


// =================================== TESTS FOR CLIPPING =======================================
// std::vector<double> Spline::lineCurveIntersection(std::vector<double> P1, std::vector<double> P2)
TEST(ClippingOperations,LineClip) {
    
}
// std::vector<std::vector<double>> Spline::getParameterLoop(std::vector<std::vector<double>> square)
TEST(ClippingOperations,ParamLoop) {
    
}
// double Spline::integrateSplineSquare(std::vector<std::vector<double>> square)
TEST(ClippingOperations,ClippedArea) {

}
// double Spline::integratedSpline(double u)
TEST(ClippingOperations,analyticIntegration) {
    
}

// =================================== TESTS FOR VISUALIZATION =======================================
// std::vector<std::vector<double>> Spline::makeRationalQuadCurve(std::vector<double> uset)
TEST(Visualization,curveMaker) {

}
