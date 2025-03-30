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

// Remaining Tests:
// KnotVectorOperations,SpanFinding
// KnotVectorOperations,BasisCoeffs
// KnotVectorOperations,BasisBounds
// InterpolationMethods,Interpolator
// SplineProperties,Curvature
// ClippingOperations,ParamLoop
// ClippingOperations,analyticIntegration
// Visualization,curveMaker (this one is making me sad)

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
    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    XCoeffsExact = {{-2.848208747679525,-0.685084207398420,1.0},
                    {-4.719066767034173,-3.742257271215488,1.771759842297852},
                    {4.719066767034200,-10.034346293927747,2.820441346083229},
                    {11.392834990718072,-14.028873960039150,3.277511770787770},
                    {11.392834990718107,-11.288537130445498,1.755102421013508},
                    {4.719066767034199,-0.452468743926028,-2.502823959473282},
                    {-4.719066767034192,14.229072309069238,-8.212312146749234},
                    {-11.392834990718100,24.155838396233065,-11.763003405514930}};
    

    std::vector<std::vector<double>> XCoeffs = circle.getXCoeffs();

    ASSERT_EQ(XCoeffsExact.size(),XCoeffs.size()) << "Different Number of XCoeffs";
    ASSERT_EQ(XCoeffsExact[0].size(),XCoeffs[0].size())<< "Different Size of XCoeffs";
    // Might need to change this for approximates
    for(int i = 0; i < XCoeffsExact.size(); i++) {
        for(int j = 0; j <XCoeffsExact[0].size(); j++) {
            EXPECT_EQ(XCoeffsExact[i][j],XCoeffs[i][j]) << "XCoeffs Differ in Span " << i << " Index " << j;
        }
    }
}
// std::vector<std::vector<double>> Spline::getYCoeffs()
TEST(GetterTests,YCoeffs) {
    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    YCoeffsExact = {{-1.179766691758548,3.444150891285807,0.0},
                    {-11.392834990718072,8.965391741942241,-0.722593359456511},
                    {-11.392834990718150,6.225054912348590,0.190852250408038},
                    {-4.719066767034176,-1.644894263644746,2.370332235060101},
                    {4.719066767034207,-12.131709301498518,5.283336412241699},
                    {11.392834990718100,-19.092356178136100,6.957648452807256},
                    {11.392834990718121,-16.352019348542427,4.826275363123294},
                    {4.719066767034185,-2.549831751496768,-2.169235015537424}};

    std::vector<std::vector<double>> YCoeffs = circle.getYCoeffs();

    ASSERT_EQ(YCoeffsExact.size(),YCoeffs.size()) << "Different Number of YCoeffs";
    ASSERT_EQ(YCoeffsExact[0].size(),YCoeffs[0].size())<< "Different Size of YCoeffs";
    // Might need to change this for approximates
    for(int i = 0; i < YCoeffsExact.size(); i++) {
        for(int j = 0; j <YCoeffsExact[0].size(); j++) {
            EXPECT_EQ(YCoeffsExact[i][j],YCoeffs[i][j]) << "YCoeffs Differ in Span " << i << " Index " << j;
        }
    }
}
// std::vector<std::vector<double>> Spline::getDCoeffs()
TEST(GetterTests,DCoeffs) {
    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    DCoeffsExact = {{3.082878933292889,-0.685084207398420,1.0},
                    {12.331515733171585,-6.850842073984197,1.913445609864562},
                    {12.331515733171500,-9.591178903577841,2.826891219729109},
                    {12.331515733171557,-12.331515733171543,4.044818699548534},
                    {12.331515733171557,-15.071852562765244,5.567228049322811},
                    {12.331515733171500,-17.812189392358860,7.394119269051892},
                    {12.331515733171557,-20.552526221952620,9.525492358735889},
                    {12.331515733171528,-23.292863051546192,11.961347318374699}};

    std::vector<std::vector<double>> DCoeffs = circle.getDCoeffs();

    ASSERT_EQ(DCoeffsExact.size(),DCoeffs.size()) << "Different Number of DCoeffs";
    ASSERT_EQ(DCoeffsExact[0].size(),DCoeffs[0].size())<< "Different Size of DCoeffs";
    // Might need to change this for approximates
    for(int i = 0; i < DCoeffsExact.size(); i++) {
        for(int j = 0; j <DCoeffsExact[0].size(); j++) {
            EXPECT_EQ(DCoeffsExact[i][j],DCoeffs[i][j]) << "DCoeffs Differ in Span " << i << " Index " << j;
        }
    }
}
// =================================== TESTS FOR KNOT VECTOR OPERATIONS =======================================
// int Spline::findSpan(double u) 
TEST(KnotVectorOperations,SpanFinding) {

}
// std::vector<std::vector<double>> Spline::BasisCoefficients(int i)
TEST(KnotVectorOperations,BasisCoeffs) {
    // Examples in Book
}
// std::vector<std::vector<double>> Spline::BasisCoefficientBounds(int i)
TEST(KnotVectorOperations,BasisBounds) {
    // Examples in Book
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
// TEST(SplineProperties,CurveCoefficients) { // Redundant with getXCoeffs, getYCoeffs, getDCoeffs
    
// }

// double Spline::getArcLength()
TEST(SplineProperties,ArcLength) {
    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    double ALExact = 6.283185307900089;
    double AL = circle.getArcLength();
    
    EXPECT_EQ(AL,ALExact) << "Circle Arc Length Fail";

    // Blob Case

}

// double Spline::getCurvature(double u)
TEST(SplineProperties,Curvature) {
    // Figure out Later
}

// double Spline::getSurfaceEnergy()
TEST(SplineProperties,SurfaceEnergy) {
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    double EkExact = 6.283185307900101;
    double Ek = circle.getArcLength();
    
    EXPECT_EQ(Ek,EkExact) << "Circle Surface Energy Fail";

    // Blob Case
}
// double Spline::getArea()
TEST(SplineProperties,Area) {
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    double AExact = 3.141592257512612;
    double A = circle.getArea();
    
    EXPECT_EQ(A,AExact) << "Circle Area Fail";

    // Blob Case
}


// =================================== TESTS FOR CLIPPING =======================================
// std::vector<double> Spline::lineCurveIntersection(std::vector<double> P1, std::vector<double> P2)
TEST(ClippingOperations,LineClip) {
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    // Simple Circle
    std::vector<double> P1 = {-1.5,-1.2};
    std::vector<double> P2 = {1.5,1.2};
    std::vector<double> u1Exact = {0.650760969641606,0.190410828172100};
    std::vector<double> u1 = circle.lineCurveIntersection(P1,P2);
    EXPECT_EQ(u1[0],u1Exact[0]) << "Simple Intersection Fail 1";
    EXPECT_EQ(u1[1],u1Exact[1]) << "Simple Intersection Fail 2";

    // Through Interpolation Points
    P1 = {0,0};
    P2 = {0,1.5};
    u1Exact = {0.333333333333334};
    u1 = circle.lineCurveIntersection(P1,P2);
    EXPECT_EQ(u1[0],u1Exact[0]) << "Through Interpolation Fail"; // Something is going wrong here, with other choices (at least in matlab)

    // No intersection Points
    P1 = {0,1.5};
    P2 = {0,1.5};
    u1Exact = {};
    u1 = circle.lineCurveIntersection(P1,P2);
    EXPECT_EQ(u1.size(),u1Exact.size()) << "No Intersection Fail";

}
// std::vector<std::vector<double>> Spline::getParameterLoop(std::vector<std::vector<double>> square)
TEST(ClippingOperations,ParamLoop) {
    
}
// double Spline::integrateSplineSquare(std::vector<std::vector<double>> square)
TEST(ClippingOperations,ClippedArea) {
    // Quarter Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    std::vector<std::vector<double>> square = {{{0,0},{1.5,0},{1.5,1.5},{0,1.5},{0,0}}};

    double A = circle.integrateSplineSquare(square);
    double AExact = 0.785398111735192;

    EXPECT_EQ(A,AExact) << "Quarter Circle Area Fail";


}
// double Spline::integratedSpline(double u)
TEST(ClippingOperations,analyticIntegration) {
    // Don't yet know how to test this well

}

// =================================== TESTS FOR VISUALIZATION =======================================
// std::vector<std::vector<double>> Spline::makeRationalQuadCurve(std::vector<double> uset)
TEST(Visualization,curveMaker) {

}
