// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com> ?????????????????
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// #include "Spline.h"
#include "irl/splines/Spline.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <fstream>
#include <cmath>
#include <iomanip>

#include <complex>

#include "gtest/gtest.h" // I don't know why this is happening

// Remaining Tests:

// =================================== TESTS FOR GETTERS =======================================
// std::vector<std::vector<double>> Spline::getControlPoints()
TEST(Spline2DGetterTests,ControlPoints) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};

    Spline spline = Spline(CP,KV,W);
    std::vector<std::vector<double>> splineCP = spline.getControlPoints();
    ASSERT_EQ(CP.size(),splineCP.size()) << "Different Number of Control Points";
    ASSERT_EQ(CP[0].size(),splineCP[0].size())<< "Different Size of Control Points";

    for(int i = 0; i < CP.size(); i++) {
        for(int j = 0; j <CP[0].size(); j++) {
            EXPECT_EQ(splineCP[i][j],CP[i][j]) << "Control Point " << i << " differ in index " << j;
        }
    }


}
// std::vector<double> Spline::getKnotVector()
TEST(Spline2DGetterTests,KnotVector) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};

    Spline spline = Spline(CP,KV,W);
    std::vector<double> splineKV = spline.getKnotVector();
    ASSERT_EQ(KV.size(),splineKV.size()) << "Different Number of Knots";

    for(int i = 0; i < KV.size(); i++) {
        EXPECT_EQ(splineKV[i],KV[i]) << "Knot Vectors Differ at index " << i;
    }
}
// std::vector<double> Spline::getWeights() 
TEST(Spline2DGetterTests,Weights) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};

    Spline spline = Spline(CP,KV,W);
    std::vector<double> splineW = spline.getWeights();
    ASSERT_EQ(W.size(),splineW.size()) << "Different Number of Weights";

    for(int i = 0; i < W.size(); i++) {
        EXPECT_EQ(splineW[i],W[i]) << "Weights Differ at index " << i;
    }
}
// std::vector<double> Spline::getBreakpoints()
TEST(Spline2DGetterTests,Breakpoints) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};
    std::vector<double> BP = {0,0.25,0.5,1};

    Spline spline = Spline(CP,KV,W);
    std::vector<double> BPs = spline.getBreakpoints();

    ASSERT_EQ(BP.size(),BPs.size()) << "Different Breakpoints Size";

    for(int i = 0; i < BPs.size(); i++) {
        EXPECT_EQ(BPs[i],BP[i]) << "Breakpoints Differ at Index" << i;
    }
}
// std::vector<std::vector<double>> Spline::getSpans() 
TEST(Spline2DGetterTests,Spans) {
    std::vector<std::vector<double>> CP = {{1,2},{2,3},{5,2},{4,3},{4,2}};
    std::vector<double> KV = {0,0,0,0.25,0.5,1,1,1};
    std::vector<double> W = {1,0.25,1,0.5,1};
    std::vector<std::vector<double>> S = {{0,0.25},{0.25,0.5},{0.5,1}};

    Spline spline = Spline(CP,KV,W);
    std::vector<std::vector<double>> Spans = spline.getSpans();

    ASSERT_EQ(S.size(),Spans.size()) << "Different Number of Control Points";
    ASSERT_EQ(S[0].size(),Spans[0].size())<< "Different Size of Control Points";

    for(int i = 0; i < S.size(); i++) {
        for(int j = 0; j <S[0].size(); j++) {
            EXPECT_EQ(Spans[i][j],S[i][j]) << "Spans " << i << " differ in index " << j;
        }
    }
}
// std::vector<std::vector<double>> Spline::getXCoeffs()
TEST(Spline2DGetterTests,XCoeffs) {
    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    std::vector<std::vector<double>> XCoeffsExact = {{-2.848208747679525,-0.685084207398420,1.0},
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
            EXPECT_NEAR(XCoeffsExact[i][j],XCoeffs[i][j],1e-12) << "XCoeffs Differ in Span " << i << " Index " << j;
        }
    }
}
// std::vector<std::vector<double>> Spline::getYCoeffs()
TEST(Spline2DGetterTests,YCoeffs) {
    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    std::vector<std::vector<double>> YCoeffsExact = {{-1.179766691758548,3.444150891285807,0.0},
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
            EXPECT_NEAR(YCoeffsExact[i][j],YCoeffs[i][j],1e-12) << "YCoeffs Differ in Span " << i << " Index " << j;
        }
    }
}
// std::vector<std::vector<double>> Spline::getDCoeffs()
TEST(Spline2DGetterTests,DCoeffs) {
    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    std::vector<std::vector<double>> DCoeffsExact = {{3.082878933292889,-0.685084207398420,1.0},
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
            EXPECT_NEAR(DCoeffsExact[i][j],DCoeffs[i][j],1e-12) << "DCoeffs Differ in Span " << i << " Index " << j;
        }
    }
}
// =================================== TESTS FOR KNOT VECTOR OPERATIONS =======================================
// int Spline::findSpan(double u) 
TEST(Spline2DKnotVectorOperations,SpanFinding) {
    std::vector<double> U = {0,0,0,0.2,0.5,0.6,0.6,0.9,1,1,1};
    std::vector<double> W = {1,1,1,1,1,1,1,1};
    std::vector<std::vector<double>> CP = {{1,2},{1,2},{1,2},{1,2},{1,2},{1,2},{1,2},{1,2}};

    Spline SpanTester = Spline(CP,U,W);

    EXPECT_EQ(SpanTester.findSpan(-1.24),-2) << "Negative Value Fail";
    EXPECT_EQ(SpanTester.findSpan(24),-1) << "Large Value Fail";
    EXPECT_EQ(SpanTester.findSpan(0),0) << "Start Value Fail";
    EXPECT_EQ(SpanTester.findSpan(1),4) << "End Value Fail";
    EXPECT_EQ(SpanTester.findSpan(0.2),1) << "Breakpoint 1 Fail";
    EXPECT_EQ(SpanTester.findSpan(0.5),2) << "Breakpoint 2 Fail";
    EXPECT_EQ(SpanTester.findSpan(0.1),0) << "In Span 0 Fail";
    EXPECT_EQ(SpanTester.findSpan(0.4),1) << "In Span 1 Fail";
    EXPECT_EQ(SpanTester.findSpan(0.5959),2) << "In Span 2 Fail";
    EXPECT_EQ(SpanTester.findSpan(0.999),4) << "In Span 4 Fail";
}

// std::vector<std::vector<double>> Spline::BasisCoefficients(int i)
TEST(Spline2DKnotVectorOperations,BasisCoeffs) {
    std::vector<double> U = {0,0,0,1,2,3,4,4,5,5,5};
    std::vector<double> W = {1,1,1,1,1,1,1,1};
    std::vector<std::vector<double>> CP = {{1,2},{1,2},{1,2},{1,2},{1,2},{1,2},{1,2},{1,2}};

    Spline BasisCoeffChecker = Spline(CP,U,W);

    // Exact Values
    std::vector<std::vector<double>> N0 = {{0,0,0}      ,{0,0,0}        ,{1,-2,1}};
    std::vector<std::vector<double>> N1 = {{0,0,0}      ,{-1.5,2,0}     ,{0.5,-2,2}};
    std::vector<std::vector<double>> N2 = {{0.5,0,0}    ,{-1,3,-1.5}    ,{0.5,-3,4.5}};
    std::vector<std::vector<double>> N3 = {{0.5,-1,0.5} ,{-1,5,-5.5}    ,{0.5,-4,8}};
    std::vector<std::vector<double>> N4 = {{0.5,-2,2}   ,{-1.5,10,-16}  ,{0,0,0}};
    std::vector<std::vector<double>> N5 = {{1,-6,9}     ,{0,0,0}        ,{1,-10,25}};
    std::vector<std::vector<double>> N6 = {{0,0,0}      ,{-2,18,-40}    ,{0,0,0}};
    std::vector<std::vector<double>> N7 = {{1,-8,16}    ,{0,0,0}        ,{0,0,0}};

    // Combine into one vector
    std::vector<std::vector<std::vector<double>>> Exacts = {N0,N1,N2,N3,N4,N5,N6,N7};

    // Checking 
    std::vector<std::vector<double>> NCurr;
    std::vector<std::vector<double>> res;
    for(int i = 0; i < Exacts.size(); i++) {
        NCurr = Exacts[i];
        res = BasisCoeffChecker.BasisCoefficients(i);
        
        // Compare Sizes
        ASSERT_EQ(res.size(),NCurr.size()) << "Different Number of Basis Coeffs";
        ASSERT_EQ(res[0].size(),NCurr[0].size())<< "Different Size of Basis Coeffs";

        for(int j = 0; j < res.size(); j++) {
            for(int k = 0; k < res[0].size(); k++){
                EXPECT_NEAR(res[j][k],NCurr[j][k],1e-12) << "Basis Function " << i << " Differs in span " << (j+1) << " Index " << k;
            }
        }
    }
}
// std::vector<std::vector<double>> Spline::BasisCoefficientBounds(int i)
TEST(Spline2DKnotVectorOperations,BasisBounds) { // Basically the same as above, but with the bounds now
    std::vector<double> U = {0,0,0,1,2,3,4,4,5,5,5};
    std::vector<double> W = {1,1,1,1,1,1,1,1};
    std::vector<std::vector<double>> CP = {{1,2},{1,2},{1,2},{1,2},{1,2},{1,2},{1,2},{1,2}};

    Spline BasisCoeffChecker = Spline(CP,U,W);

    // Exact Values
    std::vector<std::vector<double>> N0 = {{0,0},{0,0},{0,1}};
    std::vector<std::vector<double>> N1 = {{0,0},{0,1},{1,2}};
    std::vector<std::vector<double>> N2 = {{0,1},{1,2},{2,3}};
    std::vector<std::vector<double>> N3 = {{1,2},{2,3},{3,4}};
    std::vector<std::vector<double>> N4 = {{2,3},{3,4},{4,4}};
    std::vector<std::vector<double>> N5 = {{3,4},{4,4},{4,5}};
    std::vector<std::vector<double>> N6 = {{4,4},{4,5},{5,5}};
    std::vector<std::vector<double>> N7 = {{4,5},{5,5},{5,5}};

    // Combine into one vector
    std::vector<std::vector<std::vector<double>>> Exacts = {N0,N1,N2,N3,N4,N5,N6,N7};

    // Checking 
    std::vector<std::vector<double>> NCurr;
    std::vector<std::vector<double>> res;
    for(int i = 0; i < Exacts.size(); i++) {
        NCurr = Exacts[i];
        res = BasisCoeffChecker.BasisCoefficientBounds(i);

        // Compare Sizes
        ASSERT_EQ(res.size(),NCurr.size()) << "Different Number of Basis Bounds";
        ASSERT_EQ(res[0].size(),NCurr[0].size())<< "Different Size of Basis Bounds";
        for(int j = 0; j < res.size(); j++) {
            for(int k = 0; k < res[0].size(); k++){
                EXPECT_EQ(res[j][k],NCurr[j][k]) << "Basis Bound " << i << " Differs in span " << (j+1) << " Index " << k;
            }
        }
    }
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

TEST(Spline2DInterpolationMethods,PTInter) {
    // Case 1: Test Two Regular Sloped Lines
    std::vector<double> Q1 = {1.0,2.0};
    std::vector<double> Q2 = {5.0,1.0};
    std::vector<double> T1 = {1.0,1.0};
    std::vector<double> T2 = {-1.0,1.0};
    std::vector<double> ExactIntersection = {2.5,3.5};

    std::vector<std::vector<double>> Rg1 = Spline::solvePointTangentIntersection(Q1,Q2,T1,T2);
    std::vector<double> R1 = Rg1[0];

    for(int i = 0; i < ExactIntersection.size();i++) {
        EXPECT_NEAR(ExactIntersection[i],R1[i],1e-12) << "Case 1 Intersection differs at index " << i;
    }

    // Case 2: One Vertical Line
    T1 = {0.0,1.0};
    ExactIntersection = {1.0,5.0};
    Rg1 = Spline::solvePointTangentIntersection(Q1,Q2,T1,T2);
    R1 = Rg1[0];

    for(int i = 0; i < ExactIntersection.size();i++) {
        EXPECT_NEAR(ExactIntersection[i],R1[i],1e-12) << "Case 2 Intersection differs at index " << i;
    }

    // Case 3: One Horizontal Line
    T1 = {1.0,0};
    ExactIntersection = {4.0,2.0};
    Rg1 = Spline::solvePointTangentIntersection(Q1,Q2,T1,T2);
    R1 = Rg1[0];

    for(int i = 0; i < ExactIntersection.size();i++) {
        EXPECT_NEAR(ExactIntersection[i],R1[i],1e-12) << "Case 3 Intersection differs at index " << i;
    }

    // Case 4: No Intersections
    // Not in Method yet

    // Case 5: Infinite Intersections
    // Not in Method yet
}

// double Spline::makeWeight(std::vector<double> Qkm, std::vector<double> Rk, std::vector<double> Qk)
TEST(Spline2DInterpolationMethods,WeightCalculator) {
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
    EXPECT_NEAR(w,wExact,1e-12) << "Isosceles Weight Fail";

    // General Case
    Rk = {3,5};
    wExact = 0.597492359431871;
    w = Spline::makeWeight(Qkm,Rk,Qk);
    EXPECT_NEAR(w,wExact,1e-12) << "General Weight Fail";
}


// Spline Spline::LocalRQuadInterp(std::vector<std::vector<double>> Q,std::vector<std::vector<double>> T)
TEST(Spline2DInterpolationMethods,Interpolator) {
    // Compare Control Points, Knot Vector, and Weights for Circle and Blob Cases

    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    std::vector<std::vector<double>> CPExact = {{1,0},
                                                {1,0.414213562373095},
                                                {0.707106781186548,0.707106781186548},
                                                {0.414213562373095,1},
                                                {0,1},
                                                {-0.414213562373095,1},
                                                {-0.707106781186548,0.707106781186548},
                                                {-1,0.414213562373095},
                                                {-1,1.224646799147353e-16},
                                                {-1,-0.414213562373095},
                                                {-0.707106781186548,-0.707106781186548},
                                                {-0.414213562373095,-1},
                                                {-1.836970198721030e-16,-1},
                                                {0.414213562373095,-1},
                                                {0.707106781186548,-0.707106781186548},
                                                {1,-0.414213562373095},
                                                {1,0}};
    std::vector<double> UExact = {0,0,0,
                                  2.0/9.0,2.0/9.0,
                                  3.0/9.0,3.0/9.0,
                                  4.0/9.0,4.0/9.0,
                                  5.0/9.0,5.0/9.0,
                                  6.0/9.0,6.0/9.0,
                                  7.0/9.0,7.0/9.0,
                                  8.0/9.0,8.0/9.0,
                                  1,1,1};
    std::vector<double> WExact =  {1,0.9238795325112871,
                                   1,0.9238795325112871,
                                   1,0.9238795325112871,
                                   1,0.9238795325112871,
                                   1,0.9238795325112871,
                                   1,0.9238795325112871,
                                   1,0.9238795325112871,
                                   1,0.9238795325112871,1};

    std::vector<std::vector<double>> CP = circle.getControlPoints();
    std::vector<double> U = circle.getKnotVector();
    std::vector<double> W = circle.getWeights();
    
    ASSERT_EQ(CP.size(),CPExact.size()) << "Different Number of Control Points";
    ASSERT_EQ(CP[0].size(),CPExact[0].size()) << "Different Size of Control Points";
    ASSERT_EQ(U.size(),UExact.size()) << "Different Number of Knots";
    ASSERT_EQ(W.size(),WExact.size()) << "Different Number of Weights";
    
    // Loop over Control Points
    for(int i = 0; i < CP.size(); i++) {
        for(int j = 0; j < CP[0].size(); j++) {
            EXPECT_NEAR(CP[i][j],CPExact[i][j],1e-12) << "Control Point Error in Index " << i << "," << j;
        }
    }
    // Loop over Knots
    for(int i = 0; i < U.size(); i++) {
        EXPECT_NEAR(U[i],UExact[i],1e-12);
    }
    // Loop Over Weights
    for(int i = 0; i < W.size(); i++) {
        EXPECT_NEAR(W[i],WExact[i],1e-12);
    }


    // Blob Case
    std::vector<std::vector<double>> Blob = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{-3.0,0.0},{1.0,0.0}};
    std::vector<std::vector<double>> BlobT = {{0.0,1.0},{-1.0,0.0},{-1.0,0.0},{-1.0,0.0},{0.0,1.0}};

    Spline blob = Spline::LocalRQuadInterp(Blob,BlobT);

    CPExact = {{1.0,0},
               {1.0,0.414213562373095},
               {0.707106781186548,0.707106781186548},
               {0.414213562373095,1.0},
               {0,1.0},
               {-0.7071067811865481,},
               {-0.5,0.5,},
               {-0.292893218813452,0},
               {-1.0,0},
               {-2.0,0},
               {-3.0,0},
               {-6.0,0},
               {-3.666666666666666,-2.0},
               {1.0,-5.999999999999999},
               {1.0,0}};
    UExact = {0,0,0,
              1.0/4.0,1.0/4.0,
              3.0/8.0,3.0/8.0,
              1.0/2.0,1.0/2.0,
              5.0/8.0,5.0/8.0,
              3.0/4.0,3.0/4.0,
              7.0/8.0,7.0/8.0,
              1.0,1.0,1.0};
    WExact =  {1,0.923879532511287,
               1,0.923879532511287,
               1,0.572915223375334,
               1,0.572915223375333,
               1,0,
               1,0.347167685997458,
               1,0.418043000328400 ,1};

    CP = blob.getControlPoints();
    U = blob.getKnotVector();
    W = blob.getWeights();
    
    ASSERT_EQ(CP.size(),CPExact.size()) << "Different Number of Control Points";
    ASSERT_EQ(CP[0].size(),CPExact[0].size()) << "Different Size of Control Points";
    ASSERT_EQ(U.size(),UExact.size()) << "Different Number of Knots";
    ASSERT_EQ(W.size(),WExact.size()) << "Different Number of Weights";
    
    // Loop over Control Points
    for(int i = 0; i < CP.size(); i++) {
        for(int j = 0; j < CP[0].size(); j++) {
            EXPECT_NEAR(CP[i][j],CPExact[i][j],1e-12);
        }
    }
    // Loop over Knots
    for(int i = 0; i < U.size(); i++) {
        EXPECT_NEAR(U[i],UExact[i],1e-12);
    }
    // Loop Over Weights
    for(int i = 0; i < W.size(); i++) {
        EXPECT_NEAR(W[i],WExact[i],1e-12);
    }
}

// =================================== TESTS FOR SPLINE PROPERTIES =======================================
// std::vector<std::vector<double>> Spline::CurveCoefficients()
// TEST(SplineProperties,CurveCoefficients) { // Redundant with getXCoeffs, getYCoeffs, getDCoeffs
    
// }

// double Spline::getArcLength()
TEST(Spline2DSplineProperties,ArcLength) {
    // Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    double ALExact = 6.283185307900089;
    double AL = circle.getArcLength();
    
    EXPECT_NEAR(AL,ALExact,1e-12) << "Circle Arc Length Fail";

    // Blob Case
    std::vector<std::vector<double>> Blob = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{-3.0,0.0},{1.0,0.0}};
    std::vector<std::vector<double>> BlobT = {{0.0,1.0},{-1.0,0.0},{-1.0,0.0},{-1.0,0.0},{0.0,1.0}};

    Spline blob = Spline::LocalRQuadInterp(Blob,BlobT);
    ALExact = 14.340009149434596;
    AL = blob.getArcLength();
    
    EXPECT_NEAR(AL,ALExact,1e-12) << "Circle Arc Length Fail";
}

// double Spline::getCurvature(double u)
TEST(Spline2DSplineProperties,Curvature) {
    // Unit Circle
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    std::vector<double> uset = {0,0.25,0.5,0.75,1};

    double k;
    for(int i = 0; i <uset.size(); i++) {
        k = circle.getCurvature(uset[i]);
        EXPECT_NEAR(k,1.0,1e-12) << "Unit Circle Curvature Fail at u = " << uset[i];
    }

    // Non-Unit Circle
    Circle = {{2.0,0.0},{0.0,2.0},{-2.0,0.0},{0.0,-2.0},{2.0,0.0}};
    CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    circle = Spline::LocalRQuadInterp(Circle,CircleT);
    uset = {0,0.25,0.5,0.75,1};

    for(int i = 0; i <uset.size(); i++) {
        k = circle.getCurvature(uset[i]);
        EXPECT_NEAR(k,0.5,1e-12) << "R=2 Circle Curvature Fail at u = " << uset[i];
    }
}

// double Spline::getSurfaceEnergy()
TEST(Spline2DSplineProperties,SurfaceEnergy) {
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    double EkExact = 6.283185307900101;
    double Ek = circle.getSurfaceEnergy();
    
    EXPECT_NEAR(Ek,EkExact,1e-12) << "Circle Surface Energy Fail";

    // Blob Case
    std::vector<std::vector<double>> Blob = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{-3.0,0.0},{1.0,0.0}};
    std::vector<std::vector<double>> BlobT = {{0.0,1.0},{-1.0,0.0},{-1.0,0.0},{-1.0,0.0},{0.0,1.0}};

    Spline blob = Spline::LocalRQuadInterp(Blob,BlobT);
    EkExact = 10.210042233157300;
    Ek = blob.getSurfaceEnergy();
    
    EXPECT_NEAR(Ek,EkExact,1e-12) << "Circle Arc Length Fail";
}
// double Spline::getArea()
TEST(Spline2DSplineProperties,Area) {
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    double AExact = 3.141592257512612;
    double A = circle.getArea();
    
    EXPECT_NEAR(A,AExact,1e-12) << "Circle Area Fail";

    // Blob Case
    std::vector<std::vector<double>> Blob = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{-3.0,0.0},{1.0,0.0}};
    std::vector<std::vector<double>> BlobT = {{0.0,1.0},{-1.0,0.0},{-1.0,0.0},{-1.0,0.0},{0.0,1.0}};

    Spline blob = Spline::LocalRQuadInterp(Blob,BlobT);
    AExact = 12.341364192401300;
    A = blob.getArea();
    
    EXPECT_NEAR(A,AExact,1e-12) << "Blob Area Fail";
}


// =================================== TESTS FOR CLIPPING =======================================
// std::vector<double> Spline::lineCurveIntersection(std::vector<double> P1, std::vector<double> P2)
TEST(Spline2DClippingOperations,LineClip) {
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    // Simple Circle
    std::vector<double> P1 = {-1.5,-1.2};
    std::vector<double> P2 = {1.5,1.2};
    std::vector<double> u1Exact = {0.650760969641606,0.190410828172100};
    std::vector<double> u1 = circle.lineCurveIntersection(P1,P2);
    EXPECT_NEAR(u1[0],u1Exact[0],1e-12) << "Simple Intersection Fail 1";
    EXPECT_NEAR(u1[1],u1Exact[1],1e-12) << "Simple Intersection Fail 2";

    // Through Interpolation Points
    P1 = {0,0};
    P2 = {0,1.5};
    u1Exact = {0.333333333333334};
    u1 = circle.lineCurveIntersection(P1,P2);
    EXPECT_NEAR(u1[0],u1Exact[0],1e-12) << "Through Interpolation Fail"; // Something is going wrong here, with other choices (at least in matlab)

    // No intersection Points
    P1 = {0,1.5};
    P2 = {0,1.5};
    u1Exact = {};
    u1 = circle.lineCurveIntersection(P1,P2);
    EXPECT_EQ(u1.size(),u1Exact.size()) << "No Intersection Fail";

}
// std::vector<std::vector<double>> Spline::getParameterLoop(std::vector<std::vector<double>> square)
TEST(Spline2DClippingOperations,ParamLoop) {
    // Quarter Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};
    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    std::vector<std::vector<double>> square = {{{0,0},{1.5,0},{1.5,1.5},{0,1.5},{0,0}}};

    std::vector<double> Ploop = {0,0,1.0/3.0,0};
    std::vector<double> Indicators = {2,4,3,2};

    std::vector<std::vector<double>> res = circle.getParameterLoop(square);

    std::vector<double> resP = res[0];
    std::vector<double> resI = res[1];

    ASSERT_EQ(resP.size(),Ploop.size()) << "Different Number of Parameters";
    ASSERT_EQ(resI.size(),Indicators.size()) << "Different Number of Indicators";
    ASSERT_EQ(resI.size(),resP.size()) << "Different Number of Parameters and Indicators";

    for(int i = 0; i < resP.size(); i++) {
        EXPECT_NEAR(resP[i],Ploop[i],1e-12) << "Parameters differ at index " << i;
        EXPECT_EQ(resI[i],Indicators[i]) << "Indicators differ at index " << i;
    }

    // Blob Case

}
// double Spline::integrateSplineSquare(std::vector<std::vector<double>> square)
TEST(Spline2DClippingOperations,ClippedArea) {
    // Quarter Circle Case
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    std::vector<std::vector<double>> square = {{{0,0},{1.5,0},{1.5,1.5},{0,1.5},{0,0}}};

    double A = circle.integrateSplineSquare(square);
    double AExact = 0.785398111735192;

    EXPECT_NEAR(A,AExact,1e-12) << "Quarter Circle Area Fail";


}
// double Spline::integratedSpline(double u)
// TEST(ClippingOperations,analyticIntegration) { // Redundant with area clipping operations
//

// }

// =================================== TESTS FOR VISUALIZATION =======================================
// std::vector<std::vector<double>> Spline::makeRationalQuadCurve(std::vector<double> uset)
TEST(Spline2DVisualization,curveMaker) {
    // Circle 
    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);
    std::vector<double> uset = {0};
    for(int i = 1; i <= 20; i++) {
        uset.insert(uset.end(),i/20);
    }
    std::vector<std::vector<double>> curve = circle.makeRationalQuadCurve(uset);

    // Check all distances are near 1
    double dist;
    for(int i =0; i < curve.size(); i++) {
        dist = 0;
        for(int j =0; j < curve[0].size(); j++) {
            dist += pow(curve[i][j],2.0);
        }
        dist = sqrt(dist);
        EXPECT_NEAR(dist,1.0,1e-12) << "Circle Radius Wrong at Point "<< i;
    }
    
}
