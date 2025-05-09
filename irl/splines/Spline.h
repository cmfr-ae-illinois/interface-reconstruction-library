#ifndef IRL_SPLINE_SPLINE_H_
#define IRL_SPLINE_SPLINE_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <fstream>
#include <cmath>
#include <iomanip>

template <class ScalarType>
class Spline;

template <class ScalarType>
class Spline { // : public Expr<Spline<ScalarType>>
    private:
        // Variables
        std::vector<std::vector<ScalarType>> ControlPoints;
        std::vector<ScalarType> KnotVector;
        std::vector<ScalarType> Weights;

        std::vector<std::vector<ScalarType>> numerCoeffsX;
        std::vector<std::vector<ScalarType>> numerCoeffsY;
        std::vector<std::vector<ScalarType>> denomCoeffs;
        std::vector<std::vector<ScalarType>> spans;
        std::vector<ScalarType> breakpoints;

        
    public:
        // Variables
        
        // For Taylor Series
        static constexpr double ASeries[41][2] = {
            {-1.666666666666667e-1, 3.333333333333333e-1},
            {1.333333333333333e-1, 1.333333333333333e-1},
            {-8.57142857142857e-2, -8.57142857142857e-2},
            {5.079365079365079e-2, 5.079365079365079e-2},
            {-2.886002886002886e-2, -2.886002886002886e-2},
            {1.598401598401598e-2, 1.598401598401598e-2},
            {-8.7024087024087e-3, -8.7024087024087e-3},
            {4.680287033228209e-3, 4.680287033228209e-3},
            {-2.494100326917664e-3, -2.494100326917664e-3},
            {1.319629802601939e-3, 1.319629802601939e-3},
            {-6.942400265862374e-4, -6.942400265862374e-4},
            {3.635293230124298e-4, 3.635293230124298e-4},
            {-1.896186900898167e-4, -1.896186900898167e-4},
            {9.8581600152796e-5, 9.8581600152796e-5},
            {-5.110797242944491e-5, -5.110797242944491e-5},
            {2.643159786250081e-5, 2.643159786250081e-5},
            {-1.364059246832631e-5, -1.364059246832631e-5},
            {7.026314721363631e-6, 7.026314721363631e-6},
            {-3.613247313977594e-6, -3.613247313977594e-6},
            {1.855325963531499e-6, 1.855325963531499e-6},
            {-9.5139389525278e-7, -9.5139389525278e-7},
            {4.872747569336991e-7, 4.872747569336991e-7},
            {-2.492924046595037e-7, -2.492924046595037e-7},
            {1.274112023814322e-7, 1.274112023814322e-7},
            {-6.505882474542089e-8, -6.505882474542089e-8},
            {3.319227587011661e-8, 3.319227587011661e-8},
            {-1.692109727924127e-8, -1.692109727924127e-8},
            {8.61997418253746e-9, 8.61997418253746e-9},
            {-4.388255621981843e-9, -4.388255621981843e-9},
            {2.232577761324849e-9, 2.232577761324849e-9},
            {-1.135189009858826e-9, -1.135189009858826e-9},
            {5.768900973178349e-10, 5.768900973178349e-10},
            {-2.930192705126503e-10, -2.930192705126503e-10},
            {1.48761649851833e-10, 1.48761649851833e-10},
            {-7.549006672265758e-11, -7.549006672265758e-11},
            {3.82916346272267e-11, 3.82916346272267e-11},
            {-1.941527696469384e-11, -1.941527696469384e-11},
            {9.84052647841976e-12, 9.84052647841976e-12},
            {-4.985823042530466e-12, -4.985823042530466e-12},
            {2.525266498274373e-12, 2.525266498274373e-12},
            {-1.278606320361211e-12, -1.278606320361211e-12}};
        
        // Methods
        // Constructors
        Spline(std::vector<std::vector<ScalarType>> CP,std::vector<ScalarType> KV,std::vector<ScalarType> W);
        static float add(float x, float y);

        // Static Methods ***********************
        // Package for Optimization
        static std::vector<ScalarType> pack(std::vector<std::vector<ScalarType>> Q, std::vector<std::vector<ScalarType>> T);
        // Unpack from Optimization
        static std::vector<std::vector<std::vector<ScalarType>>> unpack(std::vector<ScalarType> V);

        // Tangent Finding with Akima Method
        static std::vector<std::vector<ScalarType>> AkimaTangents(std::vector<std::vector<ScalarType>> Q); 
        // Tangent Finding with Bessel Method
        static std::vector<std::vector<ScalarType>> BesselTangents(std::vector<std::vector<ScalarType>> Q);
        // A helper function for the BesselTangents method to create the u vector which is used in generating the tangents
        static std::vector<ScalarType> BesselTangentUVec(std::vector<std::vector<ScalarType>> Q);
        // Solve the intersection problem of two points (Q1,Q2) with tagent directions (T1,T2). 
        static std::vector<std::vector<ScalarType>> solvePointTangentIntersection(std::vector<ScalarType> Q1, 
                                                                              std::vector<ScalarType> Q2,
                                                                              std::vector<ScalarType> T1,
                                                                              std::vector<ScalarType> T2);
        // A global interpolation Method which returns. The points are enforced at parameter values given by uset
        static Spline GlobalPointInterp(std::vector<std::vector<ScalarType>> Q,int p,std::vector<ScalarType> uset);

        // A local non-rational quadratic interpolation method, returning a spline object using points Q and tangents T.
        static Spline LocalNRQuadInterp(std::vector<std::vector<ScalarType>> Q,std::vector<std::vector<ScalarType>> T);

        // A local rational quadratic interpolation, returning a spline object using points Q and tangents T.
        static Spline LocalRQuadInterp(std::vector<std::vector<ScalarType>> Q,std::vector<std::vector<ScalarType>> T);
        
        // A method to make knot vectors of given continuity levels
        static std::vector<ScalarType> makeKnotVector(std::vector<ScalarType> breakpoints, std::vector<int> continuity, int p);

        // A helper function in the rational quadratic interpolation method for generating weights.
        static ScalarType makeWeight(std::vector<ScalarType> Qkm, std::vector<ScalarType> Rk, std::vector<ScalarType> Qk);

        static std::array<ScalarType, 2> coeffsAreaExact(const ScalarType w);

        static std::array<ScalarType, 2> coeffsAreaSeries(const ScalarType w);

        // Dynamic Methods ********************************
        // Returns the value of a B-Spline Basis Method (N_{i,p}(u)) for a given knot vector U
        ScalarType BBasisFunction(int i, int p,ScalarType u);

        // Determines the quadratic coefficients of N_{i,2} for  a given knot vector U
        std::vector<std::vector<ScalarType>> BasisCoefficients(int i);
        
        // A function that returns the bounds of the BasisCoefficients method.
        std::vector<std::vector<ScalarType>> BasisCoefficientBounds(int i);

        // Returns the numerator and denominator of x,y spline on all spans, the spans, 
        // and the breakpoints (May need to break this up into multiple Method. Currently
        // only returns coefficients) ** Add a getSpans and getBreakpoints Method.
        std::vector<std::vector<ScalarType>> CurveCoefficients();
        
        // A function to remove repeated knots and create the breakpoint vector, returning it to the user and storing it in the variable.
        std::vector<ScalarType> makeBreakpoints(); // Privatize 

        // A function to generate spans for the spline
        std::vector<std::vector<ScalarType>> makeSpans(); // Privatize

        // Finds the span in which paramter u lies in, for a degree p basis function.
        int findSpan(int p, ScalarType u);

        // Specifically for quadratic
        int findSpan(ScalarType u);

        // Finds the arc length of the spline
        ScalarType getArcLength();

        // Gets the total area of the closed curve of the spline (Assumed closed)
        ScalarType getArea();

        // gets the signed curvature at the point given by parameter value u
        // This also is supposed to return first and second derivatives, but that may need to be a new Method
        ScalarType getCurvature(ScalarType u);
        
        // Given a parameter u, returns the tangent vector at the point on the spline
        std::vector<ScalarType> getTangent(ScalarType u);
        
        // Clips the intersection between a grid cell and the spline, returjning parameter values and indicators.
        std::vector<std::vector<ScalarType>> getParameterLoop(std::vector<std::vector<ScalarType>> square);

        // returns the total unsigned surface energy of the spline. Used to return signed value also, but this may not be good for this method.
        ScalarType getSurfaceEnergy();

        // Helper Method to return integral of xy' for area finding, evaluated at parameter value u 
        ScalarType integratedSpline(ScalarType u);

        // Helper for calculating curve integral
        std::vector<std::vector<ScalarType>> clippedBezier(ScalarType u1, ScalarType u2);

        // Method to find the area of intersection between the spline and square
        ScalarType integrateSplineSquare(std::vector<std::vector<ScalarType>> square);

        // Method to find the intersection points between a given line and the spline. Line goes from P1 to P2. Returns in order that goes from P1 to P2.
        std::vector<ScalarType> lineCurveIntersection(std::vector<ScalarType> P1, std::vector<ScalarType> P2);
        
        // A method for getting points along non-rational parabolic spline.
        std::vector<ScalarType> makeParabolicCurve(std::vector<ScalarType> uset);

        // A method for getting points along a rational quadartic spline.
        std::vector<std::vector<ScalarType>> makeRationalQuadCurve(std::vector<ScalarType> uset);

        void saveToVTK(const std::string& filename, const int nsamples = 100);
        
        
        // Getters ****************************
        std::vector<std::vector<ScalarType>> getControlPoints();
        std::vector<ScalarType> getKnotVector();
        std::vector<ScalarType> getWeights();
        std::vector<ScalarType> getBreakpoints();
        std::vector<std::vector<ScalarType>> getSpans();
        std::vector<std::vector<ScalarType>> getXCoeffs();
        std::vector<std::vector<ScalarType>> getYCoeffs();
        std::vector<std::vector<ScalarType>> getDCoeffs();

        // Setters *********************************** 
        void setControlPoints(std::vector<std::vector<ScalarType>> input);
        void setKnotVector(std::vector<ScalarType> input);
        void setWeights(std::vector<ScalarType> input);

        // Printers
        void printControlPoints();
        void printKnotVector();
        void printWeights();
        void printXCoeffs();
        void printYCoeffs();
        void printDCoeffs();
};

#include "irl/splines/Spline.tpp"

#endif