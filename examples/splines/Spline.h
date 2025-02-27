#include <iostream>
#include <vector>

class Spline {
    private:
        // Variables
        std::vector<std::vector<double>> ControlPoints;
        std::vector<double> KnotVector;
        std::vector<double> Weights;

    public: 
        // Methods
        // Constructors
        Spline(std::vector<std::vector<double>> CP,std::vector<double> KV,std::vector<double> W);
        static float add(float x, float y);

        // Static Methods ***********************

        // Tangent Finding with Akima Method
        static std::vector<std::vector<double>> AkimaTangents(std::vector<std::vector<double>> Q); 
        // Tangent Finding with Bessel Method
        static std::vector<std::vector<double>> BesselTangents(std::vector<std::vector<double>> Q);
        // Solve the intersection problem of two points (Q1,Q2) with tagent directions (T1,T2). 
        static std::vector<double> solvePointTangentIntersection(std::vector<std::vector<double>> Q1, 
                                                                 std::vector<std::vector<double>> Q2,
                                                                 std::vector<std::vector<double>> T1,
                                                                 std::vector<std::vector<double>> T2);
        // A global interpolation Method which returns. The points are enforced at parameter values given by uset
        static Spline GlobalPointInterp(std::vector<std::vector<double>> Q,int p,std::vector<double> uset);

        // A local non-rational quadratic interpolation method, returning a spline object using points Q and tangents T.
        static Spline LocalNRQuadInterp(std::vector<std::vector<double>> Q,std::vector<std::vector<double>> T);

        // A local rational quadratic interpolation, returning a spline object using points Q and tangents T.
        static Spline LocalRQuadInterp(std::vector<std::vector<double>> Q,std::vector<std::vector<double>> T);
        
        // A method to make knot vectors of given continuity levels
        static std::vector<double> makeKnotVector(std::vector<double> breakpoints, std::vector<int> continuity, int p);

        // A helper function in the rational quadratic interpolation method for generating weights.
        static double makeWeight(std::vector<double> Qkm, std::vector<double> Rk, std::vector<double> Qk);



        // Dynamic Methods ********************************
        // Returns the value of a B-Spline Basis Method (N_{i,p}(u)) for a given knot vector U
        double BBasisFunction(int i, int p,double u);

        // Determines the quadratic coefficients of N_{i,2} for  a given knot vector U
        std::vector<std::vector<double>> BasisCoefficients(int i);
        
        // A function that returns the bounds of the BasisCoefficients method.
        std::vector<std::vector<double>> BasisCoefficientBounds(int i);

        // Returns the numerator and denominator of x,y spline on all spans, the spans, 
        // and the breakpoints (May need to break this up into multiple Method. Currently
        // only returns coefficients) ** Add a getSpans and getBreakpoints Method.
        std::vector<std::vector<double>> CurveCoefficients(std::vector<double> U,std::vector<double>, std::vector<std::vector<double>> CP);

        // Finds the span in which paramter u lies in, for a degree p basis function.
        int findSpan(int p, double u, std::vector<double> U);

        // Finds the arc length of the spline
        double getArcLength(std::vector<double> U,std::vector<double> W,std::vector<std::vector<double>> CP);

        // Gets the total area of the closed curve of the spline (Assumed closed)
        double getArea(std::vector<double> U,std::vector<double> W,std::vector<std::vector<double>> CP);

        // gets the signed curvature at the point given by parameter value u
        // This also is supposed to return first and second derivatives, but that may need to be a new Method
        double getCurvature(std::vector<double> U,std::vector<double> W,std::vector<std::vector<double>> CP,double u);

        // Clips the intersection between a grid cell and the spline, returjning parameter values and indicators.
        double getParamterLoop(std::vector<double> U,std::vector<double> W,std::vector<std::vector<double>> CP,std::vector<std::vector<double>> square);

        // returns the total unsigned surface energy of the spline. Used to return signed value also, but this may not be good for this method.
        double getSurfaceEnergy(std::vector<double> U,std::vector<double> W,std::vector<std::vector<double>> CP);

        // Helper Method to return integral of xy' for area finding, evaluated at parameter value u 
        double integratedSpline(std::vector<double> U,std::vector<double> W,std::vector<std::vector<double>> CP,double u);

        // Method to find the area of intersection between the spline and square
        double integrateSplineSquare(std::vector<double> U,std::vector<double> W,std::vector<std::vector<double>> CP, std::vector<std::vector<double>> square);

        // Method to find the intersection points between a given line and the spline. Line goes from P1 to P2. Returns in order that goes from P1 to P2.
        std::vector<double> lineCurveIntersection(std::vector<double> U,std::vector<double> W,std::vector<std::vector<double>> CP, std::vector<double> P1, std::vector<double> P2);
        
        // A method for getting points along non-rational parabolic spline.
        std::vector<double> makeParabolicCurve(std::vector<std::vector<double>> CP, std::vector<double> U, std::vector<double> W, std::vector<double> uset);

        // A method for getting points along a rational quadartic spline.
        std::vector<double> makeRationalQuadCurve(std::vector<std::vector<double>> CP, std::vector<double> U, std::vector<double> W, std::vector<double> uset);



        // Getters ****************************
        std::vector<std::vector<double>> getControlPoints();
        std::vector<double> getKnotVector();
        std::vector<double> getWeights();
        // Setters *********************************** 
        void setControlPoints(std::vector<std::vector<double>> input);
        void setKnotVector(std::vector<double> input);
        void setWeights(std::vector<double> input);
};