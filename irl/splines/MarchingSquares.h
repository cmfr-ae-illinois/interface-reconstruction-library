#ifndef EXAMPLES_MARCHING_CUBES_H_
#define EXAMPLES_MARCHING_CUBES_H_
#include <iostream>
#include <array>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>

#include "examples/spline_optimization/basic_mesh.h"
#include "examples/spline_optimization/irl2d.h"
#include "examples/spline_optimization/data.h"
#include "irl/splines/Spline.h"

template <class ScalarType>
class MarchingSquares;

template <class ScalarType>
class MarchingSquares {
    private: 
        BasicMesh mesh;
        ScalarType tresh;
        int interpolationType; // 0 for Midpoint, 1 for Linear, 2 for Manson
        std::vector<std::vector<ScalarType>> result; // Keeps track of previous result from run
        /*
            Here we have the cases. These have the index of the edges for a square with the following convenction
                            * -- 2 -- *
                            |         |
                            3         1
                            |         |
                            * -- 0 -- *
            For most cases, only 2 indices are used. The second two are filled with -1 to indicate they are not used. 
            For cases 5 and 10, we use all 4 indicies
            For cases 0 and 15, no indices are used
        */
        static constexpr int caseEdges[16][4] = { // Each index has the index of points for 
            {-1,-1,-1,-1}, // Case 0
            {0 ,3 ,-1,-1}, // Case 1
            {1 ,0 ,-1,-1}, // Case 2
            {1 ,3 ,-1,-1}, // Case 3
            {2 ,1 ,-1,-1}, // Case 4
            {0 ,1 ,2 , 3}, // Case 5
            {2 ,0 ,-1,-1}, // Case 6
            {2 ,3 ,-1,-1}, // Case 7
            {3 ,2 ,-1,-1}, // Case 8
            {0 ,2 ,-1,-1}, // Case 9
            {1 ,2 ,3 , 0}, // Case 10
            {1 ,2 ,-1,-1}, // Case 11
            {3 ,1 ,-1,-1}, // Case 12
            {0 ,1 ,-1,-1}, // Case 13
            {3 ,0 ,-1,-1}, // Case 14
            {-1,-1,-1,-1}, // Case 15
        };

        // Midpoint Interpolation Function 
        // No function needed, mu = 0.5
        // Linear Interpolation Function
        ScalarType LinearInterpolation(ScalarType VOF0,ScalarType VOF1);
        // Manson Interpolation Function
        static ScalarType MansonInterpolation(ScalarType VOF0,ScalarType VOF1); // Makes needed Assumptions
        

        // Perpendicular Distance Calculator
        static ScalarType PerpendicularDistance(std::vector<ScalarType> P1, std::vector<ScalarType> P2, std::vector<ScalarType> P);
        
    public:
        // Methods
        // Constructors
        MarchingSquares();
        MarchingSquares(BasicMesh* m,ScalarType tre,int inTy);

        // Static Methods
        // Decimation Algorithms
        static std::vector<std::vector<ScalarType>> RDPDecimation(std::vector<std::vector<ScalarType>> polyline,double epsilon, bool closed);
        
        // Dynamic Methods
        
        // Function which takes in VOF Field and Returns all edges
        std::vector<std::vector<ScalarType>> run(std::vector<std::vector<ScalarType>> VOF);

        // Function which Takes in result and returns the closed contour points
        // Currently assumes only one closed contour.
        std::vector<std::vector<ScalarType>>  vertexPoints(std::vector<std::vector<ScalarType>> VOF);

        
        // Getters 
        BasicMesh getMesh(); // Done
        ScalarType getTresh();// Done
        int getInterpolationType();// Done
        std::vector<std::vector<ScalarType>> getPreviousResult();
        // Setters
        void setMesh(BasicMesh* m);// Done
        void setTresh(ScalarType t);// Done
        void setInterpolationType(int i);// Done

        // Testing
        void printMesh();
};

#include "MarchingSquares.tpp"

#endif