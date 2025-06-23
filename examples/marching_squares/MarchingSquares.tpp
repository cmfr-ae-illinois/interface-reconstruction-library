#include "MarchingSquares.h"
#include <cmath>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

#include "examples/spline_optimization/basic_mesh.h"
#include "examples/spline_optimization/irl2d.h"
#include "examples/spline_optimization/data.h"
#include "irl/splines/Spline.h"
// Constructor
template<class ScalarType>
MarchingSquares<ScalarType>::MarchingSquares() = default;

template<class ScalarType>
MarchingSquares<ScalarType>::MarchingSquares(BasicMesh* m,ScalarType tre,int inTy) {
    this->mesh = *m;
    this->tresh = tre;
    this->interpolationType = inTy;
    this->result = {{-1}};
}

// Private Functions
template<class ScalarType>
ScalarType MarchingSquares<ScalarType>::LinearInterpolation(ScalarType VOF0,ScalarType VOF1) {
    // Here we do Linear Interpolation 
    ScalarType mu = (tresh - VOF0)/(VOF1-VOF0);
    return  mu;
}

template<class ScalarType>
ScalarType MarchingSquares<ScalarType>::MansonInterpolation(ScalarType a1,ScalarType a2) { 
    // a1 and a2 are the VOF values at the adjacent points
    // Here we do Manson Interpolation - Assumes a1  < a2
    ScalarType mu = 0.5;
    if(a1 > 0.5 || a2 < 0.5) {
        mu = 0.5;
    } else { // Here we are within the 4 cases to check
        if(a2 <= 3*a1 || 3*a2 <= a1+2) { // Either or both True
            if(a2 <= 3*a1 && 3*a2 <= a1+2) { // Both True, Case 1
                mu = (a1-0.5)/(a1-a2);
            } else { // Only 1 true, check to see which case. Do case 2 by default
                // Case 2 Default
                mu = 1.5-a1-a2; 
                // Case 3 Check
                if((2*a1+2*a2-1)*(2*a1+2*a2-1) < 4*a1*(a1+a2)) {
                    mu = 1 - (2*a2-1)/(8*a1+4*a2-8*std::sqrt(a1*(a1+a2)));
                }
                // Case 4 Check
                double a1b = 1-a2;
                double a2b = 1-a1; 
                double a1f = 1-a1;
                double a2f = 1-a2;
                if((2*a1b+2*a2b-1)*(2*a1b+2*a2b-1) < 4*a1b*(a1b+a2b)) {
                    // Manson Method
                    mu = (2-a2b-1)/(8*a1b+4*a2b-8*std::sqrt(a1b*(a1b+a2b)));
                }
            }
        } else { // Neither True, Case 2
            mu = 1.5-a1-a2; 
        }
    }
    return mu;
}

// Dynamic Methods
template<class ScalarType>
std::vector<std::vector<ScalarType>> MarchingSquares<ScalarType>::run(std::vector<std::vector<ScalarType>> VOF) {
    const int a_nx = mesh.imax() +1;
    // Treshhold Meeting Matrix
    std::vector<std::vector<bool>> meetsTreshold(a_nx, std::vector<bool>(a_nx,false));
    for(int i = 0; i < meetsTreshold.size(); i++) { // Note that i corresponds to x, j corresponds to y
        for(int j = 0; j < meetsTreshold.size(); j++) {
            meetsTreshold[i][j] = VOF[i][j] >= tresh;
        }
    }

    // Case Matrix
    std::vector<std::vector<int>> cases(a_nx-1,std::vector<int>(a_nx-1,-1));
    for(int i = 0; i < cases.size(); i++){
        for(int j = 0; j < cases.size();j++) {
            cases[i][j] = meetsTreshold[i][j] + 2*meetsTreshold[i+1][j]+4*meetsTreshold[i+1][j+1]+8*meetsTreshold[i][j+1];
        }
    }

    // Edge Values
    std::vector<std::vector<ScalarType>> edges = {};
    for(int i = 0; i < cases.size(); i++) {
        for(int j = 0; j < cases.size(); j++) {
            std::vector<std::vector<ScalarType>> Points = {{0,0},{0,0},{0,0},{0,0}};
            std::vector<std::vector<ScalarType>> Shift = {{0,0},{1,0},{1,1},{0,1},{0,0}};
            
            for(int k = 0; k < 4; k++) {
                ScalarType gamma0 = VOF[i+Shift[k][0]][j+Shift[k][1]];
                ScalarType gamma1 = VOF[i+Shift[k+1][0]][j+Shift[k+1][1]];
                ScalarType x0 = ScalarType(mesh.xm(i+Shift[k][0]));
                ScalarType x1 = ScalarType(mesh.xm(i+Shift[k+1][0]));
                ScalarType y0 = ScalarType(mesh.ym(j+Shift[k][1]));
                ScalarType y1 = ScalarType(mesh.ym(j+Shift[k+1][1]));


                // Ensure Proper ordering (a1 < a2)
                ScalarType a1 = 0;
                ScalarType a2 = 0;
                std::vector<ScalarType> P1 = {0,0};
                std::vector<ScalarType> P2 = {0,0};
                if(gamma0 < gamma1) {
                    a1 = gamma0;
                    a2 = gamma1;
                    P1 = {x0,y0};
                    P2 = {x1,y1};
                } else {
                    a1 = gamma1;
                    a2 = gamma0;
                    P1 = {x1,y1};
                    P2 = {x0,y0};
                }

                // Now do Interpolation
                double mu = 0.5; // Default Value
                switch(interpolationType) {
                    case 0:// Midpoint Interpolation
                        break;
                    case 1: // Linear Interpolation
                        mu = this->LinearInterpolation(a1,a2);
                        break;
                    case 2: // Manson Interpolation
                        mu = this->MansonInterpolation(a1,a2);
                        break;
                    default: // Use Midpoint by Default
                        break;
                }
                
                // Now add the Point in 
                Points[k][0] = P1[0]*(1-mu)+P2[0]*mu;
                Points[k][1] = P1[1]*(1-mu)+P2[1]*mu;
            }

            // Now tha we have the case and the edge Points, we add the specific values to the edge list
            int inds[4] = {0,0,0,0};
            for(int k = 0; k < 4; k++) {
                inds[k] = caseEdges[cases[i][j]][k];
            }
            if(inds[0] != -1) { // If not cases 0,15, add in index 0,1
                edges.push_back({Points[inds[0]][0],Points[inds[0]][1],Points[inds[1]][0],Points[inds[1]][1]});
            }
            if(inds[2] != -1) { // Add in second edge for cases 5,10
                edges.push_back({Points[inds[2]][0],Points[inds[2]][1],Points[inds[3]][0],Points[inds[3]][1]});
            }
        }
    }
    result = edges;
    return edges;
}

template<class ScalarType>
std::vector<std::vector<ScalarType>> MarchingSquares<ScalarType>::vertexPoints(std::vector<std::vector<ScalarType>> VOF) {
    // Get Edges
    std::vector<std::vector<ScalarType>> edges = this->run(VOF);

    // Loop over edges and get ordered Results
    std::vector<std::vector<ScalarType>> ret = {{edges[0][0],edges[0][1]}};
    for(int i = 1; i <edges.size(); i++) { // This does not repeat the last point
        std::vector<ScalarType> Pcurr = ret[i-1];
        // Loop over edges and find next point
        for(int j = 0; j < edges.size(); j++) {
            std::vector<ScalarType> Pfirst = edges[j];
            if(fabs(Pfirst[0]-Pcurr[0]) < 1e-9 && fabs(Pfirst[1]-Pcurr[1]) < 1e-9) {
                ret.push_back({Pfirst[2],Pfirst[3]});
                break;
            }
        }
    }
    // Now we just need to return the ordered set of points
    return ret;
}


// Getters 
template<class ScalarType>
BasicMesh MarchingSquares<ScalarType>::getMesh() {
    return mesh;
}
template<class ScalarType>
ScalarType MarchingSquares<ScalarType>::getTresh() {
    return tresh;
}
template<class ScalarType>
int MarchingSquares<ScalarType>::getInterpolationType() {
    return interpolationType;
}
template<class ScalarType>
std::vector<std::vector<ScalarType>> MarchingSquares<ScalarType>::getPreviousResult() {
    return result;
}


// Setters
template<class ScalarType>
void MarchingSquares<ScalarType>::setMesh(BasicMesh* m) {
    mesh = *m;
}

template<class ScalarType>
void MarchingSquares<ScalarType>::setTresh(ScalarType t) {
    tresh = t;
}
template<class ScalarType>
void MarchingSquares<ScalarType>::setInterpolationType(int i) {
    interpolationType = i;
}

// Testing
template<class ScalarType>
void MarchingSquares<ScalarType>::printMesh() {
    std::cout << "Mesh Properties:\n";
    std::cout << mesh.lx() << "\n";
    std::cout << mesh.ly() << "\n";
    std::cout << mesh.imax() << "\n";
    std::cout << mesh.jmax() << "\n";
}