#include <chrono>
#include <iostream>
#include <string>

// #include "examples/PUSurfaceTension/pu_neighborhood.h"
// #include "examples/PUSurfaceTension/pu_solve.h"

#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/interface_reconstruction_methods/pu_solve.h"

#include "irl/geometry/polyhedrons/rectangular_cuboid.h"

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/moments/general_moments.h"

#include "irl/variant_reconstruction/separator_variant.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/normal.h"

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"

#include "irl/variant_reconstruction/separator_variant.h"

int main(int argc, char* argv[]) {
    const auto xy_plane = IRL::Plane(IRL::Normal(0.0, 0.0, 1.0), 0.0);
    const auto centroid = IRL::Pt(0.0,1.0,0.0);
    IRL::PUSTNeighborhood<IRL::RectangularCuboid> test;
    IRL::RectangularCuboid test2;
    const IRL::SeparatorVariant xy_plane_sep = IRL::PlanarSeparator::fromOnePlane(xy_plane);
    std::cout << "====== BEGIN NEIGHBORDHOOD TESTS ====\n";
    std::cout << "Test Initial Size: " << ((test.size() == 0) ? "Pass" : "Fail") << std::endl;
    
    test.addMember(&centroid,&xy_plane_sep);
    std::cout << "After Adding Member: " << ((test.size() == 1) ? "Pass" : "Fail") << std::endl;

    test.emptyNeighborhood();
    std::cout << "After Empty: " << ((test.size() == 0) ? "Pass" : "Fail") << std::endl;

    test.resize(25);
    std::cout << "After Resize: " << ((test.size() == 25) ? "Pass" : "Fail") << std::endl;

    // test.setCenterOfStencil(11);
    // std::cout << "Set Center: " << ((test.getCenterOfStencilIndex() == 11) ? "Pass" : "Fail") << std::endl;

    test.setMember(1,&centroid,&xy_plane_sep);

    // test.getCenterCell();
    // test.getCenterCellStoredMoments();
    // test.getCell(1);
    // test.getStoredMoments(0);
    // test.begin();
    // test.end();
    // test.cbegin();
    // test.cend();
    std::cout << "All other functions ran without problems :)\n";


    std::cout << "\n====== BEGIN WENDLAND TESTS ====\n";
    IRL::Pt xi; double delta; IRL::Pt x_eval;
    x_eval = IRL::Pt(0.0,1.0,0);

    delta = 50;
    std::cout << "Wendland Origin: " << (IRL::Wendland::phi(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland x  Origin: " << (IRL::Wendland::dphidx(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland y Origin: " << (IRL::Wendland::dphidy(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland xx Origin: " << (IRL::Wendland::ddphidxx(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland yy Origin: " << (IRL::Wendland::ddphidyy(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland xy Origin: " << (IRL::Wendland::ddphidxy(x_eval,delta, x_eval)) << std::endl;

    // Implicit Surface Tests
    std::cout << "\n====== BEGIN IMPLICIT SURFACE TESTS ====\n";
    std::vector<IRL::Pt> centroids;
    std::vector<IRL::SeparatorVariant> variantSeps;
    std::vector<IRL::PlanarSeparator> planarSeps;
    std::vector<IRL::Paraboloid> paraSeps;
    delta = 5;

    // Construct a Plane
    IRL::Normal nor(1,0,0);
    IRL::Plane plane = IRL::Plane(nor,0);
    IRL::PlanarSeparator sep = IRL::PlanarSeparator::fromOnePlane(plane);
    planarSeps.push_back(sep); planarSeps.push_back(sep); planarSeps.push_back(sep); // 3 Seps
    variantSeps.push_back(sep); variantSeps.push_back(sep); 
    // Construct a Probabola
    IRL::Normal e0(1,0,0);
    IRL::Normal e1(0,1,0);
    IRL::Normal e2(0,0,1);
    IRL::ReferenceFrame refFrame(e0,e1,e2);
    double a = 1;
    double b = 0;
    IRL::Pt datum(1,0,0);
    IRL::Paraboloid para = IRL::Paraboloid(datum,refFrame,a,b);
    paraSeps.push_back(para); paraSeps.push_back(para); paraSeps.push_back(para);
    variantSeps.push_back(para);
    std::cout << "Separator Variant Function\n";
    for(int i =0; i < variantSeps.size(); i++) {
        auto variantRet = variantSeps[i].getSignedDistanceAndGradAndHessianSep(x_eval,x_eval);
        std::cout << "F = \n" << std::get<0>(variantRet)<< "\n";
        std::cout << "Grad F = \n" << std::get<1>(variantRet) << "\n";
        std::cout << "Hess F = \n"<< std::get<2>(variantRet) << "\n";
    }   

    // Pick Centroids
    IRL::Pt p0(0,0,0);
    IRL::Pt p1(1,2,0);
    IRL::Pt p2(-1,-2,0);
    centroids.push_back(p0);
    centroids.push_back(p1);
    centroids.push_back(p2);

    // Construct Implicit Surfaces
    IRL::PUImplicitSurface planarSurface(centroids,variantSeps,delta);
    // IRL::ImplicitSurface<IRL::Paraboloid> paraSurface(centroids,paraSeps,delta);

    // Calculate Value at Point
    IRL::Pt x(0,1,1);
    std::cout << "Function Evaluation \n";
    planarSurface.F(x);

    std::cout << "X Der Function Evaluation \n";
    planarSurface.Fx(x);

    std::cout << "X Der Function Evaluation \n";
    planarSurface.Fy(x);

    std::cout << "======= Triple Function Evaluation \n";
    planarSurface.getValueAndGradAndHessian(x);
    // // PUST Solver Test
    // std::cout << "\n===== BEGIN SOLVER TESTS ==== " << std::endl;
    // IRL::Pt BL(-0.5,0.0,-0.5);
    // IRL::Pt BR(0.5,0.0,0.0);
    // IRL::Pt TL(-0.5,1.0,0.0);
    // IRL::Pt TR(0.5,1.0,0.5);
    // std::vector<IRL::Pt> poly = {BL,BR,TR,TL,BL};      
    // IRL::RectangularCuboid cell = IRL::RectangularCuboid::fromBoundingPts(BL, TR);
    // IRL::Plane plane = IRL::Plane(nor,0);
    // IRL::PlanarSeparator sep = IRL::PlanarSeparator::fromOnePlane(plane);
    // IRL::PUSTNeighborhood<IRL::RectangularCuboid,IRL::PlanarSeparator> testNeigh;
    // testNeigh.addMember(&cell,&sep);
    // testNeigh.setCenterOfStencil(0);
    
    // IRL::PUST<IRL::RectangularCuboid,IRL::PlanarSeparator> sol(testNeigh);
    // IRL::Normal solRes = sol.solve(1.0);
    // std::cout << "Force Result: " << solRes[0] << "," << solRes[1] << "," << solRes[2] << "\n";



    // std::cout << "\n======= BEGIN MAIN TEST ===== \n";
    // IRL::Pt P0,P1;
    // std::vector<IRL::Pt> inters;
    // IRL::Pt force(0.0,0.0,0.0); // Surface Tension Force
    // double STCoeff = 1;
    // IRL::Normal tangent(0.0,0.0,0.0);
    // for(int i = 0; i < poly.size()-1; i++) {
    //     P0 = poly[i];
    //     P1 = poly[i+1];
    //     std::cout << " P0 = " << P0[0] << "," << P0[1] << "," << P0[2] << "\n";
    //     std::cout << " P1 = " << P1[0] << "," << P1[1] << "," << P1[2] << "\n";
    //     inters = implicitTest.intersectEdge(P0,P1,10);

    //     std::cout << "In loop - " << inters.size() << " intersections \n";

    //     if(inters.size() > 0) {
    //         for(int j = 0; j < inters.size(); j++) {
    //             std::cout << inters[j][0] << "," << inters[j][1] << "," << inters[j][2] << "\n";
    //             double Fx = implicitTest.Fx(inters[j]);
    //             double Fy = implicitTest.Fy(inters[j]);
    //             std::cout << Fx << "," << Fy << " Gradient \n";
    //             tangent[0] = -Fy;
    //             tangent[1] = Fx;

    //             tangent.normalize();
    //             force = force + STCoeff * tangent;
    //         }
    //     }
    // }

    // std::cout << "Force = \n";
    // std::cout << force[0] << "," << force[1] << "," << force[2] << std::endl;




    return -1;
}
