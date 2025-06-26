#include <chrono>
#include <iostream>
#include <string>

#include "examples/PUSurfaceTension/pu_neighborhood.h"
#include "examples/PUSurfaceTension/pu_solve.h"


#include "irl/geometry/polyhedrons/rectangular_cuboid.h"

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/moments/general_moments.h"

#include "irl/variant_reconstruction/separator_variant.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/normal.h"

int main(int argc, char* argv[]) {
    const auto xy_plane = IRL::Plane(IRL::Normal(0.0, 0.0, 1.0), 0.0);
    IRL::PUSTNeighborhood<IRL::RectangularCuboid> test;
    IRL::RectangularCuboid test2;
    const auto xy_plane_sep = IRL::PlanarSeparator::fromOnePlane(xy_plane);
    std::cout << "====== BEGIN NEIGHBORDHOOD TESTS ====\n";
    std::cout << "Test Initial Size: " << ((test.size() == 0) ? "Pass" : "Fail") << std::endl;
    
    test.addMember(&test2,&xy_plane_sep);
    std::cout << "After Adding Member: " << ((test.size() == 1) ? "Pass" : "Fail") << std::endl;

    test.emptyNeighborhood();
    std::cout << "After Empty: " << ((test.size() == 0) ? "Pass" : "Fail") << std::endl;

    test.resize(25);
    std::cout << "After Resize: " << ((test.size() == 25) ? "Pass" : "Fail") << std::endl;

    test.setCenterOfStencil(11);
    std::cout << "Set Center: " << ((test.getCenterOfStencilIndex() == 11) ? "Pass" : "Fail") << std::endl;

    test.setMember(1,&test2,&xy_plane_sep);
    test.getCenterCell();
    test.getCenterCellStoredMoments();
    test.getCell(1);
    test.getStoredMoments(0);
    test.begin();
    test.end();
    test.cbegin();
    test.cend();
    std::cout << "All other functions ran without problems :)\n";


    std::cout << "\n====== BEGIN WENDLAND TESTS ====\n";
    IRL::Normal xi; double delta; IRL::Normal x_eval;
    x_eval = IRL::Normal(0.0,1.0,1.0);
    IRL::Normal x_cen = IRL::Normal(0.0,1.0,0.0);
    delta = 50;
    std::cout << "Wendland Origin: " << (IRL::Wendland::phi(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland x  Origin: " << (IRL::Wendland::dphidx(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland y Origin: " << (IRL::Wendland::dphidy(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland xx Origin: " << (IRL::Wendland::ddphidxx(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland yy Origin: " << (IRL::Wendland::ddphidyy(x_eval,delta, x_eval)) << std::endl;
    std::cout << "Wendland xy Origin: " << (IRL::Wendland::ddphidxy(x_eval,delta, x_eval)) << std::endl;


    std::cout << "\n====== BEGIN IMPLICIT SURFACE TESTS ====\n";
    std::vector<IRL::Normal> centroids, normals;
    double kernel_size;
    
    centroids = {x_cen};
    IRL::Normal nor(1.0,1.0,0.0);
    nor.normalize();
    normals = {nor};
    IRL::ImplicitSurface implicitTest(centroids,normals,delta);
    std::cout << implicitTest.F(x_eval) << std::endl;
    std::cout << implicitTest.Fx(x_eval) << std::endl;
    std::cout << implicitTest.Fy(x_eval) << std::endl;
    std::vector<double> res = implicitTest.HessianTerms(x_eval);
    std::cout << res[0] << ","<< res[1] << ","<< res[2] << std::endl;

    IRL::Normal x0(-1.0,0.0,0.0);
    IRL::Normal x1(1.0,0.0,0.0);

    auto result = implicitTest.intersectEdge(x0,x1,10);
    std::cout << result.size() << " Intersections Found" << std::endl;
    std::cout << result[1][0] << ","<< result[1][1] << ","<< result[1][2] << std::endl;
    return -1;
}
