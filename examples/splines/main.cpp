#include <iostream>
#include "Spline.h"

int main() {
    std::vector<std::vector<double>> V1 = {{1,2},{0.7239,3.3806},{2.1078,2},{3.6852,0.4262},{4,2}};
    std::vector<double> KV = {0,0,0,0.6667,0.6667,1,1,1};
    std::vector<double> W = {1,0.3361,1,0.5025,1};

    std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    std::vector<std::vector<double>> InterpTester = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{-3.0,0.0},{1.0,0.0}};
    std::vector<std::vector<double>> InterpTesterT = {{0.0,1.0},{-1.0,0.0},{-1.0,0.0},{-1.0,0.0},{0.0,1.0}};
    
    Spline spline = Spline(V1,KV,W);
    Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    
    std::cout << "================== Being Circle test ===========================\n";
    std::cout << "Total Area = " <<circle.getArea() << "\n";
    std::cout << "Total Arc Length = " << circle.getArcLength() << "\n";
    std::cout << "Surface Energy = " << circle.getSurfaceEnergy() << "\n";
    std::vector<std::vector<double>> square = {{{0,0},{1.5,0},{1.5,1.5},{0,1.5},{0,0}}};
    double A = circle.integrateSplineSquare(square);
    std::cout << "Clipped Area = " << A << "\n";

    circle.saveToVTK("circle_test");
    // Circle Tests Working and Matching


    std::cout << "\n================== Being Weird shape test ===========================\n";
    Spline InterpTest = Spline::LocalRQuadInterp(InterpTester,InterpTesterT);

    InterpTest.saveToVTK("InterpTester_test");
    std::cout << "Total Area = " <<InterpTest.getArea() << "\n";
    std::cout << "Total Arc Length = " << InterpTest.getArcLength() << "\n";
    std::cout << "Surface Energy = " << InterpTest.getSurfaceEnergy() << "\n";
    square = {{-4,-2},{-2,-2},{-2,0.5},{-4,0.5},{-4,-2}};
    A = InterpTest.integrateSplineSquare(square);
    std::cout << "Clipped Area = " << A << "\n";

    
    

    

    return 0;
}
