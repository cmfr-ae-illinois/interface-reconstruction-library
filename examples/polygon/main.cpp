#include <chrono>
#include <iostream>
#include <string>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"

int main(int argc, char* argv[]) {
  // Defining bounds of cuboid cell
  const auto x0 = IRL::Pt(0.0, 0.0, 0.0);
  const auto x1 = IRL::Pt(1.0, 1.0, 1.0);

  // Construct cell object
  const auto cell = IRL::RectangularCuboid::fromBoundingPts(x0, x1);

  for (int i = 0; i < 110; i++) {
    // Define plane normal and distance
    auto normal = IRL::Normal(1.0, 2.0, 3.0);
    normal.normalize();
    const auto distance = static_cast<double>(i) * 0.01;

    // Create plane object
    const auto plane = IRL::Plane(normal, distance);

    // Compute polygon
    const auto polygon = IRL::cutPlaneByHexahedron<IRL::Polygon>(cell, plane);

    // Print surface area
    std::cout << std::scientific
              << "Surface area of polygon = " << std::setprecision(3)
              << polygon.calculateVolume()
              << " and distance = " << std::setprecision(3) << distance
              << std::endl;

    // Print volume fraction
    const double vfrac =
        IRL::getVolumeFraction(cell, IRL::PlanarSeparator::fromOnePlane(plane));
    std::cout << std::scientific << "Volume fraction = " << std::setprecision(3)
              << vfrac << " and distance = " << std::setprecision(3) << distance
              << std::endl;
  }

  return 0;
}
