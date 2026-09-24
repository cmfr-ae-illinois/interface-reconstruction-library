#ifndef EXAMPLES_LEVEL_SET_RECONSTRUCTION_PU_PPIC_H_
#define EXAMPLES_LEVEL_SET_RECONSTRUCTION_PU_PPIC_H_

#include "examples/level_set_reconstruction/curvature.h"
#include "irl/generic_cutting/cut_polygon.h"
#include "irl/interface_reconstruction_methods/pu_paraboloid.h"
#include "irl/interface_reconstruction_methods/volume_fraction_matching.h"

namespace LevelSetVisualization {
inline IRL::Pt clippedPolygonCentroid(const IRL::RectangularCuboid& cell,
                                      const IRL::PlanarSeparator& interface) {
  if (interface.getNumberOfPlanes() != 1)
    throw std::runtime_error("Expected one LVIRA plane for polygon centroid");
  const auto polygon = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
      cell, interface, interface[0]);
  if (polygon.getNumberOfVertices() < 3 || polygon.calculateVolume() <= 0.0)
    throw std::runtime_error("Empty clipped LVIRA polygon in mixed cell");
  const auto centroid = polygon.calculateCentroid();
  for (int d = 0; d < 3; ++d)
    if (!std::isfinite(centroid[d]))
      throw std::runtime_error("Nonfinite LVIRA polygon centroid");
  return centroid;
}

inline IRL::Paraboloid fitPUParaboloid(const ReconstructedPU& pu,
                                       const IRL::Pt& seed,
                                       const IRL::RectangularCuboid& cell,
                                       const double fraction, const double dx) {
  const auto neighborhood = pu.fittingNeighborhood(seed);
  IRL::PUParaboloid<IRL::RectangularCuboid> solver(neighborhood, pu.radius(),
                                                   dx);
  auto paraboloid = solver.solve(seed);
  const auto& datum = paraboloid.getDatum();
  if (!std::isfinite(datum[0]) || !std::isfinite(datum[1]) ||
      !std::isfinite(datum[2]) ||
      !std::isfinite(paraboloid.getAlignedParaboloid().a()) ||
      !std::isfinite(paraboloid.getAlignedParaboloid().b()))
    throw std::runtime_error(
        "PUParaboloid::solve failed; no fallback interface was substituted");

  // Match only the position of the fitted paraboloid to the original fraction.
  IRL::SeparatorVariant matched = paraboloid;
  IRL::setDistanceToMatchVolumeFraction(cell, fraction, &matched, 1.0e-12);
  const double actual =
      static_cast<double>(IRL::getVolumeMoments<IRL::Volume>(cell, matched)) /
      cell.calculateVolume();
  if (!std::isfinite(actual) || std::abs(actual - fraction) > 1.0e-10)
    throw std::runtime_error("PU paraboloid volume-fraction matching failed");
  return std::get<IRL::Paraboloid>(matched);
}
}  // namespace LevelSetVisualization
#endif
