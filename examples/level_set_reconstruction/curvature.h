#ifndef EXAMPLES_LEVEL_SET_RECONSTRUCTION_CURVATURE_H_
#define EXAMPLES_LEVEL_SET_RECONSTRUCTION_CURVATURE_H_

#include "examples/level_set_reconstruction/pu_field.h"

namespace LevelSetVisualization {

template <class Evaluate>
bool projectToZero(IRL::Pt* point, const double dx, Evaluate evaluate) {
  for (int iteration = 0; iteration < 100; ++iteration) {
    const auto result = evaluate(*point);
    const double value = result.first;
    const auto& gradient = result.second;
    const double norm = gradient.norm();
    if (!std::isfinite(value) || !gradient.allFinite() || norm < 1.0e-12)
      return false;
    if (std::abs(value) / norm < 1.0e-10 * dx) return true;
    for (int d = 0; d < 3; ++d)
      (*point)[d] -= value * gradient[d] / (norm * norm);
  }
  return false;
}

inline IRL::Pt interfacePoint(const IRL::SeparatorVariant& interface,
                              const IRL::Pt& center, const double dx) {
  IRL::Pt point = center;
  if (!projectToZero(&point, dx, [&](const IRL::Pt& p) {
        return PU::implicitSeparatorValueandGrad(p, center, &interface);
      }))
    throw std::runtime_error("Could not project onto reconstructed interface");
  return point;
}

template <class Surface>
double referenceMeanCurvature(const Surface& surface, IRL::Pt point,
                              const double dx) {
  if (!projectToZero(&point, dx, [&](const IRL::Pt& p) {
        return std::make_pair(surface.F(p[0], p[1], p[2]),
                              surface.gradF(p[0], p[1], p[2]));
      }))
    return std::numeric_limits<double>::quiet_NaN();
  return surface.meanCurvature(point[0], point[1], point[2]);
}
}  // namespace LevelSetVisualization
#endif
