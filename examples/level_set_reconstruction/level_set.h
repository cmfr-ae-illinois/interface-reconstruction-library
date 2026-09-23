#ifndef EXAMPLES_LEVEL_SET_RECONSTRUCTION_LEVEL_SET_H_
#define EXAMPLES_LEVEL_SET_RECONSTRUCTION_LEVEL_SET_H_

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"

// Edit this surface directly, including its gradient and Hessian.
// F < 0 is the liquid; F = 0 is the interface. The domain is [-0.5, 0.5]^3.
// Keep the interface away from the boundary, or use a periodic level set:
// variant_advector's reconstruction helpers assume periodic boundaries.
struct LevelSet : IRL::GeneralImplicitSurface<double, 5> {
  double F(const double& x, const double& y, const double& z) const override {
    constexpr double radius = 0.25;
    return x * x + y * y + z * z - radius * radius;
  }

  Vec3 gradF(const double& x, const double& y, const double& z) const override {
    return Vec3(2.0 * x, 2.0 * y, 2.0 * z);
  }

  Mat3 hessF(const double&, const double&, const double&) const override {
    return 2.0 * Mat3::Identity();
  }
};

#endif
