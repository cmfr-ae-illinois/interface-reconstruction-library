#ifndef EXAMPLES_LEVEL_SET_RECONSTRUCTION_LEVEL_SET_H_
#define EXAMPLES_LEVEL_SET_RECONSTRUCTION_LEVEL_SET_H_

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"

namespace ExampleLevelSets {
using Definition = IRL::GeneralImplicitSurface<double, 5>;

struct Sphere : Definition {
  static constexpr double radius = 0.25;
  double F(const double& x, const double& y, const double& z) const override {
    return x * x + y * y + z * z - radius * radius;
  }
  Vec3 gradF(const double& x, const double& y, const double& z) const override {
    return Vec3(2 * x, 2 * y, 2 * z);
  }
  Mat3 hessF(const double&, const double&, const double&) const override {
    return 2.0 * Mat3::Identity();
  }
};

struct Ellipsoid : Definition {
  static constexpr double a = 0.30, b = 0.25, c = 0.20;
  double F(const double& x, const double& y, const double& z) const override {
    return x * x / (a * a) + y * y / (b * b) + z * z / (c * c) - 1.0;
  }
  Vec3 gradF(const double& x, const double& y, const double& z) const override {
    return Vec3(2 * x / (a * a), 2 * y / (b * b), 2 * z / (c * c));
  }
  Mat3 hessF(const double&, const double&, const double&) const override {
    Mat3 hessian = Mat3::Zero();
    hessian.diagonal() << 2 / (a * a), 2 / (b * b), 2 / (c * c);
    return hessian;
  }
  // Inherits the derivative-based, spatially varying reference curvature.
};

template <class Shape>
std::shared_ptr<const Definition> create() {
  return std::make_shared<Shape>();
}
struct Entry {
  const char* name;
  std::shared_ptr<const Definition> (*create)();
};
// Add a Definition subclass and one entry here to make a shape selectable.
inline const std::vector<Entry>& registry() {
  static const std::vector<Entry> entries = {{"sphere", &create<Sphere>},
                                             {"ellipsoid", &create<Ellipsoid>}};
  return entries;
}
}  // namespace ExampleLevelSets

// Copyable runtime-selected surface for ImplicitSurfaceCutter. The same
// selected definition supplies initialization, projection, and reference
// curvature.
struct LevelSet : IRL::GeneralImplicitSurface<double, 5> {
  explicit LevelSet(const std::string& name = "sphere") {
    for (const auto& entry : ExampleLevelSets::registry()) {
      if (name == entry.name) {
        definition_ = entry.create();
        return;
      }
    }
    throw std::invalid_argument("Unknown level set: " + name);
  }
  double F(const double& x, const double& y, const double& z) const override {
    return definition_->F(x, y, z);
  }
  Vec3 gradF(const double& x, const double& y, const double& z) const override {
    return definition_->gradF(x, y, z);
  }
  Mat3 hessF(const double& x, const double& y, const double& z) const override {
    return definition_->hessF(x, y, z);
  }

 private:
  std::shared_ptr<const ExampleLevelSets::Definition> definition_;
};
#endif
