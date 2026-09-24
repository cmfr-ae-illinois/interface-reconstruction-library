#ifndef EXAMPLES_LEVEL_SET_RECONSTRUCTION_PU_FIELD_H_
#define EXAMPLES_LEVEL_SET_RECONSTRUCTION_PU_FIELD_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "examples/variant_advector/data.h"
#include "irl/interface_reconstruction_methods/pu.h"

namespace LevelSetVisualization {
using PU = IRL::PU<IRL::RectangularCuboid>;

inline double meanCurvature(const Eigen::Vector3d& gradient,
                            const Eigen::Matrix3d& hessian) {
  const double norm = gradient.norm();
  if (!gradient.allFinite() || !hessian.allFinite() || norm < 1.0e-10)
    return std::numeric_limits<double>::quiet_NaN();
  return (gradient.squaredNorm() * hessian.trace() -
          gradient.dot(hessian * gradient)) /
         (2.0 * norm * norm * norm);
}

struct Sample {
  double value = 0.0;
  double weight = 0.0;
  double curvature = 0.0;
  double curvature_error = 0.0;
  bool supported = false;
  bool curvature_valid = false;
};

class ReconstructedPU {
 public:
  ReconstructedPU(const Data<double>& fractions,
                  const Data<IRL::SeparatorVariant>& interfaces,
                  const double radius_cells)
      : fractions_(fractions),
        interfaces_(interfaces),
        mesh_(fractions.getMesh()),
        radius_(radius_cells * mesh_.dx()) {
    if (!std::isfinite(radius_cells) || radius_cells <= 0.0 ||
        radius_cells > mesh_.getNgc())
      throw std::invalid_argument("PU radius must be > 0 and <= ghost layers");
  }

  double radius() const { return radius_; }

  IRL::PUNeighborhood<IRL::RectangularCuboid> fittingNeighborhood(
      const IRL::Pt& seed) const {
    IRL::PUNeighborhood<IRL::RectangularCuboid> neighborhood;
    // solve() accepts projections at most dx/2 from the seed. Include every
    // center whose kernel can overlap that ball, not just the starting point.
    const double reach = radius_ + 0.5 * mesh_.dx();
    const int lo[3] = {
        int(std::floor((seed[0] - reach - mesh_.x(0)) / mesh_.dx())),
        int(std::floor((seed[1] - reach - mesh_.y(0)) / mesh_.dy())),
        int(std::floor((seed[2] - reach - mesh_.z(0)) / mesh_.dz()))};
    const int hi[3] = {
        int(std::floor((seed[0] + reach - mesh_.x(0)) / mesh_.dx())),
        int(std::floor((seed[1] + reach - mesh_.y(0)) / mesh_.dy())),
        int(std::floor((seed[2] + reach - mesh_.z(0)) / mesh_.dz()))};
    const auto wrap = [](int i, int n) { return (i % n + n) % n; };
    for (int i = lo[0]; i <= hi[0]; ++i)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int k = lo[2]; k <= hi[2]; ++k) {
          const int ii = wrap(i, mesh_.getNx()), jj = wrap(j, mesh_.getNy()),
                    kk = wrap(k, mesh_.getNz());
          const double vf = fractions_(ii, jj, kk);
          if (vf < IRL::global_constants::VF_LOW ||
              vf > IRL::global_constants::VF_HIGH)
            continue;
          const IRL::Pt shift((i - ii) * mesh_.dx(), (j - jj) * mesh_.dy(),
                              (k - kk) * mesh_.dz());
          const IRL::Pt center =
              IRL::Pt(mesh_.xm(ii), mesh_.ym(jj), mesh_.zm(kk)) + shift;
          if (IRL::magnitude(center - seed) >= reach) continue;
          auto interface = interfaces_(ii, jj, kk);
          if (auto* plane = std::get_if<IRL::PlanarSeparator>(&interface)) {
            for (auto& p : *plane) p.distance() += p.normal() * shift;
          } else if (auto* paraboloid =
                         std::get_if<IRL::Paraboloid>(&interface)) {
            paraboloid->setDatum(paraboloid->getDatum() + shift);
          }
          neighborhood.addMember(&center, &interface);
        }
    neighborhood.setCenterOfStencil(0);  // Explicit-seed solve does not use it.
    return neighborhood;
  }

  Sample evaluate(const IRL::Pt& point) const {
    IRL::PUNeighborhood<IRL::RectangularCuboid> neighborhood;
    Sample sample;
    // Cell centers are at index + 1/2. Extra end cells are harmless: exact
    // radial support is checked below before adding a member.
    const int lo[3] = {
        std::max(
            mesh_.imino(),
            int(std::floor((point[0] - radius_ - mesh_.x(0)) / mesh_.dx()))),
        std::max(
            mesh_.jmino(),
            int(std::floor((point[1] - radius_ - mesh_.y(0)) / mesh_.dy()))),
        std::max(
            mesh_.kmino(),
            int(std::floor((point[2] - radius_ - mesh_.z(0)) / mesh_.dz())))};
    const int hi[3] = {
        std::min(
            mesh_.imaxo(),
            int(std::floor((point[0] + radius_ - mesh_.x(0)) / mesh_.dx()))),
        std::min(
            mesh_.jmaxo(),
            int(std::floor((point[1] + radius_ - mesh_.y(0)) / mesh_.dy()))),
        std::min(
            mesh_.kmaxo(),
            int(std::floor((point[2] + radius_ - mesh_.z(0)) / mesh_.dz())))};
    for (int i = lo[0]; i <= hi[0]; ++i)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int k = lo[2]; k <= hi[2]; ++k) {
          const double vf = fractions_(i, j, k);
          if (vf < IRL::global_constants::VF_LOW ||
              vf > IRL::global_constants::VF_HIGH)
            continue;
          const IRL::Pt center(mesh_.xm(i), mesh_.ym(j), mesh_.zm(k));
          double weight;
          IRL::Wendland::evaluate(center, radius_, point, &weight);
          if (weight <= 0.0) continue;
          neighborhood.addMember(&center, &interfaces_(i, j, k));
          sample.weight += weight;
        }
    // getPU() alone returns zero when no support exists. Never interpret
    // that as a surface. Export only cells with supported corners.
    if (sample.weight <= 1.0e-12) return sample;
    PU pu(neighborhood, radius_);
    const auto result = pu.getPUGradAndHess(point);
    sample.value = std::get<0>(result);
    sample.supported = std::isfinite(sample.value);
    const double curvature =
        meanCurvature(std::get<1>(result), std::get<2>(result));
    sample.curvature_valid = sample.supported && std::isfinite(curvature);
    if (sample.curvature_valid) sample.curvature = curvature;
    return sample;
  }

 private:
  const Data<double>& fractions_;
  const Data<IRL::SeparatorVariant>& interfaces_;
  const BasicMesh& mesh_;
  double radius_;
};
}  // namespace LevelSetVisualization
#endif
