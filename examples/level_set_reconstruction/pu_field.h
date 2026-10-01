#ifndef EXAMPLES_LEVEL_SET_RECONSTRUCTION_PU_FIELD_H_
#define EXAMPLES_LEVEL_SET_RECONSTRUCTION_PU_FIELD_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "examples/level_set_reconstruction/curvature.h"
#include "examples/variant_advector/data.h"
#include "irl/interface_reconstruction_methods/pu.h"
namespace LevelSetVisualization {

template <class WeightFunction>
using PU = IRL::PU<IRL::RectangularCuboid, WeightFunction>;

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

struct InterfaceErrors {
  std::size_t mixed_cells = 0;

  double position_error = 0.0;
  double normal_error = 0.0;
  double curvature_error = 0.0;
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

  template <class WeightFunction>
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
          WeightFunction::evaluate(center, radius_, point, &weight);
          if (weight <= 0.0) continue;
          neighborhood.addMember(&center, &interfaces_(i, j, k));
          sample.weight += weight;
        }
    // getPU() alone returns zero when no support exists. Never interpret
    // that as a surface. Export only cells with supported corners.
    if (sample.weight <= 1.0e-12) return sample;
    PU<WeightFunction> pu(neighborhood, radius_);
    const auto result = pu.getPUGradAndHess(point);
    sample.value = std::get<0>(result);
    sample.supported = std::isfinite(sample.value);
    const double curvature =
        meanCurvature(std::get<1>(result), std::get<2>(result));
    sample.curvature_valid = sample.supported && std::isfinite(curvature);
    if (sample.curvature_valid) sample.curvature = curvature;
    return sample;
  }

  template <class WeightFunction>
  InterfaceErrors computeInterfaceErrors(const LevelSet& reference) {
    InterfaceErrors errors;
    IRL::PUNeighborhood<IRL::RectangularCuboid> neighborhood;

    // Loop over computational cells
    for (int k = 0; k < mesh_.getNz(); ++k) {
      for (int j = 0; j < mesh_.getNy(); ++j) {
        for (int i = 0; i < mesh_.getNx(); ++i) {
          // Determine volume fraction
          const double vf = fractions_(i, j, k);
          if (vf > 1e-12 && vf < 1.0 - 1e-12) {
            // Mixed cell
            ++errors.mixed_cells;

            // Now that we are in a mixed cell, find the center of the cell
            const IRL::Pt cell_center(mesh_.xm(i), mesh_.ym(j), mesh_.zm(k));
            // Project onto Partition of unity
            // Create Neighborhood
            neighborhood.emptyNeighborhood();
            for (int ii = -3; ii <= 3; ++ii) {
              for (int jj = -3; jj <= 3; ++jj) {
                for (int kk = -3; kk <= 3; ++kk) {
                  const int ii_global = i + ii;
                  const int jj_global = j + jj;
                  const int kk_global = k + kk;
                  if (ii_global < 0 || ii_global >= mesh_.getNx() ||
                      jj_global < 0 || jj_global >= mesh_.getNy() ||
                      kk_global < 0 || kk_global >= mesh_.getNz())
                    continue;
                  const IRL::Pt center(mesh_.xm(ii_global), mesh_.ym(jj_global),
                                       mesh_.zm(kk_global));
                  double vf_neighbor =
                      fractions_(ii_global, jj_global, kk_global);
                  if (vf_neighbor > 1e-12 && vf_neighbor < 1.0 - 1e-12) {
                    neighborhood.addMember(
                        &center, &interfaces_(ii_global, jj_global, kk_global));
                  }
                }
              }
            }
            // Create PU object
            PU<WeightFunction> pu(neighborhood, radius_);
            // Project onto PU
            IRL::Pt projected_point = cell_center;
            if (!projectToZero(
                    &projected_point, mesh_.dx(),
                    [&](const IRL::Pt& p) { return pu.getPUAndGrad(p); })) {
              std::cout << "Warning: Could not project onto PU for cell (" << i
                        << ", " << j << ", " << k << ")\n";
              continue;
            }
            // Project Onto reference level set
            IRL::Pt reference_projected_point = cell_center;
            if (!projectToZero(&reference_projected_point, mesh_.dx(),
                               [&](const IRL::Pt& p) {
                                 return std::make_pair(
                                     reference.F(p[0], p[1], p[2]),
                                     reference.gradF(p[0], p[1], p[2]));
                               })) {
              std::cout
                  << "Warning: Could not project onto reference level set "
                     "for cell ("
                  << i << ", " << j << ", " << k << ")\n";
              continue;
            }

            // At the projected point, get the normal and mean curvature from
            // the PU
            IRL::Normal pu_normal = pu.getNormal(projected_point);
            double pu_mean_curvature = pu.getMeanCurvature(projected_point);
            // Reference
            auto reference_normal_temp = reference.gradF(
                reference_projected_point[0], reference_projected_point[1],
                reference_projected_point[2]);
            auto reference_hessian = reference.hessF(
                reference_projected_point[0], reference_projected_point[1],
                reference_projected_point[2]);
            // Cast Normal to Eigen::Vector3d and IRL:Normal
            Eigen::Vector3d reference_gradient_eigen(reference_normal_temp[0],
                                                     reference_normal_temp[1],
                                                     reference_normal_temp[2]);
            IRL::Normal reference_normal_irl = {reference_normal_temp[0],
                                                reference_normal_temp[1],
                                                reference_normal_temp[2]};
            reference_normal_irl.normalize();
            // Get Reference Mean Curvature
            double reference_mean_curvature =
                LevelSetVisualization::referenceMeanCurvature(
                    reference, reference_projected_point, mesh_.dx());
            // Compute Errors
            errors.mixed_cells++;
            // Position Error
            IRL::Pt dPOS = projected_point - reference_projected_point;
            const double distance = IRL::magnitude(dPOS);
            errors.position_error += distance * distance;
            // Normal Error - Dot porudct to get angle between normals
            pu_normal.normalize();
            reference_normal_irl.normalize();
            double DP = IRL::dotProduct(pu_normal, reference_normal_irl);
            double angle = std::acos(DP);
            errors.normal_error += angle * angle;
            // Curvature Error
            double dCURV = (std::abs(pu_mean_curvature) -
                            std::abs(reference_mean_curvature)) /
                           std::abs(reference_mean_curvature);
            errors.curvature_error += dCURV * dCURV;
          }
        }
      }
    }
    // return
    errors.position_error =
        std::sqrt(errors.position_error / errors.mixed_cells);
    errors.normal_error = std::sqrt(errors.normal_error / errors.mixed_cells);
    errors.curvature_error =
        std::sqrt(errors.curvature_error / errors.mixed_cells);
    return errors;
  }

 private:
  const Data<double>& fractions_;
  const Data<IRL::SeparatorVariant>& interfaces_;
  const BasicMesh& mesh_;
  double radius_;
};
}  // namespace LevelSetVisualization
#endif
