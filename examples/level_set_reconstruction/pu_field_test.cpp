#include <cmath>
#include <iostream>
#include <stdexcept>

#include "examples/level_set_reconstruction/level_set.h"
#include "examples/level_set_reconstruction/pu_ppic.h"

using LevelSetVisualization::PU;
using LevelSetVisualization::ReconstructedPU;

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, double tolerance,
           const char* message) {
  require(std::isfinite(actual) && std::abs(actual - expected) < tolerance,
          message);
}

int main() {
  try {
    BasicMesh mesh(8, 8, 8, 5);
    mesh.setCellBoundaries(IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, 0.5));
    Data<double> fractions(&mesh);
    Data<IRL::SeparatorVariant> interfaces(&mesh);
    const auto plane = IRL::PlanarSeparator::fromOnePlane(
        IRL::Plane(IRL::Normal(0, 0, 1), 0.03));
    const IRL::Paraboloid paraboloid(
        IRL::Pt(0, 0, 0), IRL::ReferenceFrame::fromNormal(IRL::Normal(0, 0, 1)),
        0.5, 0.5);
    const auto fill = [&](const IRL::SeparatorVariant& interface) {
      for (int i = mesh.imino(); i <= mesh.imaxo(); ++i)
        for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j)
          for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
            fractions(i, j, k) = 0.5;
            interfaces(i, j, k) = interface;
          }
    };
    fill(plane);
    ReconstructedPU field(fractions, interfaces, 2.5);
    for (const auto& p : {IRL::Pt(0.1, 0.02, 0.03), IRL::Pt(-0.5, 0.2, 0.1),
                          IRL::Pt(0.5, 0.2, 0.1), IRL::Pt(0, 0, 0)}) {
      const auto s = field.evaluate(p);
      require(s.supported && s.curvature_valid,
              "Plane support or curvature invalid");
      close(s.value, p[2] - 0.03, 1e-12, "Plane not reproduced");
      close(s.curvature, 0, 1e-10, "Plane curvature not zero");
    }
    const auto target_cell = IRL::RectangularCuboid::fromBoundingPts(
        IRL::Pt(-0.0625, -0.0625, -0.0625), IRL::Pt(0.0625, 0.0625, 0.0625));
    const auto planar_fit = LevelSetVisualization::fitPUParaboloid(
        field,
        LevelSetVisualization::clippedPolygonCentroid(target_cell, plane),
        target_cell, 0.37, mesh.dx());
    close(static_cast<double>(
              IRL::getVolumeMoments<IRL::Volume>(target_cell, planar_fit)) /
              target_cell.calculateVolume(),
          0.37, 1e-10, "Plane PU PPIC volume mismatch");
    close(planar_fit.getAlignedParaboloid().a(), 0, 1e-10,
          "Plane fit unexpectedly curved");
    close(planar_fit.getAlignedParaboloid().b(), 0, 1e-10,
          "Plane fit unexpectedly curved");
    const auto unit_cell = IRL::RectangularCuboid::fromBoundingPts(
        IRL::Pt(0, 0, 0), IRL::Pt(1, 1, 1));
    const double normalization = std::sqrt(14.0);
    const auto oblique = IRL::PlanarSeparator::fromOnePlane(IRL::Plane(
        IRL::Normal(1 / normalization, 2 / normalization, 3 / normalization),
        0.6 / normalization));
    const auto clipped_centroid =
        LevelSetVisualization::clippedPolygonCentroid(unit_cell, oblique);
    close(clipped_centroid[0], 0.2, 1e-12, "Wrong clipped triangle centroid x");
    close(clipped_centroid[1], 0.1, 1e-12, "Wrong clipped triangle centroid y");
    close(clipped_centroid[2], 0.2 / 3, 1e-12,
          "Wrong clipped triangle centroid z");
    const auto projected_center = LevelSetVisualization::interfacePoint(
        oblique, unit_cell.calculateCentroid(), 1.0);
    require(IRL::magnitude(clipped_centroid - projected_center) > 0.1,
            "Polygon centroid must not be replaced by projected cell center");
    // Default solve() remains available; explicit seed does not change weights.
    auto plane_neighborhood = field.fittingNeighborhood(IRL::Pt(0, 0, 0.03));
    plane_neighborhood.setCenterOfStencil(0);
    IRL::PUParaboloid<IRL::RectangularCuboid> seeded_solver(
        plane_neighborhood, field.radius(), mesh.dx());
    const auto seeded = seeded_solver.solve(IRL::Pt(0, 0, 0.03));
    close(seeded.getDatum()[2], 0.03, 1e-10, "Explicit projection seed failed");
    IRL::PUNeighborhood<IRL::RectangularCuboid> one;
    const IRL::Pt plane_seed(0, 0, 0.03);
    const IRL::SeparatorVariant plane_variant = plane;
    one.addMember(&plane_seed, &plane_variant);
    one.setCenterOfStencil(0);
    IRL::PUParaboloid<IRL::RectangularCuboid> default_solver(
        one, field.radius(), mesh.dx());
    close(default_solver.solve().getDatum()[2], 0.03, 1e-10,
          "Default solve regressed");
    fill(paraboloid);
    const auto apex = field.evaluate(IRL::Pt(0, 0, 0));
    close(apex.value, 0, 1e-12, "Paraboloid zero surface not reproduced");
    close(apex.curvature, 1, 1e-10, "Paraboloid apex mean curvature incorrect");

    const auto curved_fit = LevelSetVisualization::fitPUParaboloid(
        field,
        LevelSetVisualization::clippedPolygonCentroid(target_cell, plane),
        target_cell, 0.63, mesh.dx());
    close(static_cast<double>(
              IRL::getVolumeMoments<IRL::Volume>(target_cell, curved_fit)) /
              target_cell.calculateVolume(),
          0.63, 1e-10, "Curved PU PPIC volume mismatch");
    close(LevelSetVisualization::referenceMeanCurvature(
              LevelSet(), IRL::Pt(0.3, 0.1, 0.1), mesh.dx()),
          4.0, 1e-8,
          "Sphere reference curvature must be evaluated on its zero set");

    const LevelSet ellipsoid("ellipsoid");
    using Ellipsoid = ExampleLevelSets::Ellipsoid;
    close(ellipsoid.F(Ellipsoid::a, 0, 0), 0, 1e-12,
          "Selected ellipsoid function incorrect");
    close(LevelSetVisualization::referenceMeanCurvature(
              ellipsoid, IRL::Pt(0.4, 0, 0), mesh.dx()),
          0.5 * Ellipsoid::a *
              (1 / (Ellipsoid::b * Ellipsoid::b) +
               1 / (Ellipsoid::c * Ellipsoid::c)),
          1e-8, "Ellipsoid x-tip curvature incorrect");
    close(LevelSetVisualization::referenceMeanCurvature(
              ellipsoid, IRL::Pt(0, 0, 0.3), mesh.dx()),
          0.5 * Ellipsoid::c *
              (1 / (Ellipsoid::a * Ellipsoid::a) +
               1 / (Ellipsoid::b * Ellipsoid::b)),
          1e-8, "Ellipsoid z-tip curvature must vary spatially");
    bool unknown_rejected = false;
    try {
      LevelSet unknown("not-a-shape");
    } catch (const std::invalid_argument&) {
      unknown_rejected = true;
    }
    require(unknown_rejected, "Unknown level set silently accepted");

    // Mix planes and paraboloids. Compare local search with the complete
    // neighborhood, including near mesh/stencil boundaries.
    IRL::PUNeighborhood<IRL::RectangularCuboid> all;
    for (int i = mesh.imino(); i <= mesh.imaxo(); ++i)
      for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j)
        for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
          if ((i + j + k) % 2 == 0) interfaces(i, j, k) = plane;
          const IRL::Pt center(mesh.xm(i), mesh.ym(j), mesh.zm(k));
          all.addMember(&center, &interfaces(i, j, k));
        }
    PU exhaustive(all, 2.5 * mesh.dx());
    for (const auto& p : {IRL::Pt(0.12, 0.04, 0.01), IRL::Pt(-0.5, 0.04, 0.01),
                          IRL::Pt(0.125 - 1e-8, 0.04, 0.01),
                          IRL::Pt(0.125 + 1e-8, 0.04, 0.01)}) {
      const auto s = field.evaluate(p);
      close(s.value, exhaustive.getPU(p), 1e-12,
            "Support search omitted a contributor");
      close(s.curvature, exhaustive.getMeanCurvature(p), 1e-9,
            "Curvature differs from complete PU");
    }
    ReconstructedPU wide_field(fractions, interfaces, 5.0);
    PU wide_exhaustive(all, 5.0 * mesh.dx());
    for (const auto& p : {IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, 0.5),
                          IRL::Pt(0.12, 0.04, 0.01)}) {
      const auto sample = wide_field.evaluate(p);
      require(sample.supported && sample.curvature_valid,
              "Radius-five sample invalid");
      close(sample.value, wide_exhaustive.getPU(p), 1e-12,
            "Radius-five search omitted support at domain boundary");
      close(sample.curvature, wide_exhaustive.getMeanCurvature(p), 1e-9,
            "Radius-five curvature differs from exhaustive neighborhood");
    }
    // Independent derivative check from field values (not PU derivatives).
    const IRL::Pt p(0.12, 0.04, 0.01);
    constexpr double h = 1e-5;
    Eigen::Vector3d gradient;
    Eigen::Matrix3d hessian;
    const double center = field.evaluate(p).value;
    for (int a = 0; a < 3; ++a) {
      auto plus = p, minus = p;
      plus[a] += h;
      minus[a] -= h;
      const double fp = field.evaluate(plus).value,
                   fm = field.evaluate(minus).value;
      gradient[a] = (fp - fm) / (2 * h);
      hessian(a, a) = (fp - 2 * center + fm) / (h * h);
      for (int b = 0; b < a; ++b) {
        auto pp = p, pm = p, mp = p, mm = p;
        pp[a] += h;
        pp[b] += h;
        pm[a] += h;
        pm[b] -= h;
        mp[a] -= h;
        mp[b] += h;
        mm[a] -= h;
        mm[b] -= h;
        hessian(a, b) = hessian(b, a) =
            (field.evaluate(pp).value - field.evaluate(pm).value -
             field.evaluate(mp).value + field.evaluate(mm).value) /
            (4 * h * h);
      }
    }
    close(field.evaluate(p).curvature,
          LevelSetVisualization::meanCurvature(gradient, hessian), 1e-5,
          "PU curvature disagrees with finite differences");

    for (int i = mesh.imino(); i <= mesh.imaxo(); ++i)
      for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j)
        for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k)
          fractions(i, j, k) = 0;
    require(!field.evaluate(p).supported,
            "Empty support reported as zero surface");
    fractions(4, 4, 4) = 0.5;
    require(!field.evaluate(IRL::Pt(-0.5, -0.5, -0.5)).supported,
            "Outside kernel reported as zero surface");
    bool rejected = false;
    try {
      ReconstructedPU invalid(fractions, interfaces, 5.1);
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    require(rejected, "Radius exceeding ghost support accepted");
    std::cout << "Plane, paraboloid, mixed-support, derivative, PPIC matching, "
                 "reference curvature, and "
                 "unsupported-region tests passed.\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
