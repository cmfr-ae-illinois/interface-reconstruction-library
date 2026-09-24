#include "examples/level_set_reconstruction/visualization.h"
#include "examples/level_set_reconstruction/pu_ppic.h"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace LevelSetVisualization {
void storePolygonCentroids(const Data<double>& fractions,
                           const Data<IRL::SeparatorVariant>& planes,
                           Data<IRL::Pt>* centroids) {
  const auto& mesh = fractions.getMesh();
  for (int i = mesh.imin(); i <= mesh.imax(); ++i)
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j)
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double vf = fractions(i, j, k);
        if (vf < IRL::global_constants::VF_LOW ||
            vf > IRL::global_constants::VF_HIGH)
          continue;
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        (*centroids)(i, j, k) = clippedPolygonCentroid(
            cell, std::get<IRL::PlanarSeparator>(planes(i, j, k)));
      }
}

void addInterfaceDiagnostics(const Data<double>& fractions,
                             const Data<IRL::SeparatorVariant>& interfaces,
                             const LevelSet& reference,
                             const std::string& label,
                             std::vector<InterfaceScalarField>* fields) {
  const auto& mesh = fractions.getMesh();
  fields->clear();
  fields->emplace_back("mean_curvature", &mesh);
  fields->emplace_back("mean_curvature_error", &mesh);
  double sum_squared = 0.0, sum_absolute = 0.0, max_absolute = 0.0;
  std::size_t count = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i)
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j)
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double vf = fractions(i, j, k);
        if (vf < IRL::global_constants::VF_LOW ||
            vf > IRL::global_constants::VF_HIGH)
          continue;
        const auto& interface = interfaces(i, j, k);
        const auto* plane = std::get_if<IRL::PlanarSeparator>(&interface);
        if ((!plane && !std::holds_alternative<IRL::Paraboloid>(interface)) ||
            (plane && plane->getNumberOfPlanes() != 1))
          throw std::runtime_error(
              "PU visualization requires one plane or paraboloid per mixed "
              "cell");
        const IRL::Pt center(mesh.xm(i), mesh.ym(j), mesh.zm(k));
        const auto point = interfacePoint(interface, center, mesh.dx());
        const auto derivatives =
            PU::implicitSeparatorValueGradHess(point, center, &interface);
        const double curvature =
            meanCurvature(std::get<1>(derivatives), std::get<2>(derivatives));
        const double expected =
            referenceMeanCurvature(reference, point, mesh.dx());
        if (!std::isfinite(curvature) || !std::isfinite(expected))
          throw std::runtime_error(
              "Undefined interface/reference curvature in " + label);
        const double error = curvature - expected;
        const double values[] = {curvature, error};
        for (std::size_t n = 0; n < 2; ++n) {
          (*fields)[n].polygon_scalar_data(i, j, k) = values[n];
          (*fields)[n].paraboloid_scalar_data(i, j, k) = values[n];
        }
        ++count;
        sum_squared += error * error;
        sum_absolute += std::abs(error);
        max_absolute = std::max(max_absolute, std::abs(error));
      }
  if (count)
    std::cout << label
              << " mean-curvature error (one sample per mixed cell): MAE="
              << sum_absolute / count
              << ", RMSE=" << std::sqrt(sum_squared / count)
              << ", max=" << max_absolute << '\n';
}

void reconstructPUPPIC(const ReconstructedPU& pu, const Data<double>& fractions,
                       const Data<IRL::SeparatorVariant>& source,
                       const Data<IRL::Pt>& representative_points,
                       Data<IRL::SeparatorVariant>* fitted) {
  const auto& mesh = fractions.getMesh();
  std::size_t count = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i)
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j)
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        (*fitted)(i, j, k) = source(i, j, k);
        const double vf = fractions(i, j, k);
        if (vf < IRL::global_constants::VF_LOW ||
            vf > IRL::global_constants::VF_HIGH)
          continue;
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        try {
          (*fitted)(i, j, k) = fitPUParaboloid(
              pu, representative_points(i, j, k), cell, vf, mesh.dx());
        } catch (const std::exception& error) {
          throw std::runtime_error("PU PPIC cell " + std::to_string(i) + "," +
                                   std::to_string(j) + "," + std::to_string(k) +
                                   ": " + error.what());
        }
        ++count;
      }
  std::cout << "PU PPIC: " << count
            << " paraboloids; every matched fraction checked to 1e-10.\n";
}

void writePUField(const ReconstructedPU& pu, const BasicMesh& mesh,
                  const LevelSet& reference, const int sample_nx,
                  const std::string& filename) {
  const std::size_t side = std::size_t(sample_nx) + 1;
  const auto index = [side](int i, int j, int k) {
    return (std::size_t(k) * side + j) * side + i;
  };
  const auto point = [&](int i, int j, int k) {
    return IRL::Pt(mesh.x(0) + mesh.lx() * i / sample_nx,
                   mesh.y(0) + mesh.ly() * j / sample_nx,
                   mesh.z(0) + mesh.lz() * k / sample_nx);
  };
  std::vector<Sample> samples(side * side * side);
  for (int k = 0; k <= sample_nx; ++k)
    for (int j = 0; j <= sample_nx; ++j)
      for (int i = 0; i <= sample_nx; ++i) {
        auto& sample = samples[index(i, j, k)];
        const auto location = point(i, j, k);
        sample = pu.evaluate(location);
        if (sample.curvature_valid) {
          const double expected =
              referenceMeanCurvature(reference, location, mesh.dx());
          sample.curvature_valid = std::isfinite(expected);
          if (sample.curvature_valid)
            sample.curvature_error = sample.curvature - expected;
        }
      }

  std::vector<std::array<std::size_t, 8>> cells;
  for (int k = 0; k < sample_nx; ++k)
    for (int j = 0; j < sample_nx; ++j)
      for (int i = 0; i < sample_nx; ++i) {
        const std::array<std::size_t, 8> corners = {index(i, j, k),
                                                    index(i + 1, j, k),
                                                    index(i + 1, j + 1, k),
                                                    index(i, j + 1, k),
                                                    index(i, j, k + 1),
                                                    index(i + 1, j, k + 1),
                                                    index(i + 1, j + 1, k + 1),
                                                    index(i, j + 1, k + 1)};
        if (std::all_of(corners.begin(), corners.end(), [&](std::size_t n) {
              return samples[n].supported && samples[n].curvature_valid;
            }))
          cells.push_back(corners);
      }
  std::ofstream out(filename);
  out.exceptions(std::ios::badbit | std::ios::failbit);
  out << std::setprecision(17);
  out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" "
         "version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n<Piece "
         "NumberOfPoints=\""
      << samples.size() << "\" NumberOfCells=\"" << cells.size() << "\">\n";
  out << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"ascii\">\n";
  for (int k = 0; k <= sample_nx; ++k)
    for (int j = 0; j <= sample_nx; ++j)
      for (int i = 0; i <= sample_nx; ++i) {
        const auto p = point(i, j, k);
        out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
      }
  out << "</DataArray></Points>\n<Cells>\n<DataArray type=\"Int64\" "
         "Name=\"connectivity\" format=\"ascii\">\n";
  for (const auto& cell : cells) {
    for (auto n : cell) out << n << ' ';
    out << '\n';
  }
  out << "</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" "
         "format=\"ascii\">\n";
  for (std::size_t i = 1; i <= cells.size(); ++i) out << 8 * i << ' ';
  out << "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" "
         "format=\"ascii\">\n";
  for (std::size_t i = 0; i < cells.size(); ++i) out << "12 ";
  out << "\n</DataArray>\n</Cells>\n<PointData Scalars=\"pu_level_set\">\n";
  const auto write = [&](const char* name, auto value) {
    out << "<DataArray type=\"Float64\" Name=\"" << name
        << "\" format=\"ascii\">\n";
    for (const auto& sample : samples) out << value(sample) << '\n';
    out << "</DataArray>\n";
  };
  write("pu_level_set", [](const Sample& s) { return s.value; });
  write("mean_curvature", [](const Sample& s) { return s.curvature; });
  write("mean_curvature_error",
        [](const Sample& s) { return s.curvature_error; });
  out << "</PointData>\n</Piece></UnstructuredGrid></VTKFile>\n";
  out.close();
  std::cout << "PU sampling: " << cells.size() << " supported cells of "
            << std::size_t(sample_nx) * sample_nx * sample_nx << "; "
            << filename << '\n';
  if (cells.empty())
    std::cout << "No supported sampling cells; increase sampling resolution or "
                 "PU radius.\n";
}
}  // namespace LevelSetVisualization
