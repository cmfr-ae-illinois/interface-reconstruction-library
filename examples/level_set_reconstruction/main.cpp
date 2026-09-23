#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "examples/level_set_reconstruction/level_set.h"
#include "examples/variant_advector/solver.h"
#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"

namespace {
// These methods use volume fractions, not input phase centroids.
const std::vector<std::string> methods = {
    "ELVIRA",  "LVIRA",       "PLIC",    "Jibben", "iJibben",
    "Jibben2", "JibbenCubic", "JibbenM", "PU",     "MixedJibben"};

void usage(const char* executable) {
  std::cout << "Usage: " << executable
            << " [nx=16] [method=Jibben] [output=level_set_viz]\nMethods:";
  for (const auto& method : methods) std::cout << ' ' << method;
  std::cout << "\nEdit level_set.h to change F, gradF, and hessF.\n";
}

void run(const int nx, const std::string& method,
         const std::string& output_directory) {
  BasicMesh mesh(nx, nx, nx, 3);
  mesh.setCellBoundaries(IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, 0.5));
  IRL::setMinimumVolumeToTrack(10.0 * std::numeric_limits<double>::epsilon() *
                               mesh.cell_volume());
  IRL::setVolumeFractionBounds(1.0e-13);
  IRL::setVolumeFractionTolerance(1.0e-12);

  Data<double> volume_fraction(&mesh), level_set(&mesh), zero_velocity(&mesh);
  const LevelSet surface;
  double total_volume = 0.0;
  std::size_t mixed_cells = 0;
  std::cout << "Computing zeroth moments with five refinement levels..."
            << std::endl;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const IRL::ImplicitSurfaceCutter<LevelSet, IRL::Volume> cutter(surface,
                                                                       cell);
        const double fraction =
            static_cast<double>(cutter.computeVolumeMoments()) /
            mesh.cell_volume();
        if (!std::isfinite(fraction) || fraction < -1.0e-10 ||
            fraction > 1.0 + 1.0e-10) {
          throw std::runtime_error("Invalid cut volume in cell " +
                                   std::to_string(i) + "," + std::to_string(j) +
                                   "," + std::to_string(k));
        }
        volume_fraction(i, j, k) = std::clamp(fraction, 0.0, 1.0);
        total_volume += volume_fraction(i, j, k) * mesh.cell_volume();
        mixed_cells += fraction >= IRL::global_constants::VF_LOW &&
                       fraction <= IRL::global_constants::VF_HIGH;
        level_set(i, j, k) = surface.F(mesh.xm(i), mesh.ym(j), mesh.zm(k));
      }
    }
  }
  volume_fraction.updateBorder();

  // Compatibility adapter for variant_advector: only volume() is populated.
  // No first/higher volume moments or implicit surface moments are computed.
  Data<IRL::VolumeMoments> liquid(&mesh), gas(&mesh);
  Data<IRL::SeparatorVariant> interface(&mesh);
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        liquid(i, j, k) = IRL::VolumeMoments();
        gas(i, j, k) = IRL::VolumeMoments();
        liquid(i, j, k).volume() =
            volume_fraction(i, j, k) * mesh.cell_volume();
        gas(i, j, k).volume() =
            (1.0 - volume_fraction(i, j, k)) * mesh.cell_volume();
        zero_velocity(i, j, k) = 0.0;
      }
    }
  }

  std::cout << std::setprecision(16) << "Liquid volume: " << total_volume
            << "; mixed cells: " << mixed_cells << '\n'
            << "Reconstructing with " << method << "..." << std::endl;
  std::vector<InterfaceScalarField> scalar_fields;
  getReconstruction(method, liquid, gas, 0.0, zero_velocity, zero_velocity,
                    zero_velocity, &interface, &scalar_fields);

  std::filesystem::create_directories(output_directory);
  VTKOutput output(output_directory, "level_set", mesh);
  output.addData("volume_fraction", volume_fraction);
  output.addData("level_set", level_set);
  output.writeVTKFile(0.0);
  writeInterfaceToFile(liquid, interface, 0.0, &output, true);
  std::cout << "VTK output written to " << output_directory << std::endl;
}
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      usage(argv[0]);
      return 0;
    }
    if (argc > 4) throw std::invalid_argument("Too many arguments");
    const std::string nx_text = argc > 1 ? argv[1] : "16";
    std::size_t consumed = 0;
    const int nx = std::stoi(nx_text, &consumed);
    if (consumed != nx_text.size() || nx < 4)
      throw std::invalid_argument("nx must be an integer >= 4");
    const std::string method = argc > 2 ? argv[2] : "Jibben";
    if (std::find(methods.begin(), methods.end(), method) == methods.end())
      throw std::invalid_argument("Unsupported volume-fraction-only method: " +
                                  method);
    run(nx, method, argc > 3 ? argv[3] : "level_set_viz");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    usage(argv[0]);
    return 1;
  }
}
