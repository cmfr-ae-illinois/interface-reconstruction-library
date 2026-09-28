#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "examples/level_set_reconstruction/diagnostics_writer.h"
#include "examples/level_set_reconstruction/level_set.h"
#include "examples/level_set_reconstruction/visualization.h"
#include "examples/level_set_reconstruction/weights/wendland4.h"
#include "examples/level_set_reconstruction/weights/wendland6.h"
#include "examples/level_set_reconstruction/weights/wu2.h"
#include "examples/level_set_reconstruction/weights/wu4.h"
#include "examples/variant_advector/solver.h"
#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"
#include "irl/helpers/wendland.h"
namespace {
const std::vector<std::string> methods = {"LVIRA", "Jibben"};
const std::vector<std::string> weights = {"Wu2", "Wu4", "Wendland2",
                                          "Wendland4", "Wendland6"};
// error message when invalid command line arguments are provided
void usage(const char* executable) {
  std::cout
      << "Usage: " << executable
      << " [nx=16] [method=Jibben] [output=level_set_viz] "
         "[sample_nx=2*nx] [pu_radius_cells=2.5] [level_set=sphere]\nMethods:";
  for (const auto& method : methods) std::cout << ' ' << method;
  std::cout << "\nLevel sets:";
  for (const auto& entry : ExampleLevelSets::registry())
    std::cout << ' ' << entry.name;
  std::cout << "\nAdd definitions in level_set.h to extend this list.\n";
}
template <class WeightFunction>
void run(const int nx, const std::string& method,
         const std::string& output_directory, const int sample_nx,
         const double radius_cells, const LevelSet& surface) {
  // Output Variables
  int nx_out;
  int sample_nx_out;
  double radius_out;
  double position_error;
  double normal_error;
  double curvature_error;
  std::string shape;
  // Initialize Output Variables
  nx_out = nx;
  sample_nx_out = sample_nx;
  radius_out = radius_cells;
  position_error = 0.0;
  normal_error = 0.0;
  curvature_error = 0.0;
  shape = surface.name();
  // Output Assignment
  DiagnosticsWriter errors(output_directory + "/errors.txt");
  errors.add("nx", nx);
  errors.add("sample_nx", sample_nx);
  errors.add("radius_cells", radius_cells);
  errors.add("method", method);
  errors.add("shape", shape);

  errors.add("position_error", position_error);
  errors.add("normal_error", normal_error);
  errors.add("curvature_error", curvature_error);

  errors.writeHeader();
  // cutting operation for generating volume fraction field for level set of
  // choice
  const int ghost_layers =
      std::max(3, static_cast<int>(std::ceil(radius_cells)));
  BasicMesh mesh(nx, nx, nx, ghost_layers);
  mesh.setCellBoundaries(IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, 0.5));
  IRL::setMinimumVolumeToTrack(10.0 * std::numeric_limits<double>::epsilon() *
                               mesh.cell_volume());
  IRL::setVolumeFractionBounds(1.0e-13);
  IRL::setVolumeFractionTolerance(1.0e-12);

  Data<double> volume_fraction(&mesh), level_set(&mesh), zero_velocity(&mesh);
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

  // storing exact volume fractions in liquid and gas moments for reconstruction
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

  // perform reconstruction using method of choice
  std::vector<InterfaceScalarField> scalar_fields;
  LVIRA::getReconstruction(liquid, gas, 0.0, zero_velocity, zero_velocity,
                           zero_velocity, &interface, &scalar_fields);
  // Preserve the clipped LVIRA polygon centroids before Jibben replaces planes.
  Data<IRL::Pt> representative_points(&mesh);
  LevelSetVisualization::storePolygonCentroids(volume_fraction, interface,
                                               &representative_points);
  if (method == "Jibben") {
    Jibben::getReconstruction(liquid, gas, 0.0, zero_velocity, zero_velocity,
                              zero_velocity, &interface, &scalar_fields, true);
  }

  // vtk outputs
  std::filesystem::create_directories(output_directory);

  // outputting level set of choice
  VTKOutput output(output_directory, "level_set", mesh);
  output.addData("volume_fraction", volume_fraction);
  output.addData("level_set", level_set);
  output.writeVTKFile(0.0);

  // Outputting reconstructed interface and scalar fields
  interface.updateBorder();
  correctInterfaceBorders(&interface);
  // adding scalar field data on reconstructed interface
  LevelSetVisualization::addInterfaceDiagnostics<WeightFunction>(
      volume_fraction, interface, surface, method, &scalar_fields);
  writeInterfaceWithScalarToFile(liquid, interface, &scalar_fields, 0.0,
                                 &output, true);
  std::cout << "Sampling PU directly from reconstructed interfaces..."
            << std::endl;

  // outputting actual PU field
  const LevelSetVisualization::ReconstructedPU pu(volume_fraction, interface,
                                                  radius_cells);
  LevelSetVisualization::writePUField<WeightFunction>(
      pu, mesh, surface, sample_nx, output_directory + "/pu_field.vtu");
  //---------------------------------
  // Update Output Values for errors
  // --------------------------------

  // reconstructed pu paraboloid interface and scalar fields
  Data<IRL::SeparatorVariant> pu_ppic(&mesh);
  LevelSetVisualization::reconstructPUPPIC(pu, volume_fraction, interface,
                                           representative_points, &pu_ppic);
  LevelSetVisualization::addInterfaceDiagnostics<WeightFunction>(
      volume_fraction, pu_ppic, surface, "PU PPIC", &scalar_fields);
  VTKOutput ppic_output(output_directory, "pu_ppic", mesh);
  writeInterfaceWithScalarToFile(liquid, pu_ppic, &scalar_fields, 0.0,
                                 &ppic_output, true);
  std::cout << "VTK output written to " << output_directory << std::endl;
}
}  // namespace

int main(int argc, char** argv) {
  // parsing command line arguements
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      usage(argv[0]);
      return 0;
    }
    if (argc > 8) throw std::invalid_argument("Too many arguments");
    // ------------------------------------------------------------
    // Grid resolution
    // ------------------------------------------------------------
    const std::string nx_text = argc > 1 ? argv[1] : "16";
    std::size_t consumed = 0;
    const int nx = std::stoi(nx_text, &consumed);
    if (consumed != nx_text.size() || nx < 4)
      throw std::invalid_argument("nx must be an integer >= 4");

    // ------------------------------------------------------------
    // Reconstruction method
    // ------------------------------------------------------------
    const std::string method = argc > 2 ? argv[2] : "Jibben";
    if (std::find(methods.begin(), methods.end(), method) == methods.end())
      throw std::invalid_argument("Unsupported volume-fraction-only method: " +
                                  method);

    // ------------------------------------------------------------
    // Sampling resolution
    // ------------------------------------------------------------
    const std::string sample_text = argc > 4 ? argv[4] : std::to_string(2 * nx);
    const int sample_nx = std::stoi(sample_text, &consumed);
    if (consumed != sample_text.size() || sample_nx < 1)
      throw std::invalid_argument("sample_nx must be a positive integer");

    // ------------------------------------------------------------
    // PU radius
    // ------------------------------------------------------------
    const std::string radius_text = argc > 5 ? argv[5] : "2.5";
    const double radius_cells = std::stod(radius_text, &consumed);
    if (consumed != radius_text.size() || !std::isfinite(radius_cells) ||
        radius_cells <= 0.0 || radius_cells > 5.0)
      throw std::invalid_argument("pu_radius_cells must be > 0 and <= 5");

    // Periodic border filling copies from the physical domain in one pass.
    if (nx < std::max(3, static_cast<int>(std::ceil(radius_cells))))
      throw std::invalid_argument("nx must be at least ceil(pu_radius_cells)");

    // ------------------------------------------------------------
    // Shape
    // ------------------------------------------------------------
    const std::string shape = argc > 6 ? argv[6] : "sphere";

    // ------------------------------------------------------------
    // Weight function
    // ------------------------------------------------------------
    const std::string weight = argc > 7 ? argv[7] : "Wendland2";

    if (std::find(weights.begin(), weights.end(), weight) == weights.end())
      throw std::invalid_argument("Unsupported weighting function: " + weight);

    // ------------------------------------------------------------
    // Run
    // ------------------------------------------------------------
    const LevelSet surface(shape);
    std::cout << "Selected level set: " << shape << '\n';
    std::cout << "Selected weight: " << weight << '\n';
    const std::string output = argc > 3 ? argv[3] : "level_set_viz";

    if (weight == "Wu2") {
      std::cout << "Running with Wu2 weight function\n";
      run<IRL::Wu2>(nx, method, output, sample_nx, radius_cells, surface);

    } else if (weight == "Wu4") {
      std::cout << "Running with Wu4 weight function\n";
      run<IRL::Wu4>(nx, method, output, sample_nx, radius_cells, surface);

    } else if (weight == "Wendland2") {
      std::cout << "Running with Wendland2 weight function\n";
      run<IRL::Wendland>(nx, method, output, sample_nx, radius_cells, surface);

    } else if (weight == "Wendland4") {
      std::cout << "Running with Wendland4 weight function\n";
      run<IRL::Wendland4>(nx, method, output, sample_nx, radius_cells, surface);

    } else if (weight == "Wendland6") {
      std::cout << "Running with Wendland6 weight function\n";
      run<IRL::Wendland6>(nx, method, output, sample_nx, radius_cells, surface);
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    usage(argv[0]);
    return 1;
  }
}
