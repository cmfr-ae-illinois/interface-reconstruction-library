#include <fstream>
#include <functional>
#include <iostream>
#include <string>

#include <mpi.h>

#include "examples/2d_advector/reconstruction_types.h"
#include "examples/2d_advector/vtk.h"
#include "examples/pu_analysis/functions.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::string implicit_surface = "sine";
  int Nx = 64;
  int max_refine_level = 8;

  // mesh
  const IRL2D::Vec lower_domain = {-8, -8};
  const IRL2D::Vec upper_domain = {8, 8};
  BasicMesh mesh = setMesh(Nx, lower_domain, upper_domain);

  // initialization
  Data<double> velU(&mesh);
  Data<double> velV(&mesh);
  Data<double> vfrac(&mesh);
  Data<IRL2D::Moments> liquid_moments(&mesh);
  Data<IRL2D::Moments> gas_moments(&mesh);
  Data<IRL2D::Parabola> interface(&mesh);

  // implicit surface equations
  ImplicitF F;
  GradientF gradF;
  HessianF hessF;
  selectSurface(implicit_surface, F, gradF, hessF);

  // refine grid
  std::unordered_map<Point, double, PointHash> F_cache;
  auto grid = refine_grid(mesh, F, max_refine_level, gradF, hessF, F_cache);

  // compute volume fraction
  Data<double> vf(&mesh);
  int idx = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      AreaResult result;
      get_cell_area(*grid[idx], result, F, gradF, hessF);
      vf(i, j) = result.volume_fraction();
      liquid_moments(i, j).m0() = vf(i, j) * mesh.dx() * mesh.dx();
      gas_moments(i, j).m0() =
          (mesh.dx() * mesh.dx()) - liquid_moments(i, j).m0();
      ++idx;
    }
  }

  // plic reconstruction
  std::string reconstruction_method = "LVIRA";
  getReconstruction(reconstruction_method, liquid_moments, gas_moments, 0.0,
                    velU, velV, &interface);

  // interface vtk output
  if (rank == 0) {
    VTKOutput vtk_io("viz_out", "viz", mesh);
    double simulation_time = 0.0;
    vtk_io.writeVTKInterface(simulation_time, interface);
  }

  // vtk outputs
  std::vector<const Cell*> leaf_cells, mixed_leaf_cells;
  for (const auto& root : grid) {
    collect_leaf_cells(*root, leaf_cells, mixed_leaf_cells);
  }
  if (rank == 0) {
    std::string path =
        "/home/parinht2/Desktop/partition_of_unity/2d_interface/";
    write_vtr(path + "vf.vtr", vf, mesh);
    write_vtu(path + "amr.vtu", leaf_cells);
    write_vtu(path + "mixed_amr.vtu", mixed_leaf_cells);
  }

  MPI_Finalize();

  return 0;
}