#include <iostream>
#include <string>
#include <fstream>
#include <functional>

#include "examples/initialize_vf/functions.h"

int main(int argc, char* argv[]){
  // if (argc != 4){
  //   std::cout << "Incorrect amount of command line arguments supplied. \n";
  //   std::cout << "Arguments should be \n";
  //   std::cout << "Implicit Surface. Options: sphere, ellipsoid\n";
  //   std::cout << "Maximum levels of refinement \n";
  //   std::cout << "Number of cells (integer)\n";
  //   std::exit(-1);
  // }
  // std::string implicit_surface = argv[1];
  // int max_refine_level = std::atoi(argv[2]);
  // int Nx = std::atoi(argv[3]);
  std::string implicit_surface = "sphere";
  int Nx = 32;
  int max_refine_level = 0;

  // mesh
  const int GC = 5;
  const IRL::Pt lower_domain(-0.5, -0.5, -0.5);
  const IRL::Pt upper_domain(0.5, 0.5, 0.5);
  BasicMesh mesh(Nx, Nx, Nx, GC);
  mesh.setCellBoundaries(lower_domain, upper_domain);

  // implicit surface equations
  ImplicitF F;
  GradientF gradF;
  HessianF hessF;
  selectSurface(implicit_surface, F, gradF, hessF);

  // refine grid
  std::unordered_map<Eigen::Vector3d, double, Vector3dHash> F_cache;
  auto grid = refine_grid(mesh, F, gradF, hessF, max_refine_level, F_cache);
  
  // compute vf
  Data<double> vf(&mesh);
  int idx = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); i++){
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++){
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++){
        AreaResult result;
        get_cell_area(*grid[idx], result, F, gradF, hessF);
        vf(i,j,k) = result.volume_fraction();
        ++idx;
      }
    }
  }

  // vtk outputs
  std::string path = "/home/parinht2/Desktop/ppic paper/implicit surface reconstruction/initialize_vf/";
  std::vector<const Cell*> leaf_cells, mixed_leaf_cells;
  for (const auto& root : grid){
    collect_leaf_cells(*root, leaf_cells, mixed_leaf_cells);
  }

  return 0;
}