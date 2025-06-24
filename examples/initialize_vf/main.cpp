#include <iostream>
#include <string>
#include <fstream>
#include <functional>

#include "examples/initialize_vf/functions.h"

void initialize_vf(const std::string& implicit_surface, const int& max_refine_level,
                   Data<double>& vf);

int main(int argc, char* argv[]){
  // if (argc != 4){
  //   std::cout << "Incorrect amount of command line arguments supplied. \n";
  //   std::cout << "Arguments should be \n";
  //   std::cout << "Implicit Surface. Options: Circle\n";
  //   std::cout << "Maximum levels of refinement \n"
  //   std::cout << "Number of cells (integer)\n";
  //   std::exit(-1);
  // }
  // std::string implicit_surface = argv[1];
  // int max_refine_level = std::atoi(argv[2]);
  // int Nx = std::atoi(argv[3]);

  std::string implicit_surface = "circle";
  int Nx = 16;
  int max_refine_level = 5;

  // mesh
  const IRL2D::Vec lower_domain = {-0.5, -0.5};
  const IRL2D::Vec upper_domain = {0.5, 0.5};
  BasicMesh mesh = setMesh(Nx, lower_domain, upper_domain);

  // get volume fraction field
  Data<double> vf(&mesh);
  initialize_vf(implicit_surface, max_refine_level, vf); 
  std::cout << vf(6,5) << std::endl;

  // output to vtk
  write_vtr("/home/parinht2/Desktop/ppic paper/reconstruction/initialize_vf/vf.vtr", vf, mesh);


  return 0;
}

void initialize_vf(const std::string& implicit_surface, const int& max_refine_level,
                   Data<double>& vf){  

  // get mesh
  const BasicMesh& mesh = vf.getMesh();

  // implicit surface equations
  ImplicitF F;
  GradientF gradF;
  HessianF hessF;
  selectSurface(implicit_surface, F, gradF, hessF);

  // double x = mesh.x(6), y = mesh.y(5);
  // Cell cell(x, y, mesh.dx(), 0);
  // CellStatus status = get_cell_status(cell, F);
  // if (status == CellStatus::Above){
  //   std::cout << "Above" << std::endl;
  // } else if (status == CellStatus::Below){
  //   std::cout << "Below" << std::endl;
  // } else if (status == CellStatus::Mixed){
  //   std::cout << "Mixed" << std::endl;
  // }

  // auto cell_ptr = std::make_unique<Cell>(x, y, mesh.dx(), 0);
  // refine_cell(cell_ptr, F, max_refine_level);

  // refine cells
  auto grid = refine_grid(mesh, F, max_refine_level);

  // compute volume fraction
  int idx = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); i++){
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++){
      AreaResult result;
      get_cell_area(*grid[idx], result, F, gradF, hessF);
      vf(i,j) = result.volume_fraction();
      ++idx;
    }
  }
}