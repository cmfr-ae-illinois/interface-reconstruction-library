#include <iostream>
#include <string>
#include <fstream>
#include <functional>
#include <chrono>

#include "examples/initialize_vf/functions.h"

int main(int argc, char* argv[]){
  if (argc != 4){
    std::cout << "Incorrect amount of command line arguments supplied. \n";
    std::cout << "Arguments should be \n";
    std::cout << "Implicit Surface. Options: sphere, ellipsoid, torus, genus, orthocircle, wineglass (more to be added)\n";
    std::cout << "Maximum levels of refinement \n";
    std::cout << "Number of cells (integer)\n";
    std::exit(-1);
  }
  std::string implicit_surface = argv[1];
  int max_refine_level = std::atoi(argv[2]);
  int Nx = std::atoi(argv[3]);
  std::string vtk_path = "/home/parinht2/Desktop/ppic paper/implicit surface reconstruction/initialize_vf/3D/";

  // vf initializer
  vfInitializer vfi(Nx, max_refine_level, implicit_surface, vtk_path);
  vfi.run();
  double volume_amr = vfi.getTotalVolume();
  auto mixed_leaf_cells = vfi.getMixedLeafCells();

  // volume checks
  // double volume = 4.0/3.0 * M_PI * 0.15 * 0.15 * 0.15; // sphere
  double volume = 4.0/3.0 * M_PI * 0.3 * 0.15 * 0.1 ; // ellipsoid
  // double volume = 2.0 * M_PI * M_PI * 0.3 * 0.15 * 0.15; // torus
  // std::cout << "true volume = " << volume << std::endl;
  // std::cout << "amr volume = " << volume_amr << std::endl;
  // std::cout << "Percent diff = " << std::abs(volume - volume_amr) / volume * 100.0 << std::endl;
  // std::cout << "Execution time: " << elapsed.count() << " seconds\n";

  // output for convergence
  std::string conv_csv = "/home/parinht2/Desktop/ppic paper/implicit surface reconstruction/initialize_vf/3D/convergence_test/convergence.csv";
  std::ifstream test(conv_csv);
  bool exists = test.good();
  test.close();
  std::ofstream out(conv_csv, std::ios::app);
  if (out.is_open()){
    if (!exists){
      out << "Refinement Level,Refined Cell Spacing,True Volume,AMR Volume,Relative Error\n";
    }
    out << std::fixed << std::setprecision(16)
        << max_refine_level << ","
        << mixed_leaf_cells[0]->dx << ","
        << volume << ","
        << volume_amr << ","
        << std::abs(volume - volume_amr) / volume << "\n";
    out.close();
  }


  // compute surface moments
  // start = std::chrono::high_resolution_clock::now();
  // double surface_area_amr = 0.0;
  // idx = 0;
  // int N = 100;
  // for (int i = mesh.imin(); i <= mesh.imax(); i++){
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); j++){
  //     for (int k = mesh.kmin(); k <= mesh.kmax(); k++){
  //       SurfaceMoments sm;
  //       get_shape_moments(*grid[idx], sm, N, F, gradF);
  //       surface_area_amr += sm.M0;
  //       ++idx;
  //     }
  //   }
  // }
  // end = std::chrono::high_resolution_clock::now();
  // elapsed = end - start;
  // start = std::chrono::high_resolution_clock::now();
  // double surface_area_amr = 0.0;
  // int N = 20;
  // for (const Cell* c : mixed_leaf_cells){
  //   SurfaceMoments sm;
  //   get_shape_moments(*c, sm, N, F, gradF);
  //   surface_area_amr += sm.M0;
  // }
  // end = std::chrono::high_resolution_clock::now();
  // elapsed = end - start;
  // // surface area
  // double surface_area = 4.0 * M_PI * 0.15 * 0.15;
  // std::cout << "True surface area = " << surface_area << std::endl;
  // std::cout << "amr + coarea surface area = " << surface_area_amr << std::endl;
  // std::cout << "Percent diff = " << std::abs(surface_area - surface_area_amr) / surface_area * 100.0 << std::endl;
  // std::cout << "Execution time: " << elapsed.count() << " seconds\n";

  return 0;
}