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
    std::cout << "Implicit Surface. Options: sphere, ellipsoid, genus, orthocircle (more to be added)\n";
    std::cout << "Maximum levels of refinement \n";
    std::cout << "Number of cells (integer)\n";
    std::exit(-1);
  }
  std::string implicit_surface = argv[1];
  int max_refine_level = std::atoi(argv[2]);
  int Nx = std::atoi(argv[3]);
  std::string vtk_path = "/home/parinht2/Desktop/ppic paper/implicit surface reconstruction/initialize_vf/3D/";

  // initialize volume fraction field
  std::cout << "Performing amr" << std::endl;
  vfInitializer vfi(Nx, max_refine_level, implicit_surface, vtk_path);
  vfi.run();
  std::cout << "amr complete" << std::endl;
  const Surface* surface = vfi.getSurface();
  const auto& refined_grid = vfi.getRefinedGrid();

  // total mixed cells
  int total_mixed_cells = 0;
  for (const auto& cell : refined_grid) {
    if (cell->status == CellStatus::Mixed && !cell->children.empty()) {
      total_mixed_cells++;
    }
  }

  auto start = std::chrono::high_resolution_clock::now();

  // surface moments of IS using marching cubes
  // const int num_division = 10;
  // double true_area = surface->surfaceArea();
  // double amr_area = 0.0;
  // Metrics implicit_surface_metrics;
  // int current = 1;
  // for (const auto& cell : refined_grid){
  //   if (cell->status == CellStatus::Mixed && !cell->children.empty()){
  //     std::cout << "Processing cell " << current << " / " << total_mixed_cells << std::endl;
  //     SurfaceMetrics sm(*cell, 
  //                       [surface](double x, double y, double z) {
  //                         return surface->F(x, y, z); }, 
  //                       num_division);
  //     Metrics cell_metrics = sm.computeSurfaceMetrics();
  //     implicit_surface_metrics += cell_metrics;
  //     current++;
  //   }
  // }
  // amr_area = implicit_surface_metrics.sM0;
  // std::cout << "true area = " << true_area << std::endl;
  // std::cout << "amr area = " << amr_area << std::endl;
  // std::cout << "rel error = " << std::abs(true_area - amr_area) / true_area << std::endl;

  // estimate surface moment of a sphere using paraboloid surface moments
  auto& mixed_leaf_cells = vfi.getMixedLeafCells();
  double computed_area = 0.0;
  Eigen::Vector3d surface_M1(0.0, 0.0, 0.0);
  using VolumeMomentsAndSurface = IRL::AddSurfaceOutput<IRL::VolumeMoments, IRL::ParaboloidParametrizedSurfaceOutput>;
  for (const auto& c : mixed_leaf_cells){
    // build parabola
    auto cell = IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(c->x, c->y, c->z), 
                                                        IRL::Pt(c->x + c->dx, c->y + c->dx, c->z + c->dx));
    Eigen::Vector3d x_center(c->x + 0.5*c->dx, c->y + 0.5*c->dx, c->z + 0.5*c->dx);
    auto F = [surface](double x, double y, double z){ return surface->F(x,y,z); };
    auto gradF = [surface](double x, double y, double z){ return surface->gradF(x,y,z); };
    auto hessF = [surface](double x, double y, double z){ return surface->hessF(x,y,z); };
    Eigen::Vector3d x_proj = project_onto_surface(x_center, F, gradF);
    IRL::Pt x_proj_pt(x_proj(0), x_proj(1), x_proj(2));
    Eigen::Vector3d gradF_val = gradF(x_proj(0), x_proj(1), x_proj(2));
    Eigen::Matrix3d hessF_val = hessF(x_proj(0), x_proj(1), x_proj(2));
    IRL::Paraboloid paraboloid = IRL::Paraboloid::fromDerivatives(x_proj_pt, gradF_val, hessF_val);
    // compute parabola surface moments
    auto our_surface_and_moments = IRL::getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid);
    auto our_moments = our_surface_and_moments.getMoments();
    auto our_surface = our_surface_and_moments.getSurface();
    auto our_surface_area = our_surface.getSurfaceArea();
    computed_area += our_surface_area;
  }
  double true_area = surface->surfaceArea();
  std::cout << "true area = " << true_area << std::endl;
  std::cout << "computed area = " << computed_area << std::endl;
  std::cout << "rel error = " << std::abs(true_area - computed_area) / true_area << std::endl;

  // // marching cubes on paraboloid
  // auto& mixed_leaf_cells = vfi.getMixedLeafCells();
  // double computed_area = 0.0;
  // int current = 1;
  // using VolumeMomentsAndSurface = IRL::AddSurfaceOutput<IRL::VolumeMoments, IRL::ParaboloidParametrizedSurfaceOutput>;
  // for (const auto& c : mixed_leaf_cells){
  //   std::cout << "Processing cell " << current << " / " << total_mixed_cells << std::endl;
  //   // build parabola
  //   auto cell = IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(c->x, c->y, c->z), 
  //                                                       IRL::Pt(c->x + c->dx, c->y + c->dx, c->z + c->dx));
  //   Eigen::Vector3d x_center(c->x + 0.5*c->dx, c->y + 0.5*c->dx, c->z + 0.5*c->dx);
  //   auto F = [surface](double x, double y, double z){ return surface->F(x,y,z); };
  //   auto gradF = [surface](double x, double y, double z){ return surface->gradF(x,y,z); };
  //   auto hessF = [surface](double x, double y, double z){ return surface->hessF(x,y,z); };
  //   Eigen::Vector3d x_proj = project_onto_surface(x_center, F, gradF);
  //   IRL::Pt x_proj_pt(x_proj(0), x_proj(1), x_proj(2));
  //   Eigen::Vector3d gradF_val = gradF(x_proj(0), x_proj(1), x_proj(2));
  //   Eigen::Matrix3d hessF_val = hessF(x_proj(0), x_proj(1), x_proj(2));
  //   IRL::Paraboloid paraboloid = IRL::Paraboloid::fromDerivatives(x_proj_pt, gradF_val, hessF_val);
  //   // global to local coordinates (implicit function)
  //   double a = paraboloid.getAlignedParaboloid().a();
  //   double b = paraboloid.getAlignedParaboloid().b();
  //   IRL::Pt datum = paraboloid.getDatum();
  //   auto frame = paraboloid.getReferenceFrame();
  //   auto f_paraboloid = [=](const IRL::Pt& x) -> double {
  //     Eigen::Vector3d x_global(x[0], x[1], x[2]);
  //     Eigen::Vector3d x0(datum[0], datum[1], datum[2]);
  //     Eigen::Matrix3d R;
  //     for (int i = 0; i < 3; i++) R.col(i) = Eigen::Vector3d(frame[i][0], frame[i][1], frame[i][2]);
  //     Eigen::Vector3d x_local = R.transpose() * (x_global - x0);
  //     double xl = x_local[0], yl = x_local[1], zl = x_local[2];
  //     return zl + a * xl * xl + b * yl * yl;
  //   };
  //   // marching cubes
  //   IRL::MarchingCubes mc(cell, f_paraboloid);
  //   auto surface = mc.triangulate(100);
  //   auto triangles = surface.getTriangleList();
  //   // accumulating area from triangles
  //   for (const auto& tri : triangles){
  //     Eigen::Vector3d v0(tri[0][0], tri[0][1], tri[0][2]);
  //     Eigen::Vector3d v1(tri[1][0], tri[1][1], tri[1][2]);
  //     Eigen::Vector3d v2(tri[2][0], tri[2][1], tri[2][2]);
  //     Triangle triangle(v0, v1, v2);
  //     computed_area += triangle.computeM0();
  //   }
  //   current++;
  // }
  // double true_area = surface->surfaceArea();
  // std::cout << "true area = " << true_area << std::endl;
  // std::cout << "computed area = " << computed_area << std::endl;
  // std::cout << "rel error = " << std::abs(true_area - computed_area) / true_area << std::endl;

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  std::cout << "Computation time: " << duration.count() << std::endl;

  // Convergence for volume
  std::string csv_path = "/home/parinht2/Desktop/ppic paper/implicit surface reconstruction/initialize_vf/3D/convergence_test/convergence.csv";
  // performConvergence(csv_path, max_refine_level, vfi);

  // using CGAL



  // compute surface moments (coarea)
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

  return 0;
}