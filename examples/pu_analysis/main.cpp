#include <fstream>
#include <functional>
#include <iostream>
#include <string>

#include <mpi.h>

#include <sys/stat.h>

#include "examples/pu_analysis/functions.h"
#include "examples/pu_analysis/reconstruction_types.h"
#include "examples/pu_analysis/translation.h"
#include "examples/pu_analysis/vof_advection.h"
#include "examples/pu_analysis/vtk.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // implicit surface
  std::string implicit_surface = "circle";
  int Nx = 64;
  int max_refine_level = 8;
  ImplicitF F;
  GradientF gradF;
  HessianF hessF;
  selectSurface(implicit_surface, F, gradF, hessF);

  // mesh
  constexpr int GC = 5;
  BasicMesh mesh(Nx, Nx, GC);
  IRL2D::Vec my_lower_domain(-0.5, -0.5);
  IRL2D::Vec my_upper_domain(0.5, 0.5);
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  Data<double> vf(&mesh);

  // initializing variables
  Data<double> velU(&mesh);
  Data<double> velV(&mesh);
  Data<double> vfrac(&mesh);
  Data<IRL2D::Moments> liquid_moments(&mesh);
  Data<IRL2D::Moments> gas_moments(&mesh);
  Data<IRL2D::Parabola> interface(&mesh);

  // volume fraction initialization
  std::unordered_map<Point, double, PointHash> F_cache;
  auto grid = refine_grid(mesh, F, max_refine_level, gradF, hessF, F_cache);
  int idx = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      AreaResult result;
      get_cell_area(*grid[idx], result, F, gradF, hessF);
      vf(i, j) = result.volume_fraction();
      vfrac(i, j) = vf(i, j);
      liquid_moments(i, j).m0() = vfrac(i, j) * mesh.dx() * mesh.dx();
      gas_moments(i, j).m0() =
          (mesh.dx() * mesh.dx()) - liquid_moments(i, j).m0();
      ++idx;
    }
  }
  liquid_moments.updateBorder();
  gas_moments.updateBorder();

  // step function
  // ----------------------------------------------------------------------------------------------

  // common interface parameters
  // double coeff = 0.0;
  // IRL2D::ReferenceFrame frame = {IRL2D::Vec(1.0, 0.0), IRL2D::Vec(0.0, 1.0)};

  // manually assigning cells an interface
  // for (int i = mesh.imin() + 5; i <= mesh.imax() - 5; i++) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
  //     if (j == (Nx / 2) && i != 40 && i != 41) {
  //       IRL2D::Parabola parabola;
  //       parabola.frame() = frame;
  //       parabola.datum() = {mesh.xm(i), mesh.ym(j)};
  //       interface(i, j) = parabola;
  //       IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
  //           IRL2D::Vec(mesh.x(i), mesh.y(j)),
  //           IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
  //       vfrac(i, j) = IRL2D::ComputeMoments(cell, interface(i, j)).m0() /
  //                     IRL2D::ComputeArea(cell);
  //       liquid_moments(i, j) = IRL2D::ComputeMoments(cell, interface(i, j));
  //       gas_moments(i, j) = IRL2D::ComputeMoments(cell) - liquid_moments(i,
  //       j);
  //     }
  //     if (j == (Nx / 2 + 1) && (i == 40 || i == 41)) {
  //       IRL2D::Parabola parabola;
  //       parabola.frame() = frame;
  //       parabola.datum() = {mesh.xm(i), mesh.ym(j)};
  //       interface(i, j) = parabola;
  //       IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
  //           IRL2D::Vec(mesh.x(i), mesh.y(j)),
  //           IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
  //       vfrac(i, j) = IRL2D::ComputeMoments(cell, interface(i, j)).m0() /
  //                     IRL2D::ComputeArea(cell);
  //       liquid_moments(i, j) = IRL2D::ComputeMoments(cell, interface(i, j));
  //       gas_moments(i, j) = IRL2D::ComputeMoments(cell) - liquid_moments(i,
  //       j);
  //     }
  //   }
  // }

  // writing interface
  // double simulation_time = 0.0;
  // VTKOutput vtk_io("viz_out", "viz", mesh);
  // vtk_io.addData("VFrac", vfrac);
  // vtk_io.writeVTKFile(simulation_time);
  // vtk_io.writeVTKInterface(simulation_time, interface);

  // PU interface
  // const int Nx_fine = 256;
  // BasicMesh fine_mesh(Nx_fine, Nx_fine, GC);
  // fine_mesh.setCellBoundaries(IRL2D::Vec(-0.5, -0.1), IRL2D::Vec(0.5, 0.1));
  // std::vector<double> pu_F;
  // std::vector<IRL2D::Vec> centroids, normals;
  // for (int i = mesh.imin(); i <= mesh.imax(); i++) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
  //     if (vfrac(i, j) > IRL::global_constants::VF_LOW &&
  //         vfrac(i, j) < IRL::global_constants::VF_HIGH) {
  //       IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
  //           IRL2D::Vec(mesh.x(i), mesh.y(j)),
  //           IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
  //       IRL2D::BezierList clipped_interface =
  //           IRL2D::ParabolaClip(cell, interface(i, j), true);
  //       IRL2D::Vec center =
  //           (clipped_interface[0].first + clipped_interface[1].first) / 2.0;
  //       centroids.push_back(center);
  //       normals.push_back(interface(i, j).frame()[1]);
  //     }
  //   }
  // }
  // IRL2D::ImplicitSurface IS(centroids, normals, 2.5 * mesh.dx());
  // for (int j = fine_mesh.jmin(); j <= (fine_mesh.jmax() + 1); j++) {
  //   for (int i = fine_mesh.imin(); i <= (fine_mesh.imax() + 1); i++) {
  //     double Fval = IS.F(IRL2D::Vec(fine_mesh.x(i), fine_mesh.y(j)));
  //     if (std::isnan(Fval)) {
  //       Fval = 10000.;
  //     }
  //     pu_F.push_back(Fval);
  //   }
  // }
  // writeVTKStructuredPoints_PointData(
  //     "/home/parinht2/Desktop/partition_of_unity/2d_interface/step_function/"
  //     "pu_inteface.vtk",
  //     Nx_fine + 1, Nx_fine + 1, -0.5, 0.5, -0.1, 0.1, pu_F);

  // projecting points on PU interface
  // std::vector<IRL2D::Vec> pu_points;
  // for (int i = fine_mesh.imin(); i <= fine_mesh.imax(); i++) {
  //   for (int j = fine_mesh.jmin(); j <= fine_mesh.jmax(); j++) {
  //     IRL2D::Vec x0 = {fine_mesh.xm(i), fine_mesh.ym(j)};
  //     bool usePlane = false;
  //     IRL2D::Vec x_proj = IRL2D::projectToImplicitSurface(
  //         x0, centroids, normals, 2.5 * mesh.dx(), usePlane);
  //     if (std::isfinite(x_proj[0]) && std::isfinite(x_proj[1])) {
  //       pu_points.push_back(x_proj);
  //     }
  //   }
  // }
  // std::string path =
  //     "/home/parinht2/Desktop/partition_of_unity/2d_interface/step_function/";
  // std::string filepath = path + "/pu_points.vtk";
  // writeScatterVTK(pu_points, filepath);

  // // iterating using partition of unity
  // std::string reconstruction_method = "PUPLIC";
  // // std::string reconstruction_method = "PUIterate";
  // const int iter = 50;
  // for (int it = 0; it < (iter + 1); it++) {
  //   // PoU field
  //   char fname[512];
  //   std::snprintf(fname, sizeof(fname),
  //                 "/home/parinht2/Desktop/partition_of_unity/2d_interface/"
  //                 "step_function/pu_interface_%04d.vtk",
  //                 it);
  //   IRL2D::Vec minPt(-0.5, -0.1), maxPt(0.5, 0.1);
  //   dumpPoUFieldVTK(fname, interface, liquid_moments, mesh, fine_mesh, minPt,
  //                   maxPt);

  //   // interface
  //   getReconstruction(reconstruction_method, liquid_moments, gas_moments,
  //   0.0,
  //                     velU, velV, &interface);
  //   simulation_time += 1;
  //   if (it < iter) {
  //     vtk_io.writeVTKFile(simulation_time);
  //     vtk_io.writeVTKInterface(simulation_time, interface);
  //   }
  // }

  // -----------------------------------Static
  // iterations-------------------------------------------------

  // LVIRA
  std::string reconstruction_m = "LVIRA";
  getReconstruction(reconstruction_m, liquid_moments, gas_moments, 0.0, velU,
                    velV, &interface);

  // outputting PU interface
  double simulation_time = 0.0;
  VTKOutput vtk_io("viz_out", "viz", mesh);
  vtk_io.addData("VFrac", vfrac);
  vtk_io.writeVTKFile(simulation_time);
  vtk_io.writeVTKInterface(simulation_time, interface);

  // fine mesh
  const int Nx_fine = 256;
  BasicMesh fine_mesh(Nx_fine, Nx_fine, GC);
  fine_mesh.setCellBoundaries(IRL2D::Vec(-0.3, -0.3), IRL2D::Vec(0.3, 0.3));

  // pu iterations
  std::string reconstruction_method = "PUIterate";
  const int num_iter = 15;
  const int dt = 1;
  for (int it = 0; it < num_iter; it++) {
    char fname[512];
    std::snprintf(fname, sizeof(fname),
                  "/home/parinht2/Desktop/partition_of_unity/2d_interface/"
                  "static/pu_interface_%04d.vtk",
                  it);
    IRL2D::Vec minPt(-0.3, -0.3), maxPt(0.3, 0.3);
    dumpPoUFieldVTK(fname, interface, liquid_moments, mesh, fine_mesh, minPt,
                    maxPt);

    std::cout << "iteration: " << it << std::endl;
    getReconstruction(reconstruction_method, liquid_moments, gas_moments, 0.0,
                      velU, velV, &interface);
    simulation_time += dt;
    if (it < num_iter) {
      vtk_io.writeVTKFile(simulation_time);
      vtk_io.writeVTKInterface(simulation_time, interface);
    }
  }

  // -----------------------------------------------------------------------------------------------

  // //
  // -----------------------------------------------------------------------------------------------

  // // LVIRA
  // std::string reconstruction_m = "LVIRA";
  // getReconstruction(reconstruction_m, liquid_moments, gas_moments, 0.0, velU,
  //                   velV, &interface);
  // const int Nx_fine = 256;
  // BasicMesh fine_mesh(Nx_fine, Nx_fine, GC);
  // fine_mesh.setCellBoundaries(IRL2D::Vec(-0.3, -0.3), IRL2D::Vec(0.3, 0.3));
  // std::vector<double> pu_F;
  // std::vector<IRL2D::Vec> centroids, normals;
  // for (int i = mesh.imin(); i <= mesh.imax(); i++) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
  //     if (vfrac(i, j) > IRL::global_constants::VF_LOW &&
  //         vfrac(i, j) < IRL::global_constants::VF_HIGH) {
  //       IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
  //           IRL2D::Vec(mesh.x(i), mesh.y(j)),
  //           IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
  //       IRL2D::BezierList clipped_interface =
  //           IRL2D::ParabolaClip(cell, interface(i, j), true);
  //       IRL2D::Vec center =
  //           (clipped_interface[0].first + clipped_interface[1].first) / 2.0;
  //       centroids.push_back(center);
  //       normals.push_back(interface(i, j).frame()[1]);
  //     }
  //   }
  // }
  // IRL2D::ImplicitSurface IS(centroids, normals, 2.5 * mesh.dx());
  // for (int j = fine_mesh.jmin(); j <= (fine_mesh.jmax() + 1); j++) {
  //   for (int i = fine_mesh.imin(); i <= (fine_mesh.imax() + 1); i++) {
  //     double Fval = IS.F(IRL2D::Vec(fine_mesh.x(i), fine_mesh.y(j)));
  //     if (std::isnan(Fval)) {
  //       Fval = 10000.;
  //     }
  //     pu_F.push_back(Fval);
  //   }
  // }
  // writeVTKStructuredPoints_PointData(
  //     "/home/parinht2/Desktop/partition_of_unity/2d_interface/pu_inteface.vtk",
  //     Nx_fine + 1, Nx_fine + 1, -0.3, 0.3, -0.3, 0.3, pu_F);

  // // PU iteration 1
  // std::string reconstruction_method = "PU";
  // getReconstruction(reconstruction_method, liquid_moments, gas_moments, 0.0,
  //                   velU, velV, &interface);

  // // outputting parabola center points
  // std::vector<IRL2D::Vec> parabola_centroids;
  // for (int i = mesh.imin(); i <= mesh.imax(); i++) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
  //     if (vfrac(i, j) > IRL::global_constants::VF_LOW &&
  //         vfrac(i, j) < IRL::global_constants::VF_HIGH) {
  //       IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
  //           IRL2D::Vec(mesh.x(i), mesh.y(j)),
  //           IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
  //       IRL2D::BezierList clipped_interface =
  //           IRL2D::ParabolaClip(cell, interface(i, j), true);
  //       IRL2D::Vec start = clipped_interface[0].first;
  //       IRL2D::Vec end = clipped_interface[1].first;
  //       IRL2D::Vec centroid = IRL2D::centroidParabolaSegmentAnalytic(
  //           interface(i, j).datum(), interface(i, j).frame()[0],
  //           interface(i, j).frame()[1], interface(i, j).coeff(), start, end);
  //       parabola_centroids.push_back(centroid);
  //     }
  //   }
  // }
  // std::string path =
  // "/home/parinht2/Desktop/partition_of_unity/2d_interface/"; std::string
  // filepath = path + "/parabola_centroid.vtk";
  // // writeScatterVTK(parabola_centroids, filepath);

  // // outputting PU interface
  // double simulation_time = 0.0;
  // VTKOutput vtk_io("viz_out", "viz", mesh);
  // vtk_io.addData("VFrac", vfrac);
  // vtk_io.writeVTKFile(simulation_time);
  // vtk_io.writeVTKInterface(simulation_time, interface);

  // // pu iterations
  // reconstruction_method = "PUIterate";
  // const int num_iter = 10;
  // const int dt = 1;
  // for (int it = 0; it < num_iter; it++) {
  //   char fname[512];
  //   std::snprintf(fname, sizeof(fname),
  //                 "/home/parinht2/Desktop/partition_of_unity/2d_interface/"
  //                 "pu_implicit_surface/pu_interface_%04d.vtk",
  //                 it);
  //   IRL2D::Vec minPt(-0.3, -0.3), maxPt(0.3, 0.3);
  //   dumpPoUFieldVTK(fname, interface, liquid_moments, mesh, fine_mesh, minPt,
  //                   maxPt);

  //   std::cout << "iteration: " << it << std::endl;
  //   getReconstruction(reconstruction_method, liquid_moments, gas_moments,
  //   0.0,
  //                     velU, velV, &interface);
  //   simulation_time += dt;
  //   vtk_io.writeVTKFile(simulation_time);
  //   vtk_io.writeVTKInterface(simulation_time, interface);
  // }

  // //
  // -----------------------------------------------------------------------------------------------

  // amr vtk outputs
  // std::vector<const Cell*> leaf_cells, mixed_leaf_cells;
  // for (const auto& root : grid) {
  //   collect_leaf_cells(*root, leaf_cells, mixed_leaf_cells);
  // }
  // if (rank == 0) {
  //   std::string path =
  //       "/home/parinht2/Desktop/partition_of_unity/2d_interface/";
  //   write_vtr(path + "vf.vtr", vf, mesh);
  //   write_vtu(path + "amr.vtu", leaf_cells);
  //   write_vtu(path + "mixed_amr.vtu", mixed_leaf_cells);
  // }

  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------

  // implicit surface
  // std::string implicit_surface = "sine";
  // int Nx = 128;
  // int max_refine_level = 8;
  // ImplicitF F;
  // GradientF gradF;
  // HessianF hessF;
  // selectSurface(implicit_surface, F, gradF, hessF);

  // // mesh
  // BasicMesh mesh = Translation::setMesh(Nx);

  // // simulation params
  // const double U = 1.0;
  // const double L = 8.0;
  // const double T = std::abs(2.0 * L / U);
  // const double Tend = 1.0 * T;
  // // const double Tend = 0.0;
  // const double CFL = 0.4;
  // double dt = CFL * mesh.dx() / std::abs(U);
  // double simulation_time = 0.0;
  // int iteration = 0;

  // // initialization
  // Data<double> velU(&mesh);
  // Data<double> velV(&mesh);
  // Data<double> vfrac(&mesh);
  // Data<IRL2D::Moments> liquid_moments(&mesh);
  // Data<IRL2D::Moments> gas_moments(&mesh);
  // Data<IRL2D::Parabola> interface(&mesh), interface2(&mesh);
  // Translation::setVelocity(simulation_time, &velU, &velV);

  // // volume fraction initialization
  // std::unordered_map<Point, double, PointHash> F_cache;
  // auto grid = refine_grid(mesh, F, max_refine_level, gradF, hessF, F_cache);
  // int idx = 0;
  // for (int i = mesh.imin(); i <= mesh.imax(); i++) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
  //     AreaResult result;
  //     get_cell_area(*grid[idx], result, F, gradF, hessF);
  //     vfrac(i, j) = result.volume_fraction();
  //     liquid_moments(i, j).m0() = vfrac(i, j) * mesh.dx() * mesh.dx();
  //     gas_moments(i, j).m0() =
  //         (mesh.dx() * mesh.dx()) - liquid_moments(i, j).m0();
  //     ++idx;
  //   }
  // }
  // liquid_moments.updateBorder();
  // gas_moments.updateBorder();

  // // plic reconstruction
  // std::string reconstruction_method = "PU";
  // getReconstruction(reconstruction_method, liquid_moments, gas_moments, 0.0,
  //                   velU, velV, &interface);

  // // interface vtk output
  // VTKOutput vtk_io("viz_out", "viz", mesh);
  // vtk_io.addData("VelocityX", velU);
  // vtk_io.addData("VelocityY", velV);
  // vtk_io.addData("VFrac", vfrac);
  // vtk_io.writeVTKFile(simulation_time);
  // vtk_io.writeVTKInterface(simulation_time, interface);
  // std::string output_folder = "viz";
  // const int dir_err = mkdir(output_folder.c_str(), 0777);

  // // advection
  // std::string a_simulation_type = "Translation";
  // std::string a_advection_method = "FullLagL";
  // std::string a_reconstruction_method = "PU";
  // int a_visualization_frequency = 1;

  // while (simulation_time < Tend) {
  //   const double time_step_to_use = std::fmin(dt, Tend - simulation_time);
  //   Translation::setVelocity(simulation_time + 0.5 * time_step_to_use, &velU,
  //                            &velV);

  //   advectVOF(a_simulation_type, a_advection_method, a_reconstruction_method,
  //             time_step_to_use, simulation_time, velU, velV, &liquid_moments,
  //             &gas_moments, &interface);

  //   getReconstruction(a_reconstruction_method, liquid_moments, gas_moments,
  //                     time_step_to_use, velU, velV, &interface);

  //   if (a_visualization_frequency > 0 &&
  //       ((iteration + 1) % a_visualization_frequency == 0 ||
  //        time_step_to_use == Tend - simulation_time)) {
  //     for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
  //       for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
  //         vfrac(i, j) = liquid_moments(i, j).m0() / mesh.cell_volume();
  //       }
  //     }
  //     vtk_io.writeVTKFile(simulation_time);
  //     vtk_io.writeVTKInterface(simulation_time, interface);

  //     simulation_time += time_step_to_use;
  //     ++iteration;
  //     std::cout << "t = " << simulation_time << std::endl;
  //   }
  // }

  // // outputting plic data
  // std::string path =
  // "/home/parinht2/Desktop/partition_of_unity/2d_interface/"; std::string
  // filepath = path + "/plic_data.csv"; std::ofstream csvfile(filepath);
  // csvfile << "ax,ay,bx,by,tx,ty,nx,ny\n";
  // for (int ii = mesh.imin(); ii <= mesh.imax(); ii++) {
  //   for (int jj = mesh.jmin(); jj <= mesh.jmax(); jj++) {
  //     const double lvf = (liquid_moments)(ii, jj).m0() / mesh.cell_volume();
  //     if (lvf >= IRL::global_constants::VF_LOW &&
  //         lvf <= IRL::global_constants::VF_HIGH) {
  //       IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
  //           IRL2D::Vec(mesh.x(ii), mesh.y(jj)),
  //           IRL2D::Vec(mesh.x(ii + 1), mesh.y(jj + 1)));
  //       IRL2D::BezierList clipped_plic =
  //           IRL2D::ParabolaClip(cell, interface(ii, jj), true);
  //       std::pair<IRL2D::Vec, IRL2D::Vec> end_points =
  //       {clipped_plic[0].first,
  //                                                       clipped_plic[1].first};
  //       csvfile << end_points.first[0] << "," << end_points.first[1] << ","
  //               << end_points.second[0] << "," << end_points.second[1] << ","
  //               << interface(ii, jj).frame()[0][0] << ","
  //               << interface(ii, jj).frame()[0][1] << ","
  //               << interface(ii, jj).frame()[1][0] << ","
  //               << interface(ii, jj).frame()[1][1] << "\n";
  //     }
  //   }
  // }
  // csvfile.close();
  // // list of centroids and normals
  // std::vector<IRL2D::Vec> centroids;
  // std::vector<IRL2D::Vec> normals;
  // for (int i = mesh.imin(); i <= mesh.imax(); i++) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
  //     const double lvf = (liquid_moments)(i, j).m0() / mesh.cell_volume();
  //     if (lvf >= IRL::global_constants::VF_LOW &&
  //         lvf <= IRL::global_constants::VF_HIGH) {
  //       IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
  //           IRL2D::Vec(mesh.x(i), mesh.y(j)),
  //           IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
  //       IRL2D::BezierList clipped_plic =
  //           IRL2D::ParabolaClip(cell, interface(i, j), true);
  //       IRL2D::Vec center =
  //           (clipped_plic[0].first + clipped_plic[1].first) / 2.0;
  //       centroids.push_back(center);
  //       normals.push_back(interface(i, j).frame()[1]);
  //     }
  //   }
  // }
  // // PU isocontour
  // const int Nx_fine = 512;
  // const IRL2D::Vec lower_domain_fine = {-8, -1.5};
  // const IRL2D::Vec upper_domain_fine = {8, 1.5};
  // BasicMesh mesh_fine = setMesh(Nx_fine, lower_domain_fine,
  // upper_domain_fine); double delta = 2.5 * mesh.dx(); filepath = path +
  // "/pu_implicit_surface.csv"; std::ofstream pu_csvfile(filepath); csvfile <<
  // "x,y\n"; for (int i = mesh_fine.imin(); i <= mesh_fine.imax(); i++) {
  //   for (int j = mesh_fine.jmin(); j <= mesh_fine.jmax(); j++) {
  //     IRL2D::Vec x0(mesh_fine.xm(i), mesh_fine.ym(j));
  //     bool usePlane = false;
  //     IRL2D::Vec x_pu = IRL2D::projectToImplicitSurface(x0, centroids,
  //     normals,
  //                                                       delta, usePlane);
  //     // outputting projected point
  //     pu_csvfile << x_pu[0] << "," << x_pu[1] << "\n";
  //   }
  // }
  // pu_csvfile.close();

  // // partition of unity reconstruction
  // reconstruction_method = "Taubin";
  // getReconstruction(reconstruction_method, liquid_moments, gas_moments, 0.0,
  //                   velU, velV, &interface2);
  // if (rank == 0) {
  //   VTKOutput vtk_io_2("viz_out_2", "viz_2", mesh);
  //   double simulation_time = 0.0;
  //   vtk_io_2.writeVTKInterface(simulation_time, interface2);
  // }

  // // amr vtk outputs
  // std::vector<const Cell*> leaf_cells, mixed_leaf_cells;
  // for (const auto& root : grid) {
  //   collect_leaf_cells(*root, leaf_cells, mixed_leaf_cells);
  // }
  // if (rank == 0) {
  //   std::string path =
  //       "/home/parinht2/Desktop/partition_of_unity/2d_interface/";
  //   write_vtr(path + "vf.vtr", vf, mesh);
  //   write_vtu(path + "amr.vtu", leaf_cells);
  //   write_vtu(path + "mixed_amr.vtu", mixed_leaf_cells);
  // }

  MPI_Finalize();

  return 0;
}