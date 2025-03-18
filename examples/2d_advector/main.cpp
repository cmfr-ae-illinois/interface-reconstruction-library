#include <chrono>
#include <iostream>
#include <string>

#include "examples/2d_advector/deformation_2d.h"
#include "examples/2d_advector/oscillation_2d.h"
#include "examples/2d_advector/reconstruction_types.h"
#include "examples/2d_advector/rotation_2d.h"
#include "examples/2d_advector/solver.h"
#include "examples/2d_advector/vof_advection.h"

static int startSimulation(const std::string& a_simulation_type,
                           const std::string& a_advection_method,
                           const std::string& a_reconstruction_method,
                           const double a_time_step_size,
                           const double a_time_duration,
                           const int a_viz_frequency, const int a_nx);

int main(int argc, char* argv[]) {
  if (argc != 8) {
    std::cout << "Incorrect amount of command line arguments supplied. \n";
    std::cout << "Arguments should be \n";
    std::cout << "Simulation to run. Options: Rotation2D, Oscillation2D, Deformation2D\n";
    std::cout << "Advection method. Options: SemiLagQ, FullLagQ, SemiLagL, FullLagL, ReSyFullLagL\n";
    std::cout << "Reconstruction method. Options: ELVIRA, LVIRA, LVIRAQ, MOF1, MOF2, MOF2AL\n";
    std::cout << "Time step size, dt (double)\n";
    std::cout << "Simulation duration(double)\n";
    std::cout
        << "Amount of time steps between visualization output (integer)\n";
    std::cout << "Number of cells (integer)\n";
    std::exit(-1);
  }

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  std::string simulation_type = argv[1];
  std::string advection_method = argv[2];
  std::string reconstruction_method = argv[3];
  double time_step_size = std::stod(argv[4]);
  double time_duration = std::stod(argv[5]);
  int viz_frequency = atoi(argv[6]);
  int Nx = atoi(argv[7]);

  auto start = std::chrono::system_clock::now();
  startSimulation(simulation_type, advection_method, reconstruction_method,
                  time_step_size, time_duration, viz_frequency, Nx);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = end - start;

#ifdef USE_MPI
  MPI_Finalize();
  if (rank == 0) printf("Total run time: %20f \n\n", runtime.count());
#else
  printf("Total run time: %20f \n\n", runtime.count());
#endif

  return 0;
}

static int startSimulation(const std::string& a_simulation_type,
                           const std::string& a_advection_method,
                           const std::string& a_reconstruction_method,
                           const double a_time_step_size,
                           const double a_time_duration,
                           const int a_viz_frequency, const int a_nx) {
  if (a_simulation_type == "Rotation2D") {
    return runSimulation<Rotation2D>(a_simulation_type, a_advection_method,
                                     a_reconstruction_method, a_time_step_size,
                                     a_time_duration, a_viz_frequency, a_nx);
  } else if (a_simulation_type == "Oscillation2D") {
    return runSimulation<Oscillation2D>(
        a_simulation_type, a_advection_method, a_reconstruction_method,
        a_time_step_size, a_time_duration, a_viz_frequency, a_nx);
  } else if (a_simulation_type == "Deformation2D") {
    return runSimulation<Deformation2D>(
        a_simulation_type, a_advection_method, a_reconstruction_method,
        a_time_step_size, a_time_duration, a_viz_frequency, a_nx);
  } else {
    std::cout << "Unknown simulation type of : " << a_simulation_type << '\n';
    std::cout
        << "Value entries are: Rotation2D, Oscillation2D, Deformation2D. \n";
    std::exit(-1);
  }
  return -1;
}

// int main(){
//   // unit cell function test
// const auto cell = IRL2D::RectangleFromBounds(IRL2D::Vec(-0.25,-0.25), IRL2D::Vec(0.25,0.25));
// IRL2D::Print(cell);
// const auto moments = IRL2D::ComputeMoments(cell);
// std::cout << "Cell M0 = " << moments.m0() << std::endl;
// std::cout << "Cell M1 = " << moments.m1() << std::endl;
// std::cout << "Cell M2 = " << moments.m2() << std::endl;
// auto results = IRL2D::ComputeUnitSquareScaledMoments(cell, moments);
// auto unit_cell = results.first; auto scaled_moments = results.second;
// auto unit_moments = IRL2D::ComputeMoments(unit_cell);
// IRL2D::Print(unit_cell);
// // // testing recentering of moments
// // RecenterMoments(&unit_moments, unit_moments.m1()/unit_moments.m0());
// std::cout << "Unit Cell M0 = " << unit_moments.m0() << std::endl;
// std::cout << "Unit Cell M1 = " << unit_moments.m1() << std::endl;
// std::cout << "Unit Cell M2 = " << unit_moments.m2() << std::endl;
// RecenterMoments(&scaled_moments, scaled_moments.m1()/scaled_moments.m0());
// std::cout << "Scaled M0 = " << scaled_moments.m0() << std::endl;
// std::cout << "Scaled M1 = " << scaled_moments.m1() << std::endl;
// std::cout << "Scaled M2 = " << scaled_moments.m2() << std::endl;
// }