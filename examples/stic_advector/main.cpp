#include <chrono>
#include <iostream>
#include <string>

#include "examples/stic_advector/circle_rotation_2d.h"
#include "examples/stic_advector/deformation_2d.h"
#include "examples/stic_advector/reconstruction_types.h"
#include "examples/stic_advector/solver.h"
#include "examples/stic_advector/vof_advection.h"

static int startSimulation(const std::string& a_simulation_type,
                           const double a_time_step_size,
                           const double a_time_duration,
                           const int a_viz_frequency);

int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cout << "Incorrect amount of command line arguments supplied. \n";
    std::cout << "Arguments should be \n";
    std::cout
        << "Simulation to run. Options: Deformation2D, CircleRotation2D\n";
    std::cout << "Time step size, dt (double)\n";
    std::cout << "Simulation duration(double)\n";
    std::cout
        << "Amount of time steps between visualization output (integer)\n";
    std::exit(-1);
  }

  std::string simulation_type = argv[1];
  double time_step_size = std::stod(argv[2]);
  double time_duration = std::stod(argv[3]);
  int viz_frequency = atoi(argv[4]);

  auto start = std::chrono::system_clock::now();
  startSimulation(simulation_type, time_step_size, time_duration,
                  viz_frequency);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = end - start;
  printf("Total run time: %20f \n\n", runtime.count());

  return 0;
}

static int startSimulation(const std::string& a_simulation_type,
                           const double a_time_step_size,
                           const double a_time_duration,
                           const int a_viz_frequency) {
  if (a_simulation_type == "CircleRotation2D") {
    return runSimulation<CircleRotation2D>(a_time_step_size, a_time_duration,
                                           a_viz_frequency);
  } else if (a_simulation_type == "Deformation2D") {
    return runSimulation<Deformation2D>(a_time_step_size, a_time_duration,
                                        a_viz_frequency);
  } else {
    std::cout << "Unknown simulation type of : " << a_simulation_type << '\n';
    std::cout << "Value entries are: CircleRotation2D, Deformation2D. \n";
    std::exit(-1);
  }
  return -1;
}
