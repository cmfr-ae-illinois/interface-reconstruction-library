#include <chrono>
#include <iostream>
#include <string>
#include <cmath>

#include "examples/test/deformation_2d.h"
#include "examples/test/oscillation_2d.h"
#include "examples/test/reconstruction_types.h"
#include "examples/test/rotation_2d.h"
#include "examples/test/solver.h"
#include "examples/test/vof_advection.h"
#include "examples/test/basic_mesh.h"
#include "examples/test/data.h"
#include "examples/test/irl2d.h"
#include "examples/test/vtk.h"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
//#include <unsupported/Eigen/LevenbergMarquardt>

double computeVolumeFractionError(const IRL2D::BezierList& cell, 
                                  IRL2D::Parabola& interface, 
                                  const double f_star, 
                                  const double cell_area) {
    IRL2D::Moments moments = IRL2D::ComputeMoments(cell, interface);
    double f_h = moments.m0() / cell_area;
    return std::abs(f_h - f_star);
}

double computeVolumeFraction(const IRL2D::BezierList& cell, IRL2D::Parabola& interface){
  IRL2D::Moments moments = IRL2D::ComputeMoments(cell, interface);
  return moments.m0(); // assuming cell area is 1 (current case)
}

IRL2D::Vec computeCentroid(const IRL2D::BezierList& cell, IRL2D::Parabola& interface){
  IRL2D::Moments moments = IRL2D::ComputeMoments(cell, interface);
  IRL2D::Vec centroid;
  centroid = moments.m1()/moments.m0();
  return centroid;
}

IRL2D::Vec rotateVector(const IRL2D::Vec& vector, double angle){
  IRL2D::Vec newVector;
  angle = angle * M_PI / 180.0; // converting to radians

  newVector.x() = vector.x()*std::cos(angle) - vector.y()*std::sin(angle);
  newVector.y() = vector.x()*std::sin(angle) + vector.y()*std::cos(angle);

  return newVector;
}

// double bisectionMethod(const IRL2D::BezierList& cell, IRL2D::Parabola& interface,
//                        double f_star, double yMin, double yMax, double tol, double maxIter){

//   for (int iter = 0; iter < maxIter; ++iter){
//     double yMid = (yMin + yMax) / 2.0;
//     // vol frac using mid
//     interface.datum()[1] = yMid;
//     double fMid = computeVolumeFraction(cell, interface) - f_star;
//     if (std::abs(fMid) < tol) {
//       std::cout << "Converged after " << iter + 1 << " iterations.\n";
//       return yMid;
//     }
//     // lower bound
//     interface.datum()[1] = yMin;
//     double fMin = computeVolumeFraction(cell, interface) - f_star;
//     // interval update
//     if (fMin * fMid < 0){
//       yMax = yMid;
//     } else {
//       yMin = yMid;
//     }
//   }

//   std::cerr << "Max iterations reached without convergence.\n";
//   return (yMin + yMax) / 2.0; // estimate

// }

// Functor for LM algorithm
struct MOFFunctor{
  const IRL2D::BezierList& cell;
  IRL2D::Parabola& interface;
  const IRL2D::Vec centroid_star;
  const double f_star;

  typedef double Scalar;
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;

  enum{
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  MOFFunctor(const IRL2D::BezierList& cell, IRL2D::Parabola& interface,
                const IRL2D::Vec& centroid_star, const double f_star)
    : cell(cell), interface(interface), centroid_star(centroid_star), f_star(f_star) {
  }
  
  int operator()(const Eigen::VectorXd& params, Eigen::VectorXd& residuals) const {

    Eigen::Matrix2d rotationMatrix;
    double angle = params(0);  // angle in rad
    rotationMatrix << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);

    IRL2D::Parabola temp_interface = interface;

    // rotating reference frame
    Eigen::Vector2d normal (temp_interface.frame()[1][0], temp_interface.frame()[1][1]);
    Eigen::Vector2d tangent (temp_interface.frame()[0][0], temp_interface.frame()[0][1]);
    normal = rotationMatrix * normal;
    tangent = rotationMatrix * tangent;

    // new reference frame
    temp_interface.frame()[0][0] = tangent(0);
    temp_interface.frame()[0][1] = tangent(1);
    temp_interface.frame()[1][0] = normal(0);
    temp_interface.frame()[1][1] = normal(1);

    // volume fraction matching
    temp_interface = IRL2D::MatchToVolumeFraction(cell, temp_interface, f_star);

    // new centroid
    IRL2D::Moments moments = IRL2D::ComputeMoments(cell, temp_interface);
    IRL2D::Vec centroid_h = moments.m1() / moments.m0();

    // residual
    residuals(0) = centroid_h[0] - centroid_star[0];
    residuals(1) = centroid_h[1] - centroid_star[1];
    
    return 0;
  }

  int inputs() const { return 1; }
  int values() const { return 2; }

};

void RotateReferenceFrame(IRL2D::Parabola& interface, double angle){
  
  angle = angle * M_PI / 180.0; // converting to radians

  // Rotation Matrix
  auto RotationMatrix = IRL2D::Mat(IRL2D::Vec(std::cos(angle),-std::sin(angle)),
                                   IRL2D::Vec(std::sin(angle),std::cos(angle)));

  // Rotating Reference Frame
  auto tangent = IRL2D::Vec(interface.frame()[0][0], interface.frame()[0][1]);
  auto normal = IRL2D::Vec(interface.frame()[1][0], interface.frame()[1][1]);
  normal = RotationMatrix * normal;
  tangent = RotationMatrix * tangent;

  // new reference frame
  interface.frame()[0][0] = tangent[0]; 
  interface.frame()[0][1] = tangent[1];
  interface.frame()[1][0] = normal[0];
  interface.frame()[1][1] = normal[1];
}


int main(){
  const auto vertex_1 = IRL2D::Vec(0.0, 0.0); 
  const auto vertex_2 = IRL2D::Vec(1.0, 1.0);
  IRL2D::BezierList rectangle = IRL2D::RectangleFromBounds(vertex_1, vertex_2);

  IRL2D::Parabola interface;
  interface.coeff() = 0;
  interface.datum() = IRL2D::Vec(0.6,0.6);
  interface.frame() = IRL2D::Mat(IRL2D::Vec(1.0,0.0),IRL2D::Vec(0.0,1.0));
  RotateReferenceFrame(interface, -10.0);
  IRL2D::Moments moments;
  moments = IRL2D::ComputeMoments(rectangle, interface);

  std::cout << "Initial Volume Fraction = " << moments.m0() << std::endl;
  std::cout << "Initial Centroid = " << moments.m1()/moments.m0() << std::endl;


  // std::cout << "Testing rotation of reference frame by 45 degrees" << std::endl;
  // RotateReferenceFrame(interface, 45.0);
  // std::cout << interface.frame() << std::endl;

  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Finding target values" << std::endl;
  std::cout << "-----------------------------------" << std::endl;

  IRL2D::Parabola desired_interface;
  desired_interface.coeff() = 0;
  desired_interface.datum() = IRL2D::Vec(0.75, 0.75);
  desired_interface.frame() = IRL2D::Mat(IRL2D::Vec(1.0,0.0), IRL2D::Vec(0.0,1.0));
  RotateReferenceFrame(desired_interface, -45.0);

  IRL2D::Moments desired_moments;
  desired_moments = IRL2D::ComputeMoments(rectangle, desired_interface);

  double const f_star = desired_moments.m0(); // desired volume fraction
  IRL2D::Vec const centroid_star = desired_moments.m1()/desired_moments.m0();

  std::cout << "Desired Volume Fraction = " << desired_moments.m0() << std::endl;
  std::cout << "Desired Centroid = " << desired_moments.m1()/desired_moments.m0() << std::endl;
  // std::cout << "Desired Datum = " << desired_interface.datum() << std::endl;
  // std::cout << "Desired Reference Frame = " << desired_interface.frame() << std::endl;

  // LM algorithm
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Matching Centroid and Volume fraction" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  Eigen::VectorXd initialGuess(1);
  initialGuess(0) = 0.0;

  MOFFunctor errorFunction(rectangle, interface, centroid_star, f_star);
  Eigen::NumericalDiff<MOFFunctor> numericalDiffMyFunctor(errorFunction);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MOFFunctor>, double> lm(numericalDiffMyFunctor);
  lm.parameters.ftol = 1.0e-15;
  lm.parameters.xtol = 1.0e-15;
  lm.parameters.maxfev = 1000;
  //lm.parameters.iterations = 10000;
  Eigen::VectorXd angle = initialGuess;
  // lm.minimize(angle);
  // std::cout << "Angle of rotation in degrees = " << angle[0] * 180.0/M_PI << std::endl;
  // std::cout << "Angle of rotation in radians = " << angle[0] << std::endl;
  double status = lm.minimize(angle);

  //Convergence status
  // if (status == Eigen::LevenbergMarquardtSpace::ImproperInputParameters) {
  //     std::cerr << "Improper input parameters." << std::endl;
  // } else if (status == Eigen::NoConvergence) {
  //     std::cerr << "LM did not converge." << std::endl;
  // } else {
  //     std::cout << "Centroid matched successfully with angle: " << angle << " radians" << std::endl;
  // }

  // finding new reference frame
  RotateReferenceFrame(interface, angle[0] * 180.0/M_PI);
  interface = IRL2D::MatchToVolumeFraction(rectangle, interface, f_star);
  moments = IRL2D::ComputeMoments(rectangle, interface);  
  std::cout << "Volume Fraction = " << moments.m0() << std::endl;
  std::cout << "Centroid = " << moments.m1()/moments.m0() << std::endl;
  // std::cout << "Datum = " << interface.datum() << std::endl;
  // std::cout << "Reference Frame = " << interface.frame() << std::endl;
  // std::cout << "Alpha = " << interface.coeff() << std::endl;

  

}

// static int startSimulation(const std::string& a_simulation_type,
//                            const std::string& a_advection_method,
//                            const std::string& a_reconstruction_method,
//                            const double a_time_step_size,
//                            const double a_time_duration,
//                            const int a_viz_frequency, const int a_nx);

// int main(int argc, char* argv[]) {
//   if (argc != 8) {
//     std::cout << "Incorrect amount of command line arguments supplied. \n";
//     std::cout << "Arguments should be \n";
//     std::cout << "Simulation to run. Options: Rotation2D, Oscillation2D, Deformation2D\n";
//     std::cout << "Advection method. Options: SemiLagQ, FullLagQ, SemiLagL, FullLagL\n";
//     std::cout << "Reconstruction method. Options: ELVIRA, LVIRA, LVIRAQ\n";
//     std::cout << "Time step size, dt (double)\n";
//     std::cout << "Simulation duration(double)\n";
//     std::cout
//         << "Amount of time steps between visualization output (integer)\n";
//     std::cout << "Number of cells (integer)\n";
//     std::exit(-1);
//   }

// #ifdef USE_MPI
//   MPI_Init(&argc, &argv);
//   int rank, size;
//   MPI_Comm_size(MPI_COMM_WORLD, &size);
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// #endif

//   std::string simulation_type = argv[1];
//   std::string advection_method = argv[2];
//   std::string reconstruction_method = argv[3];
//   double time_step_size = std::stod(argv[4]);
//   double time_duration = std::stod(argv[5]);
//   int viz_frequency = atoi(argv[6]);
//   int Nx = atoi(argv[7]);

//   auto start = std::chrono::system_clock::now();
//   startSimulation(simulation_type, advection_method, reconstruction_method,
//                   time_step_size, time_duration, viz_frequency, Nx);
//   auto end = std::chrono::system_clock::now();
//   std::chrono::duration<double> runtime = end - start;

// #ifdef USE_MPI
//   MPI_Finalize();
//   if (rank == 0) printf("Total run time: %20f \n\n", runtime.count());
// #else
//   printf("Total run time: %20f \n\n", runtime.count());
// #endif

//   return 0;
// }

// static int startSimulation(const std::string& a_simulation_type,
//                            const std::string& a_advection_method,
//                            const std::string& a_reconstruction_method,
//                            const double a_time_step_size,
//                            const double a_time_duration,
//                            const int a_viz_frequency, const int a_nx) {
//   if (a_simulation_type == "Rotation2D") {
//     return runSimulation<Rotation2D>(a_simulation_type, a_advection_method,
//                                      a_reconstruction_method, a_time_step_size,
//                                      a_time_duration, a_viz_frequency, a_nx);
//   } else if (a_simulation_type == "Oscillation2D") {
//     return runSimulation<Oscillation2D>(
//         a_simulation_type, a_advection_method, a_reconstruction_method,
//         a_time_step_size, a_time_duration, a_viz_frequency, a_nx);
//   } else if (a_simulation_type == "Deformation2D") {
//     return runSimulation<Deformation2D>(
//         a_simulation_type, a_advection_method, a_reconstruction_method,
//         a_time_step_size, a_time_duration, a_viz_frequency, a_nx);
//   } else {
//     std::cout << "Unknown simulation type of : " << a_simulation_type << '\n';
//     std::cout
//         << "Value entries are: Rotation2D, Oscillation2D, Deformation2D. \n";
//     std::exit(-1);
//   }
//   return -1;
// }
