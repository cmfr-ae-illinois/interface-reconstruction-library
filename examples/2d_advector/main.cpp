#include <chrono>
#include <iostream>
#include <string>
#include <fstream>

#include "examples/2d_advector/deformation_2d.h"
#include "examples/2d_advector/oscillation_2d.h"
#include "examples/2d_advector/reconstruction_types.h"
#include "examples/2d_advector/rotation_2d.h"
#include "examples/2d_advector/solver.h"
#include "examples/2d_advector/vof_advection.h"

#include <Eigen/Dense>
#include <Eigen/QR>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

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

// // particle method testing -----------------------------------------------------------

// // structs for circle
// struct circle {
//   double a, b, R;
// };

// struct circleParams{
//   Eigen::Matrix4d M, B;
//   double A, D, theta;
// };

// // AF2 initialization (Pratt's approximation to GRAF)
// circleParams AF2_Initializer(const std::vector<IRL2D::Vec>& points){
//   int n = points.size();
//   Eigen::MatrixXd M(n,4);
//   for (int i = 0; i < n; i++){
//     double xi = points[i].x(), yi = points[i].y();
//     double zi = xi * xi + yi * yi;
//     M(i, 0) = zi;
//     M(i, 1) = xi;
//     M(i, 2) = yi;
//     M(i, 3) = 1.0;
//   }
//   Eigen::Matrix4d moment = M.transpose() * M;

//   Eigen::Matrix4d B; // constraint matrix
//   B << 0, 0, 0, -2,
//        0, 1, 0,  0,
//        0, 0, 1,  0,
//       -2, 0, 0,  0;
  
//   Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
//   ges.compute(moment, B);
//   // std::cout << "Generalized eigenvalues: \n " << ges.eigenvalues().transpose() << std::endl;
//   // std::cout << "Generalized eigenvectors: \n" << ges.eigenvectors() << std::endl;
  
//   double A = ges.eigenvectors().col(1)[0].real() ;
//   double B_ = ges.eigenvectors().col(1)[1].real();
//   double C = ges.eigenvectors().col(1)[2].real() ;
//   double D = ges.eigenvectors().col(1)[3].real() ;
  
//   double val = 1.0 / std::sqrt( B_ * B_ + C * C - 4.0 * A * D );
//   A = A * val; B_ = B_*val; C = C*val; D = D*val;
//   double theta = atan2(C, B_);

//   return {moment, B, A, D, theta};
// }

// // method 2: solving quartic equation
// double f(const Eigen::Matrix4d& M, const Eigen::Matrix4d& B, double eta) {
//   Eigen::Matrix4d A = M - eta * B;
//   return A.determinant();
// }

// double df(const Eigen::Matrix4d& M, const Eigen::Matrix4d& B, double eta, 
//           double h = 1e-6) {
//   return (f(M, B, eta + h) - f(M, B, eta - h)) / (2.0 * h);
// }

// double solve_eta(const Eigen::Matrix4d& M, const Eigen::Matrix4d& B, 
//                  double eta0, double tol = 1e-10, int max_iter = 100) {
//   double eta = eta0;
//   for (int i = 0; i < max_iter; ++i) {
//       double f_val = f(M, B, eta);
//       double df_val = df(M, B, eta);
//       if (std::abs(df_val) < 1e-12) {
//           std::cerr << "Derivative near zero. Stopping." << std::endl;
//           break;
//       }
//       double eta_new = eta - f_val / df_val;
//       if (std::abs(eta_new - eta) < tol) {
//           return eta_new;
//       }
//       eta = eta_new;
//   }
//   return eta;
// }

// Eigen::Vector4d compute_eigenvector(const Eigen::Matrix4d& M, const Eigen::Matrix4d& B,
//                                     const double& eta){
//   Eigen::Matrix4d A = M - eta * B;
//   Eigen::FullPivLU<Eigen::Matrix4d> lu_decomp(A);
//   Eigen::MatrixXd null_space = lu_decomp.kernel();
//   return null_space.col(0).normalized();
// }

// // functor for LMA
// struct LMAFunctor{
//   typedef double Scalar;
//   typedef Eigen::VectorXd InputType;
//   typedef Eigen::VectorXd ValueType;
//   typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
//   enum{
//     InputsAtCompileTime = Eigen::Dynamic,
//     ValuesAtCompileTime = Eigen::Dynamic
//   };

//   // variables
//   const std::vector<IRL2D::Vec>& points;

//   // constructor
//   LMAFunctor(const std::vector<IRL2D::Vec>& pts) : points(pts) {}

//   int inputs() const {return 3;} // A, D, theta
//   int values() const {return static_cast<int>(points.size()); }

//   int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
//     double A = x(0), D = x(1), theta = x(2);
//     double E = std::sqrt(1 + 4 * A * D);
//     double cos_theta = std::cos(theta);
//     double sin_theta = std::sin(theta);

//     for (int i = 0; i < values(); ++i) {
//         double xi = points[i].x(), yi = points[i].y();
//         double zi = xi * xi + yi * yi;
//         double ui = xi * cos_theta + yi * sin_theta;
//         double Pi = A * zi + E * ui + D;
//         double Qi = std::sqrt(1.0 + 4.0 * A * Pi);
//         fvec(i) = 2.0 * Pi / (1.0 + Qi);
//     }
//     return 0;
//   }

//   // analytical Jacobian
//   int df(const InputType& x, JacobianType& fjac) const {
//     double A = x(0), D = x(1), theta = x(2);
//     double E = std::sqrt(1.0 + 4.0 * A * D);
//     double cos_theta = std::cos(theta);
//     double sin_theta = std::sin(theta);

//     for (int i = 0; i < values(); ++i) {
//         double xi = points[i].x(), yi = points[i].y();
//         double zi = xi * xi + yi * yi;
//         double ui = xi * cos_theta + yi * sin_theta;
//         double Pi = A * zi + E * ui + D;
//         double Qi = std::sqrt(1 + 4 * A * Pi);
//         double di = 2.0 * Pi / (1.0 + Qi);
//         double Ri = 2.0 * (1.0 - A * di / Qi) / (Qi + 1.0);

//         fjac(i, 0) = (zi + 2.0 * D * ui / E) * Ri - (di * di / Qi); // ∂d/∂A
//         fjac(i, 1) = (2.0 * A * ui / E + 1.0) * Ri;                 // ∂d/∂D
//         fjac(i, 2) = (-xi * sin_theta + yi * cos_theta) * E * Ri;   // ∂d/∂θ
//     }
//     return 0;
//   }

// };

// circle fit_with_LMA_AF2(const std::vector<IRL2D::Vec>& points){

//   using Scalar = double;
//   circleParams init = AF2_Initializer(points);
//   Eigen::VectorXd x(3);
//   x << init.A, init.D, init.theta;

//   LMAFunctor functor(points);
//   Eigen::LevenbergMarquardt<LMAFunctor> lm(functor);
//   lm.parameters.maxfev = 200;
//   lm.parameters.xtol = 1e-12;
//   lm.minimize(x);

//   double A = x(0), D = x(1), theta = x(2);
//   double E = std::sqrt(1.0 + 4.0 * A * D);
//   double B = E * std::cos(theta);
//   double C = E * std::sin(theta);
//   double a = -B / (2.0 * A);
//   double b = -C / (2.0 * A);
//   double R = std::fabs(1.0 / (2.0 * A));

//   return {a, b, R};
// }

// // AF3 initialization (Taubin's approximation to GRAF)
// circleParams AF3_Initializer(const std::vector<IRL2D::Vec>& points){

//   int n = points.size();

//   // moment matrix
//   Eigen::MatrixXd M_(n,4);
//   for (int i = 0; i < n; i++){
//     double xi = points[i].x(), yi = points[i].y();
//     double zi = xi * xi + yi * yi;
//     M_(i, 0) = zi;
//     M_(i, 1) = xi;
//     M_(i, 2) = yi;
//     M_(i, 3) = 1.0;
//   }
//   Eigen::Matrix4d M = M_.transpose() * M_;

//   // constraint matrix
//   Eigen::Matrix4d C;
//   C.setZero();
//   C(0,0) = 4.0 * M(0,3);
//   C(0,1) = 2.0 * M(1,3);
//   C(0,2) = 2.0 * M(2,3);
//   C(1,0) = C(0,1);
//   C(1,1) = n;
//   C(2,0) = C(0,2);
//   C(2,2) = n;

//   // solving for eta (eigenvalue)
//   auto F = [&](double eta){
//     return (M - eta * C).determinant();
//   };
//   auto dF = [&](double eta){
//     const double h = 1e-6;
//     return (F(eta + h) - F(eta - h)) / (2 * h);
//   };
  
//   // Newtons iterations
//   double eta = 0.0, tol = 1e-10; int max_iter = 100;
//   for (int i = 0; i < max_iter; i++){
//     double f = F(eta);
//     double df = dF(eta);
//     if (std::abs(df) < 1e-12) {
//       std::cerr << "Derivative near zero. Stopping." << std::endl;
//       break;
//     }
//     double eta_new = eta - f / df;
//     if (std::abs(eta_new - eta) < tol) {
//       eta = eta_new;
//       break;
//     }
//     eta = eta_new;
//   }
//   std::cout << "eta (AF3) = " << eta << std::endl;

//   // results with eigenvalue solver
//   Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
//   ges.compute(M, C);

//   // computing eigenvector
//   Eigen::FullPivLU<Eigen::Matrix4d> lu_decomp(M - eta*C);
//   Eigen::MatrixXd null_space = lu_decomp.kernel();
//   auto params = null_space.col(0).normalized(); // [A B C D]
//   double A = params[0], B = params[1], C_ = params[2], D = params[3];

//   // scaling the parameters
//   double constraint = 4.0*A*A*M(0,3) + 4.0*A*B*M(1,3) +
//                       4.0*A*C_*M(2,3) + B*B*static_cast<double>(n) + 
//                       C_*C_*static_cast<double>(n);

//   double scale_factor = 1.0 / std::sqrt(constraint) * std::sqrt(80.0);
//   A *= scale_factor; B *= scale_factor; C_ *= scale_factor; D *= scale_factor;
//   double theta = std::atan2(C_,B);

//   return {M, C, A, D, theta};
// }

// double computeVfracWeight(double vfrac) {
//   const double limit_vfrac = 0.1;
//   if (vfrac < limit_vfrac) {
//       return 0.5 - 0.5 * std::cos(M_PI * vfrac / limit_vfrac);
//   } else if (vfrac > 1.0 - limit_vfrac) {
//       return 0.5 - 0.5 * std::cos(M_PI * (1.0 - vfrac) / limit_vfrac);
//   } else {
//       return 1.0;
//   }
// }

// double computeDistanceWeight(const IRL2D::Vec& pref, const IRL2D::Vec& ploc,
//                              const double& h){
//   double distance = (ploc - pref).magnitude() / h;
//   if (distance < 2.5){
//     return (1.0 + 4.0 * distance / 2.5) * pow(1.0 - distance / 2.5, 4.0);
//   } else {
//     return 0.0;
//   }
// }

// circleParams AF2W_Initializer(const std::vector<IRL2D::Vec>& points,
//                               const std::vector<double>& vfw, 
//                               const std::vector<double>& dw){

//   // Generating matrices for eigenvalue problem
//   int n = points.size();
//   Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
//   for (int i = 0; i < n; i++){
//     double xi = points[i].x(), yi = points[i].y();
//     double zi = xi * xi + yi * yi;
//     double w = vfw[i] * dw[i]; // volume fraction and distance weight
//     Eigen::Vector4d u;
//     u << zi, xi, yi, 1.0;
//     M += w * u * u.transpose();  // weighted outer product
//   }
//   Eigen::Matrix4d B; // constraint matrix
//   B << 0, 0, 0, -2,
//        0, 1, 0,  0,
//        0, 0, 1,  0,
//       -2, 0, 0,  0;
  
//   // solving the generalized eigenvalue problem
//   auto F = [&](double eta){
//     return (M - eta * B).determinant();
//   };
//   auto dF = [&](double eta){
//     const double h = 1e-6;
//     return (F(eta + h) - F(eta - h)) / (2 * h);
//   };

//   // eigenvalue calculation
//   double eta = 0.0, tol = 1e-10; int max_iter = 100;
//   for (int i = 0; i < max_iter; i++){
//     double f = F(eta);
//     double df = dF(eta);
//     if (std::abs(df) < 1e-12) {
//       std::cerr << "Derivative near zero. Stopping." << std::endl;
//       break;
//     }
//     double eta_new = eta - f / df;
//     if (std::abs(eta_new - eta) < tol) {
//       eta = eta_new;
//       break;
//     }
//     eta = eta_new;
//   }
//   std::cout << "eta (AF2W) = " << eta << std::endl;
  
//   // computing eigenvector
//   Eigen::FullPivLU<Eigen::Matrix4d> lu_decomp(M - eta*B);
//   Eigen::MatrixXd null_space = lu_decomp.kernel();
//   auto params = null_space.col(0).normalized(); // [A B C D]
//   double A = params[0], B_ = params[1], C = params[2], D = params[3];
  
//   // scaling parameters
//   double constraint = B_ * B_ + C * C - 4.0 * A * D;
//   double scale_factor = 1.0 / std::sqrt(constraint);
//   A *= scale_factor; B_ *= scale_factor; C *= scale_factor; D *= scale_factor;
//   double theta = std::atan2(C, B_);

//   return {M, B, A, D, theta};
// }

// int main(){

//   // get mesh
//   int Nx = 64;
//   BasicMesh mesh = Deformation2D::setMesh(Nx);

//   // reading in csv containing interface data
//   std::vector<IRL2D::InterfaceEndPoints> data;
//   std::string dir = "/home/parinht2/Documents/testing code/particle_var_force and circle";
//   std::string filepath = dir + "/interface.csv";
//   readCSV(filepath, data);

//   Data<IRL2D::InterfaceEndPoints> plic_data(&mesh);
//   for (auto& d : data){
//     plic_data(d.xIndex,d.yIndex) = d;
//   }

//   // parameters for particle method
//   int N = 19; double Hp = 4.0; double h = mesh.dx(); double eta = 1.0;
//   double hp = Hp * h / (static_cast<double>(N) - 1.0);
//   std::vector<IRL2D::Vec> initial_particle_positions(N), initial_particle_forces(N),
//                           final_particle_positions(N), final_particle_forces(N);
//   std::pair<IRL2D::Vec,IRL2D::Vec> target_endpoints, temp;
//   std::vector<std::pair<IRL2D::Vec, IRL2D::Vec>> line_seg_endpoints;
//   // target cell
//   int i_tar = 38, j_tar = 27;

//   // estimating curvature at a target cell using particle method
//   // for (int i = mesh.imin(); i < mesh.imax(); ++i){
//   //   for (int j = mesh.jmin(); j < mesh.jmax(); ++j){
//   //     if (plic_data(i,j).mixed == true && i == i_tar && j == j_tar){
//   //       target_endpoints = {IRL2D::Vec(plic_data(i,j).ax , plic_data(i,j).ay) ,
//   //                           IRL2D::Vec(plic_data(i,j).bx , plic_data(i,j).by)};
//   //       for (int ii = -2; ii <=2; ++ii){
//   //         for (int jj = -2; jj <=2; ++jj){
//   //           if (plic_data(i + ii,j + jj).mixed == true){
//   //             line_seg_endpoints.push_back({
//   //               IRL2D::Vec(plic_data(i+ii,j+jj).ax , plic_data(i+ii,j+jj).ay) , 
//   //               IRL2D::Vec(plic_data(i+ii,j+jj).bx , plic_data(i+ii,j+jj).by)
//   //             });
//   //           }
//   //         }
//   //       }
//   //       // initial positions
//   //       initial_particle_positions = IRL2D::InitializeParticlePositions(target_endpoints, hp, N);
//   //       // initial forces
//   //       for (int ip = 0; ip < N; ip++){
//   //         initial_particle_forces[ip] = IRL2D::ComputeParticleForce(initial_particle_positions[ip], 
//   //                                                                   line_seg_endpoints, eta);
//   //       }
//   //       // printing initial particle positions and forces
//   //       std::cout << "-------- Initial particle data -----------" << std::endl;
//   //       IRL2D::printParticleData(initial_particle_positions, initial_particle_forces);

//   //       // getting final particle and positions
//   //       IRL2D::particle_pf(initial_particle_positions, initial_particle_forces,
//   //                          final_particle_positions, final_particle_forces,
//   //                          target_endpoints, line_seg_endpoints,
//   //                          N, Hp, h, eta);
//   //       std::cout << "-------- Final particle data -----------" << std::endl;
//   //       IRL2D::printParticleData(final_particle_positions, final_particle_forces);
//   //     }
//   //   }
//   // }

//   // variable eta --------------------------------------------------------------------------
//   std::cout << "----------------Variable eta-------------------" << std::endl;
//   int Ns = 5; // stencil size
//   std::vector<std::vector<IRL2D::InterfaceEndPoints>> plicDataMat(Ns, std::vector<IRL2D::InterfaceEndPoints>(Ns));
//   for (int i = mesh.imin(); i < mesh.imax(); ++i){
//     for (int j = mesh.jmin(); j < mesh.jmax(); ++j){
//       if (plic_data(i,j).mixed == true && i == i_tar && j == j_tar){
//         target_endpoints = {IRL2D::Vec(plic_data(i,j).ax , plic_data(i,j).ay) ,
//                             IRL2D::Vec(plic_data(i,j).bx , plic_data(i,j).by)};
//         for (int ii = -2; ii <=2; ++ii){
//           for (int jj = -2; jj <=2; ++jj){
//             plicDataMat[ii+2][jj+2] = plic_data(i + ii, j + jj);
//             if (plic_data(i + ii,j + jj).mixed == true){
//               line_seg_endpoints.push_back({
//                 IRL2D::Vec(plic_data(i+ii,j+jj).ax , plic_data(i+ii,j+jj).ay) , 
//                 IRL2D::Vec(plic_data(i+ii,j+jj).bx , plic_data(i+ii,j+jj).by)
//               });
//             }
//           }
//         }
//         // initial positions
//         initial_particle_positions = IRL2D::InitializeParticlePositions(target_endpoints, hp, N);
//         // initial forces
//         for (int ip = 0; ip < N; ip++){
//           initial_particle_forces[ip] = IRL2D::ComputeParticleForce(initial_particle_positions[ip], 
//                                                                     line_seg_endpoints, eta);
//         }
//         // final particle positions and forces
//         IRL2D::curvature_vareta(initial_particle_positions, initial_particle_forces,
//                                 final_particle_positions, final_particle_forces,
//                                 target_endpoints, line_seg_endpoints,
//                                 plicDataMat, N, Hp, h);
//         // IRL2D::printParticleData(final_particle_positions, final_particle_forces);
//       }
//     }
//   }

//   // ----------------- Least squares fit of a circle ---------------------
//   std::cout << "-------------------------------------------------------------" << std::endl;
//   // generating points from interfaces
//   std::vector<IRL2D::Vec> points = IRL2D::generatePoints(line_seg_endpoints);

//   // storing data in a csv for plotting
//   dir = "/home/parinht2/Documents/testing code/particle_var_force and circle";
//   filepath = dir + "/points.csv";
//   std::ofstream csvfile(filepath);
//   csvfile << "x,y\n";
//   for (const auto& pt : points){
//     csvfile << pt.x() << "," << pt.y() << "\n";
//   }
//   csvfile.close();

//   // initial guess
//   circleParams circle1Params = AF2_Initializer(points);
//   double A = circle1Params.A, D = circle1Params.D,
//          theta = circle1Params.theta;
//   std::cout << "A = " <<  A << std::endl;
//   std::cout << "D = " << D << std::endl;
//   std::cout << "theta = " << theta << std::endl;

//   // solving quartic equation
//   double eta0 = 0.0;
//   auto M = circle1Params.M, B = circle1Params.B;
//   double eta_root = solve_eta(M, B, eta0);
//   // std::cout << "Eta (newtons method) = " << eta_root << std::endl;
//   Eigen::Vector4d eigenvector = compute_eigenvector(M, B, eta_root);
//   // std::cout << "Eigenvector(newtons method) = \n" << eigenvector << std::endl;

//   // Levenverg-Marquadt optimization
//   circle fitted_circle = fit_with_LMA_AF2(points);
//   std::cout << "Best fit circle: center = (" << fitted_circle.a << ", " << fitted_circle.b << 
//                "), R = " << fitted_circle.R << std::endl;

//   // AF3 initial guess
//   std::cout << "------------------- AF3 ---------------------" << std::endl;
//   circleParams circle2Params = AF3_Initializer(points);
//   A = circle2Params.A; D = circle2Params.D; theta = circle2Params.theta;
//   std::cout << "A = " <<  A << std::endl;
//   std::cout << "D = " << D << std::endl;
//   std::cout << "theta = " << theta << std::endl;

//   // least squares fit using weight -----------------------------------------------------------
//   std::cout << "------------------ Weights ----------------------" << std::endl;
  
//   // reference point for distance calculations
//   IRL2D::Vec pref = IRL2D::Vec( (plic_data(i_tar, j_tar).ax + plic_data(i_tar, j_tar).bx)/2.0,
//                                 (plic_data(i_tar, j_tar).ay + plic_data(i_tar, j_tar).by)/2.0);

//   // getting weight for each point
//   std::vector<double> vfracs, vfw, dw;
//   std::vector<IRL2D::Vec> ploc;
//   for (int ii = -2; ii <=2; ++ii){
//     for (int jj = -2; jj <=2; ++jj){
//       if (plic_data(i_tar + ii, j_tar + jj).mixed == true){
//         vfracs.push_back(plic_data(i_tar + ii, j_tar + jj).vf);
//         ploc.push_back(
//           IRL2D::Vec( (plic_data(i_tar+ii, j_tar+jj).ax + plic_data(i_tar+ii, j_tar+jj).bx) / 2.0,
//                       (plic_data(i_tar+ii, j_tar+jj).ay + plic_data(i_tar+ii, j_tar+jj).by) / 2.0 )
//                       );
//       }
//     }
//   }
//   for (int i = 0; i < line_seg_endpoints.size(); i++){
//     double vf_weight = computeVfracWeight(vfracs[i]);
//     double d_weight = computeDistanceWeight(pref, ploc[i], mesh.dx());
//     for (int j = 0; j < (points.size()/line_seg_endpoints.size()); j++){
//       vfw.push_back(vf_weight);
//       dw.push_back(d_weight);
//     }
//   }
//   // for (auto& i : dw){
//   //   std::cout << i << std::endl;
//   // }

//   // circle parameters
//   circleParams circle3Params = AF2W_Initializer(points, vfw, dw);
//   A = circle3Params.A; D = circle3Params.D; theta = circle3Params.theta;
//   std::cout << "A = " <<  A << std::endl;
//   std::cout << "D = " << D << std::endl;
//   std::cout << "theta = " << theta << std::endl;

//   return 0;
// }