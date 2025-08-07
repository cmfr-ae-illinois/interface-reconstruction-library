#ifndef IRL_WENDLAND_TPP_
#define IRL_WENDLAND_TPP_

#include <limits>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

namespace IRL {
// ============== Wendland Class Functions
// New
double Wendland::computeR(Pt xi, Pt x_eval) {
  Pt dx = x_eval - xi;
  return std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
}

double Wendland::eval(double r, double delta) {
  return (4 * r / delta + 1) * (1 - r / delta) * (1 - r / delta) *
         (1 - r / delta) * (1 - r / delta);
}

double Wendland::firstDer(double r, double delta) {
  return (-20 * r / (delta * delta)) * (1 - r / delta) * (1 - r / delta) *
         (1 - r / delta);
}

double Wendland::secondDer(double r, double delta) {
  return (-20 / (delta * delta)) * (1 - r / delta) * (1 - r / delta) *
         (1 - 4 * r / delta);
}

// Evaluate 1
void Wendland::evaluate(Pt xi, double delta, Pt x_eval, double* retVal) {
  // First, get r
  double r = Wendland::computeR(xi, x_eval);
  // Next Calculate F, the function value
  double F = Wendland::eval(r, delta);
  // Return
  *retVal = F;
}

// Evaluate 2
void Wendland::evaluate(Pt xi, double delta, Pt x_eval,
                        std::pair<double, Eigen::Vector3d>* retVal) {
  // First, get r
  double r = Wendland::computeR(xi, x_eval);

  // Next Calculate F, the function value
  double F = Wendland::eval(r, delta);

  // Next Calculate F',F''
  double Fp = Wendland::firstDer(r, delta);

  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);

  // Now, calculate the Gradient of r
  Eigen::Vector3d gradR = x / r;
  // Calculate Gradient
  Eigen::Vector3d gradF = Fp * gradR;
  // Return
  *retVal = std::make_pair(F, gradF);
}

// Evaluate 3
void Wendland::evaluate(
    Pt xi, double delta, Pt x_eval,
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal) {
  // First, get r
  double r = Wendland::computeR(xi, x_eval);

  // Next Calculate F, the function value
  double F = Wendland::eval(r, delta);

  // Next Calculate F',F''
  double Fp = Wendland::firstDer(r, delta);
  double Fpp = Wendland::secondDer(r, delta);

  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);

  // Now, calculate the Gradient of r
  Eigen::Vector3d gradR = x / r;

  // Finally, calculate the Hessian of r
  Eigen::Matrix3d hessR =
      (Eigen::Matrix3d::Identity() - x * x.transpose() / (r * r)) / r;

  // Calculate Return Values
  Eigen::Vector3d gradF = Fp * gradR;
  Eigen::Matrix3d hessF = Fpp * (gradR * gradR.transpose()) + Fp * hessR;
  *retVal = std::make_tuple(F, gradF, hessF);
}
}  // namespace IRL

#endif
