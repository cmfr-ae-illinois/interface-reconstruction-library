#ifndef IRL_WENDLAND4_TPP
#define IRL_WENDLAND4_TPP_

#include <limits>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

namespace IRL {
// ============== Wendland4 Class Functions
inline double Wendland4::computeR(Pt xi, Pt x_eval) {
  Pt dx = x_eval - xi;
  return std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
}

inline double Wendland4::eval(double r, double delta) {
  double rhat = r / delta;
  if (rhat <= 1.0) {
    return (35.0 * rhat * rhat + 18.0 * rhat + 3.0) * (1.0 - rhat) *
           (1.0 - rhat) * (1.0 - rhat) * (1.0 - rhat) * (1.0 - rhat) *
           (1.0 - rhat);
  } else {
    return 0.0;
  }
}

inline Eigen::Vector3d Wendland4::getGradient(const Pt& xi, const double& delta,
                                              const Pt& x_eval) {
  // First, get r
  double r = Wendland4::computeR(xi, x_eval);
  if (r > delta) {
    return Eigen::Vector3d::Zero();
  }
  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);
  // Apply Formula
  Eigen::Vector3d gradF = (-56.0 * x / (delta * delta)) * (1 - r / delta) *
                          (1 - r / delta) * (1 - r / delta) * (1 - r / delta) *
                          (1 - r / delta) * (1 + 5.0 * r / delta);
  return gradF;
}

inline Eigen::Matrix3d Wendland4::getHessian(const Pt& xi, const double& delta,
                                             const Pt& x_eval) {
  // First, get r
  double r = Wendland4::computeR(xi, x_eval);
  if (r > delta) {
    return Eigen::Matrix3d::Zero();
  }
  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);

  // Calculate Return Values
  double diagonal1 = (-56.0 / (delta * delta)) * (1.0 - r / delta) *
                     (1.0 - r / delta) * (1.0 - r / delta) * (1.0 - r / delta) *
                     (1.0 - r / delta) * (1.0 + 5.0 * r / delta);
  double otherTerms = (1680.0 / (delta * delta * delta * delta)) *
                      (1 - r / delta) * (1 - r / delta) * (1 - r / delta) *
                      (1 - r / delta);
  Eigen::Matrix3d hessF =
      diagonal1 * Eigen::Matrix3d::Identity() + otherTerms * x * x.transpose();
  return hessF;
}

// Evaluate 1
inline void Wendland4::evaluate(const Pt& xi, const double& delta,
                                const Pt& x_eval, double* retVal) {
  // First, get r
  double r = Wendland4::computeR(xi, x_eval);
  // Next Calculate F, the function value
  double F = Wendland4::eval(r, delta);
  // Return
  *retVal = F;
}

// Evaluate 2
inline void Wendland4::evaluate(const Pt& xi, const double& delta,
                                const Pt& x_eval,
                                std::pair<double, Eigen::Vector3d>* retVal) {
  // First, get r
  double r = Wendland4::computeR(xi, x_eval);
  // Next Calculate F, the function value
  double F = Wendland4::eval(r, delta);
  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);

  // Calculate Gradient
  Eigen::Vector3d gradF = getGradient(xi, delta, x_eval);
  // Return
  *retVal = std::make_pair(F, gradF);
}

// Evaluate 3
inline void Wendland4::evaluate(
    const Pt& xi, const double& delta, const Pt& x_eval,
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal) {
  // First, get r
  double r = Wendland4::computeR(xi, x_eval);

  // Next Calculate F, the function value
  double F = Wendland4::eval(r, delta);

  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);

  // Calculate Return Values
  Eigen::Vector3d gradF = getGradient(xi, delta, x_eval);
  Eigen::Matrix3d hessF = getHessian(xi, delta, x_eval);
  *retVal = std::make_tuple(F, gradF, hessF);
}
}  // namespace IRL

#endif
