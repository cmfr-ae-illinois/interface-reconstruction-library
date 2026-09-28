#ifndef IRL_WU4_TPP
#define IRL_WU4_TPP_

#include <limits>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

namespace IRL {
// ============== Wu4 Class Functions
inline double Wu4::computeR(Pt xi, Pt x_eval) {
  Pt dx = x_eval - xi;
  return std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
}

inline double Wu4::eval(double r, double delta) {
  double rhat = r / delta;
  if (rhat <= 1.0) {
    return (5.0 * rhat * rhat * rhat * rhat * rhat +
            30.0 * rhat * rhat * rhat * rhat + 72.0 * rhat * rhat * rhat +
            82.0 * rhat * rhat + 36.0 * rhat + 6.0) *
           (1.0 - rhat) * (1.0 - rhat) * (1.0 - rhat) * (1.0 - rhat) *
           (1.0 - rhat) * (1.0 - rhat);
  } else {
    return 0.0;
  }
}

inline Eigen::Vector3d Wu4::getGradient(const Pt& xi, const double& delta,
                                        const Pt& x_eval) {
  // First, get r
  double r = Wu4::computeR(xi, x_eval);
  if (r > delta) {
    return Eigen::Vector3d::Zero();
  }
  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);
  // Apply Formula
  Eigen::Vector3d gradF =
      (-11.0 * x / (delta * delta)) * (2.0 + r / delta) *
      (4.0 + 5.0 * (r / delta) * (r / delta) * (r / delta) +
       15.0 * (r / delta) * (r / delta) + 18.0 * (r / delta)) *
      (1.0 - r / delta) * (1.0 - r / delta) * (1.0 - r / delta) *
      (1.0 - r / delta) * (1.0 - r / delta);

  return gradF;
}

inline Eigen::Matrix3d Wu4::getHessian(const Pt& xi, const double& delta,
                                       const Pt& x_eval) {
  // First, get r
  double r = Wu4::computeR(xi, x_eval);
  if (r > delta) {
    return Eigen::Matrix3d::Zero();
  }
  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);

  // Calculate Return Values
  double diagonal1 = (-11.0 / (delta * delta)) * (2.0 + r / delta) *
                     (4.0 + 5.0 * (r / delta) * (r / delta) * (r / delta) +
                      15.0 * (r / delta) * (r / delta) + 18.0 * (r / delta)) *
                     (1.0 - r / delta) * (1.0 - r / delta) * (1.0 - r / delta) *
                     (1.0 - r / delta) * (1.0 - r / delta);
  double otherTerms = (99.0 / (delta * delta * delta * delta)) *
                      (16.0 + 5.0 * (r / delta) * (r / delta) * (r / delta) +
                       20.0 * (r / delta) * (r / delta) + 29.0 * (r / delta)) *
                      (1.0 - r / delta) * (1.0 - r / delta) *
                      (1.0 - r / delta) * (1.0 - r / delta);
  Eigen::Matrix3d hessF =
      diagonal1 * Eigen::Matrix3d::Identity() + otherTerms * x * x.transpose();
  return hessF;
}

// Evaluate 1
inline void Wu4::evaluate(const Pt& xi, const double& delta, const Pt& x_eval,
                          double* retVal) {
  // First, get r
  double r = Wu4::computeR(xi, x_eval);
  // Next Calculate F, the function value
  double F = Wu4::eval(r, delta);
  // Return
  *retVal = F;
}

// Evaluate 2
inline void Wu4::evaluate(const Pt& xi, const double& delta, const Pt& x_eval,
                          std::pair<double, Eigen::Vector3d>* retVal) {
  // First, get r
  double r = Wu4::computeR(xi, x_eval);
  // Next Calculate F, the function value
  double F = Wu4::eval(r, delta);
  // Now, we need to calculate the distance function derivative. To do this,
  // first make x an Eigen Vector.
  Eigen::Vector3d x(x_eval[0] - xi[0], x_eval[1] - xi[1], x_eval[2] - xi[2]);

  // Calculate Gradient
  Eigen::Vector3d gradF = getGradient(xi, delta, x_eval);
  // Return
  *retVal = std::make_pair(F, gradF);
}

// Evaluate 3
inline void Wu4::evaluate(
    const Pt& xi, const double& delta, const Pt& x_eval,
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal) {
  // First, get r
  double r = Wu4::computeR(xi, x_eval);

  // Next Calculate F, the function value
  double F = Wu4::eval(r, delta);

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
