#ifndef IRL_WU_H_
#define IRL_WU_H_

namespace IRL {
// Here we are going to make a static class for the WU function
// This is because instead of having to make different objects for different
// detla, x_eval values, we will just use the WU class functions.
// This means we will need to pass in delta,x_eval as parameters, but that
// should not be a problem maybe?
class Wu {
  // Note that for these functions, xi is the center, delta is the radius,
  // and x_eval is the location we are evaluating our function.
 public:
  // New
  // Compute the Radius
  static inline double computeR(Pt xi, Pt x_eval);  // No Update Needed
  // Compute Zeroth Derivative (Function)
  static inline double eval(double r, double delta);  // Update Needed
  // Compute WU Function
  static inline void evaluate(const Pt& xi, const double& delta,
                              const Pt& x_eval, double* retVal);
  // Compute WU Function and Gradient
  static inline void evaluate(const Pt& xi, const double& delta,
                              const Pt& x_eval,
                              std::pair<double, Eigen::Vector3d>* retVal);
  // Compute WU Function, Grad, And Hessian
  static inline void evaluate(
      const Pt& xi, const double& delta, const Pt& x_eval,
      std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal);
  // Specialized Evaluate Function for when we are evaluating near the center
  static inline Eigen::Vector3d getGradient(const Pt& xi, const double& delta,
                                            const Pt& x_eval);  // Update Needed
  static inline Eigen::Matrix3d getHessian(const Pt& xi, const double& delta,
                                           const Pt& x_eval);  // Update Needed
  // Disallow Instance Creation
  Wu() = delete;
};
}  // namespace IRL

#include "examples/level_set_reconstruction/weights/wu.tpp"

#endif