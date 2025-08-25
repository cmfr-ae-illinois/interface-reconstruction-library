#include <cmath>
#include <iostream>

#include "examples/pu_analysis/surfaces.h"

const double A = 0.8, w = 1.0;

double F_sine(const double& x, const double& y) {
  return y - A * std::sin(w * x);
}

Eigen::Vector2d gradF_sine(const double& x, const double& y) {
  double Fx = -A * w * std::cos(w * x);
  double Fy = 1.0;
  return Eigen::Vector2d(Fx, Fy);
}

Eigen::Matrix2d hessF_sine(const double& x, const double& y) {
  Eigen::Matrix2d H;
  H(0, 0) = A * w * w * std::sin(w * x);
  H(0, 1) = 0.0;
  H(1, 0) = H(0, 1);
  H(1, 1) = 0.0;
  return H;
}
