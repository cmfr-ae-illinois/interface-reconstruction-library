#include <cmath>
#include <iostream>

#include "examples/pu_analysis/surfaces.h"

const double A = 0.9, k = 20.0, L = 8.0;

double F_sine(const double& x, const double& y) {
  return y - A * std::sin(k * M_PI / L * x);
}

Eigen::Vector2d gradF_sine(const double& x, const double& y) {
  double Fx = -A * k * M_PI / L * std::cos(k * M_PI / L * x);
  double Fy = 1.0;
  return Eigen::Vector2d(Fx, Fy);
}

Eigen::Matrix2d hessF_sine(const double& x, const double& y) {
  Eigen::Matrix2d H;
  H(0, 0) = A * k * M_PI / L * k * M_PI / L * std::sin(k * M_PI / L * x);
  H(0, 1) = 0.0;
  H(1, 0) = H(0, 1);
  H(1, 1) = 0.0;
  return H;
}
