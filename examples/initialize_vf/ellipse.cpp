#include <iostream>

#include "examples/initialize_vf/surfaces.h"

const double xc = 0, yc = 0.03125;
const double a = 0.2, b = 0.05;

double F_ellipse(const double& x, const double& y){
  return std::pow((x - xc) / a, 2.0) + std::pow((y - yc) / b, 2.0) - 1.0;
}

Eigen::Vector2d gradF_ellipse(const double& x, const double& y){
  double Fx = 2.0 * (x - xc) / (a * a);
  double Fy = 2.0 * (y - yc) / (b * b);
  return Eigen::Vector2d(Fx, Fy);
}

Eigen::Matrix2d hessF_ellipse(const double& x, const double& y){
  Eigen::Matrix2d H;
  H(0,0) = 2.0 / (a * a);
  H(0,1) = 0.0;
  H(1,0) = H(0,1);
  H(1,1) = 2.0 / (b * b);
  return H;
}