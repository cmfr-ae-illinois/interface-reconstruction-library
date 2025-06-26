#include <iostream>

#include "examples/initialize_vf/surfaces.h"

const double xc = 0.0, yc = 0.0, zc = 0.0;
const double a = 1.0, b = 1.0, c = 1.0;

double F_ellipsoid(const double& x, const double& y, const double& z){
  return std::pow((x - xc) / a, 2.0) + std::pow((y - yc) / b, 2.0) +
         std::pow((z - zc) / c, 2.0) - 1.0;
}

Eigen::Vector3d gradF_ellipsoid(const double& x, const double& y, const double& z){
  double Fx = 2.0 * (x - xc) / (a * a);
  double Fy = 2.0 * (y - yc) / (b * b);
  double Fz = 2.0 * (z - zc) / (c * c);
  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Matrix3d hessF_ellipsoid(const double& x, const double& y, const double& z){
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H(0,0) = 2.0 / (a * a);
  H(1,1) = 2.0 / (b * b);
  H(2,2) = 2.0 / (c * c);
  return H;
}