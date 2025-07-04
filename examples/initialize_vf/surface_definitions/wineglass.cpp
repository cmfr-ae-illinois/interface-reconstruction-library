#include <iostream>

#include "examples/initialize_vf/surfaces.h"

double F_wineglass(const double& x, const double& y, const double& z){
  return x * x + y * y - std::log(z + 3.2) * std::log(z + 3.2) - 0.02;
}

Eigen::Vector3d gradF_wineglass(const double& x, const double& y, const double& z){
  double Fx = 2.0 * x;
  double Fy = 2.0 * y;
  double Fz = (-2*std::log(3.2 + z))/(3.2 + z);
  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Matrix3d hessF_wineglass(const double& x, const double& y, const double& z){
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H(0,0) = 2.0;
  H(1,1) = 2.0;
  H(2,2) = -2/std::pow(3.2 + z,2) + (2*std::log(3.2 + z))/std::pow(3.2 + z,2);
  return H;
}