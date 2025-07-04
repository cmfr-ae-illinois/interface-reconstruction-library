#include <iostream>

#include "examples/initialize_vf/surfaces.h"

double F_genus(const double& x, const double& y, const double& z){
  return 2.0 * y * (y * y - 3.0 * x * x) * (1.0 - z * z)
         + (x * x + y * y) * (x * x + y * y)
         - (9.0 * z * z - 1.0) * (1.0 - z * z);
}

Eigen::Vector3d gradF_genus(const double& x, const double& y, const double& z){
  double Fx = 4*x*(std::pow(x,2) + std::pow(y,2)) - 12*x*y*(1 - std::pow(z,2));
  double Fy = 4*y*(std::pow(x,2) + std::pow(y,2)) + 4*std::pow(y,2)*(1 - std::pow(z,2)) + 2*(-3*std::pow(x,2) + std::pow(y,2))*(1 - std::pow(z,2));
  double Fz = -4*y*(-3*std::pow(x,2) + std::pow(y,2))*z - 18*z*(1 - std::pow(z,2)) + 2*z*(-1 + 9*std::pow(z,2));
  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Matrix3d hessF_genus(const double& x, const double& y, const double& z){
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H(0,0) = 8*std::pow(x,2) + 4*(std::pow(x,2) + std::pow(y,2)) - 12*y*(1 - std::pow(z,2));
  H(1,1) = 8*std::pow(y,2) + 4*(std::pow(x,2) + std::pow(y,2)) + 12*y*(1 - std::pow(z,2));
  H(2,2) = -4*y*(-3*std::pow(x,2) + std::pow(y,2)) + 72*std::pow(z,2) - 18*(1 - std::pow(z,2)) + 2*(-1 + 9*std::pow(z,2));

  H(0,1) = 8*x*y - 12*x*(1 - std::pow(z,2));
  H(1,0) = H(0,1);
  H(0,2) = 24*x*y*z;
  H(2,0) = H(0,2);
  H(1,2) = -8*std::pow(y,2)*z - 4*(-3*std::pow(x,2) + std::pow(y,2))*z;
  H(2,1) = H(1,2);
  return H;
}