#include <iostream>

#include "examples/initialize_vf/surfaces.h"

const double xc = 0.0, yc = 0.0, zc = 0.0;
const double R = 0.15;

double F_sphere(const double& x, const double& y, const double& z){
  return (x - xc) * (x - xc) + (y - yc) * (y - yc) + 
         (z - zc) * (z - zc) - R * R;
}

Eigen::Vector3d gradF_sphere(const double& x, const double& y, const double& z){
  double Fx = 2.0 * (x - xc);
  double Fy = 2.0 * (y - yc);
  double Fz = 2.0 * (z - zc);
  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Matrix3d hessF_sphere(const double& x, const double& y, const double& z){
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H(0,0) = 2.0;
  H(1,1) = 2.0;
  H(2,2) = 2.0;
  return H;
}