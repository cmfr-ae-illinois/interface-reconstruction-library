#include <iostream>

#include "examples/initialize_vf/surfaces.h"

// const double xc = 0.0, yc = 0.0, zc = 0.0;
const double R = 0.3, r = 0.15;

double F_torus(const double& x, const double& y, const double& z){
  return (x*x + y*y + z*z + R*R - r*r) * (x*x + y*y + z*z + R*R - r*r) - 
         4.0 * R*R * (x*x + y*y);
}

Eigen::Vector3d gradF_torus(const double& x, const double& y, const double& z){
  double Fx = 4.0 * z * (x*x + y*y + z*z + R*R - r*r) - 8.0 * x * R*R;
  double Fy = 4.0 * y * (x*x + y*y + z*z + R*R - r*r) - 8.0 * y * R*R;
  double Fz = 4.0 * z * (x*x + y*y + z*z + R*R - r*r);
  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Matrix3d hessF_torus(const double& x, const double& y, const double& z){
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H(0,0) = 4.0 * (x*x + y*y + z*z + R*R - r*r) * 8.0 * x*x - 8.0 * R*R;
  H(1,1) = 4.0 * (x*x + y*y + z*z + R*R - r*r) * 8.0 * y*y - 8.0 * R*R;
  H(2,2) = 4.0 * (x*x + y*y + z*z + R*R - r*r) * 8.0 * z*z;
  H(0,1) = 8.0 * x * y;
  H(1,0) = H(0,1);
  H(0,2) = 8.0 * x * z;
  H(2,0) = H(0,2);
  H(1,2) = 8.0 * y * z;
  H(2,1) = H(1,2);
  return H;
}