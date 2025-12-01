#include <iostream>

#include "examples/initialize_vf/surfaces.h"

const double xc = 0, yc = 0;
const double R = 0.15;

double F_circle(const double& x, const double& y){
  return (x - xc)*(x - xc) + (y - yc)*(y - yc) - R*R;
}

Eigen::Vector2d gradF_circle(const double& x, const double& y){
  double Fx = 2*(x - xc);
  double Fy = 2*(y - yc);
  return Eigen::Vector2d(Fx, Fy);
}

Eigen::Matrix2d hessF_circle(const double& x, const double& y){
  Eigen::Matrix2d H;
  H(0,0) = 2.0;
  H(0,1) = 0.0;
  H(1,0) = H(0,1);
  H(1,1) = 2.0;
  return H;
}
