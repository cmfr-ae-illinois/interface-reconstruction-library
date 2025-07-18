#include <iostream>

#include "examples/initialize_vf/surfaces.h"


double Sphere::F(const double& x, const double& y, const double& z) const {
  return (x - xc) * (x - xc) + (y - yc) * (y - yc) + 
         (z - zc) * (z - zc) - R * R;
}

Eigen::Vector3d Sphere::gradF(const double& x, const double& y, const double& z) const {
  double Fx = 2.0 * (x - xc);
  double Fy = 2.0 * (y - yc);
  double Fz = 2.0 * (z - zc);
  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Matrix3d Sphere::hessF(const double& x, const double& y, const double& z) const {
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H(0,0) = 2.0;
  H(1,1) = 2.0;
  H(2,2) = 2.0;
  return H;
}

std::pair<IRL::Pt, IRL::Pt> Sphere::domainBounds() const {
  return {{-0.5,-0.5,-0.5},{0.5,0.5,0.5}};
}

double Sphere::volume() const {
  return 4.0 / 3.0 * M_PI * R * R * R;
}

double Sphere::surfaceArea() const {
  return 4.0 * M_PI * R * R;
}