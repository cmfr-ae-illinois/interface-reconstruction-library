#include <iostream>

#include "examples/initialize_vf/surfaces.h"


double Ellipsoid::F(const double& x, const double& y, const double& z) const {
  return std::pow((x - xc) / a, 2.0) + std::pow((y - yc) / b, 2.0) +
         std::pow((z - zc) / c, 2.0) - 1.0;
}

Eigen::Vector3d Ellipsoid::gradF(const double& x, const double& y, const double& z) const {
  double Fx = 2.0 * (x - xc) / (a * a);
  double Fy = 2.0 * (y - yc) / (b * b);
  double Fz = 2.0 * (z - zc) / (c * c);
  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Matrix3d Ellipsoid::hessF(const double& x, const double& y, const double& z) const {
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H(0,0) = 2.0 / (a * a);
  H(1,1) = 2.0 / (b * b);
  H(2,2) = 2.0 / (c * c);
  return H;
}

std::pair<IRL::Pt, IRL::Pt> Ellipsoid::domainBounds() const {
  return {{-0.5,-0.5,-0.5},{0.5,0.5,0.5}};
}

double Ellipsoid::volume() const {
  return 4.0 / 3.0 * M_PI * a * b * c;
}