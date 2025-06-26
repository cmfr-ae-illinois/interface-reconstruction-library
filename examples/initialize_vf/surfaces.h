#include <iostream>

#ifndef EXAMPLES_INITIALIZE_VF_SURFACES_H_
#define EXAMPLES_INITIALIZE_VF_SURFACES_H_

#include <Eigen/Dense>

double F_sphere(const double& x, const double& y, const double& z);
Eigen::Vector3d gradF_sphere(const double& x, const double& y, const double& z);
Eigen::Matrix3d hessF_sphere(const double& x, const double& y, const double& z);

double F_ellipsoid(const double& x, const double& y, const double& z);
Eigen::Vector3d gradF_ellipsoid(const double& x, const double& y, const double& z);
Eigen::Matrix3d hessF_ellipsoid(const double& x, const double& y, const double& z);

#endif