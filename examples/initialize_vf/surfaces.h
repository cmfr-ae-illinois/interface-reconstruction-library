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

double F_torus(const double& x, const double& y, const double& z);
Eigen::Vector3d gradF_torus(const double& x, const double& y, const double& z);
Eigen::Matrix3d hessF_torus(const double& x, const double& y, const double& z);

double F_genus(const double& x, const double& y, const double& z);
Eigen::Vector3d gradF_genus(const double& x, const double& y, const double& z);
Eigen::Matrix3d hessF_genus(const double& x, const double& y, const double& z);

double F_orthocircle(const double& x, const double& y, const double& z);
Eigen::Vector3d gradF_orthocircle(const double& x, const double& y, const double& z);
Eigen::Matrix3d hessF_orthocircle(const double& x, const double& y, const double& z);

double F_wineglass(const double& x, const double& y, const double& z);
Eigen::Vector3d gradF_wineglass(const double& x, const double& y, const double& z);
Eigen::Matrix3d hessF_wineglass(const double& x, const double& y, const double& z);

#endif