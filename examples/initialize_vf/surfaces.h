#include <iostream>

#ifndef EXAMPLES_INITIALIZE_VF_SURFACES_H_
#define EXAMPLES_INITIALIZE_VF_SURFACES_H_

#include <Eigen/Dense>

double F_circle(const double& x, const double& y);
Eigen::Vector2d gradF_circle(const double& x, const double& y);
Eigen::Matrix2d hessF_circle(const double& x, const double& y);

double F_ellipse(const double& x, const double& y);
Eigen::Vector2d gradF_ellipse(const double& x, const double& y);
Eigen::Matrix2d hessF_ellipse(const double& x, const double& y);


#endif