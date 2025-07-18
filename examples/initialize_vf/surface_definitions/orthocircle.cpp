#include <iostream>

#include "examples/initialize_vf/surfaces.h"


double Orthocircle::F(const double& x, const double& y, const double& z) const {
  return ((x * x + y * y - 1) * (x * x + y * y - 1) + z * z)
         * ((y * y + z * z - 1) * (y * y + z * z - 1) + x * x)
         * ((z * z + x * x - 1) * (z * z + x * x - 1) + y * y)
         - c1 * c1 * (1 + c2 * (x * x + y * y + z * z));
}

Eigen::Vector3d Orthocircle::gradF(const double& x, const double& y, const double& z) const {
  double Fx = -2 * std::pow(c1, 2) * c2 * x
            + 2 * x * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
                  * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            + 4 * x * (-1 + std::pow(x, 2) + std::pow(z, 2))
                  * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
                  * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2))
            + 4 * x * (-1 + std::pow(x, 2) + std::pow(y, 2))
                  * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
                  * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2));

  double Fy = -2 * std::pow(c1, 2) * c2 * y
            + 4 * y * (-1 + std::pow(y, 2) + std::pow(z, 2))
                  * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
                  * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            + 2 * y * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
                  * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2))
            + 4 * y * (-1 + std::pow(x, 2) + std::pow(y, 2))
                  * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
                  * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2));

  double Fz = -2 * std::pow(c1, 2) * c2 * z
            + 4 * z * (-1 + std::pow(y, 2) + std::pow(z, 2))
                  * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
                  * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            + 4 * z * (-1 + std::pow(x, 2) + std::pow(z, 2))
                  * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
                  * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2))
            + 2 * z * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
                  * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2));

  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Matrix3d Orthocircle::hessF(const double& x, const double& y, const double& z) const {
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();

  H(0,0) = -2 * std::pow(c1, 2) * c2
          + 2 * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
              * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
          + 4 * x * (
              4 * x * (-1 + std::pow(x, 2) + std::pow(z, 2))
                * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
              + 4 * x * (-1 + std::pow(x, 2) + std::pow(y, 2))
                * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            )
          + (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2)) * (
              32 * std::pow(x, 2) * (-1 + std::pow(x, 2) + std::pow(y, 2)) * (-1 + std::pow(x, 2) + std::pow(z, 2))
              + (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2)) * (8 * std::pow(x, 2) + 4 * (-1 + std::pow(x, 2) + std::pow(z, 2)))
              + (8 * std::pow(x, 2) + 4 * (-1 + std::pow(x, 2) + std::pow(y, 2))) * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            );

  H(1,1) = -2 * std::pow(c1, 2) * c2
          + (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
              * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
              * (8 * std::pow(y, 2) + 4 * (-1 + std::pow(y, 2) + std::pow(z, 2)))
          + 8 * y * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (
              2 * y * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
              + 4 * y * (-1 + std::pow(x, 2) + std::pow(y, 2)) * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            )
          + (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2)) * (
              16 * std::pow(y, 2) * (-1 + std::pow(x, 2) + std::pow(y, 2))
              + 2 * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
              + (8 * std::pow(y, 2) + 4 * (-1 + std::pow(x, 2) + std::pow(y, 2)))
                  * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            );

  H(2,2) = -2 * std::pow(c1, 2) * c2
          + (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
              * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
              * (8 * std::pow(z, 2) + 4 * (-1 + std::pow(y, 2) + std::pow(z, 2)))
          + (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2)) * (
              16 * std::pow(z, 2) * (-1 + std::pow(x, 2) + std::pow(z, 2))
              + (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
                  * (8 * std::pow(z, 2) + 4 * (-1 + std::pow(x, 2) + std::pow(z, 2)))
              + 2 * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            )
          + 8 * z * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (
              4 * z * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
              + 2 * z * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
            );

  H(0,1) = 4 * x * y * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
          + 16 * x * y * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
          + 8 * x * y * (-1 + std::pow(x, 2) + std::pow(y, 2)) * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
          + 16 * x * y * (-1 + std::pow(x, 2) + std::pow(y, 2)) * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
          + 8 * x * y * (-1 + std::pow(x, 2) + std::pow(y, 2)) * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2))
          + 16 * x * y * (-1 + std::pow(x, 2) + std::pow(y, 2)) * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2))
          + 8 * x * y * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2)) * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2));
  H(1,0) = H(0,1);

  H(0,2) = 8 * x * z * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
          + 16 * x * z * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
          + 4 * x * z * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
          + 16 * x * (-1 + std::pow(x, 2) + std::pow(y, 2)) * z * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
          + 8 * x * z * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2))
          + 16 * x * (-1 + std::pow(x, 2) + std::pow(y, 2)) * z * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2))
          + 8 * x * z * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2)) * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2));
  H(2,0) = H(0,2);

  H(1,2) = 8 * y * z * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
          + 16 * y * z * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2))
          + 8 * y * z * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
          + 16 * y * (-1 + std::pow(x, 2) + std::pow(y, 2)) * z * (-1 + std::pow(y, 2) + std::pow(z, 2)) * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
          + 8 * y * z * (std::pow(-1 + std::pow(x, 2) + std::pow(y, 2), 2) + std::pow(z, 2)) * (std::pow(y, 2) + std::pow(-1 + std::pow(x, 2) + std::pow(z, 2), 2))
          + 4 * y * z * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2))
          + 16 * y * (-1 + std::pow(x, 2) + std::pow(y, 2)) * z * (-1 + std::pow(x, 2) + std::pow(z, 2)) * (std::pow(x, 2) + std::pow(-1 + std::pow(y, 2) + std::pow(z, 2), 2));
  H(2,1) = H(1,2);

  return H;
}

std::pair<IRL::Pt, IRL::Pt> Orthocircle::domainBounds() const {
  return {{-2.5,-2.5,-2.5},{2.5,2.5,2.5}};
}