#include <iostream>

#ifndef EXAMPLES_INITIALIZE_VF_SURFACES_H_
#define EXAMPLES_INITIALIZE_VF_SURFACES_H_

#include <Eigen/Dense>
#include <utility>
#include <stdexcept>
#include <cmath>

#include "irl/generic_cutting/generic_cutting.h"


// Base class
class Surface {
public:
  virtual ~Surface(void) = default;

  virtual double F(const double& x, const double& y, const double& z) const = 0;
  virtual Eigen::Vector3d gradF(const double& x, const double& y, const double& z) const = 0;
  virtual Eigen::Matrix3d hessF(const double& x, const double& y, const double& z) const = 0;
  virtual std::pair<IRL::Pt, IRL::Pt> domainBounds() const = 0;

  virtual bool hasVolume() const { return false; }
  virtual double volume() const {
    throw std::runtime_error("volume() not implemented for this shape.");
  }

  virtual bool hasSurfaceArea() const { return false; }
  virtual double surfaceArea() const {
    throw std::runtime_error("surfaceArea() not implemented for this shape.");
  }
};

// Implicit surfaces -------------------------------------------------------------------

// sphere
class Sphere : public Surface {
public:
  Sphere(void) = default;

  double F(const double& x, const double& y, const double& z) const override;
  Eigen::Vector3d gradF(const double& x, const double& y, const double& z) const override;
  Eigen::Matrix3d hessF(const double& x, const double& y, const double& z) const override;
  std::pair<IRL::Pt, IRL::Pt> domainBounds() const override;

  bool hasVolume() const override { return true; }
  double volume() const override;
  bool hasSurfaceArea() const override { return true; }
  double surfaceArea() const override;

private:
  const double xc = 0.0, yc = 0.0, zc = 0.0;
  const double R = 0.15;
};

// ellipsoid
class Ellipsoid : public Surface {
public:
  Ellipsoid(void) = default;

  double F(const double& x, const double& y, const double& z) const override;
  Eigen::Vector3d gradF(const double& x, const double& y, const double& z) const override;
  Eigen::Matrix3d hessF(const double& x, const double& y, const double& z) const override;
  std::pair<IRL::Pt, IRL::Pt> domainBounds() const override;

  bool hasVolume() const override { return true; }
  double volume() const override;

private:
  const double xc = 0.0, yc = 0.0, zc = 0.0;
  const double a = 0.3, b = 0.15, c = 0.1;
};

// genus
class Genus : public Surface {
public:
  Genus(void) = default;

  double F(const double& x, const double& y, const double& z) const override;
  Eigen::Vector3d gradF(const double& x, const double& y, const double& z) const override;
  Eigen::Matrix3d hessF(const double& x, const double& y, const double& z) const override;
  std::pair<IRL::Pt, IRL::Pt> domainBounds() const override;

};

// orthocircle
class Orthocircle : public Surface {
public:
  Orthocircle(void) = default;

  double F(const double& x, const double& y, const double& z) const override;
  Eigen::Vector3d gradF(const double& x, const double& y, const double& z) const override;
  Eigen::Matrix3d hessF(const double& x, const double& y, const double& z) const override;
  std::pair<IRL::Pt, IRL::Pt> domainBounds() const override;

private:
  const double c1 = 0.075, c2 = 3;
};

#endif