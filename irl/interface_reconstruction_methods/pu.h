#ifndef IRL_PARTITION_OF_UNITY_H_
#define IRL_PARTITION_OF_UNITY_H_

#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "irl/interface_reconstruction_methods/pu_neighborhood.h"

#include "irl/moments/cell_collection.h"
#include "irl/moments/cell_grouped_moments.h"

#include "irl/geometry/general/normal.h"
#include "irl/variant_reconstruction/separator_variant.h"

namespace IRL {
// This file contains all the functions for creating and using a partition of
// unity based on wendland functions.
template <class CellType>
class PU {
 private:
  PUNeighborhod<CellType> neighborhood_m;
  const double kernel_size_m;

 public:
  // Constructors
  PU(const double a_kernel_size) : kernel_size_m(a_kernel_size) {}

  // Get Value
  double getPU(const Pt& x);
  // Get Value and Gradient
  std::pair<double, Eigen::Vector3d> getPUAndGrad(const Pt& x);
  // Get Value, Gradient, and Hessian
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> getPUGradAndHess(
      const Pt& x);
  // Get Total Weight
  double getTotalWeight(const Pt& x);
  // Find intersection between implicit curve and a provided line.
  std::vector<Pt> intersectEdge(const Pt& x0, const Pt& x1,
                                const int& Npartitions, const double& tresh,
                                bool& blocked);
  // get Normal
  Normal getNormal(const Pt& x);
  // get Mean Curvature
  double getMeanCurvature(const Pt& x);
  // Project point onto implicit surface
  const Pt projectOntoPU(const Pt& a_pt, const double dx, bool& success);
  // Signed Distance of Separator
  static double implicitSeparator(const Pt& a_pt, const Pt& a_centroid,
                                  const SeparatorVariant* a_sepPtr);
  // Signed Distance and Gradient of Separator
  static std::pair<double, Eigen::Vector3d> implicitSeparator(
      const Pt& a_pt, const Pt& a_centroid, const SeparatorVariant* a_sepPtr);
  // Signed Distance, Gradient, and Hessian of Separator
  static std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> implicitSeparator(
      const Pt& a_pt, const Pt& a_centroid, const SeparatorVariant* a_sepPtr);
}

}  // namespace IRL

#endif