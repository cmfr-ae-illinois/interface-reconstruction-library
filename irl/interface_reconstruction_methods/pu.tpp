#ifndef IRL_PARTITION_OF_UNITY_TPP_
#define IRL_PARTITION_OF_UNITY_TPP_

#include <limits>
#include <tuple>
#include <vector>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"
#include "irl/geometry/spline/rational_cubic_bezier_arc.h"
#include "irl/helpers/wendland.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/moments/general_moments.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

namespace IRL {

template <class CellType>
double PU<CellType>::getPU(const Pt& x) {
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();
  Eigen::Vector3d grad_product_sum = Eigen::Vector3d::Zero();
  for (int i = 0; i < neighborhood_m.size(); ++i) {  // Loop over separators
    // Distance Weights
    double weight = 0.0;
    Wendland::evaluate(neighborhood_m.getCentroid(i), kernel_size_m, x,
                       &weight);
    weight *= neighborhood_m.getWeight(i);
    // Add results to sums
    weight_sum += weight;

    // How get val,grad,hess of separator
    double F = PU<CellType>::implicitSeparatorValue(
        x, neighborhood_m.getCentroid(i), &neighborhood_m.getSeparator(i));
    // Now calculate F_sum, grad_product_sum, and hess_product_sum
    F_sum += weight * F;
  }
  // Put Everything Together and return.
  const double inv_weight_sum = 1.0 / safelyEpsilon(weight_sum);
  return F_sum * inv_weight_sum;
}

template <class CellType>
std::pair<double, Eigen::Vector3d> PU<CellType>::getPUAndGrad(const Pt& x) {
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();
  Eigen::Vector3d grad_product_sum = Eigen::Vector3d::Zero();
  for (int i = 0; i < neighborhood_m.size(); ++i) {  // Loop over separators
    // Distance Weights
    std::pair<double, Eigen::Vector3d> weightRet;
    Wendland::evaluate(neighborhood_m.getCentroid(i), kernel_size_m, x,
                       &weightRet);
    // Get Results
    const double weight =
        neighborhood_m.getWeight(i) * std::get<0>(weightRet);  // weights[i] *
    const Eigen::Vector3d grad_weight =
        neighborhood_m.getWeight(i) * std::get<1>(weightRet);  // weights[i] *

    // Add results to sums
    weight_sum += weight;
    grad_weight_sum += grad_weight;

    // How get val,grad,hess of separator
    std::pair<double, Eigen::Vector3d> separatorRet;
    separatorRet = PU<CellType>::implicitSeparatorValueandGrad(
        x, neighborhood_m.getCentroid(i), &neighborhood_m.getSeparator(i));

    // Get Values
    double F = std::get<0>(separatorRet);
    Eigen::Vector3d gradF = std::get<1>(separatorRet);

    // Now calculate F_sum, grad_product_sum, and hess_product_sum
    F_sum += weight * F;
    grad_product_sum += grad_weight * F + weight * gradF;
  }
  // Put Everything Together into PU Results
  const double inv_weight_sum = 1.0 / safelyEpsilon(weight_sum);
  const double PU_F = F_sum * inv_weight_sum;
  const Eigen::Vector3d PU_gradF =
      (grad_product_sum - F_sum * grad_weight_sum * inv_weight_sum) *
      inv_weight_sum;
  // Return
  return std::make_pair(PU_F, PU_gradF);
}

template <class CellType>
std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>
PU<CellType>::getPUGradAndHess(const Pt& x) {
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hess_weight_sum = Eigen::Matrix3d::Zero();
  Eigen::Vector3d grad_product_sum = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hess_product_sum = Eigen::Matrix3d::Zero();
  for (int i = 0; i < neighborhood_m.size(); ++i) {  // Loop over separators
    // Distance Weights
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> weightRet;
    Wendland::evaluate(neighborhood_m.getCentroid(i), kernel_size_m, x,
                       &weightRet);
    // Get Results
    const double weight = neighborhood_m.getWeight(i) * std::get<0>(weightRet);
    const Eigen::Vector3d grad_weight =
        neighborhood_m.getWeight(i) * std::get<1>(weightRet);
    const Eigen::Matrix3d hess_weight =
        neighborhood_m.getWeight(i) * std::get<2>(weightRet);

    // Add results to sums
    weight_sum += weight;
    grad_weight_sum += grad_weight;
    hess_weight_sum += hess_weight;

    // How get val,grad,hess of separator
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> separatorRet;
    separatorRet = PU<CellType>::implicitSeparatorValueGradHess(
        x, neighborhood_m.getCentroid(i), &neighborhood_m.getSeparator(i));

    // Get Values
    double F = std::get<0>(separatorRet);
    Eigen::Vector3d gradF = std::get<1>(separatorRet);
    Eigen::Matrix3d hessF = std::get<2>(separatorRet);

    // Now calculate F_sum, grad_product_sum, and hess_product_sum
    F_sum += weight * F;
    grad_product_sum += grad_weight * F + weight * gradF;
    hess_product_sum += hess_weight * F + grad_weight * gradF.transpose() +
                        gradF * grad_weight.transpose() + weight * hessF;
  }
  // Put Everything Together into PU Results
  const double inv_weight_sum = 1.0 / safelyEpsilon(weight_sum);
  const double PU_F = F_sum * inv_weight_sum;
  const Eigen::Vector3d PU_gradF =
      (grad_product_sum - F_sum * grad_weight_sum * inv_weight_sum) *
      inv_weight_sum;
  const Eigen::Matrix3d PU_hessF =
      (hess_product_sum + (grad_product_sum * grad_weight_sum.transpose() -
                           grad_weight_sum * grad_product_sum.transpose() -
                           F_sum * hess_weight_sum) *
                              inv_weight_sum) *
          inv_weight_sum -
      2. * (grad_product_sum - F_sum * grad_weight_sum * inv_weight_sum) *
          grad_weight_sum.transpose() * inv_weight_sum * inv_weight_sum;
  // Return
  return std::make_tuple(PU_F, PU_gradF, PU_hessF);
}

template <class CellType>
double PU<CellType>::getTotalWeight(
    const Pt& x) {  // CHECK ON TOTAL WEIGHT VALUES
  double weight_sum = 0.0;
  for (int i = 0; i < neighborhood_m.size(); ++i) {
    double weight = 0.0;
    Wendland::evaluate(neighborhood_m.getCentroid(i), kernel_size_m, x,
                       &weight);
    weight *= neighborhood_m.getWeight(i);
    weight_sum += weight;
  }
  return weight_sum;
}

template <class CellType>
std::vector<Pt> PU<CellType>::intersectEdge(const Pt& x0, const Pt& x1,
                                            const int& Npartitions,
                                            const double& thresh,
                                            bool& blocked) {
  // Split the domain into segments
  blocked = false;
  std::vector<Pt> sampleLocations = {};
  // At these locations, calculate the function value
  std::vector<double> values = {};
  // Also get the sign of these values
  std::vector<double> signs = {};
  for (int i = 0; i < Npartitions + 1; i++) {
    Pt temp =
        (1 - static_cast<double>(i) / static_cast<double>(Npartitions)) * x0 +
        (static_cast<double>(i) / static_cast<double>(Npartitions)) * x1;
    sampleLocations.push_back(temp);

    double val = this->getPU(temp);
    values.push_back(val);

    double sgn = (0.0 < val) - (val < 0.0);
    signs.push_back(sgn);
  }

  // Loop over all the partitions. If the signs are different, do a bisection
  // method to find the root
  std::vector<Pt> intersections = {};
  // std::cout << intersections.size() << std::endl;
  Pt upperX;
  Pt lowerX;
  Pt midX;
  double upperVal;
  double lowerVal;
  double midVal;

  double tol = 1e-12;
  double max_iters = 2000;
  double weight;
  for (int i = 0; i < Npartitions; i++) {
    if (signs[i] == 0 || signs[i + 1] == 0) {  // At least one root
      if (signs[i] == 0) {
        // std::cout << "Zero Intersection Found" << std::endl;
        intersections.push_back(sampleLocations[i]);
      }
      if (signs[i + 1] == 0 && i + 1 == Npartitions) {
        // std::cout << "Zero Intersection Found" << std::endl;
        intersections.push_back(sampleLocations[i + 1]);
      }
    } else if (signs[i] != signs[i + 1]) {  // Different signs, root somewhere
      // Decide which Side is upper and lower
      if (signs[i] == 1) {  // Left  is upper value
        upperX = sampleLocations[i];
        upperVal = values[i];

        lowerX = sampleLocations[i + 1];
        lowerVal = values[i + 1];
      } else {
        upperX = sampleLocations[i + 1];
        upperVal = values[i + 1];

        lowerX = sampleLocations[i];
        lowerVal = values[i];
      }
      // Apply Bisection Method
      for (int j = 0; j < max_iters; j++) {  // Do until you reach max iters
        midX = 0.5 * (upperX + lowerX);
        midVal = this->getPU(midX);
        if (std::abs(midVal) < tol) {
          break;
        } else if (midVal > 0.0) {
          upperX = midX;
        } else {
          lowerX = midX;
        }
      }  // End Bisection
      // Get weight at intersection
      weight = this->getTotalWeight(midX);
      // Add intersection if total weight is greater than threshold
      if (weight >= thresh) {
        intersections.push_back(midX);
      } else {
        // If the weight is too low, we consider the intersection to be blocked.
        std::cout << "Blocked = " << weight << "," << thresh << "\n";
        blocked = true;
      }
    }
  }  // End loop over partitions
  return intersections;
}

template <class CellType>
Normal PU<CellType>::getNormal(const Pt& x) {
  std::pair<double, Eigen::Vector3d> holdsGrad = this->getPUAndGrad(x);
  auto gradF = std::get<1>(holdsGrad);
  double Fx = gradF(0);
  double Fy = gradF(1);
  double Fz = gradF(2);
  Normal ret = Normal(Fx, Fy, Fz);
  ret.normalize();

  return ret;
}

template <class CellType>
double PU<CellType>::getMeanCurvature(const Pt& x) {
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> holdsGradAndHessian =
      this->getPUGradAndHess(x);
  auto gradF = std::get<1>(holdsGradAndHessian);
  auto hessF = std::get<2>(holdsGradAndHessian);

  double magGradF = gradF.norm();

  // grad(F)^T Hess(F) grad(F)
  double gradHessGrad = gradF.transpose() * hessF * gradF;

  // |grad(F)|^2 * Trace(Hess(F))
  double traceHess = hessF.trace();
  double gradSquaredTrace = magGradF * magGradF * traceHess;

  double numer = gradHessGrad - gradSquaredTrace;
  double denom = 2.0 * magGradF * magGradF * magGradF;

  return -numer / safelyEpsilon(denom);
}

template <class CellType>
Pt PU<CellType>::projectOntoPU(const Pt& a_pt, const double dx, bool& success) {
  success = true;
  std::pair<double, Eigen::Vector3d> holdsGrad;
  Pt projected_pt = a_pt;
  const int itmax = 50;
  for (int i = 0; i < itmax; i++) {
    holdsGrad = this->getPUAndGrad(projected_pt);
    const auto F = std::get<0>(holdsGrad);
    // Check Finite Value
    if (!std::isfinite(F)) {
      success = false;
      return a_pt;
    }
    // Return is near zero
    if (std::abs(F) < kernel_size_m * 1e-6) {
      break;
    }
    const auto gradF = std::get<1>(holdsGrad);
    // Check Finite Gradient
    if (!gradF.allFinite() ||
        gradF.squaredNorm() <= std::numeric_limits<double>::epsilon()) {
      success = false;
      return a_pt;
    }
    const double grad_norm_inv = 1.0 / safelyEpsilon(gradF.squaredNorm());

    auto step = -F * gradF * grad_norm_inv;
    for (int d = 0; d < 3; d++) {
      projected_pt[d] += step(d);
    }
    // return the point itself if max iterations reached
    if (i == (itmax - 1)) {
      success = false;
      return a_pt;
    }
  }

  double weight = this->getTotalWeight(projected_pt);
  // Make sure to find zero on actual important levelset, not on zero weight
  // edge.
  if (weight < 1e-6) {
    success = false;
    return a_pt;
  }
  // Make sure Point is good quality
  if (IRL::magnitude(projected_pt - a_pt) > 0.5 * dx) {
    success = false;
    return a_pt;
  }
  return projected_pt;
}

// Implicit Separators
template <class CellType>
double PU<CellType>::implicitSeparatorValue(const Pt& a_pt,
                                            const Pt& a_centroid,
                                            const SeparatorVariant* a_sepPtr) {
  const Pt x = a_pt - a_centroid;
  double F;
  if (const auto sepPtr = std::get_if<PlanarSeparator>(a_sepPtr)) {
    // std::cout << "Variant Plane Detected\n";
    if (sepPtr->getNumberOfPlanes() > 0) {
      const Plane plane = (*sepPtr)[0];
      const Normal n = plane.normal();
      const double d = plane.distance();

      if (n[0] != 0 || n[1] != 0 || n[2] != 0) {  // If plane exists
        F = n * a_pt - d;                         // Distance
      } else {                                    // IF plane doesn't exist
        F = 0;                                    // Zero
      }
    }
  } else if (const auto sepPtr = std::get_if<Paraboloid>(a_sepPtr)) {
    // std::cout << "Variant Paraboloid Detected\n";
    const ReferenceFrame frame = sepPtr->getReferenceFrame();
    const double a = sepPtr->getAlignedParaboloid().a();
    const double b = sepPtr->getAlignedParaboloid().b();
    // Move to local frame
    const Pt tmp = a_pt - sepPtr->getDatum();
    Pt xloc;
    for (int d = 0; d < 3; ++d) {
      xloc[d] = frame[d] * tmp;
    }

    // Taubin Distance Norm,grad,hessian
    const double dist_norm =
        1. / std::sqrt(1. + 4. * a * a * xloc[0] * xloc[0] +
                       4. * b * b * xloc[1] * xloc[1]);

    // Compute Algebraic Distance, grad, hessian
    const double F_alg =
        xloc[2] + a * xloc[0] * xloc[0] + b * xloc[1] * xloc[1];

    // Compute Signed Distance,grad,hessian
    F = F_alg * dist_norm;
  } else {
    throw std::runtime_error("No signed distance for Variant Type");
  }
  return F;
}

// Signed Distance and Gradient of Separator
template <class CellType>
std::pair<double, Eigen::Vector3d> PU<CellType>::implicitSeparatorValueandGrad(
    const Pt& a_pt, const Pt& a_centroid, const SeparatorVariant* a_sepPtr) {
  const Pt x = a_pt - a_centroid;
  double F;
  Eigen::Vector3d gradF;
  if (const auto sepPtr = std::get_if<PlanarSeparator>(a_sepPtr)) {
    // std::cout << "Variant Plane Detected\n";
    if (sepPtr->getNumberOfPlanes() > 0) {
      const Plane plane = (*sepPtr)[0];
      const Normal n = plane.normal();
      const double d = plane.distance();
      if (n[0] != 0 || n[1] != 0 || n[2] != 0) {  // If plane exists
        F = n * a_pt - d;                         // Distance
      } else {                                    // IF plane doesn't exist
        F = 0;                                    // Zero
      }
      gradF = Eigen::Vector3d(n[0], n[1], n[2]);
    }
  } else if (const auto sepPtr = std::get_if<Paraboloid>(a_sepPtr)) {
    // std::cout << "Variant Paraboloid Detected\n";
    const ReferenceFrame frame = sepPtr->getReferenceFrame();
    const double a = sepPtr->getAlignedParaboloid().a();
    const double b = sepPtr->getAlignedParaboloid().b();
    // Move to local frame
    const Pt tmp = a_pt - sepPtr->getDatum();
    Pt xloc;
    for (int d = 0; d < 3; ++d) {
      xloc[d] = frame[d] * tmp;
    }
    const Eigen::Vector3d e0(frame[0][0], frame[0][1], frame[0][2]);
    const Eigen::Vector3d e1(frame[1][0], frame[1][1], frame[1][2]);
    const Eigen::Vector3d e2(frame[2][0], frame[2][1], frame[2][2]);

    // Taubin Distance Norm,grad,hessian
    const double dist_norm =
        1. / std::sqrt(1. + 4. * a * a * xloc[0] * xloc[0] +
                       4. * b * b * xloc[1] * xloc[1]);
    const Eigen::Vector3d grad_dist_norm =
        -4. * dist_norm * dist_norm * dist_norm *
        (a * a * e0 * xloc[0] + b * b * e1 * xloc[1]);

    // Compute Algebraic Distance, grad, hessian
    const double F_alg =
        xloc[2] + a * xloc[0] * xloc[0] + b * xloc[1] * xloc[1];
    const Eigen::Vector3d grad_F_alg =
        e2 + 2. * (a * e0 * xloc[0] + b * e1 * xloc[1]);

    // Compute Signed Distance,grad,hessian
    F = F_alg * dist_norm;
    gradF = F_alg * grad_dist_norm + grad_F_alg * dist_norm;
  } else {
    throw std::runtime_error("No signed distance for Variant Type");
  }
  return std::make_pair(F, gradF);
}
// Signed Distance, Gradient, and Hessian of Separator
template <class CellType>
std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>
PU<CellType>::implicitSeparatorValueGradHess(const Pt& a_pt,
                                             const Pt& a_centroid,
                                             const SeparatorVariant* a_sepPtr) {
  const Pt x = a_pt - a_centroid;
  double F;
  Eigen::Vector3d gradF;
  Eigen::Matrix3d hessF;
  if (const auto sepPtr = std::get_if<PlanarSeparator>(a_sepPtr)) {
    // std::cout << "Variant Plane Detected\n";
    if (sepPtr->getNumberOfPlanes() > 0) {
      const Plane plane = (*sepPtr)[0];
      const Normal n = plane.normal();
      const double d = plane.distance();
      if (n[0] != 0 || n[1] != 0 || n[2] != 0) {  // If plane exists
        F = n * a_pt - d;                         // Distance
      } else {                                    // IF plane doesn't exist
        F = 0;                                    // Zero
      }
      gradF = Eigen::Vector3d(n[0], n[1], n[2]);
      hessF = Eigen::Matrix3d::Zero();
    }
  } else if (const auto sepPtr = std::get_if<Paraboloid>(a_sepPtr)) {
    // std::cout << "Variant Paraboloid Detected\n";
    const ReferenceFrame frame = sepPtr->getReferenceFrame();
    const double a = sepPtr->getAlignedParaboloid().a();
    const double b = sepPtr->getAlignedParaboloid().b();
    // Move to local frame
    const Pt tmp = a_pt - sepPtr->getDatum();
    Pt xloc;
    for (int d = 0; d < 3; ++d) {
      xloc[d] = frame[d] * tmp;
    }
    const Eigen::Vector3d e0(frame[0][0], frame[0][1], frame[0][2]);
    const Eigen::Vector3d e1(frame[1][0], frame[1][1], frame[1][2]);
    const Eigen::Vector3d e2(frame[2][0], frame[2][1], frame[2][2]);

    // Taubin Distance Norm,grad,hessian
    const double dist_norm =
        1. / std::sqrt(1. + 4. * a * a * xloc[0] * xloc[0] +
                       4. * b * b * xloc[1] * xloc[1]);
    const Eigen::Vector3d grad_dist_norm =
        -4. * dist_norm * dist_norm * dist_norm *
        (a * a * e0 * xloc[0] + b * b * e1 * xloc[1]);
    const Eigen::Matrix3d hess_dist_norm =
        4. * dist_norm * dist_norm * dist_norm * dist_norm * dist_norm *
        (a * a *
             (8. * a * a * xloc[0] * xloc[0] - 4. * b * b * xloc[1] * xloc[1] -
              1.) *
             e0 * e0.transpose() +
         b * b *
             (8. * b * b * xloc[1] * xloc[1] - 4. * a * a * xloc[0] * xloc[0] -
              1.) *
             e1 * e1.transpose() +
         12. * a * a * b * b * xloc[0] * xloc[1] *
             (e1 * e0.transpose() + e0 * e1.transpose()));

    // Compute Algebraic Distance, grad, hessian
    const double F_alg =
        xloc[2] + a * xloc[0] * xloc[0] + b * xloc[1] * xloc[1];
    const Eigen::Vector3d grad_F_alg =
        e2 + 2. * (a * e0 * xloc[0] + b * e1 * xloc[1]);
    const Eigen::Matrix3d hess_F_alg =
        2. * (a * e0 * e0.transpose() + b * e1 * e1.transpose());

    // Compute Signed Distance,grad,hessian
    F = F_alg * dist_norm;
    gradF = F_alg * grad_dist_norm + grad_F_alg * dist_norm;
    hessF = F_alg * hess_dist_norm + grad_F_alg * grad_dist_norm.transpose() +
            grad_dist_norm * grad_F_alg.transpose() + dist_norm * hess_F_alg;
  } else {
    throw std::runtime_error("No signed distance for Variant Type");
  }
  return std::make_tuple(F, gradF, hessF);
}

// Set Neighborhood
template <class CellType>
void PU<CellType>::setNeighborhood(
    const PUNeighborhood<CellType>& a_neighborhood) {
  neighborhood_m = a_neighborhood;
}

template <class CellType>
void PU<CellType>::setKernelSize(const double a_kernel_size) {
  kernel_size_m = a_kernel_size;
}

// Print Surface
template <class CellType>
void PU<CellType>::printSurface() {
  std::cout << "> Kernel Size = " << kernel_size_m << "\n";
  std::cout << "> Neighborhood Size = " << neighborhood_m.size() << "\n";
  // Loop over separators and centroids and print
  for (int i = 0; i < neighborhood_m.size(); i++) {
    std::cout << neighborhood_m.getSeparator(i) << " at "
              << neighborhood_m.getCentroid(i) << "\n";
  }
}

}  // End Namespace IRL

#endif