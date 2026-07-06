#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_

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

// =================== Implicit Surface Class Functions
// template <class SeparatorType>
// PUImplicitSurface::PUImplicitSurface(
//     const std::vector<Pt>& centroids_,
//     const std::vector<SeparatorVariant>& separators_,
//     const double& kernel_size_)
//     : centroids(centroids_),
//       separators(separators_),
//       kernel_size(kernel_size_) {}

// Separator Stuff
inline void PUImplicitSurface::printSurface() {
  // Print Kernel Size
  std::cout << "> Kernel Size = " << this->kernel_size << "\n";
  // Loop over separators and centroids and print
  for (int i = 0; i < centroids.size(); i++) {
    std::cout << separators[i] << " at " << centroids[i] << "\n";
  }
}
// Signed Distance of Separator
inline void PUImplicitSurface::implicitSeparator(
    const Pt& a_pt, const Pt& a_centroid, const SeparatorVariant* a_sepPtr,
    double* retVal) {
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
  *retVal = F;
}
// Signed Distance and Gradient of Separator
inline void PUImplicitSurface::implicitSeparator(
    const Pt& a_pt, const Pt& a_centroid, const SeparatorVariant* a_sepPtr,
    std::pair<double, Eigen::Vector3d>* retVal) {
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
  *retVal = std::make_pair(F, gradF);
}
// Signed Distance, Gradient, and Hessian of Separator
inline void PUImplicitSurface::implicitSeparator(
    const Pt& a_pt, const Pt& a_centroid, const SeparatorVariant* a_sepPtr,
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal) {
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
  *retVal = std::make_tuple(F, gradF, hessF);
}

// Evaluate Function Hess =====================================
inline void PUImplicitSurface::evaluate(
    const Pt& x, std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal) {
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hess_weight_sum = Eigen::Matrix3d::Zero();
  Eigen::Vector3d grad_product_sum = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hess_product_sum = Eigen::Matrix3d::Zero();
  for (int i = 0; i < centroids.size(); ++i) {  // Loop over separators
    // Distance Weights
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> weightRet;
    Wendland::evaluate(centroids[i], kernel_size, x, &weightRet);
    // Get Results
    const double weight = weights[i] * std::get<0>(weightRet);
    const Eigen::Vector3d grad_weight = weights[i] * std::get<1>(weightRet);
    const Eigen::Matrix3d hess_weight = weights[i] * std::get<2>(weightRet);

    // Add results to sums
    weight_sum += weight;
    grad_weight_sum += grad_weight;
    hess_weight_sum += hess_weight;

    // How get val,grad,hess of separator
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> separatorRet;
    PUImplicitSurface::implicitSeparator(x, centroids[i], &separators[i],
                                         &separatorRet);

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
  *retVal = std::make_tuple(PU_F, PU_gradF, PU_hessF);
}

// Evaluate Function Grad =================================
inline void PUImplicitSurface::evaluate(
    const Pt& x, std::pair<double, Eigen::Vector3d>* retVal) {
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();
  Eigen::Vector3d grad_product_sum = Eigen::Vector3d::Zero();
  for (int i = 0; i < centroids.size(); ++i) {  // Loop over separators
    // Distance Weights
    std::pair<double, Eigen::Vector3d> weightRet;
    Wendland::evaluate(centroids[i], kernel_size, x, &weightRet);
    // Get Results
    const double weight = std::get<0>(weightRet);                // weights[i] *
    const Eigen::Vector3d grad_weight = std::get<1>(weightRet);  // weights[i] *

    // Add results to sums
    weight_sum += weight;
    grad_weight_sum += grad_weight;

    // How get val,grad,hess of separator
    std::pair<double, Eigen::Vector3d> separatorRet;
    PUImplicitSurface::implicitSeparator(x, centroids[i], &separators[i],
                                         &separatorRet);

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
  *retVal = std::make_pair(PU_F, PU_gradF);
}

// Evaluate Function Value ======================
inline void PUImplicitSurface::evaluate(const Pt& x, double* retVal) {
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();
  Eigen::Vector3d grad_product_sum = Eigen::Vector3d::Zero();
  for (int i = 0; i < centroids.size(); ++i) {  // Loop over separators
    // Distance Weights
    double weight = 0.0;
    Wendland::evaluate(centroids[i], kernel_size, x, &weight);
    weight *= weights[i];
    // Add results to sums
    weight_sum += weight;

    // How get val,grad,hess of separator
    double F = 0.0;
    PUImplicitSurface::implicitSeparator(x, centroids[i], &separators[i], &F);
    // Now calculate F_sum, grad_product_sum, and hess_product_sum
    F_sum += weight * F;
  }
  // Put Everything Together and return.
  const double inv_weight_sum = 1.0 / safelyEpsilon(weight_sum);
  *retVal = F_sum * inv_weight_sum;
}

// Evaluate Function Hess for Circle =====================================
inline void PUImplicitSurface::evaluateCylinder(
    Pt& x, double radius, Pt& center,
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal) {
  Pt offset = x - center;
  double F = offset[0] * offset[0] + offset[1] * offset[1] - radius * radius;

  Eigen::Vector3d gradF(2 * offset[0], 2 * offset[1], 0);

  Eigen::Matrix3d hessF(3, 3);
  hessF(0, 0) = 2;
  hessF(1, 1) = 2;

  // Return
  *retVal = std::make_tuple(F, gradF, hessF);
}

// Evaluate Function Grad for Circle =====================================
inline void PUImplicitSurface::evaluateCylinder(
    Pt& x, double radius, Pt& center,
    std::pair<double, Eigen::Vector3d>* retVal) {
  Pt offset = x - center;
  double F = offset[0] * offset[0] + offset[1] * offset[1] - radius * radius;

  Eigen::Vector3d gradF(2 * offset[0], 2 * offset[1], 0);

  // Return
  *retVal = std::make_pair(F, gradF);
}

// Evaluate Function Value for Circle =====================================
inline void PUImplicitSurface::evaluateCylinder(Pt& x, double radius,
                                                Pt& center, double* retVal) {
  Pt offset = x - center;
  double F =
      (offset[0] * offset[0]) + (offset[1] * offset[1]) - (radius * radius);

  // Return
  *retVal = F;
}

inline std::vector<Pt> PUImplicitSurface::intersectEdgeCylinder(
    const Pt& x0, const Pt& x1, double radius, Pt& center,
    const int& Npartitions) {
  // Split the domain into segments
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
    double val = 0.0;
    this->evaluateCylinder(temp, radius, center, &val);
    values.push_back(val);

    double sgn = (0.0 < val) - (val < 0.0);
    signs.push_back(sgn);
    // std::cout << "===================== IN INTERSECT EDGE" <<
    // this->kernel_size
    //           << std::endl;
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
  double max_iters = 200;
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
        this->evaluateCylinder(midX, radius, center, &midVal);
        if (std::abs(midVal) < tol) {
          break;
        } else if (midVal > 0.0) {
          upperX = midX;
        } else {
          lowerX = midX;
        }
      }  // End Bisection
      intersections.push_back(midX);
    }
  }  // End loop over partitions
  return intersections;
}

// Get Total Weight at a Point
inline void PUImplicitSurface::getTotalWeight(Pt& x, double* retVal) {
  double weight_sum = 0.0;
  for (int i = 0; i < centroids.size(); ++i) {  // Loop over separators
    // Distance Weights
    double weight = 0.0;
    Wendland::evaluate(centroids[i], kernel_size, x, &weight);
    // Add results to sums
    weight_sum += weight;
  }
  // Put Everything Together and return.
  *retVal = weight_sum;
}

inline Normal PUImplicitSurface::getTangentCylinder(Pt& x, double radius,
                                                    Pt& center) {
  std::pair<double, Eigen::Vector3d> holdsGrad;
  this->evaluateCylinder(x, radius, center, &holdsGrad);
  auto gradF = std::get<1>(holdsGrad);
  double Fx = gradF(0);
  double Fy = gradF(1);

  return Normal(-Fy, Fx, 0.0);
}

// template <class SeparatorType>
inline std::vector<Pt> PUImplicitSurface::intersectEdge(const Pt& x0,
                                                        const Pt& x1,
                                                        const int& Npartitions,
                                                        const double& thresh) {
  // Split the domain into segments
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
    double val = 0.0;
    this->evaluate(temp, &val);
    values.push_back(val);

    double sgn = (0.0 < val) - (val < 0.0);
    signs.push_back(sgn);
    // std::cout << "===================== IN INTERSECT EDGE" <<
    // this->kernel_size
    //           << std::endl;
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
  double max_iters = 200;
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
        this->evaluate(midX, &midVal);
        if (std::abs(midVal) < tol) {
          break;
        } else if (midVal > 0.0) {
          upperX = midX;
        } else {
          lowerX = midX;
        }
      }  // End Bisection
      // Get weight at intersection
      this->getTotalWeight(midX, &weight);
      // Add intersection if total weight is greater than threshold
      if (weight >= thresh) {
        intersections.push_back(midX);
      } else {
        std::cout << "Blocked = " << weight << "," << thresh << "\n";
      }
    }
  }  // End loop over partitions
  return intersections;
}

inline Normal PUImplicitSurface::getTangent(Pt& x) {
  std::pair<double, Eigen::Vector3d> holdsGrad;
  this->evaluate(x, &holdsGrad);
  auto gradF = std::get<1>(holdsGrad);
  double Fx = gradF(0);
  double Fy = gradF(1);

  return Normal(-Fy, Fx, 0.0);
}

inline Normal PUImplicitSurface::getNormal(Pt& x) {
  std::pair<double, Eigen::Vector3d> holdsGrad;
  this->evaluate(x, &holdsGrad);
  auto gradF = std::get<1>(holdsGrad);
  double Fx = gradF(0);
  double Fy = gradF(1);
  double Fz = gradF(2);

  return Normal(Fx, Fy, Fz);
}

inline double PUImplicitSurface::getCurvature(Pt& x) {
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> holdsGradAndHessian;
  this->evaluate(x, &holdsGradAndHessian);
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

  return numer / denom;
}

inline const Pt PUImplicitSurface::projectOntoPU(const Pt& a_pt) {
  std::pair<double, Eigen::Vector3d> holdsGrad;
  Pt projected_pt = a_pt;
  const int itmax = 5;
  for (int i = 0; i < itmax; i++) {
    this->evaluate(projected_pt, &holdsGrad);
    const auto F = std::get<0>(holdsGrad);
    const auto gradF = std::get<1>(holdsGrad);
    const double grad_norm_inv = 1.0 / safelyEpsilon(gradF.squaredNorm());

    auto step = -F * gradF * grad_norm_inv;
    double alpha = 1.0;
    double weight_new;
    double weight_old;
    for (int d = 0; d < 3; d++) {
      projected_pt[d] += alpha * step(d);
    }
    // while (alpha > 1e-6) {  // Line search
    //   Pt new_pt = projected_pt;

    //   this->getTotalWeight(new_pt, &weight_new);
    //   this->getTotalWeight(projected_pt, &weight_old);
    //   if (weight_new <= 1e-6 || weight_old <= 1e-6) {
    //     std::cout << "Alpha = " << alpha << ", Weight New = " << weight_new
    //               << ", Weight Old = " << weight_old << "\n";
    //   }

    //   if (weight_new >= 0.1) {
    //     projected_pt = new_pt;
    //     break;
    //   }

    //   alpha *= 0.5;
    // }
  }
  return projected_pt;
}

// ============== Solver Methods
template <class CellType>
PU<CellType>::PU(void) {}

template <class CellType>
PU<CellType>::PU(PUNeighborhood<CellType> stencil_) {
  this->stencil_m = stencil_;
  // this->surface_m = this->neighborhoodToImplicitSurface(5.0);
}

template <class CellType>
PUImplicitSurface PU<CellType>::neighborhoodToImplicitSurface(double delta) {
  // const auto centroids = stencil_m.getCentroids();
  // const auto separators = stencil_m.getSeparators();
  // const auto weights = stencil_m.getWeights();
  return PUImplicitSurface(stencil_m.getCentroids(), stencil_m.getSeparators(),
                           stencil_m.getWeights(), delta);
}

template <class CellType>
Normal PU<CellType>::solveEdge(const double STin, const Pt& P0, const Pt& P1,
                               const double delta, const double Pressure,
                               const Normal& Marangoni) {
  // The Marangoni normal object holds the Xgradient, then the Y gradient, then
  // the temperature gradient of ST Marangoni = [Gx,Gy,sigma_T] Make Implicit
  // Surface std::cout << "In Solve Edge\n";

  // Something is up with paraboloids because they take like 8x the time to run.
  double STCoeff = STin;
  // The Pressure Option tells us if we should include the pressure terms or
  // not.
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(delta);
  // Find Intersection with edge
  std::vector<Pt> intersections =
      s.intersectEdge(P0, P1, 10, intersection_threshold_m);

  // Calculate some edge properties
  Pt dP = P1 - P0;
  Normal OutwardsNormal = {dP[1], -dP[0], 0};
  OutwardsNormal.normalize();
  double D = std::sqrt(dP[0] * dP[0] + dP[1] * dP[1]);
  double denom = 1.0 / (safelyEpsilon(D));
  // Give space for working variables
  Normal tangent;
  Normal gradient;
  Normal total = {0.0, 0.0, 0.0};  // Total force

  // Loop over intersections and add tangents to force
  if (intersections.size() > 0) {
    for (int j = 0; j < intersections.size(); j++) {
      tangent = s.getTangent(intersections[j]);
      tangent.normalize();
      // Here we apply the Marangoni surface tension. To do this, we assume
      // STCoeff = STCoeff + gamma_T(T-T_0). Letting gamma_T=-0.002 (Ratio for
      // water). We also pick T-T0=Gx to be the form we use. This will give us
      // th surface tension at this point We pick G = 10 for the stake of
      // simplicity
      STCoeff = STCoeff + Marangoni[2] * (Marangoni[0] * intersections[j][0] +
                                          Marangoni[1] * intersections[j][1]);
      if (Marangoni[2] == 0.0 &&
          (Marangoni[0] != 0.0 ||
           Marangoni[1] != 0.0)) {  // Force droplet breakup
        double gamma0 = 1.0;        // Base surface tension coefficient
        double R = 1.0;             // Characteristic length scale

        // Linear decrease in surface tension with x-coordinate - Al-Saud
        // Style
        STCoeff =
            gamma0 * std::max(1 - 1.25 * std::abs(intersections[j][0]), 0.1);
      }
      // If facing inside, multiply by negative
      // If facing outside, leave the same
      double Scale = OutwardsNormal * tangent;
      Scale = Scale / safelyEpsilon(std::abs(Scale));
      // Add
      total = total + Scale * STCoeff * denom * tangent;
      // Next, we need to calculate pressure contribution, if wanted
      // To determine if an intersection is going into or out of the fluid, we
      // can use the dot product of dP with the gradient (which is outwards
      // pointing). This product will contain orientation information also.
      if (Pressure >= 0.5) {
        gradient =
            Normal(tangent[1], -tangent[0],
                   0.0);  // Use this to get gradient so we don't have torecalc
        Normal dPN = {dP[0], dP[1], 0};  // Make dP into a vector
        dPN.normalize();
        Scale = gradient * dPN;                          // Dot Product
        Scale = Scale / safelyEpsilon(std::abs(Scale));  // Sign Only
        // This scale tells me if, going along the edge from P0 to P1, we are
        // going into the fluid (Scale >0) or out of the fluid (Scale <0).

        // Now we  have the sign. The domain that each intersection influences
        // is from the centerpoint out. This means that if the distance from P0
        // to the intersection if < 0.5D, then we look at the distance
        // between the intersection and P0. If it is >0.5D, then look look
        // at the distance between it and P1, which will be the same as
        // D-[the initial distance].

        Pt dPI = intersections[j] - P0;
        double L = std::sqrt(dPI[0] * dPI[0] +
                             dPI[1] * dPI[1]);  // Calculated Distance

        double Sx = L * denom;  // Also definitely positive
        if (L > 0.5 * D) {      // If More than 0.5D
          L = L - D;            // Switch to other direction (Center Out)
        }

        // Now that we have the length of section, we just need to calculat
        // the curvature, then calculat the pressure jump, then multiply it
        // all together and add

        double curv = s.getCurvature(intersections[j]);
        // This also contains some direction information, but I am unsure how
        // this needs to be positive or not. I am going to leave it for now,
        // and we will see if that is a problem. If there seems to be issues,
        // try to make this stritly positive.

        // Now we just add this all together. Note, the pressure force
        // alwaysacts inwards, so we multiply by the negative of the outwards
        // pointingnormal
        if (std::abs(Scale) >=
            1e-10) {  // Ensure that in the undefined areas, we
          // // do not add pressure terms.
          total = total + Scale * STCoeff * curv * L * denom * OutwardsNormal;
          // std::cout << "Surface Tension Coeff: " << STCoeff << "\n"
          //           << "Curvature: " << curv << "\n"
          //           << "Length Segment: " << L << "\n"
          //           << "Pressure Contribution Scale: " << Scale << "\n"
          //           << "Intersection Point: " << intersections[j] << "\n"
          //           << "Point 1: " << P0 << "\n"
          //           << "Point 2: " << P1 << "\n"
          //           << "--------------------------------------\n";
        }
      }
    }
  }
  // if (std::abs(total[0]) >= 1e-6 || std::abs(total[1]) >= 1e-6) {
  //   std::cout << "total: " << total << "\n";
  // }
  return total;
}

template <class CellType>
Normal PU<CellType>::solveFace(const double STin, const Pt& P0, const Pt& P1,
                               const Pt& P2, const Pt& P3, const double delta,
                               const double Pressure, const Normal& Marangoni) {
  const double EPSILON = machine_epsilon<double>();
  // The Marangoni normal object holds the Xgradient, then the Y gradient, then
  // the temperature gradient of ST Marangoni = [Gx,Gy,sigma_T] Make Implicit
  // Surface std::cout << "In Solve Edge\n";

  // Something is up with paraboloids because they take like 8x the time to run.
  double STCoeff = STin;
  // The Pressure Option tells us if we should include the pressure terms or
  // not.
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(delta);
  double val0, val1, val2, val3, val4;
  // Copy Points for now to get around constant constraints
  Pt P0_copy = P0;
  Pt P1_copy = P1;
  Pt P2_copy = P2;
  Pt P3_copy = P3;
  // Calculate Corner Values
  s.evaluate(P0, &val0);
  s.evaluate(P1, &val1);
  s.evaluate(P2, &val2);
  s.evaluate(P3, &val3);
  // Calculate Corner Weights
  double weight0, weight1, weight2, weight3, weight4;
  s.getTotalWeight(P0_copy, &weight0);
  s.getTotalWeight(P1_copy, &weight1);
  s.getTotalWeight(P2_copy, &weight2);
  s.getTotalWeight(P3_copy, &weight3);
  // If any of the weights are below A threshold, return 0 force since that
  // means we are not near the interface
  double factor = 100.0;
  if (weight0 < factor * EPSILON || weight1 < factor * EPSILON ||
      weight2 < factor * EPSILON || weight3 < factor * EPSILON) {
    // std::cout << "Low Weight Detected: " << weight0 << "," << weight1 << ","
    //           << weight2 << "," << weight3 << "\n";
    return Normal(0.0, 0.0, 0.0);
  }
  Pt dP = P1 - P0;
  double D = std::sqrt(dP[0] * dP[0] + dP[1] * dP[1] + dP[2] * dP[2]);
  double denom = 1.0 / (safelyEpsilon(D));
  s.evaluate(0.25 * (P0 + P1 + P2 + P3), &val4);  // Center Point Value
  // Calculate case:
  int caseValue =
      1 * (val0 > 0) + 2 * (val1 > 0) + 4 * (val2 > 0) + 8 * (val3 > 0);
  // If all corners are outside or inside, return 0 force. Otehrwise, pair
  // intersections
  std::vector<std::vector<int>> pairs = {{-1, -1}, {-1, -1}};
  switch (caseValue) {
    case 0:
      return Normal(0.0, 0.0, 0.0);
      break;
    case 1:
      pairs[0] = {0, 3};
      break;
    case 2:
      pairs[0] = {1, 0};
      break;
    case 3:
      pairs[0] = {1, 3};
      break;
    case 4:
      pairs[0] = {2, 1};
      break;
    case 5:
      if (val4 > 0) {
        pairs[0] = {0, 1};
        pairs[1] = {2, 3};
      } else {
        pairs[0] = {0, 3};
        pairs[1] = {2, 1};
      }
      break;
    case 6:
      pairs[0] = {2, 0};
      break;
    case 7:
      pairs[0] = {2, 3};
      break;
    case 8:
      pairs[0] = {3, 2};
      break;
    case 9:
      pairs[0] = {0, 2};
      break;
    case 10:
      if (val4 > 0) {
        pairs[0] = {1, 2};
        pairs[1] = {3, 0};
      } else {
        pairs[0] = {1, 0};
        pairs[1] = {3, 2};
      }
      break;
    case 11:
      pairs[0] = {1, 2};
      break;
    case 12:
      pairs[0] = {3, 1};
      break;
    case 13:
      pairs[0] = {0, 1};
      break;
    case 14:
      pairs[0] = {3, 0};
      break;
    case 15:
      return Normal(0.0, 0.0, 0.0);
      break;
    default:
      break;
  }
  std::vector<std::vector<Pt>> intersectionsSet =
      std::vector<std::vector<Pt>>(4, std::vector<Pt>{});
  // Find all Intersections with face edges
  std::vector<Pt> intersections01 =
      s.intersectEdge(P0, P1, 1, intersection_threshold_m);
  std::vector<Pt> intersections12 =
      s.intersectEdge(P1, P2, 1, intersection_threshold_m);
  std::vector<Pt> intersections23 =
      s.intersectEdge(P2, P3, 1, intersection_threshold_m);
  std::vector<Pt> intersections30 =
      s.intersectEdge(P3, P0, 1, intersection_threshold_m);

  intersectionsSet[0] = intersections01;
  intersectionsSet[1] = intersections12;
  intersectionsSet[2] = intersections23;
  intersectionsSet[3] = intersections30;

  // Calculate some edge and plane properties
  Pt dP1 = P1 - P0;
  Pt dP2 = P2 - P1;
  Normal edge1 = {dP1[0], dP1[1], dP1[2]};
  Normal edge2 = {dP2[0], dP2[1], dP2[2]};
  Normal faceNormal = IRL::crossProduct(edge1, edge2);
  faceNormal.normalize();

  // Give space for working variables
  Normal normal1, normal2, tangentStart, tangentEnd;
  Pt startPoint, endPoint, controlPoint;
  double weight = 1;  // Make Nonrational for now, but can add later if we want
  Normal total = {0.0, 0.0, 0.0};  // Total force

  // Loop over pairs, make splines, integrate, add forces to total.
  for (int j = 0; j < 2; j++) {
    if (pairs[j][0] != -1) {  // If there is a pair
      startPoint = intersectionsSet[pairs[j][0]][0];
      endPoint = intersectionsSet[pairs[j][1]][0];
      normal1 = s.getNormal(startPoint);
      normal1.normalize();
      normal2 = s.getNormal(endPoint);
      normal2.normalize();
      // Make tangent by crossing normal with face  normal
      tangentStart = IRL::crossProduct(faceNormal, normal1);
      tangentStart.normalize();
      tangentEnd = IRL::crossProduct(faceNormal, normal2);
      tangentEnd.normalize();

      RationalCubicBezierArc arc = RationalCubicBezierArc(
          startPoint, tangentStart, endPoint, tangentEnd);
      int QuadRuleOrder = 50;
      const auto& abscissea = AbscissaeGauss<double, 50>();
      const auto& weights = WeightsGauss<double, 50>();
      for (int j = 0; j < QuadRuleOrder; ++j) {
        const double t = 0.5 * (1.0 + abscissea[j]);
        const double w = weights[j];

        // Evaluation Points
        Pt pt = arc.point(t);
        // get normal,tangent, and surface tension at each point.
        Pt tangent = arc.derivative(t);
        Normal tan = Normal(tangent[0], tangent[1], tangent[2]);
        tan.normalize();

        Normal normal = s.getNormal(pt);
        normal.normalize();

        // here is where we get the surface tension coefficient
        // double ST = stencil_m.getScalar(pt);
        double ST = STCoeff + Marangoni[2] *
                                  (Marangoni[0] * pt[0] + Marangoni[1] * pt[1]);
        if (Marangoni[2] == -1.0 &&
            (Marangoni[0] != 0.0 ||
             Marangoni[1] != 0.0)) {  // Force droplet breakup
          double gamma0 = 1.0;        // Base surface tension coefficient
          double R = 1.0;             // Characteristic length scale

          // Linear decrease in surface tension with x-coordinate - Al-Saud
          // Style
          ST = gamma0 * std::max(1 - 1.25 * std::abs(pt[0]), 0.1);
        }
        Normal f = IRL::crossProduct(tan, normal);
        f.normalize();
        f = f * ST;
        total = total + 0.5 * w * f * denom;
      }

      // Debug Prints
      if (false) {  // False added as a toggle for the debug. Make true to
                    // turn on debug.
        // Arc Setup Info
        std::cout << "Start Point: " << startPoint << "\n";
        std::cout << "End Point: " << endPoint << "\n";
        std::cout << "Control Point: " << controlPoint << "\n";
        std::cout << "pairs[j][0]: " << pairs[j][0] << "\n";
        std::cout << "pairs[j][1]: " << pairs[j][1] << "\n";
        for (int i = 0; i < intersectionsSet.size(); i++) {
          std::cout << "intersectionsSet[" << i << "]: ";
          if (!intersectionsSet[i].empty()) {
            std::cout << "(" << intersectionsSet[i][0] << ")";
          } else {
            std::cout << "(empty)";
          }
          std::cout << "\n";
        }
        std::cout << "P0: " << P0 << "\n";
        std::cout << "P1: " << P1 << "\n";
        std::cout << "P2: " << P2 << "\n";
        std::cout << "P3: " << P3 << "\n\n";
        // Output Info
        std::cout << "Total: " << total << "\n";
        std::cout << "===================================================== \n";
      }
    }
  }
  return -total;
}

template <class CellType>
Normal PU<CellType>::solveEdgeCylinder(double STCoeff, Pt& P0, Pt& P1,
                                       double radius, Pt& center,
                                       double delta) {
  // Make Implicit Surface
  // std::cout << "In Solve Edge\n";
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(delta);
  // Find Intersection with edge
  std::vector<Pt> intersections =
      s.intersectEdgeCylinder(P0, P1, radius, center, 10);

  // Calculate some edge properties
  Pt dP = P1 - P0;
  Normal OutwardsNormal = {dP[1], -dP[0], 0};
  double D = std::sqrt(dP[0] * dP[0] + dP[1] * dP[1]);
  double denom = 1.0 / (safelyEpsilon(D));
  // Give space for working variables
  Normal tangent;
  Normal total = {0.0, 0.0, 0.0};  // Total force

  // Loop over intersections and add tangents to force
  if (intersections.size() > 0) {
    for (int j = 0; j < intersections.size(); j++) {
      tangent = s.getTangentCylinder(intersections[j], radius, center);
      tangent.normalize();
      // If facing inside, multiply by negative
      // If facing outside, leave the same
      double Scale = OutwardsNormal * tangent;

      tangent *= Scale / safelyEpsilon(std::abs(Scale));
      // Add
      total = total + STCoeff * denom * tangent;
    }
  }
  return total;
}

template <class CellType>
double PU<CellType>::getValue(double x, double y, double z, double delta) {
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(delta);
  Pt in = {x, y, z};
  double retVal = 0;
  s.evaluate(in, &retVal);
  // std::cout << "==== Return Value: " << retVal << std::endl;
  return retVal;
}

template <class CellType>
Normal PU<CellType>::getTangent(double x, double y, double z, double delta) {
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(delta);
  Pt in = {x, y, z};
  Normal retVal = s.getTangent(in);
  // std::cout << "==== Return Value: " << retVal << std::endl;
  return retVal;
}

template <class CellType>
double PU<CellType>::getWeight(double x, double y, double z, double delta) {
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(delta);
  Pt in = {x, y, z};
  double retVal;
  s.getTotalWeight(in, &retVal);
  // std::cout << "==== Return Value: " << retVal << std::endl;
  return retVal;
}

template <class CellType>
double PU<CellType>::getWeight(Pt& in, double delta) {
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(delta);
  double retVal;
  s.getTotalWeight(in, &retVal);
  // std::cout << "==== Return Value: " << retVal << std::endl;
  return retVal;
}

template <class CellType>
double PU<CellType>::getValueCylinder(double x, double y, double z,
                                      double radius, Pt center) {
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(radius);
  Pt in = {x, y, z};
  double retVal = 0;
  s.evaluateCylinder(in, radius, center, &retVal);
  // std::cout << "==== Return Value: " << retVal << std::endl;
  return retVal;
}

template <class CellType>
inline Normal PU<CellType>::getTangentCylinder(double x, double y, double z,
                                               double radius, Pt center) {
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(radius);
  Pt in = {x, y, z};
  Normal retVal = s.getTangentCylinder(in, radius, center);
  // std::cout << "==== Return Value: " << retVal << std::endl;
  return retVal;
}

template <class CellType>
double PU<CellType>::getWeightCylinder(double x, double y, double z,
                                       double radius, Pt center) {
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(radius);
  Pt in = {x, y, z};
  double retVal = 1.0;
  // std::cout << "==== Return Value: " << retVal << std::endl;
  return retVal;
}

template <class CellType>
void PU<CellType>::printSolver() {
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(1);
  // First Print Solver
  s.printSurface();
  std::cout << "> Threshold = " << intersection_threshold_m << "\n";
}

template <class CellType>
void PU<CellType>::setNeighborhood(PUNeighborhood<CellType> stencil_) {
  stencil_m = stencil_;
}

template <class CellType>
void PU<CellType>::setThreshold(double thresh_) {
  intersection_threshold_m = thresh_;
}

template <class CellType>
Paraboloid PU<CellType>::solve(
    const PUNeighborhood<CellType>* a_neighborhood_pointer,
    const Pt& a_centroid, const double a_delta) {
  double delta = a_delta;
  if (delta < 0.0) {
    const auto cell = a_neighborhood_pointer->getCenterCell();
    delta = 2.5 * std::pow(cell.calculateVolume(), 1.0 / 3.0);
  }

  this->setNeighborhood(*a_neighborhood_pointer);
  auto PUSurface = this->neighborhoodToImplicitSurface(delta);

  // Project provided point onto the PU surface
  const auto pt_on_PU = PUSurface.projectOntoPU(a_centroid);

  // Compute local gradient and hessian of PU approximation
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> holdsGradAndHessian;
  PUSurface.evaluate(pt_on_PU, &holdsGradAndHessian);
  const auto gradF = std::get<1>(holdsGradAndHessian);
  const auto hessF = std::get<2>(holdsGradAndHessian);

  // Return paraboloid computed from derivatives
  return IRL::Paraboloid::fromDerivatives(pt_on_PU, gradF, hessF);
}

template <class CellType>
IRL::Pt PU<CellType>::projectOntoPU(const Pt& a_pt, const double a_delta) {
  double delta = a_delta;
  if (delta < 0.0) {
    const auto cell = stencil_m.getCenterCell();
    delta = 2.5 * std::pow(cell.calculateVolume(), 1.0 / 3.0);
  }

  this->setNeighborhood(stencil_m);
  auto PUSurface = this->neighborhoodToImplicitSurface(delta);

  // Project provided point onto the PU surface
  return PUSurface.projectOntoPU(a_pt);
}

// Getting Curvature
template <class CellType>
double PU<CellType>::getCurvature(Pt& a_pt, double a_delta) {
  double delta = a_delta;
  if (delta < 0.0) {
    const auto cell = stencil_m.getCenterCell();
    delta = 2.5 * std::pow(cell.calculateVolume(), 1.0 / 3.0);
  }
  this->setNeighborhood(stencil_m);
  auto PUSurface = this->neighborhoodToImplicitSurface(delta);
  return PUSurface.getCurvature(a_pt);
}

template <class CellType>
double PU<CellType>::getCurvature(double x, double y, double z, double delta) {
  double delta_ = delta;
  if (delta_ < 0.0) {
    const auto cell = stencil_m.getCenterCell();
    delta_ = 2.5 * std::pow(cell.calculateVolume(), 1.0 / 3.0);
  }
  this->setNeighborhood(stencil_m);
  auto PUSurface = this->neighborhoodToImplicitSurface(delta_);
  Pt in = {x, y, z};
  return PUSurface.getCurvature(in);
}

}  // End Namespace IRL

#endif
