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
// Evalute Function for Ellipsoid
// ===============================================
inline void PUImplicitSurface::evaluateEllipsoid(
    Pt& x, const Normal& column1, const Normal& column2, const Normal& column3,
    const Pt& center,
    std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d>* retVal) {
  Eigen::Matrix3d A;
  // Assign columns from input to columns of A
  A(0, 0) = column1[0];
  A(1, 0) = column1[1];
  A(2, 0) = column1[2];
  A(0, 1) = column2[0];
  A(1, 1) = column2[1];
  A(2, 1) = column2[2];
  A(0, 2) = column3[0];
  A(1, 2) = column3[1];
  A(2, 2) = column3[2];
  // We know that the value is (x - center)^T * A * (x - center) - 1
  Eigen::Vector3d offset(x[0] - center[0], x[1] - center[1], x[2] - center[2]);
  double F = 1.0 - offset.transpose() * A * offset;
  // We know gradient is (A+A^T) * (x - center)
  Eigen::Vector3d gradF = -(A + A.transpose()) * offset;
  // We know hessian is (A + A^T)
  Eigen::Matrix3d hessF = -(A + A.transpose());
  *retVal = std::make_tuple(F, gradF, hessF);
}

inline void PUImplicitSurface::evaluateEllipsoid(
    Pt& x, const Normal& column1, const Normal& column2, const Normal& column3,
    const Pt& center, std::pair<double, Eigen::Vector3d>* retVal) {
  // Just like above, fill value and gradient.
  Eigen::Matrix3d A;
  // Assign columns from input to columns of A
  A(0, 0) = column1[0];
  A(1, 0) = column1[1];
  A(2, 0) = column1[2];
  A(0, 1) = column2[0];
  A(1, 1) = column2[1];
  A(2, 1) = column2[2];
  A(0, 2) = column3[0];
  A(1, 2) = column3[1];
  A(2, 2) = column3[2];
  // We know that the value is (x - center)^T * A * (x - center) - 1
  Eigen::Vector3d offset(x[0] - center[0], x[1] - center[1], x[2] - center[2]);
  double F = 1.0 - offset.transpose() * A * offset;
  // We know gradient is (A+A^T) * (x - center)
  Eigen::Vector3d gradF = -(A + A.transpose()) * offset;
  // Return
  *retVal = std::make_pair(F, gradF);
}

inline void PUImplicitSurface::evaluateEllipsoid(Pt& x, const Normal& column1,
                                                 const Normal& column2,
                                                 const Normal& column3,
                                                 const Pt& center,
                                                 double* retVal) {
  // Just like above, fill value.
  Eigen::Matrix3d A;
  // Assign columns from input to columns of A
  A(0, 0) = column1[0];
  A(1, 0) = column1[1];
  A(2, 0) = column1[2];
  A(0, 1) = column2[0];
  A(1, 1) = column2[1];
  A(2, 1) = column2[2];
  A(0, 2) = column3[0];
  A(1, 2) = column3[1];
  A(2, 2) = column3[2];
  // We know that the value is (x - center)^T * A * (x - center) - 1
  Eigen::Vector3d offset(x[0] - center[0], x[1] - center[1], x[2] - center[2]);
  double F = 1.0 - offset.transpose() * A * offset;
  // Return
  *retVal = F;
}

// Other Cylinder Functions
// =======================================================
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
        blocked = true;
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

inline Normal PUImplicitSurface::getNormalEllipsoid(Pt& x,
                                                    const Normal& column1,
                                                    const Normal& column2,
                                                    const Normal& column3,
                                                    const Pt& center) {
  std::pair<double, Eigen::Vector3d> holdsGrad;
  this->evaluateEllipsoid(x, column1, column2, column3, center, &holdsGrad);
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

  return -numer / safelyEpsilon(denom);
}

inline double PUImplicitSurface::getPlaneCurvature(
    Pt& x, Normal& planeNormal) {  // Verified
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> holdsGradAndHessian;
  this->evaluate(x, &holdsGradAndHessian);
  auto gradF = std::get<1>(holdsGradAndHessian);
  auto hessF = std::get<2>(holdsGradAndHessian);
  // Print Gradient and Hessian
  // std::cout << "=============================== In Plane Curvature Ellipsoid
  // "
  //              "===============================\n";
  // std::cout << "Point: " << x << "\n";
  // std::cout << "Gradient: " << gradF << "\n";
  // std::cout << "Hessian: \n" << hessF << "\n";

  double magGradF = gradF.norm();
  Normal surfaceNormal = this->getNormal(x);
  surfaceNormal.normalize();
  planeNormal.normalize();
  // Get Surface Tangent in Plane
  Normal gradFVec(gradF(0), gradF(1), gradF(2));
  Normal tangent = IRL::crossProduct(gradFVec, planeNormal);
  tangent.normalize();
  // std::cout << "Tangent: " << tangent << "\n";
  // Convert to Eigen::Vector3d
  Eigen::Vector3d tangentVec(tangent[0], tangent[1], tangent[2]);

  // tangent^T Hess(F) tangent
  double tangentHessTangent = tangentVec.transpose() * hessF * tangentVec;
  // std::cout << "Tangent Hess Tangent: " << tangentHessTangent << "\n";
  double numer = tangentHessTangent;
  // Denominator calc
  Normal m = IRL::crossProduct(planeNormal, tangent);
  // std::cout << "Plane Normal: " << planeNormal << "\n";
  // std::cout << "M: " << m << "\n";
  double denom = gradFVec * m;  // gradFVec dot (planeNormal cross tangent)
  // std::cout << "=============END ELIPSOID PLANE CURVATURE================\n";
  return std::abs(numer) / std::abs(safelyEpsilon(denom));
}

inline double PUImplicitSurface::getMeanCurvatureEllipsoid(
    Pt& x, const Normal& column1, const Normal& column2, const Normal& column3,
    const Pt& center) {
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> holdsGradAndHessian;
  this->evaluateEllipsoid(x, column1, column2, column3, center,
                          &holdsGradAndHessian);
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

inline double PUImplicitSurface::getPlaneCurvatureEllipsoid(
    Pt& x, const Normal& column1, const Normal& column2, const Normal& column3,
    const Pt& center, Normal& planeNormal) {  // Verified
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> holdsGradAndHessian;
  this->evaluateEllipsoid(x, column1, column2, column3, center,
                          &holdsGradAndHessian);
  auto gradF = std::get<1>(holdsGradAndHessian);
  auto hessF = std::get<2>(holdsGradAndHessian);
  // Print Gradient and Hessian
  // std::cout << "=============================== In Plane Curvature Ellipsoid
  // "
  //              "===============================\n";
  // std::cout << "Point: " << x << "\n";
  // std::cout << "Gradient: " << gradF << "\n";
  // std::cout << "Hessian: \n" << hessF << "\n";

  double magGradF = gradF.norm();
  Normal surfaceNormal =
      -this->getNormalEllipsoid(x, column1, column2, column3, center);
  surfaceNormal.normalize();
  planeNormal.normalize();
  // Get Surface Tangent in Plane
  Normal gradFVec(gradF(0), gradF(1), gradF(2));
  Normal tangent = IRL::crossProduct(gradFVec, planeNormal);
  tangent.normalize();
  // std::cout << "Tangent: " << tangent << "\n";
  // Convert to Eigen::Vector3d
  Eigen::Vector3d tangentVec(tangent[0], tangent[1], tangent[2]);

  // tangent^T Hess(F) tangent
  double tangentHessTangent = tangentVec.transpose() * hessF * tangentVec;
  // std::cout << "Tangent Hess Tangent: " << tangentHessTangent << "\n";
  double numer = tangentHessTangent;
  // Denominator calc
  Normal m = IRL::crossProduct(planeNormal, tangent);
  // std::cout << "Plane Normal: " << planeNormal << "\n";
  // std::cout << "M: " << m << "\n";
  double denom = gradFVec * m;  // gradFVec dot (planeNormal cross tangent)
  // std::cout << "=============END ELIPSOID PLANE CURVATURE================\n";
  return std::abs(numer) / std::abs(safelyEpsilon(denom));
}

inline const Pt PUImplicitSurface::projectOntoPU(const Pt& a_pt) {
  std::pair<double, Eigen::Vector3d> holdsGrad;
  Pt projected_pt = a_pt;
  const int itmax = 50;
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
  }
  return projected_pt;
}

inline const Pt PUImplicitSurface::projectOntoEllipsoid(const Pt& a_pt,
                                                        const Normal& column1,
                                                        const Normal& column2,
                                                        const Normal& column3,
                                                        const Pt& center) {
  std::pair<double, Eigen::Vector3d> holdsGrad;
  Pt projected_pt = a_pt;
  const int itmax = 50;
  for (int i = 0; i < itmax; i++) {
    this->evaluateEllipsoid(projected_pt, column1, column2, column3, center,
                            &holdsGrad);
    const auto F = std::get<0>(holdsGrad);
    const auto gradF = std::get<1>(holdsGrad);
    const double grad_norm_inv = 1.0 / safelyEpsilon(gradF.squaredNorm());

    auto step = -F * gradF * grad_norm_inv;
    double alpha = 1.0;
    for (int d = 0; d < 3; d++) {
      projected_pt[d] += alpha * step(d);
    }

    // Break if F less than a tolernace
    if (std::abs(F) < 1e-12) {
      break;
    }
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
  bool blocked = false;
  std::vector<Pt> intersections =
      s.intersectEdge(P0, P1, 10, intersection_threshold_m, blocked);

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
  // std::cout << "=========================== Solve Face \n" << std::endl;
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
  // std::cout << "Corner Values: " << val0 << "," << val1 << "," << val2 << ","
  // << val3 << "\n";
  // std::cout << "Corner Weights: " << weight0 << "," << weight1 << "," <<
  // weight2
  //           << "," << weight3 << "\n";
  Pt dP = P1 - P0;
  double D = std::sqrt(dP[0] * dP[0] + dP[1] * dP[1] + dP[2] * dP[2]);
  double denom = 1.0 / (safelyEpsilon(D));
  s.evaluate(0.25 * (P0 + P1 + P2 + P3), &val4);  // Center Point Value
  // Calculate case:
  int caseValue =
      1 * (-val0 > 0) + 2 * (-val1 > 0) + 4 * (-val2 > 0) + 8 * (-val3 > 0);
  // If all corners are outside or inside, return 0 force. Otehrwise, pair
  // intersections
  // std::cout << "Case Value: " << caseValue << "\n";
  std::vector<std::vector<int>> pairs = {{-1, -1}, {-1, -1}};
  // std::cout << "================== Case Value: " << caseValue << "\n";
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
  bool blocked01, blocked12, blocked23, blocked30;
  // std::cout << "Intersections01" << std::endl;
  std::vector<Pt> intersections01 =
      s.intersectEdge(P0, P1, 1, intersection_threshold_m, blocked01);
  // std::cout << "Intersections12" << std::endl;
  std::vector<Pt> intersections12 =
      s.intersectEdge(P1, P2, 1, intersection_threshold_m, blocked12);
  // std::cout << "Intersections23" << std::endl;
  std::vector<Pt> intersections23 =
      s.intersectEdge(P2, P3, 1, intersection_threshold_m, blocked23);
  // std::cout << "Intersections30" << std::endl;
  std::vector<Pt> intersections30 =
      s.intersectEdge(P3, P0, 1, intersection_threshold_m, blocked30);
  if (blocked01 || blocked12 || blocked23 || blocked30) {
    // std::cout << "Blocked Intersection Detected" << std::endl;
    return Normal(0.0, 0.0, 0.0);
  }
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
      // std::cout
      //     << "===================== In Face Solver =====================\n";
      // std::cout << "Start Point: " << startPoint << "\n";
      // std::cout << "End Point: " << endPoint << "\n";
      // std::cout << "Tangent Start: " << tangentStart << "\n";
      // std::cout << "Tangent End: " << tangentEnd << "\n";
      // std::cout << "Control Point 1: " << arc.control_point_1() << "\n";
      // std::cout << "Control Point 2: " << arc.control_point_2() << "\n";
      // std::cout << "weight 1: " << arc.weight(1) << "\n";
      // std::cout << "weight 1: " << arc.weight(2) << "\n";
      int QuadRuleOrder = 100;
      const auto& abscissea = AbscissaeGauss<double, 100>();
      const auto& weights = WeightsGauss<double, 100>();
      for (int j = 0; j < QuadRuleOrder; ++j) {
        const double t = 0.5 * (1.0 + abscissea[j]);
        const double w = weights[j];

        // Evaluation Points
        Pt pt = arc.point(t);
        // get normal,tangent, and surface tension at each point.
        Pt tangent = arc.derivative(t);
        Normal tan = Normal(tangent[0], tangent[1], tangent[2]);
        double speed =
            std::sqrt(tan[0] * tan[0] + tan[1] * tan[1] + tan[2] * tan[2]);
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
        total = total + 0.5 * w * f * denom * speed;
      }

      // Debug Prints
      if (false) {  // False added as a toggle for the debug. Make true to
                    // turn on debug.
        // Arc Setup Info
        std::cout << "Start Point: " << startPoint << "\n";
        std::cout << "End Point: " << endPoint << "\n";
        std::cout << "Control Point 1: " << arc.control_point_1() << "\n";
        std::cout << "Control Point 2: " << arc.control_point_2() << "\n";
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

// Solve Edge for Ellipse  ===================================
template <class CellType>
Normal PU<CellType>::solveFaceEllipsoid(
    const double STin, const Pt& P0, const Pt& P1, const Pt& P2, const Pt& P3,
    const Normal& column1, const Normal& column2, const Normal& column3,
    const Pt& center, const double Pressure, const Normal& Marangoni) {
  // The Marangoni normal object holds the Xgradient, then the Y gradient, then
  // the temperature gradient of ST Marangoni = [Gx,Gy,sigma_T] Make Implicit
  // Surface std::cout << "In Solve Edge\n";

  // Something is up with paraboloids because they take like 8x the time to run.
  double STCoeff = STin;
  // The Pressure Option tells us if we should include the pressure terms or
  // not.
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(1.0);
  std::vector<Pt> polygon = {P0, P1, P2, P3, P0};
  std::vector<Pt> intersection = {};
  // Reconstruct A
  Eigen::Matrix3d A;
  // Assign columns from input to columns of A
  A(0, 0) = column1[0];
  A(1, 0) = column1[1];
  A(2, 0) = column1[2];
  A(0, 1) = column2[0];
  A(1, 1) = column2[1];
  A(2, 1) = column2[2];
  A(0, 2) = column3[0];
  A(1, 2) = column3[1];
  A(2, 2) = column3[2];
  // Loop over edges and find intersection of ellipsoid with each edge.
  int firstIntersectionIndex = -1;
  for (int j = 0; j < polygon.size() - 1; j++) {
    // Extract Values
    Pt P0temp = polygon[j];
    Pt P1temp = polygon[j + 1];
    // Offsets
    Pt D0 = P0temp - center;
    Pt DP = P1temp - P0temp;
    // Make D0 and DP into Eigen::Vector3d
    Eigen::Vector3d D0_eigen(D0[0], D0[1], D0[2]);
    Eigen::Vector3d DP_eigen(DP[0], DP[1], DP[2]);

    // Coefficients for quadratic equation
    double a = (DP_eigen.transpose() * A * DP_eigen).value();
    double b = (D0_eigen.transpose() * A * DP_eigen).value() +
               (DP_eigen.transpose() * A * D0_eigen).value();
    double c = (D0_eigen.transpose() * A * D0_eigen).value() - 1.0;
    // Check Determinant
    double discriminant = b * b - 4 * a * c;

    // std::cout << "\n\n\n=============================\nROOTS FINDING DEBUG "
    //              "INFO - "
    //           << j << ":" << "\n";
    // std::cout << "a: " << a << ", b: " << b << ", c: " << c
    //           << ", discriminant: " << discriminant << "\n";
    // std::cout << "P0temp: " << P0temp << ", P1temp: " << P1temp
    //           << ", center: " << center << "\n";
    // std::cout << "D0: " << D0_eigen << "\n DP: " << DP_eigen << "\n";
    // std::cout << "A: \n" << A << "\n";
    // std::cout << "============================ ";
    if (2 * a <= 1e-12) {
      if (false) {  // Debug
        std::cout << "Warning: a is very small. PU_solve line 1315\n";
        std::cout << "a: " << a << ", b: " << b << ", c: " << c
                  << ", discriminant: " << discriminant << "\n";
        std::cout << "P0temp: " << P0temp << ", P1temp: " << P1temp
                  << ", center: " << center << "\n";
        std::cout << "D0: " << D0_eigen << ", DP: " << DP_eigen << "\n";
        std::cout << "A: \n" << A << "\n";
        std::cout << "============================ ";
      }

      // Solve Linear Case
      if (std::abs(b) > 1e-12) {
        double t = -c / b;
        if (t >= 0 && t <= 1) {
          Pt intersectionPoint = P0temp + t * DP;
          intersection.push_back(intersectionPoint);
          if (firstIntersectionIndex == -1) {
            firstIntersectionIndex = j;
          }
        }
      }
    } else if (discriminant >= 0.0) {
      // Calculate roots
      double sqrt_discriminant = std::sqrt(discriminant);
      double t1 = (-b - sqrt_discriminant) / safelyEpsilon(2 * a);
      double t2 = (-b + sqrt_discriminant) / safelyEpsilon(2 * a);
      // std::cout << "Roots: t1 = " << t1 << ", t2 = " << t2 << "\n";

      // Check if roots are within [0, 1]
      if (t1 >= 0 && t1 <= 1) {
        Pt intersectionPoint = P0temp + t1 * DP;
        intersection.push_back(intersectionPoint);
      }
      if (t2 >= 0 && t2 <= 1) {
        Pt intersectionPoint = P0temp + t2 * DP;
        intersection.push_back(intersectionPoint);
      }

      // If this is the first intersecting edge, store the index of the first
      // intersection point
      if (firstIntersectionIndex == -1 && !intersection.empty()) {
        firstIntersectionIndex = j;
      }
    }
  }

  // Calculate some edge and plane properties
  Pt dP1 = P1 - P0;
  Pt dP2 = P2 - P1;
  Normal edge1 = {dP1[0], dP1[1], dP1[2]};
  Normal edge2 = {dP2[0], dP2[1], dP2[2]};

  double D = std::sqrt(dP1[0] * dP1[0] + dP1[1] * dP1[1] + dP1[2] * dP1[2]);
  double denom = 1.0 / (safelyEpsilon(D));
  edge1.normalize();
  edge2.normalize();
  Normal faceNormal = IRL::crossProduct(edge1, edge2);
  faceNormal.normalize();

  // After getting all intersection points, we can pair them one after the
  // other. However, we need to ensure that we have an properly oriented points.
  // Check the outward pointing normal at the first point. If it is the same
  // direction at the first edge, then put it at the back of the list.
  // Otherwise, leave it at the front
  if (intersection.size() >= 2) {
    Normal normalAtFirstIntersection = -s.getNormalEllipsoid(
        intersection[0], column1, column2, column3, center);
    // std::cout << "First Intersection: " << intersection[0] << "\n";
    // std::cout << "First Intersection Index: " << firstIntersectionIndex <<
    // "\n"; std::cout << "First Intersection Edge: " <<
    // polygon[firstIntersectionIndex]
    // << " to " << polygon[firstIntersectionIndex + 1] << "\n";
    // std::cout << "Normal at First Intersection: " <<
    // normalAtFirstIntersection
    // << "\n";
    normalAtFirstIntersection.normalize();
    // std::cout << "Normal at First Intersection Normalized: "
    // << normalAtFirstIntersection << "\n";
    Pt dP_edge =
        polygon[firstIntersectionIndex + 1] - polygon[firstIntersectionIndex];
    Normal edgeDirection = {dP_edge[0], dP_edge[1], dP_edge[2]};
    edgeDirection.normalize();
    // std::cout << "Edge Direction: " << edgeDirection << "\n";
    double dotProduct = normalAtFirstIntersection * edgeDirection;
    // std::cout << "Dot Product: " << dotProduct << "\n";
    if (dotProduct < 0) {
      Pt firstIntersection = intersection[0];
      intersection.erase(intersection.begin());   // Remove the first point
      intersection.push_back(firstIntersection);  // Add it to the back
    }
    // std::cout << "face normal: " << faceNormal << "\n";
  }

  // Give space for working variables
  Normal normal1, normal2, tangentStart, tangentEnd;
  Pt startPoint, endPoint, controlPoint;
  double weight = 1;  // Make Nonrational for now, but can add later if we want
  Normal total = {0.0, 0.0, 0.0};  // Total force

  // Loop over pairs, make splines, integrate, add forces to total.
  // std::cout << "Intersection Size: " << intersection.size() << "\n";
  std::vector<Normal> tangents = {};
  for (int j = 0; j + 1 < intersection.size(); j += 2) {
    startPoint = intersection[j];
    endPoint = intersection[j + 1];

    // std::cout << "Start Point: " << startPoint << "\n";
    // std::cout << "End Point: " << endPoint << "\n";
    normal1 =
        -s.getNormalEllipsoid(startPoint, column1, column2, column3, center);
    normal1.normalize();

    normal2 =
        -s.getNormalEllipsoid(endPoint, column1, column2, column3, center);
    normal2.normalize();

    // Make tangent by crossing normal with face  normal
    tangentStart = IRL::crossProduct(faceNormal, normal1);
    tangentStart.normalize();
    tangentEnd = IRL::crossProduct(faceNormal, normal2);
    tangentEnd.normalize();

    tangents.push_back(tangentStart);
    tangents.push_back(tangentEnd);
    // Compute Start Curvature
    double curvature_start = s.getPlaneCurvatureEllipsoid(
        startPoint, column1, column2, column3, center, faceNormal);
    // std::cout << "Curvature Start: " << curvature_start << "\n";
    RationalBezierArc arc =
        RationalBezierArc(startPoint, tangentStart, endPoint, tangentEnd,
                          faceNormal, curvature_start);
    controlPoint = arc.control_point();
    weight = arc.weight();
    int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, 100>();
    const auto& weights = WeightsGauss<double, 100>();
    // Verified Up to here with Mathematica.
    for (int k = 0; k < QuadRuleOrder; k++) {
      const double t = 0.5 * (1.0 + abscissea[k]);
      const double w = weights[k];

      // Evaluation Points
      Pt pt = arc.point(t);
      // get normal,tangent, and surface tension at each point.
      Pt tangent = arc.derivative(t);
      Normal tan = Normal(tangent[0], tangent[1], tangent[2]);
      double speed =
          std::sqrt(tan[0] * tan[0] + tan[1] * tan[1] + tan[2] * tan[2]);
      tan.normalize();

      Normal normal =
          -s.getNormalEllipsoid(pt, column1, column2, column3, center);

      normal.normalize();

      // here is where we get the surface tension coefficient
      // double ST = stencil_m.getScalar(pt);
      double ST = STCoeff +
                  Marangoni[2] * (Marangoni[0] * pt[0] + Marangoni[1] * pt[1]);
      if (Marangoni[2] == -1.0 &&
          (Marangoni[0] != 0.0 ||
           Marangoni[1] != 0.0)) {  // Force droplet breakup
        double gamma0 = 1.0;        // Base surface tension coefficient
        double R = 1.0;             // Characteristic length scale

        // Linear decrease in surface tension with x-coordinate - Al-Saud
        // Style
        ST = gamma0 * std::max(1 - 1.25 * std::abs(pt[0]), 0.1);
      }
      ST = STCoeff;
      Normal f = IRL::crossProduct(tan, normal);
      f.normalize();
      f = f * ST;
      total = total + 0.5 * w * f * denom * speed;
      // std::cout << " ========================== \n";
      // std::cout << "Normal: " << normal << "\n";
      // std::cout << "Tangent: " << tan << "\n";
      // std::cout << "ST: " << ST << "\n";
      // std::cout << "Force: " << f << "\n";
      // std::cout << "Total: " << total << "\n";
      // std::cout << "Denom: " << denom << "\n";
      // std::cout << " ========================== \n";
    }

    // Debug Prints
    if (false) {  // False added as a toggle for the debug. Make true to
                  // turn on debug.
      // Arc Setup Info
      std::cout << "Start Point: " << startPoint << "\n";
      std::cout << "End Point: " << endPoint << "\n";
      std::cout << "Control Point: " << controlPoint << "\n";
      std::cout << "weight: " << weight << "\n";
      for (int i = 0; i < intersection.size(); i++) {
        std::cout << "intersection[" << i << "]: ";
        std::cout << "(" << intersection[i] << ")";
        std::cout << "\n";
      }

      for (int i = 0; i < tangents.size(); i++) {
        std::cout << "Tangent[" << i << "]: ";
        std::cout << "(" << tangents[i] << ")";
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
  return -total;
}

// Other
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

template <class CellType>
IRL::Pt PU<CellType>::projectOntoEllipsoid(const Pt& a_pt,
                                           const Normal& column1,
                                           const Normal& column2,
                                           const Normal& column3,
                                           const Pt& center) {
  // Implementation for projecting onto ellipsoid
  this->setNeighborhood(stencil_m);
  auto PUSurface = this->neighborhoodToImplicitSurface(
      1.0);  // Use a default delta for ellipsoid projection

  // Project provided point onto the PU surface
  return PUSurface.projectOntoEllipsoid(a_pt, column1, column2, column3,
                                        center);
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

template <class CellType>
double PU<CellType>::getMeanCurvatureEllipsoid(double x, double y, double z,
                                               const Normal& column1,
                                               const Normal& column2,
                                               const Normal& column3,
                                               const Pt& center) {
  double delta_ = 1.0;

  this->setNeighborhood(stencil_m);
  auto PUSurface = this->neighborhoodToImplicitSurface(delta_);
  Pt in = {x, y, z};
  return PUSurface.getMeanCurvatureEllipsoid(in, column1, column2, column3,
                                             center);
}

template <class CellType>
double PU<CellType>::getMeanCurvatureEllipsoid(Pt& a_pt, const Normal& column1,
                                               const Normal& column2,
                                               const Normal& column3,
                                               const Pt& center) {
  double delta = 1.0;
  this->setNeighborhood(stencil_m);
  auto PUSurface = this->neighborhoodToImplicitSurface(delta);
  return PUSurface.getMeanCurvatureEllipsoid(a_pt, column1, column2, column3,
                                             center);
}
// Getting Normal
template <class CellType>
Normal PU<CellType>::getNormal(double x, double y, double z, double delta) {
  double delta_ = delta;
  if (delta_ < 0.0) {
    const auto cell = stencil_m.getCenterCell();
    delta_ = 2.5 * std::pow(cell.calculateVolume(), 1.0 / 3.0);
  }
  this->setNeighborhood(stencil_m);
  auto PUSurface = this->neighborhoodToImplicitSurface(delta_);
  Pt in = {x, y, z};
  return PUSurface.getNormal(in);
}

template <class CellType>
Normal PU<CellType>::getNormalEllipsoid(double x, double y, double z,
                                        const Normal& column1,
                                        const Normal& column2,
                                        const Normal& column3,
                                        const Pt& center) {
  double delta_ = 1.0;

  this->setNeighborhood(stencil_m);
  auto PUSurface = this->neighborhoodToImplicitSurface(delta_);
  Pt in = {x, y, z};
  return PUSurface.getNormalEllipsoid(in, column1, column2, column3, center);
}

}  // End Namespace IRL

#endif
