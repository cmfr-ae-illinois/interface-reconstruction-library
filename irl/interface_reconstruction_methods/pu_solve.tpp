#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_SOLVE_TPP_

#include <limits>
#include <tuple>
#include <vector>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
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
      F = n * x;
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
      F = n * x;
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
      F = n * x;
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
    double weight = std::get<0>(weightRet);
    Eigen::Vector3d grad_weight = std::get<1>(weightRet);
    Eigen::Matrix3d hess_weight = std::get<2>(weightRet);

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
    double weight = std::get<0>(weightRet);
    Eigen::Vector3d grad_weight = std::get<1>(weightRet);

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
    // Add results to sums
    weight_sum += weight;

    // How get val,grad,hess of separator
    double F = 0.0;
    PUImplicitSurface::implicitSeparator(x, centroids[i], &separators[i], &F);
    // Now calculate F_sum, grad_product_sum, and hess_product_sum
    F_sum += weight * F;
  }
  // Put Everything Together and return.
  *retVal = F_sum / weight_sum;
}

// template <class SeparatorType>
inline std::vector<Pt> PUImplicitSurface::intersectEdge(
    const Pt& x0, const Pt& x1, const int& Npartitions) {
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
      // Add Intersection
      intersections.push_back(midX);
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

inline double PUImplicitSurface::getCurvature(Pt& x) {
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> holdsGradAndHessian;
  this->evaluate(x, &holdsGradAndHessian);
  auto gradF = std::get<1>(holdsGradAndHessian);
  auto hessF = std::get<2>(holdsGradAndHessian);

  double Fxx = hessF(0, 0);
  double Fyy = hessF(1, 1);
  double Fxy = hessF(0, 1);
  double Fx = gradF(0);
  double Fy = gradF(1);

  double numer = Fxx * Fy * Fy - 2 * Fxy * Fx * Fy + Fx * Fx * Fyy;
  double magGradF = std::sqrt(Fx * Fx + Fy * Fy);
  double denom = magGradF * magGradF * magGradF;

  double kz = -numer / denom;

  return kz;
}

inline const Pt PUImplicitSurface::projectOntoPU(const Pt& a_pt) {
  std::pair<double, Eigen::Vector3d> holdsGrad;
  Pt projected_pt = a_pt;
  const int itmax = 5;
  for (int i = 0; i < itmax; i++) {
    this->evaluate(a_pt, &holdsGrad);
    const auto F = std::get<0>(holdsGrad);
    const auto gradF = std::get<1>(holdsGrad);
    const double grad_norm_inv = 1.0 / safelyEpsilon(gradF.squaredNorm());
    for (int d = 0; d < 3; d++) {
      projected_pt[d] -= F * gradF(d) * grad_norm_inv;
    }
  }
  return projected_pt;
}

// ============== Solver Methods
template <class CellType>
PUST<CellType>::PUST(void) {}

template <class CellType>
PUST<CellType>::PUST(PUSTNeighborhood<CellType> stencil_) {
  this->stencil_m = stencil_;
  // this->surface_m = this->neighborhoodToImplicitSurface(5.0);
}

template <class CellType>
PUImplicitSurface PUST<CellType>::neighborhoodToImplicitSurface(double delta) {
  const auto centroids = stencil_m.getCentroids();
  const auto separators = stencil_m.getSeparators();

  return PUImplicitSurface(centroids, separators, delta);
}

template <class CellType>
Normal PUST<CellType>::solveEdge(double STCoeff, Pt& P0, Pt& P1) {
  // Make Implicit Surface
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(5.0);
  // Find Intersection with edge
  std::vector<Pt> intersections = s.intersectEdge(P0, P1, 10);

  // Calculate some edge properties
  Pt dP = P1 - P0;
  double D = std::sqrt(dP[0] * dP[0] + dP[1] * dP[1]);
  double denom = 1.0 / (safelyEpsilon(D));
  // Give space for working variables
  Normal tangent;
  Normal total = {0.0, 0.0, 0.0};  // Total force

  // Loop over intersections and add tangents to force
  if (intersections.size() > 0) {
    for (int j = 0; j < intersections.size(); j++) {
      tangent = s.getTangent(intersections[j]);
      tangent.normalize();
      total = total + STCoeff * denom * tangent;
    }
  }
  return total;
}

template <class CellType>
void PUST<CellType>::setNeighborhood(PUSTNeighborhood<CellType> stencil_) {
  stencil_m = stencil_;
}

template <class CellType>
Paraboloid PUST<CellType>::solve(
    const PUSTNeighborhood<CellType>* a_neighborhood_pointer,
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

}  // End Namespace IRL

#endif
