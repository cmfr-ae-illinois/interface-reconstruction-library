// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/pu.h"

namespace IRL {
Paraboloid PU::solve(const PUNeighborhood* a_neighborhood_pointer,
                     const double a_delta) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_m = a_neighborhood_pointer;
  delta_m = a_delta;
  return this->solve();
}

std::pair<double, Eigen::Vector3d> PU::getPUAndGrad(const Pt& a_pt) {
  // number of interfaces
  const UnsignedIndex_t ninterfaces = neighborhood_m->size();

  // information from neighborhood
  std::vector<SeparatorVariant> separators(ninterfaces);
  std::vector<Pt> centroids(ninterfaces);
  std::vector<double> weights(ninterfaces);
  for (UnsignedIndex_t i = 0; i < ninterfaces; ++i) {
    separators[i] = neighborhood_m->getSeparator(i);
    centroids[i] = neighborhood_m->getCentroid(i);
    weights[i] = neighborhood_m->getWeight(i);
  }

  const double inv_delta = 1. / delta_m;
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d gradwxF_plus_wxgradF_sum = Eigen::Vector3d::Zero();
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();

  // loop over interfaces
  for (UnsignedIndex_t i = 0; i < ninterfaces; ++i) {
    double vfrac_and_area_weight = weights[i];
    const Pt x = a_pt - centroids[i];
    const double r = IRL::magnitude(x);
    const double rhat = r * inv_delta;
    if (rhat < 1.0) {  // Then weight is non zero;
      const double weight = vfrac_and_area_weight * (1. + 4. * rhat) *
                            (1. - rhat) * (1. - rhat) * (1. - rhat) *
                            (1. - rhat);
      const Eigen::Vector3d grad_weight =
          Eigen::Vector3d(x[0], x[1], x[2]) * vfrac_and_area_weight * 20. *
          (rhat - 1.) * (rhat - 1.) * (rhat - 1.) * inv_delta * inv_delta;
      weight_sum += weight;
      grad_weight_sum += grad_weight;

      // computing PU value and gradient
      if (const PlanarSeparator* separator = std::get_if<IRL::PlanarSeparator>(
              &(separators[i]))) {  // If plane
        if (separator->getNumberOfPlanes() > 0) {
          const Plane plane = (*separator)[0];
          // Compute plane normal
          const Normal n = plane.normal();
          const double F = n * x;
          const Eigen::Vector3d gradF = Eigen::Vector3d(n[0], n[1], n[2]);
          F_sum += weight * F;
          gradwxF_plus_wxgradF_sum += grad_weight * F + weight * gradF;
        }
      } else if (const Paraboloid* paraboloid = std::get_if<IRL::Paraboloid>(
                     &(separators[i]))) {  // If paraboloid
        // Get paraboloid properties
        const ReferenceFrame& frame = paraboloid->getReferenceFrame();
        const double a = paraboloid->getAlignedParaboloid().a();
        const double b = paraboloid->getAlignedParaboloid().b();
        // Move sample point to local frame
        const Pt tmp = a_pt - paraboloid->getDatum();
        Pt xloc;
        for (int d = 0; d < 3; ++d) {
          xloc[d] = frame[d] * tmp;
        }
        const Eigen::Vector3d e0(frame[0][0], frame[0][1], frame[0][2]);
        const Eigen::Vector3d e1(frame[1][0], frame[1][1], frame[1][2]);
        const Eigen::Vector3d e2(frame[2][0], frame[2][1], frame[2][2]);
        const double dist_norm =
            1. / std::sqrt(1. + 4. * a * a * xloc[0] * xloc[0] +
                           4. * b * b * xloc[1] * xloc[1]);
        const Eigen::Vector3d grad_dist_norm =
            -4. * dist_norm * dist_norm * dist_norm *
            (a * a * e0 * xloc[0] + b * b * e1 * xloc[1]);
        const double F_alg =
            xloc[2] + a * xloc[0] * xloc[0] + b * xloc[1] * xloc[1];
        const Eigen::Vector3d grad_F_alg =
            e2 + 2. * (a * e0 * xloc[0] + b * e1 * xloc[1]);
        const double F = F_alg * dist_norm;
        const Eigen::Vector3d gradF =
            F_alg * grad_dist_norm + grad_F_alg * dist_norm;
        F_sum += weight * F;
        gradwxF_plus_wxgradF_sum += grad_weight * F + weight * gradF;
      }
    }
  }

  const double inv_weight_sum = 1.0 / weight_sum;
  const double PU_F = F_sum * inv_weight_sum;
  const Eigen::Vector3d PU_gradF =
      (gradwxF_plus_wxgradF_sum - F_sum * grad_weight_sum * inv_weight_sum) *
      inv_weight_sum;
  return std::make_pair(PU_F, PU_gradF);
}

std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> PU::getPUAndGradAndHessian(
    const Pt& a_pt) {
  // number of interfaces
  const UnsignedIndex_t ninterfaces = neighborhood_m->size();

  // information from neighborhood
  std::vector<SeparatorVariant> separators(ninterfaces);
  std::vector<Pt> centroids(ninterfaces);
  std::vector<double> weights(ninterfaces);
  for (UnsignedIndex_t i = 0; i < ninterfaces; ++i) {
    separators[i] = neighborhood_m->getSeparator(i);
    centroids[i] = neighborhood_m->getCentroid(i);
    weights[i] = neighborhood_m->getWeight(i);
  }

  const double inv_delta = 1. / delta_m;
  double weight_sum = 0.0;
  double F_sum = 0.0;
  Eigen::Vector3d grad_weight_sum = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hess_weight_sum = Eigen::Matrix3d::Zero();
  Eigen::Vector3d grad_product_sum = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hess_product_sum = Eigen::Matrix3d::Zero();

  // loop over interfaces
  for (UnsignedIndex_t i = 0; i < ninterfaces; ++i) {
    double vfrac_and_area_weight = weights[i];
    const Pt x = a_pt - centroids[i];
    const double r = IRL::magnitude(x);
    const double rhat = r * inv_delta;
    if (rhat < 1.0) {  // Then weight is non zero;
      const double weight = vfrac_and_area_weight * (1. + 4. * rhat) *
                            (1. - rhat) * (1. - rhat) * (1. - rhat) *
                            (1. - rhat);
      const Eigen::Vector3d x_vector = Eigen::Vector3d(x[0], x[1], x[2]);
      const Eigen::Vector3d grad_weight = x_vector * vfrac_and_area_weight *
                                          20. * (rhat - 1.) * (rhat - 1.) *
                                          (rhat - 1.) * inv_delta * inv_delta;
      const double hess_factor = vfrac_and_area_weight * 20. * (rhat - 1.) *
                                 (rhat - 1.) * inv_delta * inv_delta *
                                 inv_delta * inv_delta;
      Eigen::Matrix3d hess_weight = hess_factor * 3. * x_vector *
                                    x_vector.transpose() *
                                    (rhat > DBL_EPSILON ? 1.0 / rhat : 1.0);
      hess_weight += hess_factor * (delta_m * delta_m * (rhat - 1.)) *
                     Eigen::Matrix3d::Identity();

      // Increment weight sums
      weight_sum += weight;
      grad_weight_sum += grad_weight;
      hess_weight_sum += hess_weight;

      // computing PU value and gradient
      if (const PlanarSeparator* separator = std::get_if<IRL::PlanarSeparator>(
              &(separators[i]))) {  // If plane
        if (separator->getNumberOfPlanes() > 0) {
          const IRL::Plane plane = (*separator)[0];
          // Compute plane normal
          const IRL::Normal n = plane.normal();
          const double F = n * x;
          const Eigen::Vector3d gradF = Eigen::Vector3d(n[0], n[1], n[2]);
          F_sum += weight * F;
          grad_product_sum += grad_weight * F + weight * gradF;
          hess_product_sum += hess_weight * F +
                              grad_weight * gradF.transpose() +
                              gradF * grad_weight.transpose();
        }
      } else if (const Paraboloid* paraboloid = std::get_if<IRL::Paraboloid>(
                     &(separators[i]))) {  // If paraboloid
                                           // Get paraboloid properties
        const IRL::ReferenceFrame frame = paraboloid->getReferenceFrame();
        const double a = paraboloid->getAlignedParaboloid().a();
        const double b = paraboloid->getAlignedParaboloid().b();
        // Move sample point to local frame
        const IRL::Pt tmp = a_pt - paraboloid->getDatum();
        IRL::Pt xloc;
        for (int d = 0; d < 3; ++d) {
          xloc[d] = frame[d] * tmp;
        }
        const Eigen::Vector3d e0(frame[0][0], frame[0][1], frame[0][2]);
        const Eigen::Vector3d e1(frame[1][0], frame[1][1], frame[1][2]);
        const Eigen::Vector3d e2(frame[2][0], frame[2][1], frame[2][2]);
        // Compute Taubin distance norm term and grad and hessian
        const double dist_norm =
            1. / std::sqrt(1. + 4. * a * a * xloc[0] * xloc[0] +
                           4. * b * b * xloc[1] * xloc[1]);
        const Eigen::Vector3d grad_dist_norm =
            -4. * dist_norm * dist_norm * dist_norm *
            (a * a * e0 * xloc[0] + b * b * e1 * xloc[1]);
        const Eigen::Matrix3d hess_dist_norm =
            4. * dist_norm * dist_norm * dist_norm * dist_norm * dist_norm *
            (a * a *
                 (8. * a * a * xloc[0] * xloc[0] -
                  4. * b * b * xloc[1] * xloc[1] - 1.) *
                 e0 * e0.transpose() +
             b * b *
                 (8. * b * b * xloc[1] * xloc[1] -
                  4. * a * a * xloc[0] * xloc[0] - 1.) *
                 e1 * e1.transpose() +
             12. * a * a * b * b * xloc[0] * xloc[1] *
                 (e1 * e0.transpose() + e0 * e1.transpose()));
        // Compute algebraic distance term and grad and hessian
        const double F_alg =
            xloc[2] + a * xloc[0] * xloc[0] + b * xloc[1] * xloc[1];
        const Eigen::Vector3d grad_F_alg =
            e2 + 2. * (a * e0 * xloc[0] + b * e1 * xloc[1]);
        const Eigen::Matrix3d hess_F_alg =
            2. * (a * e0 * e0.transpose() + b * e1 * e1.transpose());
        // Compute signed distance and grad and hessian
        const double F = F_alg * dist_norm;
        const Eigen::Vector3d gradF =
            F_alg * grad_dist_norm + grad_F_alg * dist_norm;
        const Eigen::Matrix3d hessF =
            F_alg * hess_dist_norm + grad_F_alg * grad_dist_norm.transpose() +
            grad_dist_norm * grad_F_alg.transpose() + dist_norm * hess_F_alg;
        // Add to sums
        F_sum += weight * F;
        grad_product_sum += grad_weight * F + weight * gradF;
        hess_product_sum += hess_weight * F + grad_weight * gradF.transpose() +
                            gradF * grad_weight.transpose() + weight * hessF;
      }
    }
  }

  const double inv_weight_sum = 1.0 / weight_sum;
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
  return std::make_tuple(PU_F, PU_gradF, PU_hessF);
}

Pt PU::projectPointonPU(const Pt& a_pt) {
  Pt projected_pt = a_pt;
  const int itmax = 50;
  for (int i = 0; i < itmax; i++) {
    const auto F_and_gradF = getPUAndGrad(projected_pt);
    const double F = std::get<double>(F_and_gradF);
    if (F < delta_m * 1.e-6) {
      break;
    }
    // if ((i + 1) == itmax)
    //   std::cout << "projection incomplete F = " << F << std::endl;
    const Eigen::Vector3d gradF = std::get<Eigen::Vector3d>(F_and_gradF);
    const double grad_norm_inv = 1.0 / gradF.squaredNorm();
    for (int d = 0; d < 3; d++) {
      projected_pt[d] -= F * gradF(d) * grad_norm_inv;
    }
  }
  return projected_pt;
}

Paraboloid PU::solve(void) {
  // centroid of center of stencil
  const Pt centroid =
      neighborhood_m->getCentroid(neighborhood_m->getCenterOfStencil());

  // project centroid on PU surface
  const Pt projected_centroid = projectPointonPU(centroid);

  // compute gradient and hessian at projected centroid
  const auto F_gradF_hessF = getPUAndGradAndHessian(projected_centroid);
  const Eigen::Vector3d gradF = std::get<Eigen::Vector3d>(F_gradF_hessF);
  const Eigen::Matrix3d hessF = std::get<Eigen::Matrix3d>(F_gradF_hessF);
  auto new_normal = IRL::Normal(gradF(0), gradF(1), gradF(2));
  new_normal.normalize();

  auto paraboloid =
      Paraboloid::fromDerivatives(projected_centroid, gradF, hessF);

  if (IRL::magnitude(new_normal) < 0.9) {
    // setting datum to infinity to mark as invalid paraboloid
    const double inf = std::numeric_limits<double>::infinity();
    paraboloid.setDatum(IRL::Pt(inf, inf, inf));
  }

  return paraboloid;
}

}  // namespace IRL
