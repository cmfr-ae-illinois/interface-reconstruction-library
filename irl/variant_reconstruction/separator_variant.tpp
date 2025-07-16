// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_TPP_
#define IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_TPP_


namespace IRL {

inline void SeparatorVariant::serialize(ByteBuffer* a_buffer) const {
  const std::size_t index = this->index();
  a_buffer->pack(&index, 1);
  if (const auto separator = std::get_if<PlanarSeparator>(this)) {
    separator->serialize(a_buffer);
  } else if (const auto separator = std::get_if<Paraboloid>(this)) {
    separator->serialize(a_buffer);
  } else if (const auto separator = std::get_if<Cylinder>(this)) {
    separator->serialize(a_buffer);
  } else {
    throw std::runtime_error("Variant type cannot be serialized");
  }
}

inline void SeparatorVariant::unpackSerialized(ByteBuffer* a_buffer) {
  std::size_t index;
  a_buffer->unpack(&index, 1);
  if (index == 0) {
    PlanarSeparator separator;
    separator.unpackSerialized(a_buffer);
    (*this) = separator;
  } else if (index == 1) {
    Paraboloid separator;
    separator.unpackSerialized(a_buffer);
    (*this) = separator;
  } else if (index == 2) {
    Cylinder separator;
    separator.unpackSerialized(a_buffer);
    (*this) = separator;
  } else {
    throw std::runtime_error("Variant type cannot be unpacked");
  }
}

inline std::tuple<double,Eigen::Vector3d,Eigen::Matrix3d>
  SeparatorVariant::getSignedDistanceAndGradAndHessianSep(const Pt& a_pt, const Pt& a_centroid) const {
  const Pt x = a_pt - a_centroid;
  double F;
  Eigen::Vector3d gradF;
  Eigen::Matrix3d hessF;
  if(const auto sepPtr = std::get_if<PlanarSeparator>(this)) {
    // std::cout << "Variant Plane Detected\n";
    if(sepPtr->getNumberOfPlanes() > 0){
      const Plane plane = (*sepPtr)[0];
      const Normal n = plane.normal();
      F = n * x;
      gradF = Eigen::Vector3d(n[0],n[1],n[2]);
      hessF = Eigen::Matrix3d::Zero();
    }
  } else if (const auto sepPtr = std::get_if<Paraboloid>(this)) {
    // std::cout << "Variant Paraboloid Detected\n";
    const ReferenceFrame frame = sepPtr->getReferenceFrame();
    const double a = sepPtr->getAlignedParaboloid().a();
    const double b = sepPtr->getAlignedParaboloid().b();
    // Move to local frame
    const Pt tmp = a_pt - sepPtr->getDatum();
    Pt xloc;
    for(int d = 0; d < 3; ++d) {
      xloc[d] = frame[d]*tmp;
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
              (8. * a * a * xloc[0] * xloc[0] -
              4. * b * b * xloc[1] * xloc[1] - 1.) *
              e0 * e0.transpose() +
          b * b *
              (8. * b * b * xloc[1] * xloc[1] -
              4. * a * a * xloc[0] * xloc[0] - 1.) *
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
    gradF =
        F_alg * grad_dist_norm + grad_F_alg * dist_norm;
    hessF =
        F_alg * hess_dist_norm +
        grad_F_alg * grad_dist_norm.transpose() +
        grad_dist_norm * grad_F_alg.transpose() +
        dist_norm * hess_F_alg;
  } else {
    throw std::runtime_error("No signed distance for Variant Type");
  }
  return std::make_tuple(F,gradF,hessF);
}

}  // namespace IRL

#endif  // IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_TPP_
