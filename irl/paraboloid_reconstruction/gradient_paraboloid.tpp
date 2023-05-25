// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_TPP_
#define IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_TPP_

namespace IRL {

inline ParaboloidGradientLocalZ::ParaboloidGradientLocalZ(void) {
  gradient_m = Eigen::Matrix<double, 1, 1>::Zero();
}

inline ParaboloidGradientLocalZ::ParaboloidGradientLocalZ(
    const double a_value) {
  gradient_m = Eigen::Matrix<double, 1, 1>::Constant(a_value);
}

inline ParaboloidGradientLocalZ::ParaboloidGradientLocalZ(
    const ParaboloidGradientLocalZ& a_gradient) {
  gradient_m = a_gradient.getGrad();
}

inline ParaboloidGradientLocalZ& ParaboloidGradientLocalZ::operator+=(
    const ParaboloidGradientLocalZ& a_gradient) {
  this->getGrad() += a_gradient.getGrad();
  return (*this);
}

inline ParaboloidGradientLocalZ& ParaboloidGradientLocalZ::operator=(
    const ParaboloidGradientLocalZ& a_gradient) {
  this->getGrad() = a_gradient.getGrad();
  return (*this);
}

inline ParaboloidGradientLocalZ operator*(
    const ParaboloidGradientLocalZ& a_gradient, const double a_rhs) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = a_gradient.getGrad() * a_rhs;
  return gradient;
}

inline ParaboloidGradientLocalZ operator*(
    const double a_rhs, const ParaboloidGradientLocalZ& a_gradient) {
  return a_gradient * a_rhs;
}

inline ParaboloidGradientLocalZ operator/(
    const ParaboloidGradientLocalZ& a_gradient, const double a_rhs) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = a_gradient.getGrad() / a_rhs;
  return gradient;
}

inline ParaboloidGradientLocalZ operator+(
    const ParaboloidGradientLocalZ& a_gradient1,
    const ParaboloidGradientLocalZ& a_gradient2) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = a_gradient1.getGrad() + a_gradient2.getGrad();
  return gradient;
}

inline ParaboloidGradientLocalZ operator-(
    const ParaboloidGradientLocalZ& a_gradient1,
    const ParaboloidGradientLocalZ& a_gradient2) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = a_gradient1.getGrad() - a_gradient2.getGrad();
  return gradient;
}
inline ParaboloidGradientLocalZ operator-(
    const ParaboloidGradientLocalZ& a_gradient) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = -a_gradient.getGrad();
  return gradient;
}

inline Eigen::Matrix<double, 1, 1>& ParaboloidGradientLocalZ::getGrad(void) {
  return gradient_m;
}

inline const Eigen::Matrix<double, 1, 1>& ParaboloidGradientLocalZ::getGrad(
    void) const {
  return gradient_m;
}

inline double ParaboloidGradientLocalZ::getGradA(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradB(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradTx(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradTy(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradTz(void) const {
  return gradient_m(0);
}
inline double ParaboloidGradientLocalZ::getGradRx(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradRy(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradRz(void) const { return 0.0; }
inline void ParaboloidGradientLocalZ::setGrad(
    const ParaboloidGradientLocalZ& a_rhs) {
  gradient_m = a_rhs.getGrad();
}
inline void ParaboloidGradientLocalZ::setGradA(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradB(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTx(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTy(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTz(const double a_rhs) {
  gradient_m(0) = a_rhs;
}
inline void ParaboloidGradientLocalZ::setGradRx(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradRy(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradRz(const double a_rhs) {}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType>::ParaboloidGradientLocalBase(
    void) {
  gradient_m.setZero();
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType>::ParaboloidGradientLocalBase(
    const ScalarType a_value) {
  gradient_m.setConstant(a_value);
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType>::ParaboloidGradientLocalBase(
    const ParaboloidGradientLocalBase<double>& a_gradient) {
  for (Eigen::Index i{0}; i < gradient_m.rows(); ++i) {
    gradient_m(i) = static_cast<ScalarType>(a_gradient.getGrad()(i));
  }
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType>::ParaboloidGradientLocalBase(
    const ParaboloidGradientLocalBase<Quad_t>& a_gradient) {
  for (Eigen::Index i{0}; i < gradient_m.rows(); ++i) {
    gradient_m(i) = static_cast<ScalarType>(a_gradient.getGrad()(i));
  }
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType>&
ParaboloidGradientLocalBase<ScalarType>::operator+=(
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient) {
  for (Eigen::Index i{0}; i < gradient_m.rows(); ++i) {
    this->getGrad()(i) += static_cast<ScalarType>(a_gradient.getGrad()(i));
  }
  return (*this);
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType>&
ParaboloidGradientLocalBase<ScalarType>::operator=(
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient) {
  for (Eigen::Index i{0}; i < gradient_m.rows(); ++i) {
    this->getGrad()(i) = static_cast<ScalarType>(a_gradient.getGrad()(i));
  }
  return (*this);
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType> operator*(
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient,
    const ScalarType a_rhs) {
  ParaboloidGradientLocalBase<ScalarType> gradient(static_cast<ScalarType>(0));
  gradient.getGrad() = a_gradient.getGrad() * a_rhs;
  return gradient;
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType> operator*(
    const ScalarType a_rhs,
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient) {
  return a_gradient * a_rhs;
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType> operator/(
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient,
    const ScalarType a_rhs) {
  ParaboloidGradientLocalBase<ScalarType> gradient(static_cast<ScalarType>(0));
  gradient.getGrad() = a_gradient.getGrad() / a_rhs;
  return gradient;
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType> operator+(
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient1,
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient2) {
  ParaboloidGradientLocalBase<ScalarType> gradient(static_cast<ScalarType>(0));
  gradient.getGrad() = a_gradient1.getGrad() + a_gradient2.getGrad();
  return gradient;
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType> operator-(
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient1,
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient2) {
  ParaboloidGradientLocalBase<ScalarType> gradient(static_cast<ScalarType>(0));
  gradient.getGrad() = a_gradient1.getGrad() - a_gradient2.getGrad();
  return gradient;
}

template <class ScalarType>
inline ParaboloidGradientLocalBase<ScalarType> operator-(
    const ParaboloidGradientLocalBase<ScalarType>& a_gradient) {
  ParaboloidGradientLocalBase<ScalarType> gradient(static_cast<ScalarType>(0));
  gradient.getGrad() = -a_gradient.getGrad();
  return gradient;
}

template <class ScalarType>
inline Eigen::Matrix<ScalarType, 8, 1>&
ParaboloidGradientLocalBase<ScalarType>::getGrad(void) {
  return gradient_m;
}

template <class ScalarType>
inline const Eigen::Matrix<ScalarType, 8, 1>&
ParaboloidGradientLocalBase<ScalarType>::getGrad(void) const {
  return gradient_m;
}

template <class ScalarType>
inline ScalarType ParaboloidGradientLocalBase<ScalarType>::getGradA(
    void) const {
  return gradient_m(0);
}
template <class ScalarType>
inline ScalarType ParaboloidGradientLocalBase<ScalarType>::getGradB(
    void) const {
  return gradient_m(1);
}
template <class ScalarType>
inline ScalarType ParaboloidGradientLocalBase<ScalarType>::getGradTx(
    void) const {
  return gradient_m(2);
}
template <class ScalarType>
inline ScalarType ParaboloidGradientLocalBase<ScalarType>::getGradTy(
    void) const {
  return gradient_m(3);
}
template <class ScalarType>
inline ScalarType ParaboloidGradientLocalBase<ScalarType>::getGradTz(
    void) const {
  return gradient_m(4);
}
template <class ScalarType>
inline ScalarType ParaboloidGradientLocalBase<ScalarType>::getGradRx(
    void) const {
  return gradient_m(5);
}
template <class ScalarType>
inline ScalarType ParaboloidGradientLocalBase<ScalarType>::getGradRy(
    void) const {
  return gradient_m(6);
}
template <class ScalarType>
inline ScalarType ParaboloidGradientLocalBase<ScalarType>::getGradRz(
    void) const {
  return gradient_m(7);
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGrad(
    const ParaboloidGradientLocalBase<ScalarType>& a_rhs) {
  gradient_m = a_rhs.getGrad();
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGradA(
    const ScalarType a_rhs) {
  gradient_m(0) = a_rhs;
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGradB(
    const ScalarType a_rhs) {
  gradient_m(1) = a_rhs;
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGradTx(
    const ScalarType a_rhs) {
  gradient_m(2) = a_rhs;
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGradTy(
    const ScalarType a_rhs) {
  gradient_m(3) = a_rhs;
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGradTz(
    const ScalarType a_rhs) {
  gradient_m(4) = a_rhs;
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGradRx(
    const ScalarType a_rhs) {
  gradient_m(5) = a_rhs;
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGradRy(
    const ScalarType a_rhs) {
  gradient_m(6) = a_rhs;
}
template <class ScalarType>
inline void ParaboloidGradientLocalBase<ScalarType>::setGradRz(
    const ScalarType a_rhs) {
  gradient_m(7) = a_rhs;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_TPP_
