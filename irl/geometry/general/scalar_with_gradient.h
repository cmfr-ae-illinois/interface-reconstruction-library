// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_SCALAR_WITH_GRADIENT_H_
#define IRL_GEOMETRY_GENERAL_SCALAR_WITH_GRADIENT_H_

#include <cassert>
#include <iomanip>
#include <ostream>
#include "irl/helpers/SFINAE_boiler_plate.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class FloatType, class GradientType>
class ScalarWithGradientBase {
 public:
  using value_type = FloatType;
  using gradient_type = GradientType;
  using converted_to_double =
      ScalarWithGradientBase<double,
                             typename GradientType::converted_to_double>;
  using converted_to_quad =
      ScalarWithGradientBase<Quad_t, typename GradientType::converted_to_quad>;

  ScalarWithGradientBase(void)
      : scalar_m(static_cast<FloatType>(0)),
        gradient_m(GradientType(static_cast<FloatType>(0))) {}

  constexpr ScalarWithGradientBase(const FloatType a_value)
      : scalar_m(static_cast<FloatType>(a_value)),
        gradient_m(GradientType(static_cast<FloatType>(0))) {}

  constexpr ScalarWithGradientBase(const FloatType a_value,
                                   const GradientType& a_gradient) {
    scalar_m = a_value;
    gradient_m = a_gradient;
  }
  constexpr ScalarWithGradientBase(
      const ScalarWithGradientBase<FloatType, GradientType>& a_scalar) {
    scalar_m = a_scalar.value();
    gradient_m = a_scalar.gradient();
  }
  template <class OtherFloatType, class OtherGradientType>
  constexpr ScalarWithGradientBase(
      const ScalarWithGradientBase<OtherFloatType, OtherGradientType>&
          a_scalar) {
    scalar_m = static_cast<FloatType>(a_scalar.value());
    auto& new_grad = gradient_m.getGrad();
    auto& other_grad = a_scalar.gradient().getGrad();
    assert(other_grad.rows() == other_grad.rows() &&
           new_grad.cols() == new_grad.cols());
    for (UnsignedIndex_t i = 0; i < new_grad.rows(); ++i) {
      for (UnsignedIndex_t j = 0; j < new_grad.cols(); ++j) {
        new_grad(i, j) = static_cast<FloatType>(other_grad(i, j));
      }
    }
  }

  FloatType& value(void) { return scalar_m; }
  const FloatType& value(void) const { return scalar_m; }
  GradientType& gradient(void) { return gradient_m; }
  const GradientType& gradient(void) const { return gradient_m; }
  // operator double() const { return scalar_m; }
  ScalarWithGradientBase<FloatType, GradientType>& operator+=(
      const ScalarWithGradientBase<FloatType, GradientType>& a_rhs) {
    scalar_m += a_rhs.scalar_m;
    gradient_m = gradient_m + a_rhs.gradient_m;
    return (*this);
  }
  ScalarWithGradientBase<FloatType, GradientType>& operator*=(
      const ScalarWithGradientBase<FloatType, GradientType>& a_rhs) {
    scalar_m *= a_rhs.scalar_m;
    gradient_m = gradient_m * a_rhs.scalar_m + scalar_m * a_rhs.gradient_m;
    return (*this);
  }
  ScalarWithGradientBase<FloatType, GradientType>& operator*=(
      const FloatType a_rhs) {
    scalar_m *= a_rhs;
    gradient_m = gradient_m * a_rhs;
    return (*this);
  }
  ScalarWithGradientBase<FloatType, GradientType>& operator/=(
      const ScalarWithGradientBase<FloatType, GradientType>& a_rhs) {
    assert(a_rhs.scalar_m != FloatType(0));
    scalar_m /= a_rhs.scalar_m;
    gradient_m =
        gradient_m / a_rhs.scalar_m -
        scalar_m * a_rhs.gradient_m / (a_rhs.scalar_m * a_rhs.scalar_m);
    return (*this);
  }
  ScalarWithGradientBase<FloatType, GradientType>& operator/=(
      const FloatType a_rhs) {
    scalar_m /= a_rhs;
    gradient_m = gradient_m / a_rhs;
    return (*this);
  }
  ~ScalarWithGradientBase(void) = default;

 private:
  FloatType scalar_m;
  GradientType gradient_m;
};

template <class C>
struct has_embedded_gradient : std::false_type {};

template <class C>
struct has_embedded_gradient<const C> : has_embedded_gradient<C> {};

template <class GradientType>
struct has_embedded_gradient<ScalarWithGradientBase<double, GradientType>>
    : std::true_type {};

template <class GradientType>
struct has_embedded_gradient<ScalarWithGradientBase<Quad_t, GradientType>>
    : std::true_type {};

template <bool HasGradient, class ScalarType>
struct float_type_of_scalar_container {
  using type = ScalarType;
};

template <class ScalarType>
struct float_type_of_scalar_container<true, ScalarType> {
  using type = typename ScalarType::value_type;
};

template <class ScalarType>
using float_type = typename float_type_of_scalar_container<
    has_embedded_gradient<ScalarType>::value, ScalarType>::type;

template <bool HasGradient, class ScalarType>
struct convert_scalar_to_double_container {
  using type = double;
};

template <class ScalarType>
struct convert_scalar_to_double_container<true, ScalarType> {
  using type = typename ScalarType::converted_to_double;
};

template <class ScalarType>
using convert_to_double = typename convert_scalar_to_double_container<
    has_embedded_gradient<ScalarType>::value, ScalarType>::type;

template <bool HasGradient, class ScalarType>
struct convert_scalar_to_quad_container {
  using type = Quad_t;
};

template <class ScalarType>
struct convert_scalar_to_quad_container<true, ScalarType> {
  using type = typename ScalarType::converted_to_quad;
};

template <class ScalarType>
using convert_to_quad = typename convert_scalar_to_quad_container<
    has_embedded_gradient<ScalarType>::value, ScalarType>::type;

template <class FloatType, class GradientType>
inline bool operator==(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  return a_scalar1.value() == a_scalar2.value();
}

template <class FloatType, class GradientType>
inline bool operator!=(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  return a_scalar1.value() != a_scalar2.value();
}

template <class FloatType, class GradientType>
inline bool operator>(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  return a_scalar1.value() > a_scalar2.value();
}

template <class FloatType, class GradientType>
inline bool operator>=(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  return a_scalar1.value() >= a_scalar2.value();
}

template <class FloatType, class GradientType>
inline bool operator<(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  return a_scalar1.value() < a_scalar2.value();
}

template <class FloatType, class GradientType>
inline bool operator<=(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  return a_scalar1.value() <= a_scalar2.value();
}

template <class FloatType, class GradientType>
ScalarWithGradientBase<FloatType, GradientType> operator*(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar,
    const FloatType a_rhs) {
  ScalarWithGradientBase<FloatType, GradientType> new_scalar(a_scalar);
  new_scalar.value() *= a_rhs;
  new_scalar.gradient() = new_scalar.gradient() * a_rhs;
  return new_scalar;
}

template <class FloatType, class GradientType>
ScalarWithGradientBase<FloatType, GradientType> operator*(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  ScalarWithGradientBase<FloatType, GradientType> new_scalar(
      static_cast<FloatType>(0));
  new_scalar.value() = a_scalar1.value() * a_scalar2.value();
  if constexpr (GradientType::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar1.value() * a_scalar2.gradient().getGrad() +
        a_scalar1.gradient().getGrad() * a_scalar2.value();
    new_scalar.gradient().getHessian() =
        a_scalar1.value() * a_scalar2.gradient().getHessian() +
        a_scalar1.gradient().getHessian() * a_scalar2.value() +
        a_scalar1.gradient().getGrad() *
            a_scalar2.gradient().getGrad().transpose() +
        a_scalar2.gradient().getGrad() *
            a_scalar1.gradient().getGrad().transpose();
  } else {
    new_scalar.gradient() = a_scalar1.value() * a_scalar2.gradient() +
                            a_scalar1.gradient() * a_scalar2.value();
  }
  return new_scalar;
}

template <class FloatType, class GradientType>
ScalarWithGradientBase<FloatType, GradientType> operator*(
    const FloatType a_rhs,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar) {
  return ScalarWithGradientBase<FloatType, GradientType>(a_scalar * a_rhs);
}

template <class FloatType, class GradientType>
ScalarWithGradientBase<FloatType, GradientType> operator/(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar,
    const FloatType a_rhs) {
  ScalarWithGradientBase<FloatType, GradientType> new_scalar(a_scalar);
  new_scalar.value() /= a_rhs;
  new_scalar.gradient() = new_scalar.gradient() / a_rhs;
  return new_scalar;
}

template <class FloatType, class GradientType>
ScalarWithGradientBase<FloatType, GradientType> operator/(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  ScalarWithGradientBase<FloatType, GradientType> new_scalar(
      static_cast<FloatType>(0));
  new_scalar.value() = a_scalar1.value() / a_scalar2.value();
  if constexpr (GradientType::has_hessian) {
    ScalarWithGradientBase<FloatType, GradientType> inv_scalar2(
        static_cast<FloatType>(0));
    inv_scalar2.value() = static_cast<FloatType>(1) / a_scalar2.value();
    inv_scalar2.gradient().getGrad() = -a_scalar2.gradient().getGrad() /
                                       (a_scalar2.value() * a_scalar2.value());
    inv_scalar2.gradient().getHessian() =
        -a_scalar2.gradient().getHessian() /
            (a_scalar2.value() * a_scalar2.value()) +
        static_cast<FloatType>(2) * a_scalar2.gradient().getGrad() *
            a_scalar2.gradient().getGrad().transpose() /
            (a_scalar2.value() * a_scalar2.value() * a_scalar2.value());
    new_scalar.gradient().getGrad() =
        a_scalar1.value() * inv_scalar2.gradient().getGrad() +
        a_scalar1.gradient().getGrad() * inv_scalar2.value();
    new_scalar.gradient().getHessian() =
        a_scalar1.value() * inv_scalar2.gradient().getHessian() +
        a_scalar1.gradient().getHessian() * inv_scalar2.value() +
        a_scalar1.gradient().getGrad() *
            inv_scalar2.gradient().getGrad().transpose() +
        inv_scalar2.gradient().getGrad() *
            a_scalar1.gradient().getGrad().transpose();
  } else {
    new_scalar.gradient() = a_scalar1.gradient() / a_scalar2.value() -
                            a_scalar1.value() * a_scalar2.gradient() /
                                (a_scalar2.value() * a_scalar2.value());
  }
  return new_scalar;
}

template <class FloatType, class GradientType>
ScalarWithGradientBase<FloatType, GradientType> operator+(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  ScalarWithGradientBase<FloatType, GradientType> new_scalar(
      static_cast<FloatType>(0));
  new_scalar.value() = a_scalar1.value() + a_scalar2.value();
  new_scalar.gradient() = a_scalar1.gradient() + a_scalar2.gradient();
  return new_scalar;
}

template <class FloatType, class GradientType>
ScalarWithGradientBase<FloatType, GradientType> operator-(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar2) {
  ScalarWithGradientBase<FloatType, GradientType> new_scalar(
      static_cast<FloatType>(0));
  new_scalar.value() = a_scalar1.value() - a_scalar2.value();
  new_scalar.gradient() = a_scalar1.gradient() - a_scalar2.gradient();
  return new_scalar;
}

template <class FloatType, class GradientType>
ScalarWithGradientBase<FloatType, GradientType> operator-(
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar) {
  ScalarWithGradientBase<FloatType, GradientType> new_scalar(
      static_cast<FloatType>(0));
  new_scalar.value() = -a_scalar.value();
  new_scalar.gradient() = -a_scalar.gradient();
  return new_scalar;
}

// /************************************ SQRT
// ************************************/

// template <class ScalarWithGradientType>
// enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
//             ScalarWithGradientType>
// SqrtMoments(const ScalarWithGradientType& a_scalar) {
//   return sqrt(a_scalar);
// }

// // template <>
// // inline enable_if_t<!has_embedded_gradient<double>::value, double>
// // SqrtMoments(
// //     const double& a_scalar) {
// //   return std::sqrt(a_scalar);
// // }

// // template <>
// // inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t>
// // SqrtMoments(
// //     const double& a_scalar) {
// //   return sqrtq(a_scalar);
// // }

// template <class ScalarWithGradientType>
// inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
//                    ScalarWithGradientType>
// SqrtMoments(const ScalarWithGradientType& a_scalar) {
//   using FloatType = typename ScalarWithGradientType::value_type;
//   ScalarWithGradientType new_scalar(static_cast<FloatType>(0));
//   new_scalar.value() = sqrt(a_scalar.value());
//   if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
//     new_scalar.gradient().getGrad() =
//         a_scalar.gradient().getGrad() /
//         (static_cast<FloatType>(2) * new_scalar.value());
//     new_scalar.gradient().getHessian() =
//         a_scalar.gradient().getHessian() /
//             (static_cast<FloatType>(2) * new_scalar.value()) -
//         a_scalar.gradient().getGrad() *
//             new_scalar.gradient().getGrad().transpose() /
//             (static_cast<FloatType>(2) * new_scalar.value() *
//              new_scalar.value());
//   } else {
//     new_scalar.gradient() =
//         a_scalar.gradient() / (static_cast<FloatType>(2) *
//         new_scalar.value());
//   }
//   return new_scalar;
// }

// /************************************ POW
// ************************************/

// template <class ScalarWithGradientType>
// enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
//             ScalarWithGradientType>
// PowMoments(const ScalarWithGradientType& a_scalar,
//            const typename ScalarWithGradientType::value_type a_power) {
//   return pow(a_scalar, a_power);
// }

// // template <>
// // inline enable_if_t<!has_embedded_gradient<double>::value, double>
// PowMoments(
// //     const double& a_scalar, const double a_power) {
// //   return std::pow(a_scalar, a_power);
// // }
// // template <>
// // inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t>
// PowMoments(
// //     const Quad_t& a_scalar, const Quad_t a_power) {
// //   return powq(a_scalar, a_power);
// // }

// template <class ScalarWithGradientType>
// inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
//                    ScalarWithGradientType>
// PowMoments(const ScalarWithGradientType& a_scalar,
//            const typename ScalarWithGradientType::value_type a_power) {
//   using FloatType = typename ScalarWithGradientType::value_type;
//   ScalarWithGradientType new_scalar(static_cast<FloatType>(0));
//   new_scalar.value() = pow(a_scalar.value(), a_power);
//   if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
//     new_scalar.gradient().getGrad() =
//         a_power * a_scalar.gradient().getGrad() *
//         pow(a_scalar.value(), a_power - static_cast<FloatType>(1));
//     new_scalar.gradient().getHessian() =
//         a_power * a_scalar.gradient().getHessian() *
//             pow(a_scalar.value(), a_power - static_cast<FloatType>(1)) +
//         a_power * a_scalar.gradient().getGrad() *
//             (a_power - static_cast<FloatType>(1)) *
//             a_scalar.gradient().getGrad().transpose() *
//             pow(a_scalar.value(), a_power - static_cast<FloatType>(2));
//   } else {
//     new_scalar.gradient() =
//         a_power * a_scalar.gradient() *
//         pow(a_scalar.value(), a_power - static_cast<FloatType>(1));
//   }
//   return new_scalar;
// }

// /************************************ LOG
// ************************************/

// template <class ScalarWithGradientType>
// enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
//             ScalarWithGradientType>
// LogMoments(const ScalarWithGradientType& a_scalar) {
//   return log(a_scalar);
// }

// // template <>
// // inline enable_if_t<!has_embedded_gradient<double>::value, double>
// LogMoments(
// //     const double& a_scalar) {
// //   return std::log(a_scalar);
// // }
// // template <>
// // inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t>
// LogMoments(
// //     const Quad_t& a_scalar) {
// //   return logq(a_scalar);
// // }

// template <class ScalarWithGradientType>
// inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
//                    ScalarWithGradientType>
// LogMoments(const ScalarWithGradientType& a_scalar) {
//   using FloatType = typename ScalarWithGradientType::value_type;
//   ScalarWithGradientType new_scalar(static_cast<FloatType>(0));
//   new_scalar.value() = log(a_scalar.value());
//   if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
//     new_scalar.gradient().getGrad() =
//         a_scalar.gradient().getGrad() / a_scalar.value();
//     new_scalar.gradient().getHessian() =
//         a_scalar.gradient().getHessian() / a_scalar.value() -
//         a_scalar.gradient().getGrad() *
//             a_scalar.gradient().getGrad().transpose() /
//             (a_scalar.value() * a_scalar.value());
//   } else {
//     new_scalar.gradient() = a_scalar.gradient() / a_scalar.value();
//   }
//   return new_scalar;
// }

// /************************************ ATAN
// ************************************/

// template <class ScalarWithGradientType>
// enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
//             ScalarWithGradientType>
// ArctanMoments(const ScalarWithGradientType& a_scalar) {
//   return atan(a_scalar);
// }

// // template <>
// // inline enable_if_t<!has_embedded_gradient<double>::value, double>
// // ArctanMoments(
// //     const double& a_scalar) {
// //   return std::atan(a_scalar);
// // }
// // template <>
// // inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t>
// // ArctanMoments(
// //     const Quad_t& a_scalar) {
// //   return atanq(a_scalar);
// // }

// template <class ScalarWithGradientType>
// inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
//                    ScalarWithGradientType>
// ArctanMoments(const ScalarWithGradientType& a_scalar) {
//   using FloatType = typename ScalarWithGradientType::value_type;
//   ScalarWithGradientType new_scalar(static_cast<FloatType>(0));
//   new_scalar.value() = atan(a_scalar.value());
//   if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
//     new_scalar.gradient().getGrad() =
//         a_scalar.gradient().getGrad() /
//         (static_cast<FloatType>(1) + a_scalar.value() * a_scalar.value());
//     new_scalar.gradient().getHessian() =
//         a_scalar.gradient().getHessian() /
//             (static_cast<FloatType>(1) + a_scalar.value() * a_scalar.value())
//             -
//         a_scalar.gradient().getGrad() *
//             (static_cast<FloatType>(2) *
//              a_scalar.gradient().getGrad().transpose() * a_scalar.value()) /
//             ((static_cast<FloatType>(1) + a_scalar.value() *
//             a_scalar.value()) *
//              (static_cast<FloatType>(1) + a_scalar.value() *
//              a_scalar.value()));
//   } else {
//     new_scalar.gradient() =
//         a_scalar.gradient() /
//         (static_cast<FloatType>(1) + a_scalar.value() * a_scalar.value());
//   }
//   return new_scalar;
// }

// /*********************************** ATANH
// ***********************************/

// template <class ScalarWithGradientType>
// enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
//             ScalarWithGradientType>
// ArctanhMoments(const ScalarWithGradientType& a_scalar) {
//   return atanh(a_scalar);
// }

// // template <>
// // inline enable_if_t<!has_embedded_gradient<double>::value, double>
// // ArctanhMoments(const double& a_scalar) {
// //   return std::atanh(a_scalar);
// // }
// // template <>
// // inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t>
// // ArctanhMoments(const Quad_t& a_scalar) {
// //   return atanhq(a_scalar);
// // }

// template <class ScalarWithGradientType>
// inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
//                    ScalarWithGradientType>
// ArctanhMoments(const ScalarWithGradientType& a_scalar) {
//   using FloatType = typename ScalarWithGradientType::value_type;
//   ScalarWithGradientType new_scalar(static_cast<FloatType>(0));
//   new_scalar.value() = atanh(a_scalar.value());
//   if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
//     new_scalar.gradient().getGrad() =
//         a_scalar.gradient().getGrad() /
//         (static_cast<FloatType>(1) - a_scalar.value() * a_scalar.value());
//     new_scalar.gradient().getHessian() =
//         a_scalar.gradient().getHessian() /
//             (static_cast<FloatType>(1) - a_scalar.value() * a_scalar.value())
//             +
//         a_scalar.gradient().getGrad() *
//             (static_cast<FloatType>(2) *
//              a_scalar.gradient().getGrad().transpose() * a_scalar.value()) /
//             ((static_cast<FloatType>(1) - a_scalar.value() *
//             a_scalar.value()) *
//              (static_cast<FloatType>(1) - a_scalar.value() *
//              a_scalar.value()));
//   } else {
//     new_scalar.gradient() =
//         a_scalar.gradient() /
//         (static_cast<FloatType>(1) - a_scalar.value() * a_scalar.value());
//   }
//   return new_scalar;
// }

template <class FloatType, class GradientType>
inline std::ostream& operator<<(
    std::ostream& out,
    const ScalarWithGradientBase<FloatType, GradientType>& a_scalar) {
  out << std::scientific << std::setprecision(6) << std::showpos
      << a_scalar.value() << std::setprecision(3)
      << " (grad = " << a_scalar.gradient() << ")";
  return out;
}

template <class GradientType>
using ScalarWithGradient = ScalarWithGradientBase<double, GradientType>;

}  // namespace IRL

#include "irl/geometry/general/scalar_with_gradient.tpp"

#endif