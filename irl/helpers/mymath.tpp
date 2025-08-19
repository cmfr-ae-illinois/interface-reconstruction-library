// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_HELPERS_MYMATH_TPP_
#define IRL_HELPERS_MYMATH_TPP_

namespace IRL {

inline constexpr double deg2Rad(const double a_degree) {
  return a_degree * M_PI / 180.0;
}

inline constexpr double rad2Deg(const double a_radian) {
  return a_radian * 180.0 / M_PI;
}

inline double angleNormalize(const double a_radian) {
  return std::max(0.0,
                  a_radian - (std::floor(0.5 * a_radian / M_PI)) * 2.0 * M_PI);
}
inline double signedAngleNormalize(const double a_radian) {
  return a_radian < 0.0 ? -(2.0 * M_PI - angleNormalize(a_radian))
                        : angleNormalize(a_radian);
}

template <class DataType>
inline typename DataType::value_type magnitude(const DataType& a_vector) {
  return sqrt(a_vector[0] * a_vector[0] + a_vector[1] * a_vector[1] +
              a_vector[2] * a_vector[2]);
}

template <class DataType>
inline typename DataType::value_type squaredMagnitude(
    const DataType& a_vector) {
  return a_vector[0] * a_vector[0] + a_vector[1] * a_vector[1] +
         a_vector[2] * a_vector[2];
}

// template <class DataType>
// inline DataType crossProduct(const DataType& a_vector_0,
//                              const DataType& a_vector_1) {
//   return DataType(
//       a_vector_0[1] * a_vector_1[2] - a_vector_0[2] * a_vector_1[1],
//       a_vector_0[2] * a_vector_1[0] - a_vector_0[0] * a_vector_1[2],
//       a_vector_0[0] * a_vector_1[1] - a_vector_0[1] * a_vector_1[0]);
// }

template <class DataType>
inline DataType crossProductNormalized(const DataType& a_vector_0,
                                       const DataType& a_vector_1) {
  DataType object_to_return = crossProduct(a_vector_0, a_vector_1);
  object_to_return /= safelyTiny(magnitude(object_to_return));
  return object_to_return;
}

template <class T1, class T2>
inline typename T1::value_type dotProduct(const T1& a_vector_0,
                                          const T2& a_vector_1) {
  static_assert(
      std::is_same<typename T1::value_type, typename T2::value_type>::value,
      "Trying to dot two vectors of incompatible types (e.g., double and "
      "quad).");
  return a_vector_0[0] * a_vector_1[0] + a_vector_0[1] * a_vector_1[1] +
         a_vector_0[2] * a_vector_1[2];
}

template <class DataType>
inline typename DataType::value_type scalarTripleProduct(
    const DataType& a_vector_0, const DataType& a_vector_1,
    const DataType& a_vector_2) {
  return dotProduct(a_vector_0, crossProduct(a_vector_1, a_vector_2));
}

template <>
inline bool isnan(double a_scalar) {
  return std::isnan(a_scalar);
}

template <>
inline bool isnan(Quad_t a_scalar) {
  return isnanq(a_scalar) == 1;
}

template <class ScalarType>
inline enable_if_t<has_embedded_gradient<ScalarType>::value, bool> isnan(
    ScalarType a_scalar) {
  using FloatType = float_type<ScalarType>;
  return isnan<FloatType>(a_scalar.value);
}

template <class ScalarType>
inline ScalarType machine_epsilon(void) {
  using FloatType = float_type<ScalarType>;
  if constexpr (std::is_same_v<FloatType, double>) {
    return ScalarType(DBL_EPSILON);
  } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
    return ScalarType(1.0e5q * FLT128_EPSILON);
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain machine epsilon for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType machine_pi(void) {
  using FloatType = float_type<ScalarType>;
  if constexpr (std::is_same_v<FloatType, double>) {
    return ScalarType(M_PI);
  } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
    return ScalarType(M_PIq);
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain machine pi for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType abs(const ScalarType a_scalar) {
  using FloatType = float_type<ScalarType>;
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::fabs(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return fabsq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    if (a_scalar.value() < FloatType(0)) {
      return -a_scalar;
    } else {
      return a_scalar;
    }
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain absolute value for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType fabs(const ScalarType a_scalar) {
  using FloatType = float_type<ScalarType>;
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::fabs(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return fabsq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    if (a_scalar.value() < FloatType(0)) {
      return -a_scalar;
    } else {
      return a_scalar;
    }
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain absolute value for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType sqrt(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::sqrt(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return sqrtq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::sqrt(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = sqrtq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      new_scalar.gradient().getGrad() =
          a_scalar.gradient() /
          (FloatType(2) * safelyEpsilon(new_scalar.value()));
      new_scalar.gradient().getHessian() =
          a_scalar.gradient().getHessian() /
              (FloatType(2) * new_scalar.value()) -
          a_scalar.gradient().getGrad() *
              new_scalar.gradient().getGrad().transpose() /
              (FloatType(2) * new_scalar.value() * new_scalar.value());
    } else {
      new_scalar.gradient() =
          a_scalar.gradient() /
          (FloatType(2) * safelyEpsilon(new_scalar.value()));
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain sqrt for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType approxsqrt(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::sqrt(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return sqrtq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::sqrt(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = sqrtq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      new_scalar.gradient().getGrad() =
          a_scalar.gradient().getGrad() / (FloatType(2) * new_scalar.value());
      new_scalar.gradient().getHessian() =
          a_scalar.gradient().getHessian() /
              (FloatType(2) * new_scalar.value()) -
          a_scalar.gradient().getGrad() *
              new_scalar.gradient().getGrad().transpose() /
              (FloatType(2) * new_scalar.value() * new_scalar.value());
    } else {
      new_scalar.gradient() =
          a_scalar.gradient() /
          (FloatType(2) * safelyEpsilon(sqrt(a_scalar.value())));
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain approxsqrt for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType invsqrt(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return 1.0 / std::sqrt(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return 1.0q / sqrtq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = 1.0 / std::sqrt(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = 1.0q / sqrtq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      new_scalar.gradient() =
          -a_scalar.gradient() /
          (FloatType(2) * safelyEpsilon(pow(a_scalar.value(), FloatType(1.5))));
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain approxsqrt for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType approxinvsqrt(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return 1.0 / std::sqrt(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return 1.0q / sqrtq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = 1.0 / std::sqrt(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = 1.0q / sqrtq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      new_scalar.gradient() =
          -a_scalar.gradient() /
          (FloatType(2) * safelyEpsilon(pow(a_scalar.value(), FloatType(1.5))));
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain approxsqrt for unknown float type");
  }
}

template <class ScalarType, class PowerType>
inline ScalarType pow(const ScalarType a_scalar, const PowerType a_power) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::pow(a_scalar, static_cast<double>(a_power));
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return powq(a_scalar, static_cast<Quad_t>(a_power));
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    FloatType new_power = FloatType(0);
    if constexpr (has_embedded_gradient<PowerType>::value) {
      new_power = static_cast<FloatType>(a_power.value());
    } else {
      new_power = static_cast<FloatType>(a_power);
    }
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::pow(a_scalar.value(), new_power);
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = powq(a_scalar.value(), new_power);
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      new_scalar.gradient() = a_scalar.gradient() * new_power *
                              pow(a_scalar.value(), new_power - FloatType(1));
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain pow for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType log(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::log(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return logq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::log(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = logq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      new_scalar.gradient() =
          a_scalar.gradient() / safelyEpsilon(a_scalar.value());
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain log for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType atan2(const ScalarType a_scalar_y,
                        const ScalarType a_scalar_x) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::atan2(a_scalar_y, a_scalar_x);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return atan2q(a_scalar_y, a_scalar_x);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::atan2(a_scalar_y.value(), a_scalar_x.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = atan2q(a_scalar_y.value(), a_scalar_x.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      // new_scalar.gradient() =
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain atan2 for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType atan(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::atan(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return atanq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::atan(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = atanq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      new_scalar.gradient() =
          a_scalar.gradient() /
          (FloatType(1) + a_scalar.value() * a_scalar.value());
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain atan for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType atanh(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::atanh(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return atanhq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::atanh(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = atanhq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      new_scalar.gradient() =
          a_scalar.gradient() /
          safelyEpsilon(FloatType(1) - a_scalar.value() * a_scalar.value());
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain atanh for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType cos(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::cos(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return cosq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::cos(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = cosq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      new_scalar.gradient() = -a_scalar.gradient() * sin(a_scalar.value());
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain cos for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType sin(const ScalarType a_scalar) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::sin(a_scalar);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    return sinq(a_scalar);
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    ScalarType new_scalar(FloatType(0));
    if constexpr (std::is_same_v<FloatType, double>) {
      new_scalar.value() = std::sin(a_scalar.value());
    } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
      new_scalar.value() = sinq(a_scalar.value());
    }
    if constexpr (ScalarType::gradient_type::has_hessian) {
      // new_scalar.gradient().getGrad() =
      // new_scalar.gradient().getHessian() =
    } else {
      new_scalar.gradient() = a_scalar.gradient() * cos(a_scalar.value());
    }
    return new_scalar;
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain sin for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType copysign(const ScalarType a_scalar, const ScalarType a_sign) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::copysign(a_scalar, a_sign);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    if (a_sign >= 0.0q) {
      return fabs(a_scalar);
    } else {
      return -fabs(a_scalar);
    }
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    using FloatType = float_type<ScalarType>;
    if (a_sign.value() >= FloatType(0)) {
      return fabs(a_scalar);
    } else {
      return -fabs(a_scalar);
    }
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain copysign for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType minimum(const ScalarType a_scalar1,
                          const ScalarType a_scalar2) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::min(a_scalar1, a_scalar2);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    if (a_scalar1 < a_scalar2) {
      return a_scalar1;
    } else {
      return a_scalar2;
    }
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    if (a_scalar1.value() < a_scalar2.value()) {
      return a_scalar1;
    } else {
      return a_scalar2;
    }
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain minimum for unknown float type");
  }
}

template <class ScalarType>
inline ScalarType maximum(const ScalarType a_scalar1,
                          const ScalarType a_scalar2) {
  if constexpr (std::is_same_v<ScalarType, double>) {
    return std::max(a_scalar1, a_scalar2);
  } else if constexpr (std::is_same_v<ScalarType, Quad_t>) {
    if (a_scalar1 >= a_scalar2) {
      return a_scalar1;
    } else {
      return a_scalar2;
    }
  } else if constexpr (has_embedded_gradient<ScalarType>::value) {
    if (a_scalar1.value() >= a_scalar2.value()) {
      return a_scalar1;
    } else {
      return a_scalar2;
    }
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain maximum for unknown float type");
  }
}

}  // namespace IRL

#endif  // IRL_HELPERS_MYMATH_TPP_
