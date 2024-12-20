#ifndef IRL_QUADRATIC_RECONSTRUCTION_QUADRATIC_HELPER_TPP_
#define IRL_QUADRATIC_RECONSTRUCTION_QUADRATIC_HELPER_TPP_

namespace IRL {

  // solveQuadraticBetween0And1 Doesn't seems to be used 
  // TODO : Remove?

  // // Returns solution to quadratic equation solve.
  // // The smallest solution will always be first.
  // template <class ScalarType>
  // inline StackVector<ScalarType, 2> solveQuadraticBetween0And1(
  //     const ScalarType a, const ScalarType b, const ScalarType c) {
  //   ScalarType discriminant = b * b - static_cast<ScalarType>(4) * a * c;

  //   if (discriminant > static_cast<ScalarType>(0)) {
  //       if (a != static_cast<ScalarType>(0)) {
  //         /* First fast try in 32-bit precision */
  //         const ScalarType approx_discriminant = approxsqrt(discriminant);
  //         const ScalarType approx_q = -static_cast<ScalarType>(0.5) *
  //                                     (b + copysign(approx_discriminant, b));
  //         const ScalarType approx_sol1 = approx_q / safelyTiny(a);
  //         const ScalarType approx_sol2 = c / safelyTiny(approx_q);
  //         if ((approx_sol1 < -0.01 || approx_sol1 > 1.01) &&
  //             (approx_sol2 < -0.01 || approx_sol2 > 1.01)) {
  //             return StackVector<ScalarType, 2>();
  //         }

  //         /* Real calculation */
  //         discriminant = sqrt(discriminant);
  //         const ScalarType q =
  //             -static_cast<ScalarType>(0.5) * (b + copysign(discriminant, b));
  //         const ScalarType sol1 = q / safelyTiny(a);
  //         const ScalarType sol2 = c / safelyTiny(q);
  //         return sol1 < sol2 ? StackVector<ScalarType, 2>({sol1, sol2})
  //                             : StackVector<ScalarType, 2>({sol2, sol1});
  //       } else {
  //         return StackVector<ScalarType, 2>({-c / b});
  //       }
  //   }
  //   return StackVector<ScalarType, 2>();
  // };

  // Returns solution to quadratic equation solve.
  // The smallest solution will always be first.
  template <class ScalarType>
  inline StackVector<ScalarType, 2> solveQuadratic(const ScalarType a,
                                                  const ScalarType b,
                                                  const ScalarType c) {
    ScalarType discriminant = b * b - static_cast<ScalarType>(4) * a * c;
    // if (discriminant > static_cast<ScalarType>(0)) {
    //   discriminant = sqrt(discriminant);
    //   const ScalarType q =
    //       -(b + copysign(discriminant, b)) / static_cast<ScalarType>(2);
    //   const ScalarType sol1 = q / safelyTiny(a);
    //   const ScalarType sol2 = c / safelyTiny(q);
    //   // if (!isnan(sol1) && !isnan(sol2) &&
    //   //     !(fabs(q) < machine_epsilon<ScalarType>() &&
    //   //       fabs(a) < machine_epsilon<ScalarType>()) &&
    //   //     !(fabs(c) < machine_epsilon<ScalarType>() &&
    //   //       fabs(q) < machine_epsilon<ScalarType>())) {
    //   return sol1 < sol2 ? StackVector<ScalarType, 2>({sol1, sol2})
    //                      : StackVector<ScalarType, 2>({sol2, sol1});
    //   // }
    // }

    if (discriminant > static_cast<ScalarType>(0)) {
      if (a != static_cast<ScalarType>(0)) {
        discriminant = sqrt(discriminant);
        const ScalarType q =
            -static_cast<ScalarType>(0.5) * (b + copysign(discriminant, b));
        // if (b == static_cast<ScalarType>(0) && c ==
        // static_cast<ScalarType>(0))
        // {
        //   return StackVector<ScalarType, 2>(
        //       {static_cast<ScalarType>(0), static_cast<ScalarType>(0)});
        // } else if (q == static_cast<ScalarType>(0)) {
        //   return StackVector<ScalarType, 2>({static_cast<ScalarType>(0)});
        // }
        const ScalarType sol1 = q / safelyTiny(a);
        const ScalarType sol2 = c / safelyTiny(q);
        return sol1 < sol2 ? StackVector<ScalarType, 2>({sol1, sol2})
                            : StackVector<ScalarType, 2>({sol2, sol1});

      } else {
        return StackVector<ScalarType, 2>({-c / b});
      }
    }
    return StackVector<ScalarType, 2>();
  };

  // Returns solution to quadratic equation solve.
  // The smallest solution will always be first.
  template <class GradientType, class ScalarType>
  inline StackVector<std::pair<ScalarType, GradientType>, 2>
  solveQuadraticWithGradient(const ScalarType a, const ScalarType b,
                          const ScalarType c, const GradientType& a_grad,
                          const GradientType& b_grad,
                          const GradientType& c_grad) {
    using return_type = typename std::pair<ScalarType, GradientType>;
    const ScalarType discriminant = b * b - static_cast<ScalarType>(4) * a * c;
    if (discriminant > static_cast<ScalarType>(0)) {
      const auto d_grad = static_cast<ScalarType>(2) * b_grad * b -
                          static_cast<ScalarType>(4) * (a_grad * c + a * c_grad);
      if (a != static_cast<ScalarType>(0)) {
        const ScalarType sqrtd = sqrt(discriminant);
        const auto sqrtd_grad =
            d_grad / safelyEpsilon(sqrtd) / static_cast<ScalarType>(2);
        const ScalarType q =
            -(b + copysign(sqrtd, b)) / static_cast<ScalarType>(2);
        const auto q_grad =
            -(b_grad + sqrtd_grad * copysign(static_cast<ScalarType>(1), b)) /
            static_cast<ScalarType>(2);
        if (b == static_cast<ScalarType>(0) && c == static_cast<ScalarType>(0)) {
            const auto zero_return = std::pair<ScalarType, GradientType>(
                {static_cast<ScalarType>(0),
                GradientType(static_cast<ScalarType>(0))});
            return StackVector<return_type, 2>({zero_return, zero_return});
        } else if (q == static_cast<ScalarType>(0)) {
            const auto zero_return = std::pair<ScalarType, GradientType>(
                {static_cast<ScalarType>(0),
                GradientType(static_cast<ScalarType>(0))});
            return StackVector<return_type, 2>({zero_return});
        }
        const ScalarType sol1 = q / a;
        const ScalarType sol2 = c / q;
        const auto sol1_grad = (q_grad * a - q * a_grad) / safelyEpsilon(a * a);
        const auto sol2_grad = (c_grad * q - c * q_grad) / safelyEpsilon(q * q);
        const auto return_sol1 =
            std::pair<ScalarType, GradientType>({sol1, sol1_grad});
        const auto return_sol2 =
            std::pair<ScalarType, GradientType>({sol2, sol2_grad});
        return sol1 < sol2
                    ? StackVector<return_type, 2>({return_sol1, return_sol2})
                    : StackVector<return_type, 2>({return_sol2, return_sol1});

      } else {
        const ScalarType sol = -c / b;
        const auto sol_grad = (-c_grad * b + c * b_grad) / safelyEpsilon(b * b);
        const auto return_sol =
            std::pair<ScalarType, GradientType>({sol, sol_grad});
        return StackVector<return_type, 2>({return_sol});
      }
    }
    return StackVector<return_type, 2>();
  };  

}// namespace IRL

#endif // IRL_QUADRATIC_RECONSTRUCTION_QUADRATIC_HELPER_TPP_