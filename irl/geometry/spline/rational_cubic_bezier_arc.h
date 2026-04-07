#ifndef IRL_SPLINE_RATIONAL_CUBIC_BEZIER_ARC_H_
#define IRL_SPLINE_RATIONAL_CUBIC_BEZIER_ARC_H_

#include <cstdint>
#include <utility>

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/geometry/general/pt.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"

#include <Eigen/Dense>

namespace IRL {

/// \brief Rational Cubic Bezier arc defined by endpoints, 2 control points, and
/// 4 weights.
template <class ScalarType>
class RationalCubicBezierArc {
 public:
  using value_type = ScalarType;
  /// \brief Default Constructor;
  RationalCubicBezierArc(void);

  // Other Constrcutors here

  // const Scalar

 private:
}
