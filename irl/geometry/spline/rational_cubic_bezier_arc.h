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
class RationalCubicBezierArcBase {
 public:
  using value_type = ScalarType;
  /// \brief Default Constructor;
  RationalCubicBezierArcBase(void);

  /// \brief Constructor that initializes a rational cubic Bèzier arc
  RationalCubicBezierArcBase(const PtBase<ScalarType>& a_start_pt,
                             const PtBase<ScalarType>& a_control_pt_1,
                             const PtBase<ScalarType>& a_control_pt_2,
                             const PtBase<ScalarType>& a_end_pt,
                             const ScalarType a_weight_1,
                             const ScalarType a_weight_2);

  /// \brief Constructor that initializes a nonrational cubic Bèzier arc
  RationalCubicBezierArcBase(const PtBase<ScalarType>& a_start_pt,
                             const PtBase<ScalarType>& a_control_pt_1,
                             const PtBase<ScalarType>& a_control_pt_2,
                             const PtBase<ScalarType>& a_end_pt);

  /// \brief Constructor that initializes a nonrational cubic Bèzier arc by
  /// computing the control points
  RationalCubicBezierArcBase(const PtBase<ScalarType>& a_start_pt,
                             const NormalBase<ScalarType>& a_start_tangent,
                             const PtBase<ScalarType>& a_end_pt,
                             const NormalBase<ScalarType>& a_end_tangent);
  /// \brief Return const weight, given index
  const ScalarType& weight(const int index) const;
  /// \brief Return const reference to stored start point.
  const PtBase<ScalarType>& start_point(void) const;
  /// \brief Return const referene to stored control point 1.
  const PtBase<ScalarType>& control_point_1(void) const;
  /// \brief Return const referene to stored control point 2.
  const PtBase<ScalarType>& control_point_2(void) const;
  /// \brief Return const reference to stored end point.
  const PtBase<ScalarType>& end_point(void) const;
  /// \brief Return point evaluted at parameter a_t.
  const PtBase<ScalarType> point(const ScalarType a_t) const;
  /// \brief Return derivative of the curve with respect to a_t
  const PtBase<ScalarType> derivative(const ScalarType a_t) const;
  /// \brief return approximation of arc_length
  ScalarType arc_length(void) const;
  /// \brief Return split arcs.
  std::pair<RationalCubicBezierArcBase, RationalCubicBezierArcBase> split(
      const ScalarType a_t) const;
  /// \brief Return arc in global coordinates
  RationalCubicBezierArcBase<ScalarType> moveToReferenceFrame(
      const PtBase<ScalarType>& datum,
      const ReferenceFrameBase<ScalarType>& frame) const;
  /// \brief Return const weight.
  ScalarType& weight(int index);
  /// \brief Return const reference to stored start point.
  PtBase<ScalarType>& start_point(void);
  /// \brief Return const reference to stored control point 1.
  PtBase<ScalarType>& control_point_1(void);
  /// \brief Return const reference to stored control point 2.
  PtBase<ScalarType>& control_point_2(void);
  /// \brief Return const reference to stored end point.
  PtBase<ScalarType>& end_point(void);
  /// \Brief Export Spline to VTK File
  void saveToVTK(const std::string& filename, const int nsamples = 100);
  /// \brief overload += operator.
  RationalCubicBezierArcBase& operator+=(const PtBase<ScalarType>& a_rhs);
  /// \brief unary minus operator
  RationalCubicBezierArcBase operator-(void) const;

  /// \brief Default destructor.
  ~RationalCubicBezierArcBase(void) = default;

 private:
  PtBase<ScalarType> start_point_m;      // Start point.
  PtBase<ScalarType> control_point_1_m;  // Control point 1.
  PtBase<ScalarType> control_point_2_m;  // Control point 2.
  PtBase<ScalarType> end_point_m;        // End point.
  ScalarType weight_0_m;                 // Weight 0
  ScalarType weight_1_m;                 // Weight 1
  ScalarType weight_2_m;                 // Weight 2
  ScalarType weight_3_m;                 // Weight 3
};

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out,
    const RationalCubicBezierArcBase<ScalarType>& a_rational_cubic_bezier_arc);

using RationalCubicBezierArc = RationalCubicBezierArcBase<double>;

}  // namespace IRL

#include "irl/geometry/spline/rational_cubic_bezier_arc.tpp"

#endif
