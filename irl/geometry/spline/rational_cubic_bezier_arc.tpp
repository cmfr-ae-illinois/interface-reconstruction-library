#ifndef IRL_SPLINE_RATIONAL_CUBIC_BEZIER_ARC_TPP_
#define IRL_SPLINE_RATIONAL_BEZIER_ARC_TPP_

#include "irl/quadratic_reconstruction/rational_bezier_arc.h"

namespace IRL {
// Constructors
template <class ScalarType>
inline RationalCubicBezierArcBase<ScalarType>::RationalCubicBezierArcBase(void)
    : weight_0_m{ScalarType(0)},
      weight_1_m{ScalarType(0)},
      weight_2_m{ScalarType(0)},
      weight_3_m{ScalarType(0)},
      start_point_m{ScalarType(0), ScalarType(0), ScalarType(0)},
      control_point_1_m{ScalarType(0), ScalarType(0), ScalarType(0)},
      control_point_2_m{ScalarType(0), ScalarType(0), ScalarType(0)},
      end_point_m{ScalarType(0), ScalarType(0), ScalarType(0)} {}

template <class ScalarType>
inline RationalCubicBezierArcBase<ScalarType>::RationalCubicBezierArcBase(
    const PtBase<ScalarType>& a_start_pt,
    const PtBase<ScalarType>& a_control_pt_1,
    const PtBase<ScalarType>& a_control_pt_2,
    const PtBase<ScalarType>& a_end_pt, const ScalarType a_weight_1,
    const ScalarType a_weight_2)
    : weight_0_m(ScalarType(1)),
      weight_1_m(a_weight_1),
      weight_2_m(a_weight_2),
      weight_3_m(ScalarType(1)),
      start_point_m(a_start_pt),
      control_point_1_m(a_control_pt_1),
      control_point_2_m(a_control_pt_2),
      end_point_m(a_end_pt) {}

template <class ScalarType>
inline RationalCubicBezierArcBase<ScalarType>::RationalCubicBezierArcBase(
    const PtBase<ScalarType>& a_start_pt,
    const PtBase<ScalarType>& a_control_pt_1,
    const PtBase<ScalarType>& a_control_pt_2,
    const PtBase<ScalarType>& a_end_pt)
    : weight_0_m(ScalarType(1)),
      weight_1_m(ScalarType(1)),
      weight_2_m(ScalarType(1)),
      weight_3_m(ScalarType(1)),
      start_point_m(a_start_pt),
      control_point_1_m(a_control_pt_1),
      control_point_2_m(a_control_pt_2),
      end_point_m(a_end_pt) {}

template <class ScalarType>
inline RationalCubicBezierArcBase<ScalarType>::RationalCubicBezierArcBase(
    const PtBase<ScalarType>& a_start_pt,
    const NormalBase<ScalarType>& a_start_tangent,
    const PtBase<ScalarType>& a_end_pt,
    const NormalBase<ScalarType>& a_end_tangent)
    : weight_0_m(ScalarType(1)),
      weight_1_m(ScalarType(1)),
      weight_2_m(ScalarType(1)),
      weight_3_m(ScalarType(1)),
      start_point_m(a_start_pt),
      end_point_m(a_end_pt) {
  // Define Numbers
  ScalarType ONE = ScalarType(1);
  ScalarType TWO = ScalarType(2);
  ScalarType THREE = ScalarType(3);
  ScalarType FIVE = ScalarType(5);
  ScalarType SIX = ScalarType(6);
  ScalarType SEVEN = ScalarType(7);
  ScalarType NINE = ScalarType(9);
  ScalarType EIGHTEEN = ScalarType(18);
  ScalarType TWENTYSEVEN = ScalarType(27);
  ScalarType THIRTYTWO = ScalarType(32);
  // Define Supplementary Variables
  const PtBase<ScalarType> dP = a_end_pt - a_start_pt;
  const NormalBase<ScalarType> dPN = {dP[0], dP[1], dP[2]};

  const ScalarType a11 = THIRTYTWO / TWENTYSEVEN;
  const ScalarType a12 = a_start_tangent * a_end_tangent / EIGHTEEN;
  const ScalarType a21 = a_start_tangent * a_end_tangent / EIGHTEEN;
  const ScalarType a22 = FIVE / NINE;

  const ScalarType b1 = -SEVEN * (a_start_tangent * dPN) / EIGHTEEN;
  const ScalarType b2 = SEVEN * (a_end_tangent * dPN) / SIX;

  // Solve
  if (fabs(a11 * a22 - a12 * a21) <= ScalarType(1.0e-12) ||
      // This case should only happen if unit
      // tangents are not supplied.
      (fabs(b1) <= ScalarType(1.0e-12) && fabs(b2) <= ScalarType(1.0e-12))) {
    // This happens if the tangents are parallel, and are also perpendicular to
    // the chord. In this case, we will limit how far the control points can go.
    // The tangent direction is correct, but the bending energy is not
    // minimized.
    const ScalarType alpha = TWO * dPN.calculateMagnitude();
    control_point_1_m =
        a_start_pt + alpha * a_start_tangent /
                         (THREE * a_start_tangent.calculateMagnitude());
    control_point_2_m =
        a_end_pt -
        alpha * a_end_tangent / (THREE * a_end_tangent.calculateMagnitude());
  } else {  // This construction is bending energy minimizing. That is to say,
            // the integral of the square of the second derivative magnitude is
            // minimized.
    const ScalarType inv_det = ONE / (a11 * a22 - a12 * a21);
    const ScalarType alpha_0 = inv_det * (a22 * b1 - a12 * b2);
    const ScalarType alpha_1 = inv_det * (a11 * b2 - a21 * b1);

    control_point_1_m = a_start_pt + alpha_0 * a_start_tangent / THREE;
    control_point_2_m = a_end_pt - alpha_1 * a_end_tangent / THREE;
  }
}

// Evaluation Methods
template <class ScalarType>
inline const PtBase<ScalarType> RationalCubicBezierArcBase<ScalarType>::point(
    const ScalarType t) const {
  /* Defining constants and types */
  const ScalarType ZERO = ScalarType(0);
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);

  assert(t >= ZERO && t <= ONE);
  if (weight_1_m > ScalarType(1.0e15) || weight_2_m > ScalarType(1.0e15)) {
    // Implement later, we will always be using weight = 1 for now.
    return PtBase<ScalarType>(-1, -1, -1);
    std::cout << "Warning: weight is too large for accurate point evaluation, "
                 "returning "
                 "(-1,-1,-1)"
              << std::endl;
  } else {
    const ScalarType denominator =
        (ONE - t) * (ONE - t) * (ONE - t) +
        THREE * weight_1_m * t * (ONE - t) * (ONE - t) +
        THREE * weight_2_m * t * t * (ONE - t) + t * t * t;
    auto numerator = PtBase<ScalarType>(
        (ONE - t) * (ONE - t) * (ONE - t) * start_point_m +
        THREE * weight_1_m * t * (ONE - t) * (ONE - t) * control_point_1_m +
        THREE * weight_2_m * t * t * (ONE - t) * control_point_2_m +
        t * t * t * end_point_m);
    return numerator / denominator;
  }
}

template <class ScalarType>
inline const PtBase<ScalarType>
RationalCubicBezierArcBase<ScalarType>::derivative(const ScalarType t) const {
  /* Defining constants and types */
  const ScalarType ZERO = ScalarType(0);
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);
  const ScalarType SIX = ScalarType(6);
  assert(t >= ZERO && t <= ONE);
  if (weight_1_m > ScalarType(1.0e15) || weight_2_m > ScalarType(1.0e15)) {
    // Implement later, we will always be using weight = 1 for now.
    return PtBase<ScalarType>(-1, -1, -1);
    std::cout << "Warning: weight is too large for derivative evaluation, "
                 "returning "
                 "(-1,-1,-1)"
              << std::endl;
  } else {
    if (ZERO > t || t > ONE) {
      std::cout << "Warning: t is out of bounds for derivative evaluation, "
                   "returning "
                   "(0,0,0)"
                << std::endl;
      return PtBase<ScalarType>(ZERO, ZERO, ZERO);
    }
    const ScalarType B0 = (ONE - t) * (ONE - t) * (ONE - t);
    const ScalarType B1 = THREE * t * (ONE - t) * (ONE - t);
    const ScalarType B2 = THREE * t * t * (ONE - t);
    const ScalarType B3 = t * t * t;

    const ScalarType dB0 = -THREE * (ONE - t) * (ONE - t);
    const ScalarType dB1 = THREE * (THREE * t - ONE) * (t - ONE);
    const ScalarType dB2 = THREE * t * (TWO - THREE * t);
    const ScalarType dB3 = THREE * t * t;
    // Set up Demoninator
    ScalarType denominator =
        weight_0_m * B0 + weight_1_m * B1 + weight_2_m * B2 + weight_3_m * B3;
    ScalarType invDenominator = ONE / denominator;
    ScalarType dDenominator = weight_0_m * dB0 + weight_1_m * dB1 +
                              weight_2_m * dB2 + weight_3_m * dB3;
    // Set up Numerator
    auto numerator = PtBase<ScalarType>(
        weight_0_m * B0 * start_point_m + weight_1_m * B1 * control_point_1_m +
        weight_2_m * B2 * control_point_2_m + weight_3_m * B3 * end_point_m);
    auto dNumerator = PtBase<ScalarType>(weight_0_m * dB0 * start_point_m +
                                         weight_1_m * dB1 * control_point_1_m +
                                         weight_2_m * dB2 * control_point_2_m +
                                         weight_3_m * dB3 * end_point_m);
    return (dNumerator * denominator - numerator * dDenominator) *
           invDenominator * invDenominator;
  }
}

template <class ScalarType>
inline ScalarType RationalCubicBezierArcBase<ScalarType>::arc_length(
    void) const {
  /* Defining constants and types */
  const ScalarType ZERO = ScalarType(0);
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);
  const ScalarType FIVE = ScalarType(5);
  const ScalarType EIGHT = ScalarType(8);
  const ScalarType EIGHTEEN = ScalarType(18);

  // 3-point quadrature rule for arc_length calculation
  const ScalarType t0 = (ONE - sqrt(THREE / FIVE)) / TWO;
  const ScalarType t1 = ONE / TWO;
  const ScalarType t2 = (ONE + sqrt(THREE / FIVE)) / TWO;
  const ScalarType w0 = FIVE / EIGHTEEN;
  const ScalarType w1 = EIGHT / EIGHTEEN;
  const ScalarType w2 = FIVE / EIGHTEEN;
  PtBase<ScalarType> pt_0 = derivative(t0);
  PtBase<ScalarType> pt_1 = derivative(t1);
  PtBase<ScalarType> pt_2 = derivative(t2);
  const ScalarType norm0 =
      sqrt(pt_0[0] * pt_0[0] + pt_0[1] * pt_0[1] + pt_0[2] * pt_0[2]);
  const ScalarType norm1 =
      sqrt(pt_1[0] * pt_1[0] + pt_1[1] * pt_1[1] + pt_1[2] * pt_1[2]);
  const ScalarType norm2 =
      sqrt(pt_2[0] * pt_2[0] + pt_2[1] * pt_2[1] + pt_2[2] * pt_2[2]);
  return w0 * norm0 + w1 * norm1 + w2 * norm2;
}

template <class ScalarType>
inline RationalCubicBezierArcBase<ScalarType>
RationalCubicBezierArcBase<ScalarType>::moveToReferenceFrame(
    const PtBase<ScalarType>& datum,
    const ReferenceFrameBase<ScalarType>& frame) const {
  // for translatioon
  Eigen::Matrix<ScalarType, 3, 1> x0;
  x0 << datum[0], datum[1], datum[2];

  // rotation matrix
  Eigen::Matrix<ScalarType, 3, 3> R;
  for (int i = 0; i < 3; i++) {
    R.col(i) << frame[i][0], frame[i][1], frame[i][2];
  }

  auto transform = [&](const PtBase<ScalarType>& p_local) {
    Eigen::Matrix<ScalarType, 3, 1> x_loc;
    x_loc << p_local[0], p_local[1], p_local[2];
    Eigen::Matrix<ScalarType, 3, 1> x_global = x0 + R * x_loc;
    return PtBase<ScalarType>(x_global[0], x_global[1], x_global[2]);
  };

  // transformed bezier arc
  return RationalCubicBezierArcBase<ScalarType>(
      transform(start_point_m), transform(control_point_1_m),
      transform(control_point_2_m), transform(end_point_m), weight_1_m,
      weight_2_m);
}

// Getters
template <class ScalarType>
inline const ScalarType& RationalCubicBezierArcBase<ScalarType>::weight(
    const int index) const {
  assert(index >= 0 && index <= 3);
  if (index == 0) {
    return weight_0_m;
  } else if (index == 1) {
    return weight_1_m;
  } else if (index == 2) {
    return weight_2_m;
  } else {
    return weight_3_m;
  }
}

template <class ScalarType>
inline ScalarType& RationalCubicBezierArcBase<ScalarType>::weight(
    const int index) {
  assert(index >= 0 && index <= 3);
  if (index == 0) {
    return weight_0_m;
  } else if (index == 1) {
    return weight_1_m;
  } else if (index == 2) {
    return weight_2_m;
  } else {
    return weight_3_m;
  }
}

template <class ScalarType>
inline const PtBase<ScalarType>&
RationalCubicBezierArcBase<ScalarType>::start_point(void) const {
  return start_point_m;
}

template <class ScalarType>
inline PtBase<ScalarType>& RationalCubicBezierArcBase<ScalarType>::start_point(
    void) {
  return start_point_m;
}

template <class ScalarType>
inline const PtBase<ScalarType>&
RationalCubicBezierArcBase<ScalarType>::end_point(void) const {
  return end_point_m;
}

template <class ScalarType>
inline PtBase<ScalarType>& RationalCubicBezierArcBase<ScalarType>::end_point(
    void) {
  return end_point_m;
}

template <class ScalarType>
inline const PtBase<ScalarType>&
RationalCubicBezierArcBase<ScalarType>::control_point_1(void) const {
  return control_point_1_m;
}

template <class ScalarType>
inline PtBase<ScalarType>&
RationalCubicBezierArcBase<ScalarType>::control_point_1(void) {
  return control_point_1_m;
}

template <class ScalarType>
inline const PtBase<ScalarType>&
RationalCubicBezierArcBase<ScalarType>::control_point_2(void) const {
  return control_point_2_m;
}

template <class ScalarType>
inline PtBase<ScalarType>&
RationalCubicBezierArcBase<ScalarType>::control_point_2(void) {
  return control_point_2_m;
}

// Overloaded Operators
template <class ScalarType>
inline RationalCubicBezierArcBase<ScalarType>&
RationalCubicBezierArcBase<ScalarType>::operator+=(
    const PtBase<ScalarType>& a_rhs) {
  this->start_point() += a_rhs;
  this->control_point_1() += a_rhs;
  this->control_point_2() += a_rhs;
  this->end_point() += a_rhs;
  // This is effectively a translation of the curve, so we do not need to update
  // the weights.
  return (*this);
}

template <class ScalarType>
inline RationalCubicBezierArcBase<ScalarType>
RationalCubicBezierArcBase<ScalarType>::operator-(void) const {
  // This inverts the direction of the spline
  return RationalCubicBezierArcBase(end_point_m, control_point_2_m,
                                    control_point_1_m, start_point_m,
                                    weight_3_m, weight_2_m);
}

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out,
    const RationalCubicBezierArcBase<ScalarType>& a_rational_cubic_bezier_arc) {
  std::streamsize old_precision = out.precision();
  out.precision(16);
  out << "Rational Cubic Bezier Arc: " << std::endl;
  out << "Start Point: " << a_rational_cubic_bezier_arc.start_point()
      << std::endl;
  out << "Control Point 1: " << a_rational_cubic_bezier_arc.control_point_1()
      << std::endl;
  out << "Control Point 2: " << a_rational_cubic_bezier_arc.control_point_2()
      << std::endl;
  out << "End Point: " << a_rational_cubic_bezier_arc.end_point() << std::endl;
  out << "Weights: " << a_rational_cubic_bezier_arc.weight(0) << ", "
      << a_rational_cubic_bezier_arc.weight(1) << ", "
      << a_rational_cubic_bezier_arc.weight(2) << ", "
      << a_rational_cubic_bezier_arc.weight(3) << std::endl;
  out.precision(old_precision);
  return out;
}

template <class ScalarType>
inline void RationalCubicBezierArcBase<ScalarType>::saveToVTK(
    const std::string& filename, const int nsamples) {
  ScalarType ZERO = ScalarType(0);
  ScalarType ONE = ScalarType(1);
  std::vector<ScalarType> uset(nsamples, ZERO);
  std::vector<PtBase<ScalarType>> curve(nsamples,
                                        PtBase<ScalarType>(ZERO, ZERO, ZERO));
  for (int i = 0; i < nsamples; i++) {
    uset[i] = (ONE) * static_cast<ScalarType>(i) /
              static_cast<ScalarType>(nsamples - 1);

    curve[i] = this->point(uset[i]);
  }

  std::ofstream file;
  file.open(filename + std::string(".vtu"));
  file << "<VTKFile type=\"UnstructuredGrid\">\n<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << nsamples << "\" NumberOfCells=\"" << 1
       << "\">\n";
  file << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";
  for (int i = 0; i < nsamples; i++) {
    file << std::scientific << std::setprecision(15) << curve[i][0] << " "
         << curve[i][1] << curve[i][2] << " \n";
  }
  file << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" "
          "Name=\"connectivity\" format=\"ascii\">\n";
  int count = 0;
  for (int i = 0; i < nsamples; i++) {
    file << count++ << " ";
  }
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  file << nsamples << " ";
  file << "\n</DataArray>\n";
  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  file << "7 ";
  file << "\n</DataArray>\n";
  file << "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  file.close();
}

}  // Namespace IRL

#endif