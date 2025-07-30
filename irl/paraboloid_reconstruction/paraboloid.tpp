// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_
#define IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_

namespace IRL {

template <class ScalarType>
inline ParaboloidBase<ScalarType>::ParaboloidBase(void)
    : datum_m(),
      frame_m(),
      paraboloid_m(),
      place_infinite_shortcut_m({false, false}) {}

template <class ScalarType>
inline ParaboloidBase<ScalarType>::ParaboloidBase(
    const PtBase<ScalarType>& a_datum,
    const ReferenceFrameBase<ScalarType>& a_reference_frame,
    const ScalarType a_coef_a, const ScalarType a_coef_b)
    : datum_m(a_datum),
      frame_m(a_reference_frame),
      paraboloid_m(std::array<ScalarType, 2>({a_coef_a, a_coef_b})),
      place_infinite_shortcut_m({false, false}) {
  // assert(frame_m.isOrthonormalBasis());
}

template <class ScalarType>
inline ParaboloidBase<ScalarType> ParaboloidBase<ScalarType>::fromDerivatives(
    const PtBase<ScalarType>& a_datum,
    const Eigen::Vector<ScalarType, 3>& a_gradF,
    const Eigen::Matrix<ScalarType, 3, 3>& a_hessF) {
  using Vector2Type = Eigen::Vector<ScalarType, 2>;
  using Vector3Type = Eigen::Vector<ScalarType, 3>;
  using Matrix22Type = Eigen::Matrix<ScalarType, 2, 2>;
  using Matrix33Type = Eigen::Matrix<ScalarType, 3, 3>;
  using Matrix32Type = Eigen::Matrix<ScalarType, 3, 2>;
  using NormalType = NormalBase<ScalarType>;
  using ReferenceFrameType = ReferenceFrameBase<ScalarType>;
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType HALF = ONE / TWO;
  const Matrix33Type hessF = 0.5 * (a_hessF + a_hessF.transpose());

  // This uses the method described in
  // https://www.geometrictools.com/Documentation/PrincipalCurvature.pdf
  const ScalarType inv_gradF_norm = ONE / a_gradF.norm();
  const NormalType normal =
      NormalType(a_gradF(0), a_gradF(1), a_gradF(2)) * inv_gradF_norm;
  ReferenceFrameType frame = ReferenceFrameType::fromNormal(normal);
  Matrix32Type J;
  for (UnsignedIndex_t i = 0; i < 3; i++) {
    for (UnsignedIndex_t j = 0; j < 2; j++) {
      J(i, j) = frame[j][i];
    }
  }
  const Matrix22Type A = (J.transpose() * hessF * J) * inv_gradF_norm;
  ////////// TODO: compute eigenvalues and eigenvectors "by hand"
  Eigen::EigenSolver<Matrix22Type> eigensolver(A);
  const ScalarType eval1 = eigensolver.eigenvalues()(0).real();
  const ScalarType eval2 = eigensolver.eigenvalues()(1).real();
  const Vector2Type evec1 =
      Vector2Type(eigensolver.eigenvectors()(0, 0).real(),
                  eigensolver.eigenvectors()(1, 0).real());
  //////////
  const Vector3Type T1 = J * evec1;
  frame[0] = NormalType(T1(0), T1(1), T1(2));
  frame[0].normalize();
  frame[1] = crossProduct(frame[2], frame[0]);
  return ParaboloidBase<ScalarType>(a_datum, frame, HALF * eval1, HALF * eval2);
}

template <class ScalarType>
inline ParaboloidBase<ScalarType> ParaboloidBase<ScalarType>::createAlwaysAbove(
    void) {
  ParaboloidBase<ScalarType> par;
  par.markAsAlwaysAbove();
  return par;
}

template <class ScalarType>
inline ParaboloidBase<ScalarType> ParaboloidBase<ScalarType>::createAlwaysBelow(
    void) {
  ParaboloidBase<ScalarType> par;
  par.markAsAlwaysBelow();
  return par;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::setDatum(
    const PtBase<ScalarType>& a_datum) {
  datum_m = a_datum;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::setReferenceFrame(
    const ReferenceFrameBase<ScalarType>& a_reference_frame) {
  assert(a_reference_frame.isOrthonormalBasis());
  frame_m = a_reference_frame;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::setAlignedParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid) {
  paraboloid_m = a_aligned_paraboloid;
}

template <class ScalarType>
inline const PtBase<ScalarType>& ParaboloidBase<ScalarType>::getDatum(
    void) const {
  return datum_m;
}

template <class ScalarType>
inline const ReferenceFrameBase<ScalarType>&
ParaboloidBase<ScalarType>::getReferenceFrame(void) const {
  return frame_m;
}

template <class ScalarType>
inline const AlignedParaboloidBase<ScalarType>&
ParaboloidBase<ScalarType>::getAlignedParaboloid(void) const {
  return paraboloid_m;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::markAsRealReconstruction(void) {
  place_infinite_shortcut_m[0] = false;
  place_infinite_shortcut_m[1] = false;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::markAsAlwaysAbove(void) {
  place_infinite_shortcut_m[0] = true;
  place_infinite_shortcut_m[1] = false;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::markAsAlwaysBelow(void) {
  place_infinite_shortcut_m[0] = false;
  place_infinite_shortcut_m[1] = true;
}

template <class ScalarType>
inline bool ParaboloidBase<ScalarType>::isAlwaysAbove(void) const {
  return place_infinite_shortcut_m[0];
}

template <class ScalarType>
inline bool ParaboloidBase<ScalarType>::isAlwaysBelow(void) const {
  return place_infinite_shortcut_m[1];
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::regenerateAtLocation(
    const PtBase<ScalarType>& a_pt) {
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType HALF = ONE / TWO;

  // Get paraboloid information
  const ScalarType a = paraboloid_m.a();
  const ScalarType b = paraboloid_m.b();

  // Bring point in local frame of reference
  const Pt pt_tmp = a_pt - datum_m;
  PtBase<ScalarType> local_pt;
  for (UnsignedIndex_t n = 0; n < 3; ++n) {
    local_pt[n] = frame_m[n] * pt_tmp;
  }
  local_pt[2] = -a * local_pt[0] * local_pt[0] - b * local_pt[1] * local_pt[1];

  // Compute local derivatives
  const Eigen::Vector<ScalarType, 3> gradF(TWO * a * local_pt[0],
                                           TWO * b * local_pt[1], ONE);
  Eigen::Matrix<ScalarType, 3, 3> hessF =
      Eigen::Matrix<ScalarType, 3, 3>::Zero();
  hessF(0, 0) = TWO * a;
  hessF(1, 1) = TWO * b;
  Eigen::Matrix<ScalarType, 3, 3> adjHessF =
      Eigen::Matrix<ScalarType, 3, 3>::Zero();
  adjHessF(2, 2) = TWO * TWO * a * b;
  auto new_normal = NormalBase<ScalarType>(gradF(0), gradF(1), gradF(2));
  new_normal.normalize();

  // Based on Goldman 2005
  ScalarType H = gradF.transpose() * (hessF * gradF);
  H -= gradF.squaredNorm() * hessF.trace();
  H /= TWO * gradF.squaredNorm() * gradF.norm();
  ScalarType K = gradF.transpose() * (adjHessF * gradF);
  K /= gradF.squaredNorm() * gradF.squaredNorm();
  const ScalarType k1 = -H + sqrt(maximum(H * H - K, ZERO));
  const ScalarType k2 = -H - sqrt(maximum(H * H - K, ZERO));

  // Compute principal directions
  const ScalarType B = a - b, C = -gradF(1) * a, E = gradF(0) * b;
  const ScalarType U = TWO * gradF(0) * gradF(1) * a;
  const ScalarType V = TWO * (B - C * gradF(1) - E * gradF(0));
  const ScalarType W = -TWO * gradF(1) * gradF(0) * b;
  const ScalarType delta = V * V - TWO * TWO * U * W;
  const ScalarType X1 = -V - sqrt(maximum(0., delta));
  auto T1 =
      NormalBase<ScalarType>(X1, TWO * U, -X1 * gradF(0) - TWO * U * gradF(1));
  const ScalarType normT1 = magnitude(T1);
  ReferenceFrame new_local_frame;
  if (normT1 > DBL_EPSILON) {
    T1 *= ONE / normT1;
    new_local_frame[2] = new_normal;
    new_local_frame[0] = T1;
    new_local_frame[1] = crossProduct(new_local_frame[2], new_local_frame[0]);
  } else {  /// The point is umbilical
    new_local_frame = ReferenceFrameBase<ScalarType>::fromNormal(new_normal);
  }

  // Now move frame back to canonical frame and return paraboloid
  ReferenceFrame new_frame;
  for (int n = 0; n < 3; n++) {
    new_frame[n] = NormalBase<ScalarType>(ZERO, ZERO, ZERO);
    for (int d = 0; d < 3; d++) {
      new_frame[n] += new_local_frame[n][d] * frame_m[d];
    }
    new_frame[n].normalize();
  }

  // Bring point back to canonical frame of reference
  PtBase<ScalarType> projection(ZERO, ZERO, ZERO);
  for (UnsignedIndex_t n = 0; n < 3; ++n) {
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      projection[n] += frame_m[d][n] * local_pt[d];
    }
  }

  datum_m += projection;
  frame_m = new_frame;
  paraboloid_m.a() = HALF * k1;
  paraboloid_m.b() = HALF * k2;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::serialize(ByteBuffer* a_buffer) const {
  datum_m.serialize(a_buffer);
  for (std::size_t d = 0; d < 3; ++d) {
    frame_m[d].serialize(a_buffer);
  }
  paraboloid_m.serialize(a_buffer);
  const UnsignedIndex_t bool_to_int =
      (place_infinite_shortcut_m[0] ? 1 : 0) +
      2 * (place_infinite_shortcut_m[1] ? 1 : 0);
  a_buffer->pack(&bool_to_int, 1);
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::unpackSerialized(ByteBuffer* a_buffer) {
  datum_m.unpackSerialized(a_buffer);
  NormalBase<ScalarType> normal;
  for (std::size_t d = 0; d < 3; ++d) {
    frame_m[d].unpackSerialized(a_buffer);
  }
  paraboloid_m.unpackSerialized(a_buffer);
  UnsignedIndex_t int_to_bool = 0;
  a_buffer->unpack(&int_to_bool, 1);
  place_infinite_shortcut_m[0] = (int_to_bool % 2 == 1) ? true : false;
  place_infinite_shortcut_m[1] = (int_to_bool / 2 == 1) ? true : false;
}

template <class ScalarType>
inline PtBase<ScalarType> conicCenter(
    const PlaneBase<ScalarType>& a_plane,
    const AlignedParaboloidBase<ScalarType>& a_paraboloid) {
  const auto& face_normal = a_plane.normal();
  const auto& face_distance = a_plane.distance();
  return PtBase<ScalarType>(
      face_normal[0] / safelyTiny(static_cast<ScalarType>(2) *
                                  a_paraboloid.a() * face_normal[2]),
      face_normal[1] / safelyTiny(static_cast<ScalarType>(2) *
                                  a_paraboloid.b() * face_normal[2]),
      -face_normal[0] * face_normal[0] /
              safelyTiny(static_cast<ScalarType>(2) * a_paraboloid.a() *
                         face_normal[2] * face_normal[2]) -
          face_normal[1] * face_normal[1] /
              safelyTiny(static_cast<ScalarType>(2) * a_paraboloid.b() *
                         face_normal[2] * face_normal[2]) +
          face_distance / safelyTiny(face_normal[2]));
}

template <class ScalarType>
inline NormalBase<ScalarType> getParaboloidSurfaceNormal(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtBase<ScalarType>& a_pt) {
  return NormalBase<ScalarType>(
      static_cast<ScalarType>(2) * a_paraboloid.a() * a_pt[0],
      static_cast<ScalarType>(2) * a_paraboloid.b() * a_pt[1],
      static_cast<ScalarType>(1));
};

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient getParaboloidSurfaceNormalWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtTypeWithGradient& a_pt) {
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(static_cast<ScalarType>(0)),
       B_grad = gradient_type(static_cast<ScalarType>(0));
  const auto& pt = a_pt.getPt();
  const auto& pt_grad = a_pt.getData();
  const auto surface_normal = PtBase<ScalarType>(
      static_cast<ScalarType>(2) * A * pt[0],
      static_cast<ScalarType>(2) * B * pt[1], static_cast<ScalarType>(1));
  auto surface_normal_withgrad = PtTypeWithGradient(surface_normal);
  auto& surface_normal_grad = surface_normal_withgrad.getData();
  surface_normal_grad[0] =
      static_cast<ScalarType>(2) * (A_grad * pt[0] + A * pt_grad[0]);
  surface_normal_grad[1] =
      static_cast<ScalarType>(2) * (B_grad * pt[1] + B * pt_grad[1]);
  // surface_normal_grad[2] = 0.0;
  return surface_normal_withgrad;
};

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongLineOntoParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const ScalarType a = (a_paraboloid.a() * a_line[0] * a_line[0] +
                        a_paraboloid.b() * a_line[1] * a_line[1]);
  const ScalarType b = (a_line[2] +
                        static_cast<ScalarType>(2) * a_paraboloid.a() *
                            a_starting_pt[0] * a_line[0] +
                        static_cast<ScalarType>(2) * a_paraboloid.b() *
                            a_starting_pt[1] * a_line[1]);
  const ScalarType c = (a_starting_pt[2] +
                        a_paraboloid.a() * a_starting_pt[0] * a_starting_pt[0] +
                        a_paraboloid.b() * a_starting_pt[1] * a_starting_pt[1]);
  // check if starting point is on paraboloid (then solution = 0)
  // if (std::fabs(c) < 5.0 * DBL_EPSILON) {
  //   return a_starting_pt;
  // } else {
  const StackVector<ScalarType, 2> solutions =
      solveQuadratic<ScalarType>(a, b, c);
  if (solutions.size() == 0) {
    std::cout << "No solution found for projection on paraboloid" << a_line
              << a_starting_pt << std::endl;
  }
  assert(solutions.size() > 0);
  if (solutions.size() == 1) {
    return a_starting_pt + a_line * solutions[0];
  } else {
    if (abs(solutions[0]) < abs(solutions[1])) {
      return a_starting_pt + a_line * solutions[0];
    } else {
      return a_starting_pt + a_line * solutions[1];
    }
  }
  // }
}

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongHalfLineOntoParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const ScalarType a = (a_paraboloid.a() * a_line[0] * a_line[0] +
                        a_paraboloid.b() * a_line[1] * a_line[1]);
  const ScalarType b = (a_line[2] +
                        static_cast<ScalarType>(2) * a_paraboloid.a() *
                            a_starting_pt[0] * a_line[0] +
                        static_cast<ScalarType>(2) * a_paraboloid.b() *
                            a_starting_pt[1] * a_line[1]);
  const ScalarType c = (a_starting_pt[2] +
                        a_paraboloid.a() * a_starting_pt[0] * a_starting_pt[0] +
                        a_paraboloid.b() * a_starting_pt[1] * a_starting_pt[1]);
  // check if starting point is on paraboloid (then solution = 0)
  // if (std::fabs(c) < DBL_EPSILON) {
  //   return Pt(static_cast<ScalarType>(DBL_MAX),
  //             static_cast<ScalarType>(DBL_MAX),
  //             static_cast<ScalarType>(DBL_MAX));
  // } else {
  const StackVector<ScalarType, 2> solutions =
      solveQuadratic<ScalarType>(a, b, c);
  if (solutions.size() == 0) {
    return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                              static_cast<ScalarType>(DBL_MAX),
                              static_cast<ScalarType>(DBL_MAX));
  }
  // assert(solutions.size() > 0);
  if (solutions.size() == 1) {
    // assert(solutions[0] >= static_cast<ScalarType>(0));
    if (solutions[0] < machine_epsilon<ScalarType>()) {
      return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX));
    }
    return a_starting_pt + a_line * solutions[0];
  } else {
    // assert(maximum(solutions[0], solutions[1]) >=
    // static_cast<ScalarType>(0));
    if ((solutions[1] < static_cast<ScalarType>(0))) {
      return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX));

    } else {
      const ScalarType distance_along_line =
          solutions[0] > static_cast<ScalarType>(0)
              ? minimum(solutions[0], solutions[1])
              : maximum(solutions[0], solutions[1]);
      return a_starting_pt + a_line * distance_along_line;
    }
  }
  // }
}

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient projectPtAlongHalfLineOntoParaboloidWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtTypeWithGradient& a_line, const PtTypeWithGradient& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  ScalarType EPSILON;
  if constexpr (std::is_same<ScalarType, Quad_t>::value) {
    EPSILON = FLT128_EPSILON;
  } else {
    EPSILON = DBL_EPSILON;
  }
  const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(static_cast<ScalarType>(0)),
       B_grad = gradient_type(static_cast<ScalarType>(0));
  A_grad.setGradA(static_cast<ScalarType>(1));
  B_grad.setGradB(static_cast<ScalarType>(1));
  const PtBase<ScalarType>& line = a_line.getPt();
  const auto& line_grad = a_line.getData();
  const PtBase<ScalarType>& starting_pt = a_starting_pt.getPt();
  const auto& starting_pt_grad = a_starting_pt.getData();
  const ScalarType a = (A * line[0] * line[0] + B * line[1] * line[1]);
  const ScalarType b =
      (line[2] + static_cast<ScalarType>(2) * A * starting_pt[0] * line[0] +
       static_cast<ScalarType>(2) * B * starting_pt[1] * line[1]);
  const ScalarType c = (starting_pt[2] + A * starting_pt[0] * starting_pt[0] +
                        B * starting_pt[1] * starting_pt[1]);
  // check if starting point is on paraboloid(then solution = 0)
  if (fabs(c) < static_cast<ScalarType>(5) * EPSILON) {
    return a_starting_pt;
  } else {
    const auto a_grad =
        A_grad * line[0] * line[0] + B_grad * line[1] * line[1] +
        static_cast<ScalarType>(2) * A * line_grad[0] * line[0] +
        static_cast<ScalarType>(2) * B * line_grad[1] * line[1];
    const auto b_grad =
        line_grad[2] +
        static_cast<ScalarType>(2) * A_grad * starting_pt[0] * line[0] +
        static_cast<ScalarType>(2) * A * starting_pt_grad[0] * line[0] +
        static_cast<ScalarType>(2) * A * starting_pt[0] * line_grad[0] +
        static_cast<ScalarType>(2) * B_grad * starting_pt[1] * line[1] +
        static_cast<ScalarType>(2) * B * starting_pt_grad[1] * line[1] +
        static_cast<ScalarType>(2) * B * starting_pt[1] * line_grad[1];
    const auto c_grad =
        starting_pt_grad[2] + A_grad * starting_pt[0] * starting_pt[0] +
        static_cast<ScalarType>(2) * A * starting_pt_grad[0] * starting_pt[0] +
        B_grad * starting_pt[1] * starting_pt[1] +
        static_cast<ScalarType>(2) * B * starting_pt_grad[1] * starting_pt[1];
    const auto solutions =
        solveQuadraticWithGradient(a, b, c, a_grad, b_grad, c_grad);
    if (solutions.size() == 0) {
      std::cout << "Solution not found on half-line: " << line << starting_pt
                << std::endl;
    }
    // assert(solutions.size() > 0);
    if (solutions.size() == 1) {
      const auto sol = solutions[0].first;
      const auto sol_grad = solutions[0].second;
      // assert(sol >= static_cast<ScalarType>(0));
      auto intersection =
          PtTypeWithGradient(PtBase<ScalarType>(starting_pt + sol * line));
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        intersection.getData()[d] =
            starting_pt_grad[d] + sol_grad * line[d] + sol * line_grad[d];
      }
      return intersection;
    } else {
      const auto sol0 = solutions[0].first;
      const auto sol1 = solutions[1].first;
      // assert(maximum(sol0, sol1) >= static_cast<ScalarType>(0));
      if (sol0 >= static_cast<ScalarType>(0)) {
        const auto sol0_grad = solutions[0].second;
        auto intersection =
            PtTypeWithGradient(PtBase<ScalarType>(starting_pt + sol0 * line));
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          intersection.getData()[d] =
              starting_pt_grad[d] + sol0_grad * line[d] + sol0 * line_grad[d];
        }
        return intersection;
      } else {
        const auto sol1_grad = solutions[1].second;
        auto intersection =
            PtTypeWithGradient(PtBase<ScalarType>(starting_pt + sol1 * line));
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          intersection.getData()[d] =
              starting_pt_grad[d] + sol1_grad * line[d] + sol1 * line_grad[d];
        }
        return intersection;
      }
    }
  }
}

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out, const ParaboloidBase<ScalarType>& a_paraboloid) {
  const auto& datum = a_paraboloid.getDatum();
  const auto& frame = a_paraboloid.getReferenceFrame();
  const auto& aligned_paraboloid = a_paraboloid.getAlignedParaboloid();

  out << "Datum: " << datum << '\n';
  out << "Frame: \n"
      << frame[0] << '\n'
      << frame[1] << '\n'
      << frame[2] << '\n';
  out << "Aligned Paraboloid: " << aligned_paraboloid << '\n';
  out << "is always above? " << a_paraboloid.isAlwaysAbove() << '\n';
  out << "is always below? " << a_paraboloid.isAlwaysBelow() << '\n';
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_
