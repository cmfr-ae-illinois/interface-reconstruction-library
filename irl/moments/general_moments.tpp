// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2023 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_GENERAL_MOMENTS_TPP_
#define IRL_MOMENTS_GENERAL_MOMENTS_TPP_

#include <cassert>
#include <iostream>
#include <vector>

#include "irl/geometry/general/pt.h"

namespace IRL {

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>::GeneralMomentsBase(void)
    : moments_m() {
  static_assert(DIM == 2 || DIM == 3,
                "Dimension argument of GeneralMoments template must be 2 or 3");
}
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>::GeneralMomentsBase(
    const GeneralMomentsBase<ORDER, DIM, double>& a_moments) {
  for (UnsignedIndex_t i = 0; i < a_moments.size(); ++i) {
    moments_m[i] = ScalarType(a_moments[i]);
  }
  static_assert(DIM == 2 || DIM == 3,
                "Dimension argument of GeneralMoments template must be 2 or 3");
}
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>::fromScalarConstant(
    const ScalarType a_value) {
  GeneralMomentsBase<ORDER, DIM, ScalarType> mom;
  std::fill(mom.moments_m.begin(), mom.moments_m.end(), a_value);
  return mom;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
template <class GeometryType>
GeneralMomentsBase<ORDER, DIM, ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>::calculateMoments(
    GeometryType* a_geometry) {
  return GeneralMomentsBase<ORDER, DIM, ScalarType>(
      a_geometry->template calculateGeneralMoments<ORDER>());
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
typename GeneralMomentsBase<ORDER, DIM, ScalarType>::storage
GeneralMomentsBase<ORDER, DIM, ScalarType>::calculateCentralMoments(
    void) const {
  // Implementation of factorial requires this assertion
  static_assert(ORDER < 10,
                "Calculation of CentralMoments supported only up to order 10.");

  storage central;
  central[0] = (*this)[0];
  std::fill(central.begin() + 1, central.begin() + DIM + 1, ScalarType(0));

  auto fact = [](const int a_number) {
    return a_number < 2 ? ScalarType(1) : [a_number]() {
      int sum = 1;
      for (int i = 2; i < a_number + 1; ++i) {
        sum *= i;
      }
      return ScalarType(sum);
    }();
  };

  auto binomial_coefficient = [fact](const int a_start, const int a_end) {
    return fact(a_end) / (fact(a_start) * fact(a_end - a_start));
  };

  if constexpr (DIM == 2) {
    auto get_linear_index = [fact](const int i, const int j) {
      // TODO: There should be some kind of formula to find this much more
      // quickly.
      const int order = i + j;
      int m = (order) * (order + 1) / 3;
      for (int ii = order; ii >= 0; --ii) {
        const auto jj = order - ii;
        if (ii == i && jj == j) return m;
      }
      assert(false);
      return -1;
    };

    PtBase<ScalarType> c =
        PtBase<ScalarType>((*this)[1], (*this)[2], ScalarType(0)) / (*this)[0];

    // TODO: There is definitely a better way to build
    // this that builds on previous results
    for (int corder = 2, m = 1 + DIM; corder <= ORDER; ++corder) {
      for (int i = corder; i >= 0; --i, ++m) {
        const auto j = corder - i;
        central[m] = ScalarType(0);
        for (int ii = 0; ii < i + 1; ++ii) {
          for (int jj = 0; jj < j + 1; ++jj) {
            central[m] += binomial_coefficient(jj, j) *
                          binomial_coefficient(ii, i) * pow(-c[0], i - ii) *
                          pow(-c[1], j - jj) *
                          (*this)[get_linear_index(ii, jj)];
          }
        }
      }
    }
  } else {
    auto get_linear_index = [fact](const int i, const int j, const int k) {
      // TODO: There should be some kind of formula to find this much more
      // quickly.
      const int order = i + j + k;
      int m = (order) * (order + 1) * (order + 2) / 6;
      for (int ii = order; ii >= 0; --ii) {
        for (int jj = order - ii; jj >= 0; --jj, ++m) {
          const auto kk = order - ii - jj;
          if (ii == i && jj == j && kk == k) return m;
        }
      }
      assert(false);
      return -1;
    };

    // TODO: There is definitely a better way to build
    // this that builds on previous results
    PtBase<ScalarType> c =
        PtBase<ScalarType>((*this)[1], (*this)[2], (*this)[3]) / (*this)[0];
    for (int corder = 2, m = 1 + DIM; corder <= ORDER; ++corder) {
      for (int i = corder; i >= 0; --i) {
        for (int j = corder - i; j >= 0; --j, ++m) {
          const auto k = corder - i - j;

          central[m] = ScalarType(0);
          for (int ii = 0; ii < i + 1; ++ii) {
            for (int jj = 0; jj < j + 1; ++jj) {
              for (int kk = 0; kk < k + 1; ++kk) {
                central[m] += binomial_coefficient(kk, k) *
                              binomial_coefficient(jj, j) *
                              binomial_coefficient(ii, i) * pow(-c[0], i - ii) *
                              pow(-c[1], j - jj) * pow(-c[2], k - kk) *
                              (*this)[get_linear_index(ii, jj, kk)];
              }
            }
          }
        }
      }
    }
  }

  return central;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
typename GeneralMomentsBase<ORDER, DIM, ScalarType>::storage
GeneralMomentsBase<ORDER, DIM, ScalarType>::calculateInvariantMoments(
    void) const {
  auto central = this->calculateCentralMoments();
  storage invariants;
  if constexpr (DIM == 2) {
    std::copy(central.begin(), central.begin() + 1 + DIM, invariants.begin());
    for (int corder = 2, m = 1 + DIM; corder <= ORDER; ++corder) {
      for (int i = corder; i >= 0; --i, ++m) {
        const auto j = corder - i;
        invariants[m] =
            central[m] /
            std::pow(central[0],
                     ScalarType(1) + ScalarType(0.5) * ScalarType(i + j));
      }
    }
  } else {
    std::copy(central.begin(), central.begin() + 1 + DIM, invariants.begin());
    for (int corder = 2, m = 1 + DIM; corder <= ORDER; ++corder) {
      for (int i = corder; i >= 0; --i) {
        for (int j = corder - i; j >= 0; --j, ++m) {
          const auto k = corder - i - j;
          invariants[m] =
              central[m] /
              std::pow(central[0],
                       ScalarType(1) + ScalarType(i + j + k) / ScalarType(3));
        }
      }
    }
  }
  return invariants;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
void GeneralMomentsBase<ORDER, DIM, ScalarType>::normalizeAsInvariant(void) {
  for (std::size_t m = 1; m < 1 + DIM; ++m) {
    (*this)[m] = (*this)[m] / (*this)[0];
  }
  if constexpr (DIM == 2) {
    for (int corder = 2, m = 1 + DIM; corder <= ORDER; ++corder) {
      for (int i = corder; i >= 0; --i, ++m) {
        const auto j = corder - i;
        (*this)[m] =
            (*this)[m] /
            std::pow((*this)[0],
                     ScalarType(1) + ScalarType(0.5) * ScalarType(i + j));
      }
    }
  } else {
    for (int corder = 2, m = 1 + DIM; corder <= ORDER; ++corder) {
      for (int i = corder; i >= 0; --i) {
        for (int j = corder - i; j >= 0; --j, ++m) {
          const auto k = corder - i - j;
          (*this)[m] =
              (*this)[m] /
              std::pow((*this)[0],
                       ScalarType(1) + ScalarType(i + j + k) / ScalarType(3));
        }
      }
    }
  }
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
ScalarType& GeneralMomentsBase<ORDER, DIM, ScalarType>::volume(void) {
  return moments_m[0];
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
const ScalarType GeneralMomentsBase<ORDER, DIM, ScalarType>::volume(
    void) const {
  return moments_m[0];
}

// template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
// const VolumeBase<ScalarType> GeneralMomentsBase<ORDER, DIM,
// ScalarType>::volume(
//     void) const {
//   return moments_m[0];
// }

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
void GeneralMomentsBase<ORDER, DIM, ScalarType>::multiplyByVolume(void) {
  for (std::size_t i = 1; i < moments_m.size(); ++i) {
    moments_m[i] *= moments_m[0];
  }
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
void GeneralMomentsBase<ORDER, DIM, ScalarType>::normalizeByVolume(void) {
  if (moments_m[0] != ScalarType(0)) {
    for (std::size_t i = 1; i < moments_m.size(); ++i) {
      moments_m[i] /= moments_m[0];
    }
  }
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>&
GeneralMomentsBase<ORDER, DIM, ScalarType>::operator+=(
    const GeneralMomentsBase<ORDER, DIM, ScalarType>& a_rhs) {
  for (std::size_t i = 0; i < moments_m.size(); ++i) {
    moments_m[i] += a_rhs[i];
  }
  return *this;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>&
GeneralMomentsBase<ORDER, DIM, ScalarType>::operator*=(const ScalarType a_rhs) {
  for (std::size_t i = 0; i < moments_m.size(); ++i) {
    moments_m[i] *= a_rhs;
  }
  return *this;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>&
GeneralMomentsBase<ORDER, DIM, ScalarType>::operator/=(const ScalarType a_rhs) {
  for (std::size_t i = 0; i < moments_m.size(); ++i) {
    moments_m[i] /= a_rhs;
  }
  return *this;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>&
GeneralMomentsBase<ORDER, DIM, ScalarType>::operator=(
    const ScalarType a_value) {
  std::fill(moments_m.begin(), moments_m.end(), ScalarType(a_value));
  return *this;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
GeneralMomentsBase<ORDER, DIM, ScalarType>&
GeneralMomentsBase<ORDER, DIM, ScalarType>::operator=(
    const GeneralMomentsBase<ORDER, DIM, double>& a_rhs) {
  for (UnsignedIndex_t i = 0; i < a_rhs.size(); ++i) {
    moments_m[i] = ScalarType(a_rhs[i]);
  }
  return *this;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
UnsignedIndex_t GeneralMomentsBase<ORDER, DIM, ScalarType>::size(void) const {
  return moments_m.size();
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
ScalarType& GeneralMomentsBase<ORDER, DIM, ScalarType>::operator[](
    const UnsignedIndex_t a_index) {
  assert(a_index < this->size());
  return moments_m[a_index];
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
ScalarType GeneralMomentsBase<ORDER, DIM, ScalarType>::operator[](
    const UnsignedIndex_t a_index) const {
  assert(a_index < this->size());
  return moments_m[a_index];
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
std::ostream& operator<<(
    std::ostream& out,
    const GeneralMomentsBase<ORDER, DIM, ScalarType>& a_moments) {
  static constexpr UnsignedIndex_t PER_LINE = 10;
  out << "There are " << a_moments.size() << " total moments.\n";
  UnsignedIndex_t index = 0;

  static_assert(DIM == 3, "DIM == 2 not implemented yet.");

  std::vector<std::string> terms(a_moments.size());
  terms[0] = "Volume";
  std::size_t max_size = terms[0].size();
  if constexpr (DIM == 3) {
    for (int corder = 1, m = 1; corder <= ORDER; ++corder) {
      for (int i = corder; i >= 0; --i) {
        for (int j = corder - i; j >= 0; --j, ++m) {
          const auto k = corder - i - j;
          std::string term = "";

          if (i > 0) {
            term += "x^" + std::to_string(i);
          }
          if (j > 0) {
            term += (term != "") ? " * y^" + std::to_string(j)
                                 : "y^" + std::to_string(j);
          }
          if (k > 0) {
            term += (term != "") ? " * z^" + std::to_string(k)
                                 : "z^" + std::to_string(k);
          }

          terms[m] = term;
          max_size = std::max(max_size, terms[m].size());
        }
      }
    }
  }

  for (UnsignedIndex_t i = 0; i < a_moments.size(); ++i) {
    terms[i].resize(max_size, ' ');
    out << terms[i] << " " << a_moments[i] << '\n';
  }
  return out;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
inline GeneralMomentsBase<ORDER, DIM, ScalarType> operator+(
    const GeneralMomentsBase<ORDER, DIM, ScalarType>& a_vm1,
    const GeneralMomentsBase<ORDER, DIM, ScalarType>& a_vm2) {
  GeneralMomentsBase<ORDER, DIM, ScalarType> mom;
  for (std::size_t i = 0; i < a_vm1.size(); ++i) {
    mom[i] = a_vm1[i] + a_vm2[i];
  }
  return mom;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
inline GeneralMomentsBase<ORDER, DIM, ScalarType> operator-(
    const GeneralMomentsBase<ORDER, DIM, ScalarType>& a_vm1,
    const GeneralMomentsBase<ORDER, DIM, ScalarType>& a_vm2) {
  GeneralMomentsBase<ORDER, DIM, ScalarType> mom;
  for (std::size_t i = 0; i < a_vm1.size(); ++i) {
    mom[i] = a_vm1[i] - a_vm2[i];
  }
  return mom;
}

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
inline GeneralMomentsBase<ORDER, DIM, ScalarType> operator*(
    const ScalarType a_multiplier,
    const GeneralMomentsBase<ORDER, DIM, ScalarType>& a_vm) {
  GeneralMomentsBase<ORDER, DIM, ScalarType> mom;
  for (std::size_t i = 0; i < a_vm.size(); ++i) {
    mom[i] = a_multiplier * a_vm[i];
  }
  return mom;
}
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
inline GeneralMomentsBase<ORDER, DIM, ScalarType> operator*(
    const GeneralMomentsBase<ORDER, DIM, ScalarType>& a_vm,
    const ScalarType a_multiplier) {
  return a_multiplier * a_vm;
}

}  // namespace IRL

#endif  // IRL_MOMENTS_GENERAL_MOMENTS_TPP_
