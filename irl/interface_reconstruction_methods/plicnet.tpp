// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Andrew Cahaly <ajc428@cornell.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PLICNET_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PLICNET_TPP_

#include "irl/interface_reconstruction_methods/plicnet_weights_and_biases.h"

namespace IRL {

inline PLICNet::PLICNet(void)
    : stencil_m0_m(0.0),
      stencil_m1_m(0.0, 0.0, 0.0),
      dir_m(0),
      dir2_m(0),
      flipped_m(false),
      central_cell_defined_m(false),
      data_m{} {}

inline void PLICNet::clear(void) {
  stencil_m0_m = 0.0;
  stencil_m1_m = Pt(0.0, 0.0, 0.0);
  dir_m = 0;
  dir2_m = 0;
  flipped_m = false;
  central_cell_defined_m = false;
  data_m.fill(0.0);
}

inline void PLICNet::setMember(const Pt& a_lower_cell_pt,
                               const Pt& a_upper_cell_pt,
                               const double a_volume_fraction,
                               const Pt& a_centroid_0, const Pt& a_centroid_1,
                               const int a_i, const int a_j, const int a_k) {
  assert(a_i >= -1 && a_i <= 1);
  assert(a_j >= -1 && a_j <= 1);
  assert(a_k >= -1 && a_k <= 1);
  const std::array<double, 3> dx = {a_upper_cell_pt[0] - a_lower_cell_pt[0],
                                    a_upper_cell_pt[1] - a_lower_cell_pt[1],
                                    a_upper_cell_pt[2] - a_lower_cell_pt[2]};
  const Pt center = 0.5 * (a_lower_cell_pt + a_upper_cell_pt);
  data_m[7 * ((a_i + 1) * 9 + (a_j + 1) * 3 + (a_k + 1))] = a_volume_fraction;
  for (int n = 0; n < 3; ++n) {
    const double m1_0 = (a_centroid_0[n] - center[n]) / dx[n];
    const double m1_1 = (a_centroid_1[n] - center[n]) / dx[n];
    data_m[7 * ((a_i + 1) * 9 + (a_j + 1) * 3 + (a_k + 1)) + n + 1] =
        std::fabs(m1_0) > 0.5 ? 0.0 : m1_0;
    data_m[7 * ((a_i + 1) * 9 + (a_j + 1) * 3 + (a_k + 1)) + n + 4] =
        std::fabs(m1_1) > 0.5 ? 0.0 : m1_1;
  }
  if (a_i == 0 && a_j == 0 && a_k == 0) {
    cell_m =
        RectangularCuboid::fromBoundingPts(a_lower_cell_pt, a_upper_cell_pt);
    central_cell_defined_m = true;
  }
}

inline PlanarSeparator PLICNet::getPlanarSeparator(void) {
  if (central_cell_defined_m == false) {
    throw std::runtime_error(
        "PLICNet::getPlanarSeparator() called before central cell was "
        "defined.");
  }
  this->calculateNormal();
  const double vfrac = flipped_m ? 1.0 - data_m[7 * 13] : data_m[7 * 13];
  double distance = findDistanceOnePlane(cell_m, vfrac, normal_m);
  return PlanarSeparator::fromOnePlane(Plane(normal_m, distance));
}

inline void PLICNet::calculateNormal(void) {
  if (central_cell_defined_m == false) {
    throw std::runtime_error(
        "PLICNet::getNormal() called before central cell was defined.");
  }
  this->reflectMoments();
  std::array<double, 100> tmp1{};
  std::array<double, 100> tmp2{};
  std::array<double, 100> tmp3{};
  tmp1.fill(0.0);
  tmp2.fill(0.0);
  tmp3.fill(0.0);
  for (int j = 0; j < 100; ++j) {
    tmp1[j] = plicnet::lay1_bias[j];
    for (int i = 0; i < 189; ++i)
      tmp1[j] += data_m[i] * plicnet::lay1_weight[j][i];
    if (tmp1[j] < 0.0) tmp1[j] = 0.0;
  }

  for (int j = 0; j < 100; ++j) {
    tmp2[j] = plicnet::lay2_bias[j];
    for (int i = 0; i < 100; ++i)
      tmp2[j] += tmp1[i] * plicnet::lay2_weight[j][i];
    if (tmp2[j] < 0.0) tmp2[j] = 0.0;
  }

  for (int j = 0; j < 100; ++j) {
    tmp3[j] = plicnet::lay3_bias[j];
    for (int i = 0; i < 100; ++i)
      tmp3[j] += tmp2[i] * plicnet::lay3_weight[j][i];
    if (tmp3[j] < 0.0) tmp3[j] = 0.0;
  }

  for (int j = 0; j < 3; ++j) {
    normal_m[j] = plicnet::lay4_bias[j];
    for (int i = 0; i < 100; ++i)
      normal_m[j] += tmp3[i] * plicnet::lay4_weight[j][i];
  }

  double temp = 0.0;
  switch (dir2_m) {
    case 1:
      temp = normal_m[0];
      normal_m[0] = normal_m[1];
      normal_m[1] = temp;
      break;
    case 2:
      temp = normal_m[1];
      normal_m[1] = normal_m[2];
      normal_m[2] = temp;
      break;
    case 3:
      temp = normal_m[0];
      normal_m[0] = normal_m[2];
      normal_m[2] = temp;
      break;
    case 4:
      temp = normal_m[1];
      normal_m[1] = normal_m[2];
      normal_m[2] = temp;
      temp = normal_m[0];
      normal_m[0] = normal_m[1];
      normal_m[1] = temp;
      break;
    case 5:
      temp = normal_m[0];
      normal_m[0] = normal_m[2];
      normal_m[2] = temp;
      temp = normal_m[0];
      normal_m[0] = normal_m[1];
      normal_m[1] = temp;
      break;
  }

  switch (dir_m) {
    case 1:
      normal_m[0] = -normal_m[0];
      break;
    case 2:
      normal_m[1] = -normal_m[1];
      break;
    case 3:
      normal_m[2] = -normal_m[2];
      break;
    case 4:
      normal_m[0] = -normal_m[0];
      normal_m[1] = -normal_m[1];
      break;
    case 5:
      normal_m[0] = -normal_m[0];
      normal_m[2] = -normal_m[2];
      break;
    case 6:
      normal_m[1] = -normal_m[1];
      normal_m[2] = -normal_m[2];
      break;
    case 7:
      normal_m[0] = -normal_m[0];
      normal_m[1] = -normal_m[1];
      normal_m[2] = -normal_m[2];
      break;
  }

  if (!flipped_m) {
    normal_m[0] = -normal_m[0];
    normal_m[1] = -normal_m[1];
    normal_m[2] = -normal_m[2];
  }
  const std::array<double, 3> dx = {cell_m.calculateSideLength(0),
                                    cell_m.calculateSideLength(1),
                                    cell_m.calculateSideLength(2)};
  normal_m[0] = normal_m[0] * dx[0];
  normal_m[1] = normal_m[1] * dx[1];
  normal_m[2] = normal_m[2] * dx[2];
  normal_m.normalize();
}

inline void PLICNet::flip(void) {
  if (data_m[7 * 13] < 0.5) return;
  for (int n = 0; n < 27; ++n) {
    data_m[7 * n] = 1.0 - data_m[7 * n];
    const double m1x = data_m[7 * n + 1];
    const double m1y = data_m[7 * n + 2];
    const double m1z = data_m[7 * n + 3];
    data_m[7 * n + 1] = data_m[7 * n + 4];
    data_m[7 * n + 2] = data_m[7 * n + 5];
    data_m[7 * n + 3] = data_m[7 * n + 6];
    data_m[7 * n + 4] = m1x;
    data_m[7 * n + 5] = m1y;
    data_m[7 * n + 6] = m1z;
  }
  flipped_m = true;
}

inline void PLICNet::computeStencilMoments(void) {
  this->flip();
  stencil_m0_m = 0.0;
  stencil_m1_m = Pt(0.0, 0.0, 0.0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        const double volume_fraction = data_m[7 * (i * 9 + j * 3 + k)];
        const Pt centroid(data_m[7 * (i * 9 + j * 3 + k) + 1],
                          data_m[7 * (i * 9 + j * 3 + k) + 2],
                          data_m[7 * (i * 9 + j * 3 + k) + 3]);
        stencil_m0_m += volume_fraction;
        stencil_m1_m += volume_fraction * (centroid + Pt(i - 1, j - 1, k - 1));
      }
    }
  }
  stencil_m1_m *= 1.0 / safelyEpsilon(stencil_m0_m);
}

inline void PLICNet::reflectMomentsX(void) {
  double temp;
  for (int k = 0; k <= 2; ++k) {
    for (int j = 0; j <= 2; ++j) {
      for (int i = 0; i <= 2; ++i) {
        if (i == 0) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 1 || n == 4) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  -data_m[7 * (2 * 9 + j * 3 + k) + n];
              data_m[7 * (2 * 9 + j * 3 + k) + n] = -temp;
            } else {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (2 * 9 + j * 3 + k) + n];
              data_m[7 * (2 * 9 + j * 3 + k) + n] = temp;
            }
          }
        } else if (i == 1) {
          data_m[7 * (i * 9 + j * 3 + k) + 1] =
              -data_m[7 * (i * 9 + j * 3 + k) + 1];
          data_m[7 * (i * 9 + j * 3 + k) + 4] =
              -data_m[7 * (i * 9 + j * 3 + k) + 4];
        }
      }
    }
  }
}

inline void PLICNet::reflectMomentsY(void) {
  double temp;
  for (int k = 0; k <= 2; ++k) {
    for (int j = 0; j <= 2; ++j) {
      for (int i = 0; i <= 2; ++i) {
        if (j == 0) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 2 || n == 5) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  -data_m[7 * (i * 9 + 2 * 3 + k) + n];
              data_m[7 * (i * 9 + 2 * 3 + k) + n] = -temp;
            } else {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (i * 9 + 2 * 3 + k) + n];
              data_m[7 * (i * 9 + 2 * 3 + k) + n] = temp;
            }
          }
        } else if (j == 1) {
          data_m[7 * (i * 9 + j * 3 + k) + 2] =
              -data_m[7 * (i * 9 + j * 3 + k) + 2];
          data_m[7 * (i * 9 + j * 3 + k) + 5] =
              -data_m[7 * (i * 9 + j * 3 + k) + 5];
        }
      }
    }
  }
}

inline void PLICNet::reflectMomentsZ(void) {
  double temp;
  for (int k = 0; k <= 2; ++k) {
    for (int j = 0; j <= 2; ++j) {
      for (int i = 0; i <= 2; ++i) {
        if (k == 0) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 3 || n == 6) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  -data_m[7 * (i * 9 + j * 3 + 2) + n];
              data_m[7 * (i * 9 + j * 3 + 2) + n] = -temp;
            } else {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (i * 9 + j * 3 + 2) + n];
              data_m[7 * (i * 9 + j * 3 + 2) + n] = temp;
            }
          }
        } else if (k == 1) {
          data_m[7 * (i * 9 + j * 3 + k) + 3] =
              -data_m[7 * (i * 9 + j * 3 + k) + 3];
          data_m[7 * (i * 9 + j * 3 + k) + 6] =
              -data_m[7 * (i * 9 + j * 3 + k) + 6];
        }
      }
    }
  }
}

inline void PLICNet::reflectMomentsXY(void) {
  double temp;
  for (int k = 0; k <= 2; ++k) {
    for (int j = 0; j <= 2; ++j) {
      for (int i = 0; i <= 2; ++i) {
        if (i == j) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 1 || n == 4) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (i * 9 + j * 3 + k) + n + 1];
              data_m[7 * (i * 9 + j * 3 + k) + n + 1] = temp;
            }
          }
        } else if (i > j) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 1 || n == 4) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (j * 9 + i * 3 + k) + n + 1];
              data_m[7 * (j * 9 + i * 3 + k) + n + 1] = temp;
              temp = data_m[7 * (j * 9 + i * 3 + k) + n];
              data_m[7 * (j * 9 + i * 3 + k) + n] =
                  data_m[7 * (i * 9 + j * 3 + k) + n + 1];
              data_m[7 * (i * 9 + j * 3 + k) + n + 1] = temp;
            } else if (n == 0 || n == 3 || n == 6) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (j * 9 + i * 3 + k) + n];
              data_m[7 * (j * 9 + i * 3 + k) + n] = temp;
            }
          }
        }
      }
    }
  }
}

inline void PLICNet::reflectMomentsYZ(void) {
  double temp;
  for (int k = 0; k <= 2; ++k) {
    for (int j = 0; j <= 2; ++j) {
      for (int i = 0; i <= 2; ++i) {
        if (j == k) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 2 || n == 5) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (i * 9 + j * 3 + k) + n + 1];
              data_m[7 * (i * 9 + j * 3 + k) + n + 1] = temp;
            }
          }
        } else if (j > k) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 2 || n == 5) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (i * 9 + k * 3 + j) + n + 1];
              data_m[7 * (i * 9 + k * 3 + j) + n + 1] = temp;
              temp = data_m[7 * (i * 9 + k * 3 + j) + n];
              data_m[7 * (i * 9 + k * 3 + j) + n] =
                  data_m[7 * (i * 9 + j * 3 + k) + n + 1];
              data_m[7 * (i * 9 + j * 3 + k) + n + 1] = temp;
            } else if (n == 0 || n == 1 || n == 4) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (i * 9 + k * 3 + j) + n];
              data_m[7 * (i * 9 + k * 3 + j) + n] = temp;
            }
          }
        }
      }
    }
  }
}

inline void PLICNet::reflectMomentsXZ(void) {
  double temp;
  for (int k = 0; k <= 2; ++k) {
    for (int j = 0; j <= 2; ++j) {
      for (int i = 0; i <= 2; ++i) {
        if (i == k) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 1 || n == 4) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (i * 9 + j * 3 + k) + n + 2];
              data_m[7 * (i * 9 + j * 3 + k) + n + 2] = temp;
            }
          }
        } else if (i > k) {
          for (int n = 0; n <= 6; ++n) {
            if (n == 1 || n == 4) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (k * 9 + j * 3 + i) + n + 2];
              data_m[7 * (k * 9 + j * 3 + i) + n + 2] = temp;
              temp = data_m[7 * (k * 9 + j * 3 + i) + n];
              data_m[7 * (k * 9 + j * 3 + i) + n] =
                  data_m[7 * (i * 9 + j * 3 + k) + n + 2];
              data_m[7 * (i * 9 + j * 3 + k) + n + 2] = temp;
            } else if (n == 0 || n == 2 || n == 5) {
              temp = data_m[7 * (i * 9 + j * 3 + k) + n];
              data_m[7 * (i * 9 + j * 3 + k) + n] =
                  data_m[7 * (k * 9 + j * 3 + i) + n];
              data_m[7 * (k * 9 + j * 3 + i) + n] = temp;
            }
          }
        }
      }
    }
  }
}

inline void PLICNet::reflectMoments(void) {
  this->computeStencilMoments();
  Pt new_center = stencil_m1_m;
  double temp;

  if (std::abs(new_center[0]) <= plicnet::tol) new_center[0] = 0;
  if (std::abs(new_center[1]) <= plicnet::tol) new_center[1] = 0;
  if (std::abs(new_center[2]) <= plicnet::tol) new_center[2] = 0;

  if (new_center[0] < 0 && new_center[1] >= 0 && new_center[2] >= 0) {
    dir_m = 1;
    this->reflectMomentsX();
    new_center[0] = -new_center[0];
  } else if (new_center[0] >= 0 && new_center[1] < 0 && new_center[2] >= 0) {
    dir_m = 2;
    this->reflectMomentsY();
    new_center[1] = -new_center[1];
  } else if (new_center[0] >= 0 && new_center[1] >= 0 && new_center[2] < 0) {
    dir_m = 3;
    this->reflectMomentsZ();
    new_center[2] = -new_center[2];
  } else if (new_center[0] < 0 && new_center[1] < 0 && new_center[2] >= 0) {
    dir_m = 4;
    this->reflectMomentsX();
    this->reflectMomentsY();
    new_center[0] = -new_center[0];
    new_center[1] = -new_center[1];
  } else if (new_center[0] < 0 && new_center[1] >= 0 && new_center[2] < 0) {
    dir_m = 5;
    this->reflectMomentsX();
    this->reflectMomentsZ();
    new_center[0] = -new_center[0];
    new_center[2] = -new_center[2];
  } else if (new_center[0] >= 0 && new_center[1] < 0 && new_center[2] < 0) {
    dir_m = 6;
    this->reflectMomentsY();
    this->reflectMomentsZ();
    new_center[1] = -new_center[1];
    new_center[2] = -new_center[2];
  } else if (new_center[0] < 0 && new_center[1] < 0 && new_center[2] < 0) {
    dir_m = 7;
    this->reflectMomentsX();
    this->reflectMomentsY();
    this->reflectMomentsZ();
    new_center[0] = -new_center[0];
    new_center[1] = -new_center[1];
    new_center[2] = -new_center[2];
  }

  if (std::abs(new_center[0] - new_center[1]) <= plicnet::tol &&
      (new_center[0] - new_center[2]) > plicnet::tol) {
    dir2_m = 0;
  } else if (std::abs(new_center[1] - new_center[2]) <= plicnet::tol &&
             (new_center[0] - new_center[1]) > plicnet::tol) {
    dir2_m = 0;
  } else if (std::abs(new_center[0] - new_center[1]) <= plicnet::tol &&
             (new_center[2] - new_center[0]) > plicnet::tol) {
    dir2_m = 3;
    this->reflectMomentsXZ();
    temp = new_center[0];
    new_center[0] = new_center[2];
    new_center[2] = temp;
  } else if (std::abs(new_center[0] - new_center[2]) <= plicnet::tol &&
             (new_center[1] - new_center[0]) > plicnet::tol) {
    dir2_m = 1;
    this->reflectMomentsXY();
    temp = new_center[0];
    new_center[0] = new_center[1];
    new_center[1] = temp;
  } else if (std::abs(new_center[0] - new_center[2]) <= plicnet::tol &&
             (new_center[0] - new_center[1]) > plicnet::tol) {
    dir2_m = 2;
    this->reflectMomentsYZ();
    temp = new_center[1];
    new_center[1] = new_center[2];
    new_center[2] = temp;
  } else if (std::abs(new_center[1] - new_center[2]) <= plicnet::tol &&
             (new_center[1] - new_center[0]) > plicnet::tol) {
    dir2_m = 3;
    this->reflectMomentsXZ();
    temp = new_center[0];
    new_center[0] = new_center[2];
    new_center[2] = temp;
  } else if (new_center[1] > new_center[0] && new_center[0] >= new_center[2]) {
    dir2_m = 1;
    this->reflectMomentsXY();
    temp = new_center[0];
    new_center[0] = new_center[1];
    new_center[1] = temp;
  } else if (new_center[2] > new_center[1] && new_center[0] >= new_center[2]) {
    dir2_m = 2;
    this->reflectMomentsYZ();
    temp = new_center[1];
    new_center[1] = new_center[2];
    new_center[2] = temp;
  } else if (new_center[2] > new_center[1] && new_center[1] >= new_center[0]) {
    dir2_m = 3;
    this->reflectMomentsXZ();
    temp = new_center[0];
    new_center[0] = new_center[2];
    new_center[2] = temp;
  } else if (new_center[1] > new_center[0]) {
    dir2_m = 4;
    this->reflectMomentsXY();
    this->reflectMomentsYZ();
    temp = new_center[0];
    new_center[0] = new_center[1];
    new_center[1] = temp;
    temp = new_center[1];
    new_center[1] = new_center[2];
    new_center[2] = temp;
  } else if (new_center[2] > new_center[1]) {
    dir2_m = 5;
    this->reflectMomentsXY();
    this->reflectMomentsXZ();
    temp = new_center[0];
    new_center[0] = new_center[1];
    new_center[1] = temp;
    temp = new_center[0];
    new_center[0] = new_center[2];
    new_center[2] = temp;
  }
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PLICNET_TPP_
