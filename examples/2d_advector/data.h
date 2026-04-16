// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_DATA_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_DATA_H_

#include <algorithm>
#include <cstring>

#include "examples/2d_advector/basic_mesh.h"
#include "examples/2d_advector/irl2d.h"

/// \brief A basic multi-dimensional data container.
template <class ContainedType>
class Data {
 public:
  using value_t = ContainedType;

  /// \brief Default constructor
  Data(void) = default;

  /// \brief Construct and set initial size.
  explicit Data(const BasicMesh* a_mesh_ptr) : mesh_m(a_mesh_ptr) {
    data_m = new ContainedType[this->getMesh().size()];
  }

  // Copy constructor
  Data(const Data& other) {
    delete[] this->data_m;
    const BasicMesh& mesh = other.getMesh();
    mesh_m = other.mesh_m;
    data_m = new ContainedType[mesh.size()];
    std::memcpy(data_m, other.data_m, sizeof(ContainedType) * mesh.size());
  }

  // Move constructor
  Data(Data&& other) {
    delete[] data_m;
    mesh_m = other.mesh_m;
    data_m = other.data_m;
    other.mesh_m = nullptr;
    other.data_m = nullptr;
  }

  // Copy assignment
  Data& operator=(const Data& other) {
    if (this != &other) {
      delete[] data_m;
      const BasicMesh& mesh = other.getMesh();
      mesh_m = other.mesh_m;
      data_m = new ContainedType[mesh.size()];
      std::memcpy(data_m, other.data_m, sizeof(ContainedType) * mesh.size());
    }
    return (*this);
  }

  // Move assignment
  Data& operator=(Data&& other) {
    if (this != &other) {
      delete[] data_m;
      mesh_m = other.mesh_m;
      data_m = other.data_m;
      other.mesh_m = nullptr;
      other.data_m = nullptr;
    }
    return (*this);
  }

  /// \brief Return the pointer to the mesh.
  const BasicMesh& getMesh(void) { return *mesh_m; }
  const BasicMesh& getMesh(void) const { return *mesh_m; }

  /// \brief Provide access to data.
  ContainedType& operator()(const int i, const int j) {
    return data_m[this->calculateIndex(i, j)];
  }

  /// \brief Provide const access to data.
  const ContainedType& operator()(const int i, const int j) const {
    return data_m[this->calculateIndex(i, j)];
  }

  /// \brief Update the border for periodicity.
  void updateBorder(void);
  void updateLowerX(void);
  void updateUpperX(void);
  void updateLowerY(void);
  void updateUpperY(void);

  /// \brief Destructor to delete memory allocated during construction.
  ~Data(void) { delete[] data_m; }

 private:
  /// \brief Calculate the index, where the first real (non-ghost) cell is at 0
  /// and the fastest changing indices are k, j, i.
  std::size_t calculateIndex(const int i, const int j) {
    assert(i >= this->getMesh().imino() && i <= this->getMesh().imaxo());
    assert(j >= this->getMesh().jmino() && j <= this->getMesh().jmaxo());
    assert(static_cast<std::size_t>(
               (i + this->getMesh().getNgc()) * this->getMesh().getNyo() +
               (j + this->getMesh().getNgc())) < this->getMesh().size());
    return static_cast<std::size_t>((i + this->getMesh().getNgc()) *
                                        this->getMesh().getNyo() +
                                    (j + this->getMesh().getNgc()));
  }

  /// \brief Const getInd
  std::size_t calculateIndex(const int i, const int j) const {
    assert(i >= this->getMesh().imino() && i <= this->getMesh().imaxo());
    assert(j >= this->getMesh().jmino() && j <= this->getMesh().jmaxo());
    assert(static_cast<std::size_t>(
               (i + this->getMesh().getNgc()) * this->getMesh().getNyo() +
               (j + this->getMesh().getNgc())) < this->getMesh().size());
    return static_cast<std::size_t>((i + this->getMesh().getNgc()) *
                                        this->getMesh().getNyo() +
                                    (j + this->getMesh().getNgc()));
  }

  ContainedType* data_m = nullptr;
  const BasicMesh* mesh_m = nullptr;
};

template <class ContainedType>
void Data<ContainedType>::updateBorder(void) {
  this->updateLowerX();
  this->updateUpperX();
  this->updateLowerY();
  this->updateUpperY();
}

template <class ContainedType>
void Data<ContainedType>::updateLowerX(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      (*this)(i, j) = (*this)(mesh.imax() + 1 + i, j);
    }
  }
}
template <class ContainedType>
void Data<ContainedType>::updateUpperX(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      (*this)(i, j) = (*this)(i - mesh.getNx(), j);
    }
  }
}

template <class ContainedType>
void Data<ContainedType>::updateLowerY(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      (*this)(i, j) = (*this)(i, mesh.jmax() + 1 + j);
    }
  }
}
template <class ContainedType>
void Data<ContainedType>::updateUpperY(void) {
  const BasicMesh& mesh = this->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      (*this)(i, j) = (*this)(i, j - mesh.getNy());
    }
  }
}

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_DATA_H_
