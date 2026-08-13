// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This class is heavily inspired from the FArrayBox class in AmReX.

#ifndef IRL_SEPUARRAYBOX_TPP_
#define IRL_SEPUARRAYBOX_TPP_

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_Utility.H>
#include <AMReX_VectorIO.H>

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>

namespace amrex {

#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
bool SepUnionArrayBox::do_initval = true;
#else
bool SepUnionArrayBox::do_initval = false;
#endif

std::unique_ptr<SepUFABio> SepUnionArrayBox::sufabio;

namespace {
bool initialized = false;
}

void SepUnionArrayBox::Initialize() {
  if (initialized) {
    return;
  }

  sufabio = std::make_unique<SepUFABio>();

  amrex::ExecOnFinalize(SepUnionArrayBox::Finalize);
  initialized = true;
}

void SepUnionArrayBox::Finalize() {
  sufabio.reset();
  initialized = false;
}

SepUnionArrayBox::SepUnionArrayBox(Arena* ar) noexcept
    : BaseFab<IRL::SeparatorUnion>(ar) {}

SepUnionArrayBox::SepUnionArrayBox(const Box& b, int n, Arena* ar)
    : BaseFab<IRL::SeparatorUnion>(b, n, ar) {
#ifndef AMREX_USE_GPU
  // For debugging purposes
  if (do_initval) {
    setVal<RunOn::Host>(IRL::SeparatorUnion(
        IRL::Plane(IRL::Normal(std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max()),
                   std::numeric_limits<double>::max())));
  }
#endif
}

SepUnionArrayBox::SepUnionArrayBox(const Box& b, int n, bool alloc, bool shared,
                                   Arena* ar)
    : BaseFab<IRL::SeparatorUnion>(b, n, alloc, shared, ar) {
#ifndef AMREX_USE_GPU
  // For debugging purposes
  if (alloc && do_initval) {
    setVal<RunOn::Host>(IRL::SeparatorUnion(
        IRL::Plane(IRL::Normal(std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max()),
                   std::numeric_limits<double>::max())));
  }
#endif
}

SepUnionArrayBox::SepUnionArrayBox(const SepUnionArrayBox& rhs,
                                   MakeType make_type, int scomp, int ncomp)
    : BaseFab<IRL::SeparatorUnion>(rhs, make_type, scomp, ncomp) {}

void SepUnionArrayBox::resize(const Box& b, int N, Arena* ar) {
  BaseFab<IRL::SeparatorUnion>::resize(b, N, ar);
  // For debugging purposes
  if (do_initval) {
#if defined(AMREX_USE_GPU)
    if (Gpu::inLaunchRegion() && arena()->isDeviceAccessible()) {
      setVal<RunOn::Device>(IRL::SeparatorUnion(
          IRL::Plane(IRL::Normal(std::numeric_limits<double>::max(),
                                 std::numeric_limits<double>::max(),
                                 std::numeric_limits<double>::max()),
                     std::numeric_limits<double>::max())));
      Gpu::streamSynchronize();
    } else
#endif
    {
      setVal<RunOn::Host>(IRL::SeparatorUnion(
          IRL::Plane(IRL::Normal(std::numeric_limits<double>::max(),
                                 std::numeric_limits<double>::max(),
                                 std::numeric_limits<double>::max()),
                     std::numeric_limits<double>::max())));
    }
  }
}

std::string SepUnionArrayBox::getClassName() {
  return std::string("amrex::SepUnionArrayBox");
}

SepUFABio const& SepUnionArrayBox::getFABio() { return *sufabio; }

void SepUnionArrayBox::readFrom(std::istream& is) {
  std::string type;
  is >> type;
  if (type != "SepUFAB") {
    amrex::Error(std::string("SepUnionArrayBox::readFrom: SepUFAB is expected, "
                             "but instead we have ") +
                 type);
  }

  IntDescriptor data_descriptor;
  is >> data_descriptor;

  Box tmp_box;
  int tmp_ncomp;
  is >> tmp_box;
  is >> tmp_ncomp;
  AMREX_ASSERT(tmp_ncomp >= 0 && tmp_ncomp < std::numeric_limits<int>::max());
  is.ignore(99999, '\n');

  if (this->box() != tmp_box || this->nComp() != tmp_ncomp) {
    this->resize(tmp_box, tmp_ncomp);
  }

#ifdef AMREX_USE_GPU
  if (this->arena()->isManaged() || this->arena()->isDevice()) {
    SepUnionArrayBox hostfab(this->box(), this->nComp(), The_Pinned_Arena());
    sufabio->read(is, hostfab);
    Gpu::htod_memcpy_async(
        this->dataPtr(), hostfab.dataPtr(),
        hostfab.size() * sizeof(SepUnionArrayBox::value_type));
    Gpu::streamSynchronize();
  } else
#endif
  {
    sufabio->read(is, *this);
  }
}

void SepUFABio::write_header(std::ostream& os, const SepUnionArrayBox& fab,
                             int nvar) {
  AMREX_ASSERT(nvar <= fab.nComp());
  os << "SepUFAB " << FPC::NativeIntDescriptor();
  os << fab.box() << ' ' << nvar << '\n';
}

void readSepUnionData(IRL::SeparatorUnion* data, std::size_t size,
                      std::istream& is) {
  IRL::SeparatorUnion value;
  for (std::size_t j = 0; j < size; ++j) {
    is.read((char*)&value, sizeof(IRL::SeparatorUnion));
    data[j] = static_cast<IRL::SeparatorUnion>(value);
  }
}

void SepUFABio::read(std::istream& is, SepUnionArrayBox& fab) {
  readSepUnionData(fab.dataPtr(), fab.size(), is);
}

}  // namespace amrex

#endif  // IRL_SEPUARRAYBOX_TPP_
