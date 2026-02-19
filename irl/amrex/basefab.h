// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_AMREX_BASEFAB_H_
#define IRL_AMREX_BASEFAB_H_

#include <AMReX_BaseFab.H>

#include "irl/variant_reconstruction/separator_union.h"

namespace IRL {

class BaseFabSepUnion : public amrex::BaseFab<SeparatorUnion> {
 public:
  using amrex::BaseFab<SeparatorUnion>::BaseFab;

 private:
  static bool do_initval;
  static std::unique_ptr<int> svfabio;
};

}  // namespace IRL

#include "irl/amrex/basefab.tpp"

#endif  // IRL_AMREX_BASEFAB_H_
