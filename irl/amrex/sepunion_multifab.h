// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SEPUMULTIFAB_H_
#define IRL_SEPUMULTIFAB_H_

#include <AMReX_FabArray.H>

#include "irl/amrex/sepunion_arraybox.h"

namespace amrex {
using SepUnionMultiFab = FabArray<SepUnionArrayBox>;
}  // namespace amrex

#include "irl/amrex/sepunion_multifab.tpp"

#endif  // IRL_SEPUMULTIFAB_H_
