// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This class is heavily inspired from the FArrayBox class in AmReX. 

#ifndef IRL_SEPUARRAYBOX_H_
#define IRL_SEPUARRAYBOX_H_

#include <AMReX_Config.H>

#include <AMReX_BaseFab.H>
#include <AMReX_Box.H>
#include <AMReX_SPACE.H>

#include <iosfwd>

#include "irl/variant_reconstruction/separator_union.h"

namespace amrex {

class SepUnionArrayBox;

class SepUFABio {
 public:
  static void write_header(std::ostream& os, const SepUnionArrayBox& fab,
                           int nvar);
  static void read(std::istream& is, SepUnionArrayBox& fab);
};

/**
* \class amrex::SepUnionArrayBox
* \brief  A Fortran Array of IRL::SeparatorUnion

*  SepUnionArrayBox is derived from BaseFab<IRL::SeparatorUnion>.

*  The C pre-processor macro AMREX_SPACEDIM must be defined to use
*  this class.  The internal precision of FARRAYBOX objects is
*  set by defining either BL_USE_FLOAT or BL_USE_DOUBLE

*  This is NOT a polymorphic class.

*  This class does NOT provide a copy constructor or assignment operator.

*  This class is heavily inspired from the class amrex::IArrayBox
*/
class SepUnionArrayBox : public BaseFab<IRL::SeparatorUnion> {
 public:
  friend class SepUFABio;
  using BaseFab<IRL::SeparatorUnion>::BaseFab;

  //! Construct an invalid FAB with no memory.
  SepUnionArrayBox() noexcept = default;

  explicit SepUnionArrayBox(Arena* ar) noexcept;

  SepUnionArrayBox(const Box& b, int ncomp, Arena* ar);

  /**
   * \brief Construct an initial FAB with the data space allocated but
   * not initialized. ncomp is the number of components
   * (variables) at each data point in the Box.
   */
  explicit SepUnionArrayBox(const Box& b, int ncomp = 1, bool alloc = true,
                            bool shared = false, Arena* ar = nullptr);

  SepUnionArrayBox(const SepUnionArrayBox& rhs, MakeType make_type, int scomp,
                   int ncomp);

  explicit SepUnionArrayBox(Array4<IRL::SeparatorUnion> const& a) noexcept
      : BaseFab<IRL::SeparatorUnion>(a) {}

  explicit SepUnionArrayBox(Array4<IRL::SeparatorUnion> const& a,
                            IndexType t) noexcept
      : BaseFab<IRL::SeparatorUnion>(a, t) {}

  explicit SepUnionArrayBox(Array4<IRL::SeparatorUnion const> const& a) noexcept
      : BaseFab<IRL::SeparatorUnion>(a) {}

  explicit SepUnionArrayBox(Array4<IRL::SeparatorUnion const> const& a,
                            IndexType t) noexcept
      : BaseFab<IRL::SeparatorUnion>(a, t) {}

  //!  The destructor.
  ~SepUnionArrayBox() noexcept override = default;

  SepUnionArrayBox(SepUnionArrayBox&& rhs) noexcept = default;
  SepUnionArrayBox& operator=(SepUnionArrayBox&&) = default;

  SepUnionArrayBox(const SepUnionArrayBox&) = delete;
  SepUnionArrayBox& operator=(const SepUnionArrayBox&) = delete;

  //! Set the fab to the value r.
  template <RunOn run_on>
  SepUnionArrayBox& operator=(int v) noexcept;

  //! For debugging purposes we hide BaseFab version and do some extra work.
  void resize(const Box& b, int N = 1, Arena* ar = nullptr);

  void readFrom(std::istream& is);

  static void Initialize();
  static void Finalize();

  static SepUFABio const& getFABio();

  static std::string getClassName();

 private:
  static bool do_initval;
  static std::unique_ptr<SepUFABio> sufabio;
};

template <RunOn run_on>
SepUnionArrayBox& SepUnionArrayBox::operator=(int v) noexcept {
  BaseFab<IRL::SeparatorUnion>::operator= <run_on>(v);
  return *this;
}

}  // namespace amrex

#include "irl/amrex/sepunion_arraybox.tpp"

#endif  // IRL_SEPUARRAYBOX_H_
