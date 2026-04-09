!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module f_PUNeigh_class
  use f_SeparatorVariant_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_PUNeigh
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_PUNeigh

  type, public :: PUNeigh_type
    type(c_PUNeigh) :: c_object
  contains
    final :: PUNeigh_class_delete
  end type PUNeigh_type

  interface new
    module procedure PUNeigh_class_new
  end interface
  interface reserve
    module procedure PUNeigh_class_reserve
  end interface
  interface setSize
    module procedure PUNeigh_class_setSize
  end interface
  interface addMember
    module procedure PUNeigh_class_addMember
  end interface
  interface setMember
    module procedure PUNeigh_class_setMember
  end interface
  interface emptyNeighborhood
    module procedure PUNeigh_class_emptyNeighborhood
  end interface
  interface setCenterOfStencil
    module procedure PUNeigh_class_setCenterOfStencil
  end interface

  interface

    subroutine F_PUNeigh_new(this) &
      bind(C, name="c_PUNeigh_new")
      import
      implicit none
      type(c_PUNeigh) :: this
    end subroutine F_PUNeigh_new

    subroutine F_PUNeigh_delete(this) &
      bind(C, name="c_PUNeigh_delete")
      import
      implicit none
      type(c_PUNeigh) :: this
    end subroutine F_PUNeigh_delete

    subroutine F_PUNeigh_reserve(this, a_size) &
      bind(C, name="c_PUNeigh_reserve")
      import
      implicit none
      type(c_PUNeigh) :: this
      integer(C_INT) :: a_size
    end subroutine F_PUNeigh_reserve

    subroutine F_PUNeigh_setSize(this, a_size) &
      bind(C, name="c_PUNeigh_setSize")
      import
      implicit none
      type(c_PUNeigh) :: this
      integer(C_INT) :: a_size
    end subroutine F_PUNeigh_setSize

    subroutine F_PUNeigh_addMember(this, a_separator, a_centroid, a_weight) &
      bind(C, name="c_PUNeigh_addMember")
      import
      implicit none
      type(c_PUNeigh) :: this
      type(c_SeparatorVariant) :: a_separator
      real(C_DOUBLE), dimension(1:3) :: a_centroid
      real(C_DOUBLE) :: a_weight 
    end subroutine F_PUNeigh_addMember

    subroutine F_PUNeigh_setMember(this, a_index, a_separator, a_centroid, a_weight) &
      bind(C, name="c_PUNeigh_setMember")
      import
      implicit none
      type(c_PUNeigh) :: this
      integer(C_INT) :: a_index
      type(c_SeparatorVariant) :: a_separator
      real(C_DOUBLE), dimension(1:3) :: a_centroid
      real(C_DOUBLE) :: a_weight
    end subroutine F_PUNeigh_setMember

    subroutine F_PUNeigh_emptyNeighborhood(this) &
      bind(C, name="c_PUNeigh_emptyNeighborhood")
      import
      implicit none
      type(c_PUNeigh) :: this
    end subroutine F_PUNeigh_emptyNeighborhood

    subroutine F_PUNeigh_setCenterOfStencil(this, a_index) &
      bind(C, name="c_PUNeigh_setCenterOfStencil")
      import
      implicit none
      type(c_PUNeigh) :: this
      integer(C_INT) :: a_index
    end subroutine F_PUNeigh_setCenterOfStencil

  end interface


  contains

    subroutine PUNeigh_class_new(this)
      implicit none
      type(PUNeigh_type), intent(inout) :: this
      call F_PUNeigh_new(this%c_object)
    end subroutine PUNeigh_class_new

    impure elemental subroutine PUNeigh_class_delete(this)
      implicit none
      type(PUNeigh_type), intent(in) :: this
      call F_PUNeigh_delete(this%c_object)
    end subroutine PUNeigh_class_delete

    subroutine PUNeigh_class_reserve(this, a_size)
      implicit none
      type(PUNeigh_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_size
      call F_PUNeigh_reserve(this%c_object,a_size)
    end subroutine PUNeigh_class_reserve

    subroutine PUNeigh_class_setSize(this, a_size)
      implicit none
      type(PUNeigh_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_size
      call F_PUNeigh_setSize(this%c_object,a_size)
    end subroutine PUNeigh_class_setSize

    subroutine PUNeigh_class_addMember(this, a_separator, a_centroid, a_weight)
      implicit none
      type(PUNeigh_type), intent(in) :: this
      type(SeparatorVariant_type), intent(in) :: a_separator
      real(IRL_double), dimension(1:3), intent(in) :: a_centroid
      real(IRL_double), intent(in) :: a_weight
      call F_PUNeigh_addMember(this%c_object, a_separator%c_object, a_centroid, a_weight)
    end subroutine PUNeigh_class_addMember

    subroutine PUNeigh_class_setMember(this, a_index, a_separator, a_centroid, a_weight)
      implicit none
      type(PUNeigh_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      type(SeparatorVariant_type), intent(in) :: a_separator
      real(IRL_double), dimension(1:3), intent(in) :: a_centroid
      real(IRL_double), intent(in) :: a_weight
      call F_PUNeigh_setMember(this%c_object, a_index, a_separator%c_object, a_centroid, a_weight)
    end subroutine PUNeigh_class_setMember

    subroutine PUNeigh_class_emptyNeighborhood(this)
      implicit none
      type(PUNeigh_type), intent(in) :: this
      call F_PUNeigh_emptyNeighborhood(this%c_object)
    end subroutine PUNeigh_class_emptyNeighborhood

    subroutine PUNeigh_class_setCenterOfStencil(this, a_index)
      implicit none
      type(PUNeigh_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      call F_PUNeigh_setCenterOfStencil(this%c_object, a_index)
    end subroutine PUNeigh_class_setCenterOfStencil

end module f_PUNeigh_class
