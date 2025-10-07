!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module f_JibbenNeigh_class
  use f_Poly_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_JibbenNeigh
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_JibbenNeigh

  type, public :: JibbenNeigh_type
    type(c_JibbenNeigh) :: c_object
  contains
    final :: JibbenNeigh_class_delete
  end type JibbenNeigh_type

  interface new
    module procedure JibbenNeigh_class_new
  end interface
  interface reserve
    module procedure JibbenNeigh_class_reserve
  end interface
  interface setSize
    module procedure JibbenNeigh_class_setSize
  end interface
  interface localize
    module procedure JibbenNeigh_class_localize
  end interface
  interface setDelta
    module procedure JibbenNeigh_class_setDelta
  end interface
  interface addMember
    module procedure JibbenNeigh_class_addMember
  end interface
  interface emptyNeighborhood
    module procedure JibbenNeigh_class_emptyNeighborhood
  end interface
  interface setCenterOfStencil
    module procedure JibbenNeigh_class_setCenterOfStencil
  end interface

  interface

    subroutine F_JibbenNeigh_new(this) &
      bind(C, name="c_JibbenNeigh_new")
      import
      implicit none
      type(c_JibbenNeigh) :: this
    end subroutine F_JibbenNeigh_new

    subroutine F_JibbenNeigh_delete(this) &
      bind(C, name="c_JibbenNeigh_delete")
      import
      implicit none
      type(c_JibbenNeigh) :: this
    end subroutine F_JibbenNeigh_delete

    subroutine F_JibbenNeigh_reserve(this, a_size) &
      bind(C, name="c_JibbenNeigh_reserve")
      import
      implicit none
      type(c_JibbenNeigh) :: this
      integer(C_INT) :: a_size
    end subroutine F_JibbenNeigh_reserve

    subroutine F_JibbenNeigh_setSize(this, a_size) &
      bind(C, name="c_JibbenNeigh_setSize")
      import
      implicit none
      type(c_JibbenNeigh) :: this
      integer(C_INT) :: a_size
    end subroutine F_JibbenNeigh_setSize

    subroutine F_JibbenNeigh_localize(this) &
      bind(C, name="c_JibbenNeigh_localize")
      import
      implicit none
      type(c_JibbenNeigh) :: this
    end subroutine F_JibbenNeigh_localize

    subroutine F_JibbenNeigh_setDelta(this, a_delta) &
      bind(C, name="c_JibbenNeigh_setDelta")
      import
      implicit none
      type(c_JibbenNeigh) :: this
      real(C_DOUBLE) :: a_delta
    end subroutine F_JibbenNeigh_setDelta

    subroutine F_JibbenNeigh_addMember(this, a_polygon, a_weight) &
      bind(C, name="c_JibbenNeigh_addMember")
      import
      implicit none
      type(c_JibbenNeigh) :: this
      type(c_Poly) :: a_polygon ! Pointer to Poly
      real(C_DOUBLE) :: a_weight ! Pointer to double with weight
    end subroutine F_JibbenNeigh_addMember

    subroutine F_JibbenNeigh_emptyNeighborhood(this) &
      bind(C, name="c_JibbenNeigh_emptyNeighborhood")
      import
      implicit none
      type(c_JibbenNeigh) :: this
    end subroutine F_JibbenNeigh_emptyNeighborhood

    subroutine F_JibbenNeigh_setCenterOfStencil(this, a_center_cell_index) &
      bind(C, name="c_JibbenNeigh_setCenterOfStencil")
      import
      implicit none
      type(c_JibbenNeigh) :: this
      integer(C_INT) :: a_center_cell_index
    end subroutine F_JibbenNeigh_setCenterOfStencil

  end interface


  contains

    subroutine JibbenNeigh_class_new(this)
      implicit none
      type(JibbenNeigh_type), intent(inout) :: this
      call F_JibbenNeigh_new(this%c_object)
    end subroutine JibbenNeigh_class_new

    impure elemental subroutine JibbenNeigh_class_delete(this)
      implicit none
      type(JibbenNeigh_type), intent(in) :: this
      call F_JibbenNeigh_delete(this%c_object)
    end subroutine JibbenNeigh_class_delete

    subroutine JibbenNeigh_class_reserve(this, a_size)
      implicit none
      type(JibbenNeigh_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_size
      call F_JibbenNeigh_reserve(this%c_object,a_size)
    end subroutine JibbenNeigh_class_reserve

    subroutine JibbenNeigh_class_setSize(this, a_size)
      implicit none
      type(JibbenNeigh_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_size
      call F_JibbenNeigh_setSize(this%c_object,a_size)
    end subroutine JibbenNeigh_class_setSize

    subroutine JibbenNeigh_class_localize(this)
      implicit none
      type(JibbenNeigh_type), intent(in) :: this
      call F_JibbenNeigh_localize(this%c_object)
    end subroutine JibbenNeigh_class_localize

    subroutine JibbenNeigh_class_setDelta(this, a_delta)
      implicit none
      type(JibbenNeigh_type), intent(in) :: this
      real(IRL_double), intent(in) :: a_delta
      call F_JibbenNeigh_setDelta(this%c_object,a_delta)
    end subroutine JibbenNeigh_class_setDelta

    subroutine JibbenNeigh_class_addMember(this, a_polygon, a_weight)
      implicit none
      type(JibbenNeigh_type), intent(in) :: this
      type(Poly_type), intent(in) :: a_polygon
      real(IRL_double), intent(in) :: a_weight
      call F_JibbenNeigh_addMember(this%c_object, a_polygon%c_object, a_weight)
    end subroutine JibbenNeigh_class_addMember

    subroutine JibbenNeigh_class_emptyNeighborhood(this)
      implicit none
      type(JibbenNeigh_type), intent(in) :: this
      call F_JibbenNeigh_emptyNeighborhood(this%c_object)
    end subroutine JibbenNeigh_class_emptyNeighborhood

    subroutine JibbenNeigh_class_setCenterOfStencil(this, a_center_cell_index)
      implicit none
      type(JibbenNeigh_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_center_cell_index
      call F_JibbenNeigh_setCenterOfStencil(this%c_object, a_center_cell_index)
    end subroutine JibbenNeigh_class_setCenterOfStencil

end module f_JibbenNeigh_class
