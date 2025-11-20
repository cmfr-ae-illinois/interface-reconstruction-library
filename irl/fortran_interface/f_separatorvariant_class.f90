!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_SeparatorVariant_class.f90
!!
!! This file allows use of the IRL SeparatorVariant
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's SeparatorVariant class along with enabling
!! some of its methods.
module f_SeparatorVariant_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_ObjServer_SeparatorVariant_class
  implicit none

  type, public, bind(C) :: c_SeparatorVariant
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning  = .false.
  end type c_SeparatorVariant

  type, public :: SeparatorVariant_type
    type(c_SeparatorVariant) :: c_object
  contains
    final :: SeparatorVariant_class_delete
  end type SeparatorVariant_type

  interface new
    module procedure SeparatorVariant_class_new
    module procedure SeparatorVariant_class_newFromObjectAllocationServer
  end interface
  interface setNumberOfPlanes
    module procedure SeparatorVariant_class_setNumberOfPlanes
  end interface
  interface setPlane
    module procedure SeparatorVariant_class_setPlane
  end interface
  interface copy
    module procedure SeparatorVariant_class_copy
  end interface
  interface getNumberOfPlanes
    module procedure SeparatorVariant_class_getNumberOfPlanes
  end interface
  interface getPlane
    module procedure SeparatorVariant_class_getPlane
  end interface
  interface setPrincipalCurvatures
    module procedure SeparatorVariant_class_setPrincipalCurvatures
  end interface
  interface getPrincipalCurvatures
    module procedure SeparatorVariant_class_getPrincipalCurvatures
  end interface
  interface isFlipped
    module procedure SeparatorVariant_class_isFlipped
  end interface
  interface isPlane
    module procedure SeparatorVariant_class_isPlane
  end interface
  interface isParaboloid
    module procedure SeparatorVariant_class_isParaboloid
  end interface
  interface printToScreen
    module procedure SeparatorVariant_class_printToScreen
  end interface
  interface shiftOrigin
    module procedure SeparatorVariant_class_shift
  end interface

  interface

    subroutine F_SeparatorVariant_new(this) &
      bind(C, name="c_SeparatorVariant_new")
      import
      implicit none
      type(c_SeparatorVariant) :: this
    end subroutine F_SeparatorVariant_new

    subroutine F_SeparatorVariant_newFromObjectAllocationServer(this, a_object_allocation_server) &
      bind(C, name="c_SeparatorVariant_newFromObjectAllocationServer")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      type(c_ObjServer_SeparatorVariant) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<SeparatorVariant>
    end subroutine F_SeparatorVariant_newFromObjectAllocationServer

    subroutine F_SeparatorVariant_delete(this) &
      bind(C, name="c_SeparatorVariant_delete")
      import
      implicit none
      type(c_SeparatorVariant) :: this
    end subroutine F_SeparatorVariant_delete

    subroutine F_SeparatorVariant_setNumberOfPlanes(this, a_number_to_set) &
      bind(C, name="c_SeparatorVariant_setNumberOfPlanes")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      integer(C_INT), intent(in) :: a_number_to_set ! scalar
    end subroutine F_SeparatorVariant_setNumberOfPlanes

    subroutine F_SeparatorVariant_setPlane(this, a_plane_index_to_set,a_normal, a_distance) &
      bind(C, name="c_SeparatorVariant_setPlane")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      integer(C_INT), intent(in) :: a_plane_index_to_set ! scalar
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
      real(C_DOUBLE), intent(in) :: a_distance ! scalar
    end subroutine F_SeparatorVariant_setPlane

    subroutine F_SeparatorVariant_copy(this, a_other_SeparatorVariant) &
      bind(C, name="c_SeparatorVariant_copy")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      type(c_SeparatorVariant) :: a_other_SeparatorVariant
    end subroutine F_SeparatorVariant_copy

    function F_SeparatorVariant_getNumberOfPlanes(this) result(a_number_of_planes) &
      bind(C, name="c_SeparatorVariant_getNumberOfPlanes")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      integer(C_INT) :: a_number_of_planes
    end function F_SeparatorVariant_getNumberOfPlanes

    subroutine F_SeparatorVariant_getPlane(this, a_index, a_plane_listed) &
      bind(C, name="c_SeparatorVariant_getPlane")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_plane_listed
    end subroutine F_SeparatorVariant_getPlane

    subroutine F_SeparatorVariant_setPrincipalCurvatures(this, a_curvatures) &
      bind(C, name="c_SeparatorVariant_setPrincipalCurvatures")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_curvatures !  dimension(1:2)
    end subroutine F_SeparatorVariant_setPrincipalCurvatures

    subroutine F_SeparatorVariant_getPrincipalCurvatures(this, a_curvatures) &
      bind(C, name="c_SeparatorVariant_getPrincipalCurvatures")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_curvatures
    end subroutine F_SeparatorVariant_getPrincipalCurvatures

    function F_SeparatorVariant_isFlipped(this) result(a_flipped) &
      bind(C, name="c_SeparatorVariant_isFlipped")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      logical(C_BOOL) :: a_flipped
    end function F_SeparatorVariant_isFlipped

    function F_SeparatorVariant_isPlane(this) result(a_isplane) &
      bind(C, name="c_SeparatorVariant_isPlane")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      logical(C_BOOL) :: a_isplane
    end function F_SeparatorVariant_isPlane

    function F_SeparatorVariant_isParaboloid(this) result(a_isparaboloid) &
      bind(C, name="c_SeparatorVariant_isParaboloid")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      logical(C_BOOL) :: a_isparaboloid
    end function F_SeparatorVariant_isParaboloid

    subroutine F_SeparatorVariant_printToScreen(this) &
      bind(C, name="c_SeparatorVariant_printToScreen")
      import
      implicit none
      type(c_SeparatorVariant) :: this
    end subroutine F_SeparatorVariant_printToScreen

    subroutine F_SeparatorVariant_shift(this, a_shift) &
      bind(C, name="c_SeparatorVariant_shift")
      import
      implicit none
      type(c_SeparatorVariant) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_shift !  dimension(1:3)
    end subroutine F_SeparatorVariant_shift

  end interface

  contains

    subroutine SeparatorVariant_class_new(this)
      implicit none
      type(SeparatorVariant_type), intent(inout) :: this
      call F_SeparatorVariant_new(this%c_object)
    end subroutine SeparatorVariant_class_new

    subroutine SeparatorVariant_class_newFromObjectAllocationServer(this, a_object_allocation_server)
      implicit none
      type(SeparatorVariant_type), intent(inout) :: this
      type(ObjServer_SeparatorVariant_type), intent(in) :: a_object_allocation_server
      call F_SeparatorVariant_newFromObjectAllocationServer(this%c_object, a_object_allocation_server%c_object)
    end subroutine SeparatorVariant_class_newFromObjectAllocationServer

    impure elemental subroutine SeparatorVariant_class_delete(this)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      call F_SeparatorVariant_delete(this%c_object)
    end subroutine SeparatorVariant_class_delete

    subroutine SeparatorVariant_class_setNumberOfPlanes(this, a_number_to_set)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_number_to_set
      call F_SeparatorVariant_setNumberOfPlanes(this%c_object, a_number_to_set)
    end subroutine SeparatorVariant_class_setNumberOfPlanes

    subroutine SeparatorVariant_class_setPlane(this, a_plane_index_to_set,a_normal, a_distance)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index_to_set
      real(IRL_double), dimension(1:3), intent(in) :: a_normal
      real(IRL_double), intent(in) :: a_distance
      call F_SeparatorVariant_setPlane(this%c_object, a_plane_index_to_set, a_normal, a_distance)
    end subroutine SeparatorVariant_class_setPlane

    subroutine SeparatorVariant_class_copy(this, a_other_SeparatorVariant)
      implicit none
      type(SeparatorVariant_type), intent(inout) :: this
      type(SeparatorVariant_type), intent(in) :: a_other_SeparatorVariant
      call F_SeparatorVariant_copy(this%c_object, a_other_SeparatorVariant%c_object)
    end subroutine SeparatorVariant_class_copy

    function SeparatorVariant_class_getNumberOfPlanes(this) result(a_number_of_planes)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_planes
      a_number_of_planes = F_SeparatorVariant_getNumberOfPlanes(this%c_object)
    end function SeparatorVariant_class_getNumberOfPlanes

    function SeparatorVariant_class_getPlane(this, a_index) result(a_plane_listed)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(4) :: a_plane_listed
      call F_SeparatorVariant_getPlane(this%c_object, a_index, a_plane_listed)
    end function SeparatorVariant_class_getPlane

    subroutine SeparatorVariant_class_setPrincipalCurvatures(this, a_curvatures)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      real(IRL_double), dimension(1:2), intent(in) :: a_curvatures
      call F_SeparatorVariant_setPrincipalCurvatures(this%c_object, a_curvatures)
    end subroutine SeparatorVariant_class_setPrincipalCurvatures

    function SeparatorVariant_class_getPrincipalCurvatures(this) result(a_curvatures)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      real(IRL_double), dimension(2) :: a_curvatures
      call F_SeparatorVariant_getPrincipalCurvatures(this%c_object, a_curvatures)
    end function SeparatorVariant_class_getPrincipalCurvatures

    function SeparatorVariant_class_isFlipped(this) result(a_flipped)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      logical(1) :: a_flipped
      a_flipped = F_SeparatorVariant_isFlipped(this%c_object)
      return
    end function SeparatorVariant_class_isFlipped

    function SeparatorVariant_class_isPlane(this) result(a_isplane)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      logical(1) :: a_isplane
      a_isplane = F_SeparatorVariant_isPlane(this%c_object)
      return
    end function SeparatorVariant_class_isPlane

    function SeparatorVariant_class_isParaboloid(this) result(a_isparaboloid)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      logical(1) :: a_isparaboloid
      a_isparaboloid = F_SeparatorVariant_isParaboloid(this%c_object)
      return
    end function SeparatorVariant_class_isParaboloid

    subroutine SeparatorVariant_class_printToScreen(this)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      call F_SeparatorVariant_printToScreen(this%c_object)
    end subroutine SeparatorVariant_class_printToScreen

    subroutine SeparatorVariant_class_shift(this, a_shift)
      implicit none
      type(SeparatorVariant_type), intent(in) :: this
      real(IRL_double), dimension(1:3), intent(in) :: a_shift
      call F_SeparatorVariant_shift(this%c_object, a_shift)
    end subroutine SeparatorVariant_class_shift


end module f_SeparatorVariant_class
