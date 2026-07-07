!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_SeparatorUnion_class.f90
!!
!! This file allows use of the IRL SeparatorUnion
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's SeparatorUnion class along with enabling
!! some of its methods.
module f_SeparatorUnion_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  ! use f_ObjServer_SeparatorUnion_class
  implicit none

  type, public, bind(C) :: c_SeparatorUnion
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning  = .false.
  end type c_SeparatorUnion

  type, public :: SeparatorUnion_type
    type(c_SeparatorUnion) :: c_object
  contains
    final :: SeparatorUnion_class_delete
  end type SeparatorUnion_type

  integer, parameter :: SEPARATORUNION_SIZE = 128
  type, public, bind(C) :: SeparatorUnion_type_raw
    integer(c_int8_t) :: raw(128)
  end type SeparatorUnion_type_raw

  interface new
    module procedure SeparatorUnion_class_new
  end interface
  interface setNumberOfPlanes 
    module procedure SeparatorUnion_class_setNumberOfPlanes
  end interface
  interface setPlane
    module procedure SeparatorUnion_class_setPlane
  end interface
  interface setToOnePlane
    module procedure SeparatorUnion_class_setToOnePlane_raw
  end interface
  interface setToFull
    module procedure SeparatorUnion_class_setToFull_raw
  end interface
  interface setToEmpty
    module procedure SeparatorUnion_class_setToEmpty_raw
  end interface
  interface copy
    module procedure SeparatorUnion_class_copy
  end interface
  interface getNumberOfPlanes
    module procedure SeparatorUnion_class_getNumberOfPlanes
  end interface
  interface getPlane
    module procedure SeparatorUnion_class_getPlane
    module procedure SeparatorUnion_class_getPlane_raw
  end interface
  interface reflect
    module procedure SeparatorUnion_class_reflect_raw
  end interface
  interface isFlipped
    module procedure SeparatorUnion_class_isFlipped
  end interface
  interface isOnePlane
    module procedure SeparatorUnion_class_isOnePlane_raw
  end interface
  interface isParaboloid
    module procedure SeparatorUnion_class_isParaboloid_raw
  end interface
  interface isTypeDefined
    module procedure SeparatorUnion_class_isTypeDefined_raw
  end interface
  interface isFull
    module procedure SeparatorUnion_class_isFull_raw
  end interface
  interface isEmpty
    module procedure SeparatorUnion_class_isEmpty_raw
  end interface
  interface printToScreen
    module procedure SeparatorUnion_class_printToScreen
  end interface
  interface shiftOrigin
    module procedure SeparatorUnion_class_shift
  end interface
  interface getMeanCurvature
    module procedure SeparatorUnion_class_getMeanCurvature_raw
  end interface

  interface

    subroutine F_SeparatorUnion_new(this) &
      bind(C, name="c_SeparatorUnion_new")
      import
      implicit none
      type(c_SeparatorUnion) :: this
    end subroutine F_SeparatorUnion_new

    subroutine F_SeparatorUnion_delete(this) &
      bind(C, name="c_SeparatorUnion_delete")
      import
      implicit none
      type(c_SeparatorUnion) :: this
    end subroutine F_SeparatorUnion_delete

    subroutine F_SeparatorUnion_setNumberOfPlanes(this, a_number_to_set) &
      bind(C, name="c_SeparatorUnion_setNumberOfPlanes")
      import
      implicit none
      type(c_SeparatorUnion) :: this
      integer(C_INT), intent(in) :: a_number_to_set ! scalar
    end subroutine F_SeparatorUnion_setNumberOfPlanes

    subroutine F_SeparatorUnion_setPlane(this, a_plane_index_to_set,a_normal, a_distance) &
      bind(C, name="c_SeparatorUnion_setPlane")
      import
      implicit none
      type(c_SeparatorUnion) :: this
      integer(C_INT), intent(in) :: a_plane_index_to_set ! scalar
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
      real(C_DOUBLE), intent(in) :: a_distance ! scalar
    end subroutine F_SeparatorUnion_setPlane

    subroutine F_SeparatorUnion_copy(this, a_other_SeparatorUnion) &
      bind(C, name="c_SeparatorUnion_copy")
      import
      implicit none
      type(c_SeparatorUnion) :: this
      type(c_SeparatorUnion) :: a_other_SeparatorUnion
    end subroutine F_SeparatorUnion_copy

    function F_SeparatorUnion_getNumberOfPlanes(this) result(a_number_of_planes) &
      bind(C, name="c_SeparatorUnion_getNumberOfPlanes")
      import
      implicit none
      type(c_SeparatorUnion) :: this
      integer(C_INT) :: a_number_of_planes
    end function F_SeparatorUnion_getNumberOfPlanes

    subroutine F_SeparatorUnion_getPlane(this, a_index, a_plane_listed) &
      bind(C, name="c_SeparatorUnion_getPlane")
      import
      implicit none
      type(c_SeparatorUnion) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_plane_listed
    end subroutine F_SeparatorUnion_getPlane

    function F_SeparatorUnion_isFlipped(this) result(a_flipped) &
      bind(C, name="c_SeparatorUnion_isFlipped")
      import
      implicit none
      type(c_SeparatorUnion) :: this
      logical(C_BOOL) :: a_flipped
    end function F_SeparatorUnion_isFlipped

    subroutine F_SeparatorUnion_printToScreen(this) &
      bind(C, name="c_SeparatorUnion_printToScreen")
      import
      implicit none
      type(c_SeparatorUnion) :: this
    end subroutine F_SeparatorUnion_printToScreen

    subroutine F_SeparatorUnion_shift(this, a_shift) &
      bind(C, name="c_SeparatorUnion_shift")
      import
      implicit none
      type(c_SeparatorUnion) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_shift !  dimension(1:3)
    end subroutine F_SeparatorUnion_shift

    ! subroutine F_SeparatorUnion_setToOnePlane_raw(this, a_normal, a_distance) bind(C, name="c_SeparatorUnion_setToOnePlane_raw")
    !   import
    !   implicit none
    !   type(C_PTR) :: this
    !   real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
    !   real(C_DOUBLE), intent(in) :: a_distance ! scalar
    ! end subroutine F_SeparatorUnion_setToOnePlane_raw

    subroutine F_SeparatorUnion_setToOnePlane_raw(this, a_normal, a_distance) bind(C, name="c_SeparatorUnion_setToOnePlane_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
      real(C_DOUBLE), intent(in) :: a_distance ! scalar
    end subroutine F_SeparatorUnion_setToOnePlane_raw

    subroutine F_SeparatorUnion_setToFull_raw(this) bind(C, name="c_SeparatorUnion_setToFull_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
    end subroutine F_SeparatorUnion_setToFull_raw

    subroutine F_SeparatorUnion_setToEmpty_raw(this) bind(C, name="c_SeparatorUnion_setToEmpty_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
    end subroutine F_SeparatorUnion_setToEmpty_raw

    subroutine F_SeparatorUnion_getPlane_raw(this, a_index, a_plane_listed) &
      bind(C, name="c_SeparatorUnion_getPlane_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_plane_listed
    end subroutine F_SeparatorUnion_getPlane_raw

    subroutine F_SeparatorUnion_reflect_raw(this, a_reflected_sep, a_dir, a_loc) bind(C, name="c_SeparatorUnion_reflect_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this, a_reflected_sep
      integer(C_INT), intent(in) :: a_dir ! scalar
      real(C_DOUBLE), intent(in) :: a_loc ! scalar
    end subroutine F_SeparatorUnion_reflect_raw

    function F_SeparatorUnion_isOnePlane_raw(this) result(a_isoneplane) &
      bind(C, name="c_SeparatorUnion_isOnePlane_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      logical(C_BOOL) :: a_isoneplane
    end function F_SeparatorUnion_isOnePlane_raw

    function F_SeparatorUnion_isParaboloid_raw(this) result(a_isparaboloid) &
      bind(C, name="c_SeparatorUnion_isParaboloid_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      logical(C_BOOL) :: a_isparaboloid
    end function F_SeparatorUnion_isParaboloid_raw

    function F_SeparatorUnion_isTypeDefined_raw(this) result(a_isdef) &
      bind(C, name="c_SeparatorUnion_isTypeDefined_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      logical(C_BOOL) :: a_isdef
    end function F_SeparatorUnion_isTypeDefined_raw

    function F_SeparatorUnion_isFull_raw(this) result(a_isfull) &
      bind(C, name="c_SeparatorUnion_isFull_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      logical(C_BOOL) :: a_isfull
    end function F_SeparatorUnion_isFull_raw

    function F_SeparatorUnion_isEmpty_raw(this) result(a_isempty) &
      bind(C, name="c_SeparatorUnion_isEmpty_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      logical(C_BOOL) :: a_isempty
    end function F_SeparatorUnion_isEmpty_raw
    
    function F_SeparatorUnion_getMeanCurvature_raw(this) result(a_curvature) &
      bind(C, name="c_SeparatorUnion_getMeanCurvature_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      real(C_DOUBLE) :: a_curvature
    end function F_SeparatorUnion_getMeanCurvature_raw

  end interface

  contains

    subroutine SeparatorUnion_class_new(this)
      implicit none
      type(SeparatorUnion_type), intent(inout) :: this
      call F_SeparatorUnion_new(this%c_object)
    end subroutine SeparatorUnion_class_new

    impure elemental subroutine SeparatorUnion_class_delete(this)
      implicit none
      type(SeparatorUnion_type), intent(in) :: this
      call F_SeparatorUnion_delete(this%c_object)
    end subroutine SeparatorUnion_class_delete

    subroutine SeparatorUnion_class_setNumberOfPlanes(this, a_number_to_set)
      implicit none
      type(SeparatorUnion_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_number_to_set
      call F_SeparatorUnion_setNumberOfPlanes(this%c_object, a_number_to_set)
    end subroutine SeparatorUnion_class_setNumberOfPlanes

    subroutine SeparatorUnion_class_setPlane(this, a_plane_index_to_set,a_normal, a_distance)
      implicit none
      type(SeparatorUnion_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index_to_set
      real(IRL_double), dimension(1:3), intent(in) :: a_normal
      real(IRL_double), intent(in) :: a_distance
      call F_SeparatorUnion_setPlane(this%c_object, a_plane_index_to_set, a_normal, a_distance)
    end subroutine SeparatorUnion_class_setPlane

    subroutine SeparatorUnion_class_copy(this, a_other_SeparatorUnion)
      implicit none
      type(SeparatorUnion_type), intent(inout) :: this
      type(SeparatorUnion_type), intent(in) :: a_other_SeparatorUnion
      call F_SeparatorUnion_copy(this%c_object, a_other_SeparatorUnion%c_object)
    end subroutine SeparatorUnion_class_copy

    function SeparatorUnion_class_getNumberOfPlanes(this) result(a_number_of_planes)
      implicit none
      type(SeparatorUnion_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_planes
      a_number_of_planes = F_SeparatorUnion_getNumberOfPlanes(this%c_object)
    end function SeparatorUnion_class_getNumberOfPlanes

    function SeparatorUnion_class_getPlane(this, a_index) result(a_plane_listed)
      implicit none
      type(SeparatorUnion_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(4) :: a_plane_listed
      call F_SeparatorUnion_getPlane(this%c_object, a_index, a_plane_listed)
    end function SeparatorUnion_class_getPlane

    function SeparatorUnion_class_isFlipped(this) result(a_flipped)
      implicit none
      type(SeparatorUnion_type), intent(in) :: this
      logical(1) :: a_flipped
      a_flipped = F_SeparatorUnion_isFlipped(this%c_object)
      return
    end function SeparatorUnion_class_isFlipped

    subroutine SeparatorUnion_class_printToScreen(this)
      implicit none
      type(SeparatorUnion_type), intent(in) :: this
      call F_SeparatorUnion_printToScreen(this%c_object)
    end subroutine SeparatorUnion_class_printToScreen

    subroutine SeparatorUnion_class_shift(this, a_shift)
      implicit none
      type(SeparatorUnion_type), intent(in) :: this
      real(IRL_double), dimension(1:3), intent(in) :: a_shift
      call F_SeparatorUnion_shift(this%c_object, a_shift)
    end subroutine SeparatorUnion_class_shift

    ! Routines for raw type
    subroutine SeparatorUnion_class_setToOnePlane_raw(this, a_normal, a_distance)
      implicit none
      type(SeparatorUnion_type_raw), intent(inout) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_normal !  dimension(1:3)
      real(C_DOUBLE), intent(in) :: a_distance ! scalar
      call F_SeparatorUnion_setToOnePlane_raw(this, a_normal, a_distance)
    end subroutine SeparatorUnion_class_setToOnePlane_raw

    subroutine SeparatorUnion_class_setToFull_raw(this)
      implicit none
      type(SeparatorUnion_type_raw), intent(inout) :: this
      call F_SeparatorUnion_setToFull_raw(this)
    end subroutine SeparatorUnion_class_setToFull_raw

    subroutine SeparatorUnion_class_setToEmpty_raw(this)
      implicit none
      type(SeparatorUnion_type_raw), intent(inout) :: this
      call F_SeparatorUnion_setToEmpty_raw(this)
    end subroutine SeparatorUnion_class_setToEmpty_raw

    function SeparatorUnion_class_getPlane_raw(this, a_index) result(a_plane_listed)
      implicit none
      type(SeparatorUnion_type_raw), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(4) :: a_plane_listed
      call F_SeparatorUnion_getPlane_raw(this, a_index, a_plane_listed)
    end function SeparatorUnion_class_getPlane_raw

    subroutine SeparatorUnion_class_reflect_raw(this, a_reflected_sep, a_dir, a_loc)
      implicit none
      type(SeparatorUnion_type_raw), intent(inout) :: this
      type(SeparatorUnion_type_raw), intent(in) ::a_reflected_sep
      integer(C_INT), intent(in) :: a_dir
      real(C_DOUBLE), intent(in) :: a_loc
      call F_SeparatorUnion_reflect_raw(this, a_reflected_sep, a_dir, a_loc)
    end subroutine SeparatorUnion_class_reflect_raw

    function SeparatorUnion_class_isOnePlane_raw(this) result(a_isoneplane)
      implicit none
      type(SeparatorUnion_type_raw), intent(in) :: this
      logical(1) :: a_isoneplane
      a_isoneplane = F_SeparatorUnion_isOnePlane_raw(this)
      return
    end function SeparatorUnion_class_isOnePlane_raw

    function SeparatorUnion_class_isParaboloid_raw(this) result(a_isparaboloid)
      implicit none
      type(SeparatorUnion_type_raw), intent(in) :: this
      logical(1) :: a_isparaboloid
      a_isparaboloid = F_SeparatorUnion_isParaboloid_raw(this)
      return
    end function SeparatorUnion_class_isParaboloid_raw

    function SeparatorUnion_class_isTypeDefined_raw(this) result(a_isdef)
      implicit none
      type(SeparatorUnion_type_raw), intent(in) :: this
      logical(1) :: a_isdef
      a_isdef = F_SeparatorUnion_isTypeDefined_raw(this)
      return
    end function SeparatorUnion_class_isTypeDefined_raw

    function SeparatorUnion_class_isFull_raw(this) result(a_isfull)
      implicit none
      type(SeparatorUnion_type_raw), intent(in) :: this
      logical(1) :: a_isfull
      a_isfull = F_SeparatorUnion_isFull_raw(this)
      return
    end function SeparatorUnion_class_isFull_raw

    function SeparatorUnion_class_isEmpty_raw(this) result(a_isempty)
      implicit none
      type(SeparatorUnion_type_raw), intent(in) :: this
      logical(1) :: a_isempty
      a_isempty = F_SeparatorUnion_isEmpty_raw(this)
      return
    end function SeparatorUnion_class_isEmpty_raw

    function SeparatorUnion_class_getMeanCurvature_raw(this) result(a_curvature)
      implicit none
      type(SeparatorUnion_type_raw), intent(in) :: this
      real(IRL_double) :: a_curvature
      a_curvature = F_SeparatorUnion_getMeanCurvature_raw(this)
    end function SeparatorUnion_class_getMeanCurvature_raw


end module f_SeparatorUnion_class
