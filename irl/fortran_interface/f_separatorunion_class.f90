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

module f_SeparatorUnion_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  integer, parameter :: SEP_UNION_NUM_DOUBLES = 16
  type, public, bind(C) :: SeparatorUnion_type_raw
    real(C_DOUBLE) :: raw_data(SEP_UNION_NUM_DOUBLES) = 0.0_C_DOUBLE !!! Must add up to 128 bits with 8-byte alignment to match C++ object
  end type SeparatorUnion_type_raw

  interface setToOnePlane
    module procedure SeparatorUnion_class_setToOnePlane_raw
  end interface
  interface setToFull
    module procedure SeparatorUnion_class_setToFull_raw
  end interface
  interface setToEmpty
    module procedure SeparatorUnion_class_setToEmpty_raw
  end interface
  interface getPlane
    module procedure SeparatorUnion_class_getPlane_raw
  end interface
  interface reflect
    module procedure SeparatorUnion_class_reflect_raw
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
  interface getMeanCurvature
    module procedure SeparatorUnion_class_getMeanCurvature_raw
  end interface

  interface

    subroutine F_SeparatorUnion_setToOnePlane_raw(this, a_normal, a_distance) bind(C, name="c_SeparatorUnion_setToOnePlane_raw")
      import
      implicit none
      type(SeparatorUnion_type_raw) :: this
      real(C_DOUBLE), dimension(3), intent(in) :: a_normal !  dimension(1:3)
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
      real(C_DOUBLE), dimension(4), intent(out) :: a_plane_listed
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

    subroutine SeparatorUnion_class_setToOnePlane_raw(this, a_normal, a_distance)
      implicit none
      type(SeparatorUnion_type_raw), intent(inout) :: this
      real(C_DOUBLE), dimension(3), intent(in) :: a_normal !  dimension(1:3)
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
      type(SeparatorUnion_type_raw), intent(in) :: a_reflected_sep
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
      real(C_DOUBLE) :: a_curvature
      a_curvature = F_SeparatorUnion_getMeanCurvature_raw(this)
    end function SeparatorUnion_class_getMeanCurvature_raw

end module f_SeparatorUnion_class
