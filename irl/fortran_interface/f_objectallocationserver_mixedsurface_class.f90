!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module f_ObjServer_MixedPolygonBezierSurface_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_ObjServer_MixedPolygonBezierSurface
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_ObjServer_MixedPolygonBezierSurface

  type, public :: ObjServer_MixedPolygonBezierSurface_type
    type(c_ObjServer_MixedPolygonBezierSurface) :: c_object
  contains
    final :: ObjServer_MixedPolygonBezierSurface_class_delete
  end type ObjServer_MixedPolygonBezierSurface_type

  interface new
    module procedure ObjServer_MixedPolygonBezierSurface_class_new
  end interface

  interface

    subroutine F_ObjServer_MixedPolygonBezierSurface_new(this, a_number_to_allocate) &
      bind(C, name="c_ObjServer_MixedPolygonBezierSurface_new")
      import
      implicit none
      type(c_ObjServer_MixedPolygonBezierSurface) :: this
      integer(C_SIZE_T) :: a_number_to_allocate
    end subroutine F_ObjServer_MixedPolygonBezierSurface_new

    subroutine F_ObjServer_MixedPolygonBezierSurface_delete(this) &
      bind(C, name="c_ObjServer_MixedPolygonBezierSurface_delete")
      import
      implicit none
      type(c_ObjServer_MixedPolygonBezierSurface) :: this
    end subroutine F_ObjServer_MixedPolygonBezierSurface_delete

  end interface

  contains

    subroutine ObjServer_MixedPolygonBezierSurface_class_new(this, a_number_to_allocate)
      implicit none
      type(ObjServer_MixedPolygonBezierSurface_type), intent(inout) :: this
      integer(IRL_LargeOffsetIndex_t) :: a_number_to_allocate
      call F_ObjServer_MixedPolygonBezierSurface_new(this%c_object, a_number_to_allocate)
    end subroutine ObjServer_MixedPolygonBezierSurface_class_new

    impure elemental subroutine ObjServer_MixedPolygonBezierSurface_class_delete(this)
      implicit none
      type(ObjServer_MixedPolygonBezierSurface_type), intent(in) :: this
      call F_ObjServer_MixedPolygonBezierSurface_delete(this%c_object)
    end subroutine ObjServer_MixedPolygonBezierSurface_class_delete

end module f_ObjServer_MixedPolygonBezierSurface_class
