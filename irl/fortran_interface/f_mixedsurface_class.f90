!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module f_MixedPolygonBezierSurface_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  use f_SeparatorVariant_class
  use f_SeparatorUnion_class
  use f_RectCub_class
  use f_ObjServer_MixedPolygonBezierSurface_class
  implicit none

  type, public, bind(C) :: c_MixedPolygonBezierSurface
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning  = .false.
  end type c_MixedPolygonBezierSurface

  type, public :: MixedPolygonBezierSurface_type
    type(c_MixedPolygonBezierSurface) :: c_object
  contains
    final :: MixedPolygonBezierSurface_delete
  end type MixedPolygonBezierSurface_type

  interface new
    module procedure MixedPolygonBezierSurface_new
    module procedure MixedPolygonBezierSurface_newFromObjectAllocationServer
  end interface
  interface getNumberOfPoints
    module procedure MixedPolygonBezierSurface_getNumberOfPoints
  end interface  
  interface getNumberOfPolygons
    module procedure MixedPolygonBezierSurface_getNumberOfPolygons
  end interface  
  interface getNumberOfTriangles
    module procedure MixedPolygonBezierSurface_getNumberOfTriangles
  end interface  
  interface getPolygonSize
    module procedure MixedPolygonBezierSurface_getPolygonSize
  end interface  
  interface getTri
    module procedure MixedPolygonBezierSurface_getTri
  end interface
  interface getPt
    module procedure MixedPolygonBezierSurface_getPt
  end interface
  interface zeroMixedSurface
    module procedure MixedPolygonBezierSurface_clear
  end interface  
  interface getSurface
    module procedure MixedPolygonBezierSurface_getSurface_RectCub_Variant
    module procedure MixedPolygonBezierSurface_getSurface_RectCub_SepUnion_raw
  end interface  


  interface

    subroutine F_MixedPolygonBezierSurface_new(this) &
      bind(C, name="c_MixedPolygonBezierSurface_new")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
    end subroutine F_MixedPolygonBezierSurface_new

    subroutine F_MixedPolygonBezierSurface_newFromObjectAllocationServer(this, a_object_allocation_server) &
      bind(C, name="c_MixedPolygonBezierSurface_newFromObjectAllocationServer")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
      type(c_ObjServer_MixedPolygonBezierSurface) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<MixedPolygonBezierSurface>
    end subroutine F_MixedPolygonBezierSurface_newFromObjectAllocationServer

    subroutine F_MixedPolygonBezierSurface_delete(this) &
      bind(C, name="c_MixedPolygonBezierSurface_delete")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
    end subroutine F_MixedPolygonBezierSurface_delete

    function F_MixedPolygonBezierSurface_getNumberOfPoints(this) result(a_number_of_pts) &
      bind(C, name="c_MixedPolygonBezierSurface_getNumberOfPoints")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
      integer(C_INT) :: a_number_of_pts
    end function F_MixedPolygonBezierSurface_getNumberOfPoints

    function F_MixedPolygonBezierSurface_getNumberOfPolygons(this) result(a_number_of_polys) &
      bind(C, name="c_MixedPolygonBezierSurface_getNumberOfPolygons")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
      integer(C_INT) :: a_number_of_polys
    end function F_MixedPolygonBezierSurface_getNumberOfPolygons

    function F_MixedPolygonBezierSurface_getNumberOfTriangles(this) result(a_number_of_tris) &
      bind(C, name="c_MixedPolygonBezierSurface_getNumberOfTriangles")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
      integer(C_INT) :: a_number_of_tris
    end function F_MixedPolygonBezierSurface_getNumberOfTriangles

    function F_MixedPolygonBezierSurface_getPolygonSize(this, a_poly_index) result(a_poly_size) &
      bind(C, name="c_MixedPolygonBezierSurface_getPolygonSize")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
      integer(C_INT) :: a_poly_index
      integer(C_INT) :: a_poly_size
    end function F_MixedPolygonBezierSurface_getPolygonSize

    subroutine F_MixedPolygonBezierSurface_getTri(this, a_index_tri, a_conn_tri) &
      bind(C, name="c_MixedPolygonBezierSurface_getQuadraticTriangleConnectivity")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
      integer(C_INT) :: a_index_tri
      integer(C_INT), dimension(*), intent(out) :: a_conn_tri 
    end subroutine F_MixedPolygonBezierSurface_getTri

    subroutine F_MixedPolygonBezierSurface_getPt(this, a_index_vert, a_pt) &
      bind(C, name="c_MixedPolygonBezierSurface_getPt")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
      integer(C_INT) :: a_index_vert
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt
    end subroutine F_MixedPolygonBezierSurface_getPt

    subroutine F_MixedPolygonBezierSurface_clear(this) &
      bind(C, name="c_MixedPolygonBezierSurface_clear")
      import
      implicit none
      type(c_MixedPolygonBezierSurface) :: this
    end subroutine F_MixedPolygonBezierSurface_clear

    subroutine F_MixedPolygonBezierSurface_getSurface_RectCub_Variant(a_cuboid, a_variant, a_surface) &
    bind(C, name="c_MixedPolygonBezierSurface_getSurface_RectCub_Variant")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_cuboid 
      type(c_SeparatorVariant) :: a_variant 
      type(c_MixedPolygonBezierSurface) :: a_surface 
    end subroutine F_MixedPolygonBezierSurface_getSurface_RectCub_Variant

    subroutine F_MixedPolygonBezierSurface_getSurface_RectCub_SepUnion_raw(a_cuboid, a_sepunion, a_surface) &
    bind(C, name="c_MixedPolygonBezierSurface_getSurface_RectCub_SepUnion_raw")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_RectCub) :: a_cuboid 
      type(SeparatorUnion_type_raw) :: a_sepunion 
      type(c_MixedPolygonBezierSurface) :: a_surface 
    end subroutine F_MixedPolygonBezierSurface_getSurface_RectCub_SepUnion_raw

  end interface


  contains

    subroutine MixedPolygonBezierSurface_new(this)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(inout) :: this
      call F_MixedPolygonBezierSurface_new(this%c_object)
    end subroutine MixedPolygonBezierSurface_new

    subroutine MixedPolygonBezierSurface_newFromObjectAllocationServer(this, a_object_allocation_server)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(inout) :: this
      type(ObjServer_MixedPolygonBezierSurface_type), intent(in) :: a_object_allocation_server
      call F_MixedPolygonBezierSurface_newFromObjectAllocationServer(this%c_object, a_object_allocation_server%c_object)
    end subroutine MixedPolygonBezierSurface_newFromObjectAllocationServer

    impure elemental subroutine MixedPolygonBezierSurface_delete(this)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(in) :: this
      call F_MixedPolygonBezierSurface_delete(this%c_object)
    end subroutine MixedPolygonBezierSurface_delete

    function MixedPolygonBezierSurface_getNumberOfPoints(this) result(a_number_of_pts)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_pts
      a_number_of_pts = F_MixedPolygonBezierSurface_getNumberOfPoints(this%c_object)
      return
    end function MixedPolygonBezierSurface_getNumberOfPoints

    function MixedPolygonBezierSurface_getNumberOfPolygons(this) result(a_number_of_polys)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_polys
      a_number_of_polys = F_MixedPolygonBezierSurface_getNumberOfPolygons(this%c_object)
      return
    end function MixedPolygonBezierSurface_getNumberOfPolygons

    function MixedPolygonBezierSurface_getNumberOfTriangles(this) result(a_number_of_tris)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_number_of_tris
      a_number_of_tris = F_MixedPolygonBezierSurface_getNumberOfTriangles(this%c_object)
      return
    end function MixedPolygonBezierSurface_getNumberOfTriangles

    function MixedPolygonBezierSurface_getPolygonSize(this, a_poly_index) result(a_poly_size)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_poly_index
      integer(IRL_UnsignedIndex_t) :: a_poly_size
      a_poly_size = F_MixedPolygonBezierSurface_getPolygonSize(this%c_object, a_poly_index)
      return
    end function MixedPolygonBezierSurface_getPolygonSize

    function MixedPolygonBezierSurface_getTri(this, a_index_tri) result(a_conn)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index_tri
      integer(IRL_UnsignedIndex_t), dimension(6) :: a_conn
      call F_MixedPolygonBezierSurface_getTri(this%c_object, a_index_tri, a_conn)
      return
    end function MixedPolygonBezierSurface_getTri

    function MixedPolygonBezierSurface_getPt(this, a_index_vert) result(a_pt)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index_vert
      real(IRL_double), dimension(4) :: a_pt
      call F_MixedPolygonBezierSurface_getPt(this%c_object, a_index_vert, a_pt)
      return
    end function MixedPolygonBezierSurface_getPt

    impure elemental subroutine MixedPolygonBezierSurface_clear(this)
      implicit none
      type(MixedPolygonBezierSurface_type), intent(in) :: this
      call F_MixedPolygonBezierSurface_clear(this%c_object)
    end subroutine MixedPolygonBezierSurface_clear

    subroutine MixedPolygonBezierSurface_getSurface_RectCub_Variant(a_cuboid, a_variant, a_surface)
      use, intrinsic :: iso_c_binding
      implicit none
        type(RectCub_type), intent(in) :: a_cuboid
        type(SeparatorVariant_type), intent(in) :: a_variant
        type(MixedPolygonBezierSurface_type), intent(inout) :: a_surface
        call F_MixedPolygonBezierSurface_getSurface_RectCub_Variant(a_cuboid%c_object, a_variant%c_object, a_surface%c_object)
    end subroutine MixedPolygonBezierSurface_getSurface_RectCub_Variant 

    subroutine MixedPolygonBezierSurface_getSurface_RectCub_SepUnion_raw(a_cuboid, a_sepunion, a_surface)
      use, intrinsic :: iso_c_binding
      implicit none
        type(RectCub_type), intent(in) :: a_cuboid
        type(SeparatorUnion_type_raw), intent(in) :: a_sepunion
        type(MixedPolygonBezierSurface_type), intent(inout) :: a_surface
        call F_MixedPolygonBezierSurface_getSurface_RectCub_SepUnion_raw(a_cuboid%c_object, a_sepunion, a_surface%c_object)
    end subroutine MixedPolygonBezierSurface_getSurface_RectCub_SepUnion_raw 


end module f_MixedPolygonBezierSurface_class
