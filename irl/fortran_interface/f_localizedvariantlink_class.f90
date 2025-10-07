!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_LocalizedVariantLink_class.f90
!!
!! This file allows use of the IRL LocalizedSeparatorVariantLink
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's LocalizedSeparatorVariantLink class along with enabling
!! some of its methods.
module f_LocVariantLink_class
  use f_DefinedTypes
  use f_ObjServer_LocVariantLink_class
  use f_PlanarLoc_class
  use f_SeparatorVariant_class
  use, intrinsic :: iso_c_binding
  implicit none

  type, public, bind(C) :: c_LocVariantLink
    type(C_PTR), private :: object = C_NULL_PTR
    logical(C_BOOL), private :: is_owning = .false.
  end type c_LocVariantLink

  type, public :: LocVariantLink_type
    type(c_LocVariantLink) :: c_object
  contains
    final :: LocVariantLink_class_delete
  end type LocVariantLink_type

  interface new
    module procedure LocVariantLink_class_new
    module procedure LocVariantLink_class_newFromObjectAllocationServer
  end interface
  interface setId
    module procedure LocVariantLink_class_setId
  end interface
  interface getId
    module procedure LocVariantLink_class_getId
  end interface
  interface setEdgeConnectivity
    module procedure LocVariantLink_class_setEdgeConnectivity
  end interface
  interface setEdgeConnectivityNull
    module procedure LocVariantLink_class_setEdgeConnectivityNull
  end interface
  interface printToScreen
    module procedure LocVariantLink_class_printToScreen
  end interface

  interface

    subroutine F_LocVariantLink_new(this, a_planar_localizer, a_separator_variant) &
      bind(C, name="c_LocVariantLink_new")
      import
      implicit none
      type(c_LocVariantLink) :: this
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_SeparatorVariant) :: a_separator_variant ! Pointer to SeparatorVariant
    end subroutine F_LocVariantLink_new

    subroutine F_LocVariantLink_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                            a_planar_localizer, a_separator_variant) &
      bind(C, name="c_LocVariantLink_newFromObjectAllocationServer")
      import
      implicit none
      type(c_LocVariantLink) :: this
      type(c_ObjServer_LocVariantLink) :: a_object_allocation_server ! Pointer to ObjectAllocationServer<LocVariantLink>
      type(c_PlanarLoc) :: a_planar_localizer ! Pointer to PlanarLoc
      type(c_SeparatorVariant) :: a_separator_variant ! Pointer to SeparatorVariant
    end subroutine F_LocVariantLink_newFromObjectAllocationServer

    subroutine F_LocVariantLink_delete(this) &
      bind(C, name="c_LocVariantLink_delete")
      import
      implicit none
      type(c_LocVariantLink) :: this
    end subroutine F_LocVariantLink_delete

    subroutine F_LocVariantLink_setId(this, a_id) &
      bind(C, name="c_LocVariantLink_setId")
      import
      implicit none
      type(c_LocVariantLink) :: this
      integer(C_INT), intent(in) :: a_id ! Scalar >= 0
    end subroutine F_LocVariantLink_setId

    function F_LocVariantLink_getId(this) result(a_id) &
      bind(C, name="c_LocVariantLink_getId")
      import
      implicit none
      type(c_LocVariantLink) :: this
      integer(C_INT) :: a_id ! Scalar >= 0
    end function F_LocVariantLink_getId

    subroutine F_LocVariantLink_setEdgeConnectivity(this, a_plane_index, a_ptr_to_neighbor) &
      bind(C, name="c_LocVariantLink_setEdgeConnectivity")
      import
      implicit none
      type(c_LocVariantLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
      type(c_LocVariantLink) :: a_ptr_to_neighbor ! Ptr to neighboring LocVariantLink on other side of plane.
    end subroutine F_LocVariantLink_setEdgeConnectivity

    subroutine F_LocVariantLink_setEdgeConnectivityNull(this, a_plane_index) &
      bind(C, name="c_LocVariantLink_setEdgeConnectivityNull")
      import
      implicit none
      type(c_LocVariantLink) :: this
      integer(C_INT), intent(in) :: a_plane_index ! Index for plane of localizer
    end subroutine F_LocVariantLink_setEdgeConnectivityNull
    
    subroutine F_LocVariantLink_printToScreen(this) &
      bind(C, name="c_LocVariantLink_printToScreen")
      import
      implicit none
      type(c_LocVariantLink) :: this
    end subroutine F_LocVariantLink_printToScreen

  end interface

  contains

    subroutine LocVariantLink_class_new(this, a_planar_localizer, a_separator_variant)
      implicit none
      type(LocVariantLink_type), intent(inout) :: this
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(SeparatorVariant_type), intent(in) :: a_separator_variant
      call F_LocVariantLink_new(this%c_object, a_planar_localizer%c_object, a_separator_variant%c_object)
    end subroutine LocVariantLink_class_new

    subroutine LocVariantLink_class_newFromObjectAllocationServer(this, a_object_allocation_server, &
                                                                a_planar_localizer, a_separator_variant)
      implicit none
      type(LocVariantLink_type), intent(inout) :: this
      type(ObjServer_LocVariantLink_type), intent(in) :: a_object_allocation_server
      type(PlanarLoc_type), intent(in) :: a_planar_localizer
      type(SeparatorVariant_type), intent(in) :: a_separator_variant
      call F_LocVariantLink_newFromObjectAllocationServer(this%c_object, &
          a_object_allocation_server%c_object, a_planar_localizer%c_object, a_separator_variant%c_object)
    end subroutine LocVariantLink_class_newFromObjectAllocationServer

    impure elemental subroutine LocVariantLink_class_delete(this)
      implicit none
      type(LocVariantLink_type), intent(in) :: this
      call F_LocVariantLink_delete(this%c_object)
    end subroutine LocVariantLink_class_delete

    subroutine LocVariantLink_class_setId(this, a_id)
      implicit none
      type(LocVariantLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_id
      ! a_id must be >= 0
      call F_LocVariantLink_setId(this%c_object, a_id)
    end subroutine LocVariantLink_class_setId

    function LocVariantLink_class_getId(this) result(a_id)
      implicit none
      type(LocVariantLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t) :: a_id
      ! a_id must be >= 0
      a_id = F_LocVariantLink_getId(this%c_object)
    end function LocVariantLink_class_getId

    subroutine LocVariantLink_class_setEdgeConnectivity(this, a_plane_index, &
        a_neighboring_LocVariantLink)
      implicit none
      type(LocVariantLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      type(LocVariantLink_type) :: a_neighboring_LocVariantLink
      call F_LocVariantLink_setEdgeConnectivity(this%c_object, a_plane_index, a_neighboring_LocVariantLink%c_object)
    end subroutine LocVariantLink_class_setEdgeConnectivity

    subroutine LocVariantLink_class_setEdgeConnectivityNull(this, a_plane_index)
      implicit none
      type(LocVariantLink_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_plane_index
      call F_LocVariantLink_setEdgeConnectivityNull(this%c_object, a_plane_index)
    end subroutine LocVariantLink_class_setEdgeConnectivityNull
    
    subroutine LocVariantLink_class_printToScreen(this)
      implicit none
      type(LocVariantLink_type), intent(in) :: this
      call F_LocVariantLink_printToScreen(this%c_object)
    end subroutine LocVariantLink_class_printToScreen

end module f_LocVariantLink_class
