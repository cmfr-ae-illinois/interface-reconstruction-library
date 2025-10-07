!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_objectallocationserver_SeparatorVariantarator_class.f90
!!
!! This file allows use of the IRL ObjectAllocationServer<SeparatorVariantarator>
!! class through a fortran interface.

!> \brief A fortran type class that allows the creation of
!! IRL's ObjectAllocationServer<SeparatorVariant> class along with enabling
!! some of its methods.
module f_ObjServer_SeparatorVariant_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_ObjServer_SeparatorVariant
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_ObjServer_SeparatorVariant

  type, public :: ObjServer_SeparatorVariant_type
    type(c_ObjServer_SeparatorVariant) :: c_object
  contains
    final :: ObjServer_SeparatorVariant_class_delete
  end type ObjServer_SeparatorVariant_type

  interface new
    module procedure ObjServer_SeparatorVariant_class_new
  end interface

  interface

    subroutine F_ObjServer_SeparatorVariant_new(this, a_number_to_allocate) &
      bind(C, name="c_ObjServer_SeparatorVariant_new")
      import
      implicit none
      type(c_ObjServer_SeparatorVariant) :: this
      integer(C_SIZE_T) :: a_number_to_allocate
    end subroutine F_ObjServer_SeparatorVariant_new

    subroutine F_ObjServer_SeparatorVariant_delete(this) &
      bind(C, name="c_ObjServer_SeparatorVariant_delete")
      import
      implicit none
      type(c_ObjServer_SeparatorVariant) :: this
    end subroutine F_ObjServer_SeparatorVariant_delete

  end interface

  contains

    subroutine ObjServer_SeparatorVariant_class_new(this, a_number_to_allocate)
      implicit none
      type(ObjServer_SeparatorVariant_type), intent(inout) :: this
      integer(IRL_LargeOffsetIndex_t) :: a_number_to_allocate
      call F_ObjServer_SeparatorVariant_new(this%c_object, a_number_to_allocate)
    end subroutine ObjServer_SeparatorVariant_class_new

    impure elemental subroutine ObjServer_SeparatorVariant_class_delete(this)
      implicit none
      type(ObjServer_SeparatorVariant_type), intent(in) :: this
      call F_ObjServer_SeparatorVariant_delete(this%c_object)
    end subroutine ObjServer_SeparatorVariant_class_delete

end module f_ObjServer_SeparatorVariant_class
