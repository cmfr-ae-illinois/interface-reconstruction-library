!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module f_ReconstructionMetricsInterface
  use f_DefinedTypes
  use f_SeparatorVariant_class
  use f_JibbenNeigh_class
  use f_PUNeigh_RectCub_class
  implicit none

  ! ---------------------------------------------------------------------------------------

  interface reconstructionMetricWithJibben3D
    module procedure reconstructionMetricWithJibben3D_JibbenNeigh
  end interface reconstructionMetricWithJibben3D

  interface angularVarianceMetricWithJibben3D
    module procedure angularVarianceMetricWithJibben3D_JibbenNeigh
  end interface angularVarianceMetricWithJibben3D

  interface volumeErrorSquaredWithJibben3D
    module procedure volumeErrorSquaredWithJibben3D_JibbenNeigh
  end interface volumeErrorSquaredWithJibben3D

  ! ---------------------------------------------------------------------------------------

  interface
  function F_reconstructionMetricWithJibben3D(a_jibben_neighborhood) &
    bind(C, name="c_reconstructionMetricWithJibben3D")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_JibbenNeigh), intent(in) :: a_jibben_neighborhood
      real(C_DOUBLE) :: F_reconstructionMetricWithJibben3D
    end function F_reconstructionMetricWithJibben3D
  end interface

  interface
  function F_angularVarianceMetricWithJibben3D(a_jibben_neighborhood) &
    bind(C, name="c_angularVarianceMetricWithJibben3D")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_JibbenNeigh), intent(in) :: a_jibben_neighborhood
      real(C_DOUBLE) :: F_angularVarianceMetricWithJibben3D
    end function F_angularVarianceMetricWithJibben3D
  end interface

  interface
  function F_volumeErrorSquaredWithJibben3D(a_jibben_neighborhood, a_dx) &
    bind(C, name="c_volumeErrorSquaredWithJibben3D")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(c_JibbenNeigh), intent(in) :: a_jibben_neighborhood
      real(C_DOUBLE), intent(in) :: a_dx
      real(C_DOUBLE) :: F_volumeErrorSquaredWithJibben3D
    end function F_volumeErrorSquaredWithJibben3D
  end interface

  contains

  ! ---------------------------------------------------------------------------------------

  function reconstructionMetricWithJibben3D_JibbenNeigh(a_jibben_neighborhood) result(metric)
    use, intrinsic :: iso_c_binding
    implicit none
    type(JibbenNeigh_type), intent(in) :: a_jibben_neighborhood
    real(IRL_double) :: metric

    metric = F_reconstructionMetricWithJibben3D(a_jibben_neighborhood%c_object)
  end function reconstructionMetricWithJibben3D_JibbenNeigh

  function angularVarianceMetricWithJibben3D_JibbenNeigh(a_jibben_neighborhood) result(metric)
    use, intrinsic :: iso_c_binding
    implicit none
    type(JibbenNeigh_type), intent(in) :: a_jibben_neighborhood
    real(IRL_double) :: metric

    metric = F_angularVarianceMetricWithJibben3D(a_jibben_neighborhood%c_object)
  end function angularVarianceMetricWithJibben3D_JibbenNeigh

  function volumeErrorSquaredWithJibben3D_JibbenNeigh(a_jibben_neighborhood, a_dx) result(metric)
    use, intrinsic :: iso_c_binding
    implicit none
    type(JibbenNeigh_type), intent(in) :: a_jibben_neighborhood
    real(IRL_double), intent(in) :: a_dx
    real(IRL_double) :: metric

    metric = F_volumeErrorSquaredWithJibben3D(a_jibben_neighborhood%c_object, a_dx)
  end function volumeErrorSquaredWithJibben3D_JibbenNeigh
  
end module f_ReconstructionMetricsInterface
