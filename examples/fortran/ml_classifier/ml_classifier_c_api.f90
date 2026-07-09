module ml_classifier_c_api
  use, intrinsic :: iso_c_binding
  use ml_classifier
  implicit none
contains

  function ml_classifier_fortran(vfrac, liq_bary) result(predicted_class) &
    bind(C, name="ml_classifier_fortran")
    real(c_double), dimension(5, 5, 5), intent(inout) :: vfrac
    real(c_double), dimension(5, 5, 5, 3), intent(inout) :: liq_bary
    integer(c_int) :: predicted_class
    predicted_class = get_class(vfrac, liq_bary)-1  ! Adjust for Fortran 1-based indexing
  end function ml_classifier_fortran

end module ml_classifier_c_api