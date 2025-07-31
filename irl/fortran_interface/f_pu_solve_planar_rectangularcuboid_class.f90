

!! Interface for PUSolve Class's solve function

module f_PUSolve_RectCub_class
    use f_RectCub_class
    use f_SeparatorVariant_class
    use, intrinsic :: iso_c_binding
    use f_DefinedTypes
    implicit none

    type, public, bind(C) :: c_PUST_RectCub
        type(C_PTR), private :: object = C_NULL_PTR
    end type c_PUST_RectCub

    type, public :: PUST_RectCub_type
        type(c_PUST_RectCub) :: c_object
    contains
        final :: PUST_RectCub_class_delete
    end type PUST_RectCub_type

    interface new
        module procedure PUST_RectCub_class_new
    end interface
    interface solveEdge
        module procedure PUST_RectCub_class_solveEdge
    end interface

    interface

        subroutine F_PUST_RectCub_new(this) &
            bind(C, name = "c_PUST_RectCub_new")
            import
            implicit none
            type(c_PUST_RectCub) :: this
        end subroutine F_PUST_RectCub_new

        subroutine F_PUST_RectCub_delete(this) & 
            bind(C, name = "c_PUST_RectCub_delete")
            import
            implicit none
            type(c_PUST_RectCub) :: this
        end subroutine F_PUST_RectCub_delete

        subroutine F_PUST_RectCub_solveEdge(this,surface_tension_coefficient, start_point, end_point, a_force)
            bind(C, name = "c_PUST_RectCub_solveEdge")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: surface_tension_coefficient
            real(C_DOUBLE), dimension(*), intent(in) :: start_point ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: end_point ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(out) :: a_force
        end subroutine F_PUST_RectCub_solveEdge
    
    end interface

    contains

        subroutine PUST_RectCub_class_new(this)
            implicit none
            type(PUST_RectCub_type), intent(inout) :: this
            call F_PUST_RectCub_new(this%c_object)
        end subroutine PUST_RectCub_class_new

        impure elemental subroutine PUST_RectCub_class_delete(this)
            implicit none
            type(PUST_RectCub_type), intent(in) :: this
            call F_PUST_RectCub_delete(this%c_object)
        end subroutine PUST_RectCub_class_delete

        subroutine PUST_RectCub_class_solveEdge(this,surface_tension_coefficient, start_point, end_point, a_force)
            implicit none
            type(PUST_RectCub_type) :: this
            real(C_DOUBLE) :: surface_tension_coefficient
            real(IRL_double), dimension(1:3), intent(out) :: start_point
            real(IRL_double), dimension(1:3), intent(out) :: end_point
            real(IRL_double), dimension(1:3) :: a_force
            call F_PUST_RectCub_solveEdge(this%c_object, surface_tension_coefficient,start_point,end_point, a_force)
            return
        end subroutine PUST_RectCub_class_solveEdge

end module f_PUSolve_RectCub_class