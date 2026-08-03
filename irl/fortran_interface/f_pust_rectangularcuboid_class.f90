!! Interface for PUST Class
module f_PUST_RectCub_class
    use f_RectCub_class
    use f_PUNeigh_RectCub_class
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
    interface solveFace
        module procedure PUST_RectCub_class_solveFace
    end interface
    interface setNeighborhood
        module procedure PUST_RectCub_class_setNeighborhood
    end interface
    interface setKernelSize
        module procedure PUST_RectCub_class_setKernelSize
    end interface
    interface getValue
        module procedure PUST_RectCub_class_getValue
    end interface
    interface getWeight
        module procedure PUST_RectCub_class_getWeight
    end interface
    interface getMeanCurvaturePU
        module procedure PUST_RectCub_class_getMeanCurvature
    end interface
    interface printSolver
        module procedure PUST_RectCub_class_printSolver
    end interface
    interface projectToPU
        module procedure PUST_RectCub_class_projectToPU
    end interface
    interface getNormalPU
        module procedure PUST_RectCub_class_getNormal
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

        subroutine F_PUST_RectCub_setNeighborhood(this,neighborhood)&
            bind(C, name = "c_PUST_RectCub_setNeighborhood")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            type(c_PUNeigh_RectCub) :: neighborhood
        end subroutine F_PUST_RectCub_setNeighborhood

        subroutine F_PUST_RectCub_setKernelSize(this,kernel_size)&
            bind(C, name = "c_PUST_RectCub_setKernelSize")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: kernel_size
        end subroutine F_PUST_RectCub_setKernelSize

        subroutine F_PUST_RectCub_solveEdge(this,surface_tension_coefficient, start_point, end_point,delta,&
             Pressure, Marangoni, a_force) bind(C, name = "c_PUST_RectCub_solveEdge")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: surface_tension_coefficient
            real(C_DOUBLE) :: delta
            real(C_DOUBLE), dimension(*), intent(in) :: start_point ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: end_point ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: Marangoni ! dimension(1:3)
            real(C_DOUBLE) :: Pressure
            real(C_DOUBLE), dimension(*), intent(out) :: a_force
        end subroutine F_PUST_RectCub_solveEdge

        subroutine F_PUST_RectCub_solveFace(this,surface_tension_coefficient, P0, P1,P2,P3,delta,&
             Pressure, Marangoni, a_force) bind(C, name = "c_PUST_RectCub_solveFace")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: surface_tension_coefficient
            real(C_DOUBLE) :: delta
            real(C_DOUBLE), dimension(*), intent(in) :: P0 ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: P1 ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: P2 ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: P3 ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: Marangoni ! dimension(1:3)
            real(C_DOUBLE) :: Pressure
            real(C_DOUBLE), dimension(*), intent(out) :: a_force
        end subroutine F_PUST_RectCub_solveFace

        subroutine F_PUST_RectCub_getValue(this,x,y,z,delta,value)&
            bind(C, name="c_PUST_RectCub_getValue")
            import
            implicit none 
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE) :: value
        end subroutine F_PUST_RectCub_getValue

        subroutine F_PUST_RectCub_getWeight(this,x,y,z,delta,weight)&
            bind(C, name="c_PUST_RectCub_getWeight")
            import
            implicit none 
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE) :: weight
        end subroutine F_PUST_RectCub_getWeight

        subroutine F_PUST_RectCub_getMeanCurvature(this,x,y,z,delta,curv)&
            bind(C, name="c_PUST_RectCub_getMeanCurvature")
            import
            implicit none 
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE) :: curv
        end subroutine F_PUST_RectCub_getMeanCurvature

        subroutine F_PUST_RectCub_printSolver(this)&
            bind(C, name="c_PUST_RectCub_printSolver")
            import
            implicit none 
            type(c_PUST_RectCub) :: this
        end subroutine F_PUST_RectCub_printSolver

        subroutine F_PUST_RectCub_projectToPU(this,P0,delta,Pout)&
            bind(C,name ="c_PUST_RectCub_projectToPU")
            import 
            implicit none
            type(c_PUST_RectCub) :: this 
            real(C_DOUBLE) :: delta
            real(C_DOUBLE), dimension(*), intent(in) :: P0 ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(out) :: Pout ! dimension(1:3)
        end subroutine F_PUST_RectCub_projectToPU

        subroutine F_PUST_RectCub_getNormal(this,x,y,z,delta,normal)&
            bind(C,name ="c_PUST_RectCub_getNormal")
            import 
            implicit none
            type(c_PUST_RectCub) :: this 
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE), dimension(*), intent(out) :: normal ! dimension(1:3)
        end subroutine F_PUST_RectCub_getNormal
    end interface

    contains

        subroutine PUST_RectCub_class_new(this)
            implicit none
            type(PUST_RectCub_Type), intent(inout) :: this
            call F_PUST_RectCub_new(this%c_object)
        end subroutine PUST_RectCub_class_new

        impure elemental subroutine PUST_RectCub_class_delete(this)
            implicit none
            type(PUST_RectCub_Type), intent(in) :: this
            call F_PUST_RectCub_delete(this%c_object)
        end subroutine PUST_RectCub_class_delete
        
        subroutine PUST_RectCub_class_setNeighborhood(this,neighborhood)
            implicit none
            type(PUST_RectCub_Type), intent(in) :: this
            type(PUNeigh_RectCub_type), intent(in) :: neighborhood
            call F_PUST_RectCub_setNeighborhood(this%c_object,neighborhood%c_object)
        end subroutine PUST_RectCub_class_setNeighborhood

        subroutine PUST_RectCub_class_setKernelSize(this,kernel_size)
            implicit none
            type(PUST_RectCub_type), intent(in) :: this
            real(C_DOUBLE) :: kernel_size
            call F_PUST_RectCub_setKernelSize(this%c_object,kernel_size)
        end subroutine PUST_RectCub_class_setKernelSize

        subroutine PUST_RectCub_class_solveEdge(this,surface_tension_coefficient, start_point, end_point, &
            delta,Pressure,Marangoni, a_force)
            implicit none
            type(PUST_RectCub_Type) :: this
            real(C_DOUBLE) :: surface_tension_coefficient
            real(C_DOUBLE) :: delta
            real(IRL_double), dimension(1:3), intent(in) :: start_point
            real(IRL_double), dimension(1:3), intent(in) :: end_point
            real(C_DOUBLE), dimension(*), intent(in) :: Marangoni ! dimension(1:3)
            real(C_DOUBLE) :: Pressure
            real(IRL_double), dimension(1:3), intent(out) :: a_force
            call F_PUST_RectCub_solveEdge(this%c_object, surface_tension_coefficient,start_point,end_point,&
            delta,Pressure,Marangoni, a_force)
            return
        end subroutine PUST_RectCub_class_solveEdge

        subroutine PUST_RectCub_class_solveFace(this,surface_tension_coefficient, P0, P1, &
            P2,P3,delta,Pressure,Marangoni, a_force)
            implicit none
            type(PUST_RectCub_Type) :: this
            real(C_DOUBLE) :: surface_tension_coefficient
            real(C_DOUBLE) :: delta
            real(IRL_double), dimension(1:3), intent(in) :: P0
            real(IRL_double), dimension(1:3), intent(in) :: P1
            real(IRL_double), dimension(1:3), intent(in) :: P2
            real(IRL_double), dimension(1:3), intent(in) :: P3
            real(C_DOUBLE), dimension(*), intent(in) :: Marangoni ! dimension(1:3)
            real(C_DOUBLE) :: Pressure
            real(IRL_double), dimension(1:3), intent(out) :: a_force
            call F_PUST_RectCub_solveFace(this%c_object, surface_tension_coefficient,P0,P1,&
            P2,P3,delta,Pressure,Marangoni, a_force)
            return
        end subroutine PUST_RectCub_class_solveFace

        subroutine PUST_RectCub_class_getValue(this,x,y,z,delta,value)
            implicit none 
            type(PUST_RectCub_Type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE) :: value
            call F_PUST_RectCub_getValue(this%c_object,x,y,z,delta,value)
        end subroutine PUST_RectCub_class_getValue

        subroutine PUST_RectCub_class_getWeight(this,x,y,z,delta,weight)
            implicit none
            type(PUST_RectCub_Type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE) :: weight
            call F_PUST_RectCub_getWeight(this%c_object,x,y,z,delta,weight)
        end subroutine PUST_RectCub_class_getWeight

        subroutine PUST_RectCub_class_getMeanCurvature(this,x,y,z,delta,curv)
            implicit none
            type(PUST_RectCub_Type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE) :: curv
            call F_PUST_RectCub_getMeanCurvature(this%c_object,x,y,z,delta,curv)
        end subroutine PUST_RectCub_class_getMeanCurvature

        subroutine PUST_RectCub_class_printSolver(this)
            implicit none
            type(PUST_RectCub_Type) :: this
            call F_PUST_RectCub_printSolver(this%c_object)
        end subroutine PUST_RectCub_class_printSolver

        subroutine PUST_RectCub_class_projectToPU(this, P0,delta,Pout)
            implicit none
            type(PUST_RectCub_Type) :: this
            real(C_DOUBLE) :: delta
            real(IRL_double), dimension(1:3), intent(in) :: P0
            real(IRL_double), dimension(1:3), intent(out) :: Pout

            call F_PUST_RectCub_projectToPU(this%c_object, P0,delta,Pout)
            return
        end subroutine PUST_RectCub_class_projectToPU

        subroutine PUST_RectCub_class_getNormal(this,x,y,z,delta,normal)
            implicit none
            type(PUST_RectCub_Type) :: this 
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(IRL_double), dimension(*), intent(inout) :: normal ! dimension(1:3)
            call F_PUST_RectCub_getNormal(this%c_object,x,y,z,delta,normal)
        end subroutine PUST_RectCub_class_getNormal
end module f_PUST_RectCub_class