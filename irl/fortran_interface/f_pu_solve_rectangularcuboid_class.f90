

!! Interface for PUSolve Class's solve function

module f_PUSolve_RectCub_class
    use f_RectCub_class
    use f_PUSTNeigh_RectCub_class
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
        module procedure PUST_RectCub_class_solveEdge,PUST_rectCub_class_solveEdgeCylinder
    end interface
    interface setNeighborhood
        module procedure PUST_RectCub_class_setNeighborhood
    end interface
    interface setThreshold
        module procedure PUST_RectCub_class_setThreshold
    end interface
    interface getValue
        module procedure PUST_RectCub_class_getValue,PUST_RectCub_class_getValueCylinder
    end interface
    interface getTangent
        module procedure PUST_RectCub_class_getTangent,PUST_RectCub_class_getTangentCylinder
    end interface
    interface getWeight
        module procedure PUST_RectCub_class_getWeight,PUST_RectCub_class_getWeightCylinder
    end interface
    interface printSolver
        module procedure PUST_RectCub_class_printSolver
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
            type(c_PUSTNeigh_RectCub) :: neighborhood
        end subroutine F_PUST_RectCub_setNeighborhood

        subroutine F_PUST_RectCub_setThreshold(this,threshold)&
            bind(C,name = "c_PUST_RectCub_setThreshold")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: threshold
        end subroutine F_PUST_RectCub_setThreshold

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

        subroutine F_PUST_RectCub_solveEdgeCylinder(this,surface_tension_coefficient, start_point,&
            end_point,radius,center,delta, a_force) bind(C, name = "c_PUST_RectCub_solveEdgeCylinder")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: surface_tension_coefficient
            real(C_DOUBLE) :: radius
            real(C_DOUBLE) :: delta
            real(C_DOUBLE), dimension(*), intent(in) :: start_point ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: end_point ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: center ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(out) :: a_force
        end subroutine F_PUST_RectCub_solveEdgeCylinder

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

        subroutine F_PUST_RectCub_getValueCylinder(this,x,y,z,radius,center,value)&
            bind(C, name="c_PUST_RectCub_getValueCylinder")
            import
            implicit none 
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: radius
            real(C_DOUBLE), dimension(*), intent(in) :: center
            real(C_DOUBLE) :: value
        end subroutine F_PUST_RectCub_getValueCylinder
        
        subroutine F_PUST_RectCub_getTangent(this,x,y,z,delta,tangent)&
            bind(C, name="c_PUST_RectCub_getTangent")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE), dimension(*), intent(out) :: tangent
        end subroutine F_PUST_RectCub_getTangent

        subroutine F_PUST_RectCub_getTangentCylinder(this,x,y,z,radius,center,tangent)&
            bind(C, name="c_PUST_RectCub_getTangent")
            import
            implicit none
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: radius
            real(C_DOUBLE), dimension(*), intent(in) :: center
            real(C_DOUBLE), dimension(*), intent(out) :: tangent
        end subroutine F_PUST_RectCub_getTangentCylinder

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

        subroutine F_PUST_RectCub_getWeightCylinder(this,x,y,z,radius,center,weight)&
            bind(C, name="c_PUST_RectCub_getWeight")
            import
            implicit none 
            type(c_PUST_RectCub) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: radius
            real(C_DOUBLE), dimension(*), intent(in) :: center
            real(C_DOUBLE) :: weight
        end subroutine F_PUST_RectCub_getWeightCylinder

        subroutine F_PUST_RectCub_printSolver(this)&
            bind(C, name="c_PUST_RectCub_printSolver")
            import
            implicit none 
            type(c_PUST_RectCub) :: this
        end subroutine F_PUST_RectCub_printSolver
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
        
        subroutine PUST_RectCub_class_setNeighborhood(this,neighborhood)
            implicit none
            type(PUST_RectCub_type), intent(in) :: this
            type(PUSTNeigh_RectCub_type), intent(in) :: neighborhood
            call F_PUST_RectCub_setNeighborhood(this%c_object,neighborhood%c_object)
        end subroutine PUST_RectCub_class_setNeighborhood

        subroutine PUST_RectCub_class_setThreshold(this,threshold)
            implicit none
            type(PUST_RectCub_type), intent(in) :: this
            real(C_DOUBLE) :: threshold
            call F_PUST_RectCub_setThreshold(this%c_object,threshold)
        end subroutine PUST_RectCub_class_setThreshold

        subroutine PUST_RectCub_class_solveEdge(this,surface_tension_coefficient, start_point, end_point, &
            delta,Pressure,Marangoni, a_force)
            implicit none
            type(PUST_RectCub_type) :: this
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

        subroutine PUST_RectCub_class_solveEdgeCylinder(this,surface_tension_coefficient, start_point, end_point,&
            radius,center,delta, a_force)
            implicit none
            type(PUST_RectCub_type) :: this
            real(C_DOUBLE) :: surface_tension_coefficient
            real(C_DOUBLE) :: radius
            real(C_DOUBLE) :: delta
            real(C_DOUBLE), dimension(*), intent(in) :: start_point ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: end_point ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(in) :: center ! dimension(1:3)
            real(C_DOUBLE), dimension(*), intent(out) :: a_force
            call F_PUST_RectCub_solveEdgeCylinder(this%c_object, surface_tension_coefficient, start_point, &
            end_point,radius,center,delta, a_force)
            return
        end subroutine PUST_RectCub_class_solveEdgeCylinder

        subroutine PUST_RectCub_class_getValue(this,x,y,z,delta,value)
            implicit none 
            type(PUST_RectCub_type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE) :: value
            call F_PUST_RectCub_getValue(this%c_object,x,y,z,delta,value)
        end subroutine PUST_RectCub_class_getValue

        subroutine PUST_RectCub_class_getValueCylinder(this,x,y,z,radius,center,value)
            implicit none 
            type(PUST_RectCub_type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: radius
            real(C_DOUBLE), dimension(*), intent(in) :: center
            real(C_DOUBLE) :: value
            call F_PUST_RectCub_getValueCylinder(this%c_object,x,y,z,radius,center,value)
        end subroutine PUST_RectCub_class_getValueCylinder

        subroutine PUST_RectCub_class_getTangent(this,x,y,z,delta,tangent)
            implicit none
            type(PUST_RectCub_type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE), dimension(*), intent(out) :: tangent
            call F_PUST_RectCub_getTangent(this%c_object,x,y,z,delta,tangent)
        end subroutine PUST_RectCub_class_getTangent

        subroutine PUST_RectCub_class_getTangentCylinder(this,x,y,z,radius,center,tangent)
            implicit none
            type(PUST_RectCub_type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: radius
            real(C_DOUBLE), dimension(*), intent(in) :: center
            real(C_DOUBLE), dimension(*), intent(out) :: tangent
            call F_PUST_RectCub_getTangentCylinder(this%c_object,x,y,z,radius,center,tangent)
        end subroutine PUST_RectCub_class_getTangentCylinder

        subroutine PUST_RectCub_class_getWeight(this,x,y,z,delta,weight)
            implicit none
            type(PUST_RectCub_type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: delta
            real(C_DOUBLE) :: weight
            call F_PUST_RectCub_getWeight(this%c_object,x,y,z,delta,weight)
        end subroutine PUST_RectCub_class_getWeight

        subroutine PUST_RectCub_class_getWeightCylinder(this,x,y,z,radius,center,weight)
            implicit none
            type(PUST_RectCub_type) :: this
            real(C_DOUBLE) :: x
            real(C_DOUBLE) :: y
            real(C_DOUBLE) :: z
            real(C_DOUBLE) :: radius
            real(C_DOUBLE), dimension(*), intent(in) :: center
            real(C_DOUBLE) :: weight
            call F_PUST_RectCub_getWeightCylinder(this%c_object,x,y,z,radius,center,weight)
        end subroutine PUST_RectCub_class_getWeightCylinder

        subroutine PUST_RectCub_class_printSolver(this)
            implicit none
            type(PUST_RectCub_type) :: this
            call F_PUST_RectCub_printSolver(this%c_object)
        end subroutine PUST_RectCub_class_printSolver
end module f_PUSolve_RectCub_class