

!! Interface for PUSTNeighborhood Class

module f_PUSTNeigh_RectCub_class
    use f_RectCub_class
    use f_PlanarSep_class
    use, intrinsic :: iso_c_binding
    use f_DefinedTypes
    implicit none

    type, public, bind(C) :: c_PUSTNeigh_RectCub
        type(C_PTR), private :: object = C_NULL_PTR
    end type c_PUSTNeigh_RectCub

    type, public :: PUSTNeigh_RectCub_type
        type(c_PUSTNeigh_RectCub) :: c_object 
    contains
        final :: PUSTSNeigh_RectCub_class_delete
    end type PUSTNeigh_RectCub_type

    interface new
        module procedure PUSTNeigh_RectCub_class_new
    end interface
    interface setSize
        module procedure PUSTNeigh_RectCub_class_setSize
    end interface
    interface setMember
        module procedure PUSTNeigh_RectCub_class_setMember
    end interface
    interface addMember
        module procedure PUSTNeigh_RectCub_class_addMember
    end interface
    interface emptyNeighborhood
        module procedure PUSTNeigh_RectCub_class_emptyNeighborhood
    end interface
    interface setCenterOfStencil
        module procedure PUSTNeigh_RectCub_class_setCenterOfStencil
    end interface

  interface
    
    subroutine F_PUSTNeigh_RectCub_new(this) &
        bind(C,name = "c_PUSTNeigh_RectCub_new")
        import
        implicit none
        type(c_PUSTNeigh_RectCub) :: this
    end subroutine F_PUSTNeigh_RectCub_new

    subroutine F_PUSTNeigh_RectCub_delete(this) &
        bind(C, name ="c_PUSTNeigh_RectCub_delete")
        import
        implicit none
        type(c_PUSTNeigh_RectCub) :: this
    end subroutine F_PUSTNeigh_RectCub_delete

    subroutine F_PUSTNeigh_RectCub_setSize(this,a_size) &
        bind(C, name="c_PUSTNeigh_RectCub_setSize")
        import
        implicit none
        type(c_PUSTNeigh_RectCub) :: this
        integer(C_INT) :: a_size
    end subroutine F_PUSTNeigh_RectCub_setSize

    subroutine F_PUSTNeigh_RectCub_setMember(this,a_index, a_rectangular_cuboid, &
            a_planar_separator) &
        bind(C, name = "c_PUSTNeigh_RectCub_setMember")
        import
        implicit none
        integer(C_INT), intent(in) :: a_index
        type(c_PUSTNeigh_RectCub) :: this
        type(c_RectCub) :: a_rectangular_cuboid
        type(c_PlanarSep) :: a_planar_separator
    end subroutine F_PUSTNeigh_RectCub_setMember

    subroutine F_PUSTNeigh_RectCub_addMember(this,a_rectangular_cuboid, &
            a_planar_separator) &
        bind(C, name = "c_PUSTNeigh_RectCub_addMember")
        import
        implicit none 
        type(c_PUSTNeigh_RectCub) :: this
        type(c_RectCub) :: a_rectangular_cuboid
        type(c_PlanarSep) :: a_volume_fraction
    end subroutine F_PUSTNeigh_RectCub_addMember

    subroutine F_PUSTNeigh_RectCub_emptyNeighborhood(this) &
        bind(C, name="c_PUSTNeigh_RectCub_emptyNeighborhood")
        import 
        implicit none 
        type(c_PUSTNeigh_RectCub) :: this
    end subroutine F_PUSTNeigh_RectCub_emptyNeighborhood

    subroutine F_PUSTNeigh_RectCub_setCenterOfStencil(this,a_center_cell_index) &
        bind(C, name="c_PUSTNeigh_RectCub_setCenterOfStencil")
        import 
        implicit none
        type(c_PUSTNeigh_RectCub) :: this
        integer(C_INT) :: a_center_cell_index
    end subroutine F_PUSTNeigh_RectCub_setCenterOfStencil

  end interface

  contains


    subroutine PUSTNeigh_RectCub_class_new(this)
        implicit none
        type(PUSTNeigh_RectCub_type), intent(inout) :: this
        call F_PUSTNeigh_RectCub_new(this%c_object)
    end subroutine PUSTNeigh_RectCub_class_new

    impure elemental subroutine PUSTNeigh_RectCub_class_delete(this)
        implicit none
        type(PUSTNeigh_RectCub_type), intent(in) :: this
        call F_PUSTNeigh_RectCub_delete(this%c_object)
    end subroutine PUSTNeigh_RectCub_class_delete

    subroutine PUSTNeigh_RectCub_class_setSize(this, a_size)
        implicit none
        type(PUSTNeigh_RectCub_type), intent(in) :: this
        integer(IRL_UnsignedIndex_t), intent(in) :: a_size
        call F_PUSTNeigh_RectCub_setSize(this%c_object,a_size)
    end subroutine PUSTNeigh_RectCub_class_setSize

    subroutine PUSTNeigh_RectCub_class_setMember(this,a_index, a_rectangular_cuboid, &
            a_planar_separator)
        implicit none
        type(PUSTNeigh_RectCub_type), intent(in) :: this
        integer(IRL_SignedIndex_t), intent(in) :: a_index
        type(RectCub_type), intent(in) :: a_rectangular_cuboid
        type(PlanarSep_type), intent(in) :: a_planar_separator
        call F_PUSTNeigh_RectCub_setMember(this%c_object,a_index,a_rectangular_cuboid%c_object, &
            a_planar_separator%c_object)
    end subroutine PUSTNeigh_RectCub_class_setMember

    subroutine PUSTNeigh_RectCub_class_addMember(this, a_rectangular_cuboid, &
            a_planar_separator)
        implicit none
        type(PUSTNeigh_RectCub_type), intent(in) :: this
        type(RectCub_type), intent(in) :: a_rectangular_cuboid
        type(PlanarSep_type), intent(in) :: a_planar_separator
        call F_PUSTNeigh_RectCub_addMember(this%c_object, a_rectangular_cuboid%c_object, &
            a_planar_separator%c_object)
    end subroutine PUSTNeigh_RectCub_class_addMember 

    subroutine PUSTNeigh_RectCub_class_emptyNeighborhood(this)
        implicit none
        type(PUSTNeigh_RectCub_type), intent(in) ::this
        call F_PUSTNeigh_RectCub_emptyNeighborhood(this%c_object)
    end subroutine PUSTNeigh_RectCub_class_emptyNeighborhood

    subroutine PUSTNeigh_RectCub_class_setCenterOfStencil(this,a_center_cell_index)
        implicit none
        type(PUSTNeigh_RectCub_type), intent(in) :: this
        integer(IRL_UnsignedIndex_t), intent(in) :: a_center_cell_index
        call F_PUSTNeigh_RectCub_setCenterOfStencil(this%c_object, a_center_cell_index)
    end subroutine PUSTNeigh_RectCub_class_setCenterOfStencil

end module f_PUSTNeigh_RectCub_class