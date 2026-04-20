

!! Interface for PUNeighborhood Class

module f_PUNeigh_RectCub_class
    use f_RectCub_class
    use f_SeparatorVariant_class
    use, intrinsic :: iso_c_binding
    use f_DefinedTypes
    implicit none

    type, public, bind(C) :: c_PUNeigh_RectCub
        type(C_PTR), private :: object = C_NULL_PTR
    end type c_PUNeigh_RectCub

    type, public :: PUNeigh_RectCub_type
        type(c_PUNeigh_RectCub) :: c_object 
    contains
        final :: PUNeigh_RectCub_class_delete
    end type PUNeigh_RectCub_type

    interface new
        module procedure PUNeigh_RectCub_class_new
    end interface
    interface setSize
        module procedure PUNeigh_RectCub_class_setSize
    end interface
    interface reserve
        module procedure PUNeigh_RectCub_class_reserve
    end interface
    interface setMember
        module procedure PUNeigh_RectCub_class_setMember
    end interface
    interface addMember
        module procedure PUNeigh_RectCub_class_addMember
    end interface
    interface emptyNeighborhood
        module procedure PUNeigh_RectCub_class_emptyNeighborhood
    end interface

  interface
    
    subroutine F_PUNeigh_RectCub_new(this) &
        bind(C,name = "c_PUNeigh_RectCub_new")
        import
        implicit none
        type(c_PUNeigh_RectCub) :: this
    end subroutine F_PUNeigh_RectCub_new

    subroutine F_PUNeigh_RectCub_delete(this) &
        bind(C, name ="c_PUNeigh_RectCub_delete")
        import
        implicit none
        type(c_PUNeigh_RectCub) :: this
    end subroutine F_PUNeigh_RectCub_delete

    subroutine F_PUNeigh_RectCub_setSize(this,a_size) &
        bind(C, name="c_PUNeigh_RectCub_setSize")
        import
        implicit none
        type(c_PUNeigh_RectCub) :: this
        integer(C_INT) :: a_size
    end subroutine F_PUNeigh_RectCub_setSize

    subroutine F_PUNeigh_RectCub_reserve(this,a_size) &
        bind(C, name="c_PUNeigh_RectCub_reserve")
        import
        implicit none
        type(c_PUNeigh_RectCub) :: this
        integer(C_INT) :: a_size
    end subroutine F_PUNeigh_RectCub_reserve

    subroutine F_PUNeigh_RectCub_setMember(this,a_index, a_centroid, a_weight, &
            a_separator) &
        bind(C, name = "c_PUNeigh_RectCub_setMember")
        import
        implicit none
        integer(C_INT), intent(in) :: a_index
        type(c_PUNeigh_RectCub) :: this
        real(C_DOUBLE), dimension(*) :: a_centroid ! dimension(1:3)
        real(C_DOUBLE) :: a_weight
        type(c_SeparatorVariant) :: a_separator
    end subroutine F_PUNeigh_RectCub_setMember

    subroutine F_PUNeigh_RectCub_addMember(this,a_centroid, a_weight, &
            a_separator, a_scalar) &
        bind(C, name = "c_PUNeigh_RectCub_addMember")
        import
        implicit none 
        type(c_PUNeigh_RectCub) :: this
        real(C_DOUBLE), dimension(*) :: a_centroid ! dimension(1:3)
        real(C_DOUBLE) :: a_weight
        type(c_SeparatorVariant) :: a_separator
        real(C_DOUBLE) :: a_scalar
    end subroutine F_PUNeigh_RectCub_addMember

    subroutine F_PUNeigh_RectCub_emptyNeighborhood(this) &
        bind(C, name="c_PUNeigh_RectCub_emptyNeighborhood")
        import 
        implicit none 
        type(c_PUNeigh_RectCub) :: this
    end subroutine F_PUNeigh_RectCub_emptyNeighborhood

  end interface

  contains


    subroutine PUNeigh_RectCub_class_new(this)
        implicit none
        type(PUNeigh_RectCub_type), intent(inout) :: this
        call F_PUNeigh_RectCub_new(this%c_object)
    end subroutine PUNeigh_RectCub_class_new

    impure elemental subroutine PUNeigh_RectCub_class_delete(this)
        implicit none
        type(PUNeigh_RectCub_type), intent(in) :: this
        call F_PUNeigh_RectCub_delete(this%c_object)
    end subroutine PUNeigh_RectCub_class_delete

    subroutine PUNeigh_RectCub_class_setSize(this, a_size)
        implicit none
        type(PUNeigh_RectCub_type), intent(in) :: this
        integer(IRL_UnsignedIndex_t), intent(in) :: a_size
        call F_PUNeigh_RectCub_setSize(this%c_object,a_size)
    end subroutine PUNeigh_RectCub_class_setSize

    subroutine PUNeigh_RectCub_class_reserve(this, a_size)
        implicit none
        type(PUNeigh_RectCub_type), intent(in) :: this
        integer(IRL_UnsignedIndex_t), intent(in) :: a_size
        call F_PUNeigh_RectCub_reserve(this%c_object,a_size)
    end subroutine PUNeigh_RectCub_class_reserve

    subroutine PUNeigh_RectCub_class_setMember(this,a_index, a_centroid, a_weight, &
            a_separator)
        implicit none
        type(PUNeigh_RectCub_type), intent(in) :: this
        integer(IRL_SignedIndex_t), intent(in) :: a_index
        real(IRL_double), dimension(1:3), intent(in) :: a_centroid
        real(IRL_double), intent(in) :: a_weight
        type(SeparatorVariant_type), intent(in) :: a_separator
        call F_PUNeigh_RectCub_setMember(this%c_object,a_index,a_centroid,a_weight, &
            a_separator%c_object)
    end subroutine PUNeigh_RectCub_class_setMember

    subroutine PUNeigh_RectCub_class_addMember(this, a_centroid, a_weight, &
            a_separator,a_scalar)
        implicit none
        type(PUNeigh_RectCub_type), intent(in) :: this
        real(IRL_double), dimension(1:3), intent(in) :: a_centroid
        real(IRL_double), intent(in) :: a_weight
        type(SeparatorVariant_type), intent(in) :: a_separator
        real(IRL_double), intent(in) :: a_scalar
        call F_PUNeigh_RectCub_addMember(this%c_object, a_centroid, a_weight, &
            a_separator%c_object,a_scalar)
    end subroutine PUNeigh_RectCub_class_addMember 

    subroutine PUNeigh_RectCub_class_emptyNeighborhood(this)
        implicit none
        type(PUNeigh_RectCub_type), intent(in) ::this
        call F_PUNeigh_RectCub_emptyNeighborhood(this%c_object)
    end subroutine PUNeigh_RectCub_class_emptyNeighborhood


end module f_PUNeigh_RectCub_class