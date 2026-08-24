
module amrex_sepunionmfab_module

  use iso_c_binding
  use amrex_error_module
  use amrex_fort_module, only : ndims => amrex_spacedim, amrex_real
  use amrex_box_module
  use amrex_boxarray_module
  use amrex_distromap_module
  use amrex_geometry_module
  use amrex_string_module
  use amrex_omp_module
  use amrex_multifab_module
  use irl_fortran_interface
  
  implicit none

  private

  public :: amrex_sepunionmfab_destroy ! List first to avoid XL compiler bug
  public :: amrex_sepunionmfab_build, amrex_sepunionmfab_rebuild
  public :: amrex_sepunionmfab_write, amrex_sepunionmfab_read
  public :: amrex_sepunionmfab_build_alias

  type, public   :: amrex_sepunionmfab
     logical               :: owner = .false.
     type   (c_ptr)        :: p     =  c_null_ptr
     integer(c_int)        :: nc    =  0
     integer(c_int)        :: ng(3) =  0
     type(amrex_boxarray)  :: ba
     type(amrex_distromap) :: dm
   contains
     generic   :: assignment(=) => amrex_sepunionmfab_assign, amrex_sepunionmfab_install  ! shallow copy
     procedure :: move          => amrex_sepunionmfab_move     ! transfer ownership
     procedure :: ncomp         => amrex_sepunionmfab_ncomp
     procedure :: nghost        => amrex_sepunionmfab_nghost
     procedure :: nghostvect    => amrex_sepunionmfab_nghost_vect
     procedure :: nodal_type    => amrex_sepunionmfab_nodal_type   ! get index type
     generic   :: dataPtr       => amrex_sepunionmfab_dataptr_iter, amrex_sepunionmfab_dataptr_int
     generic   :: copy          => amrex_sepunionmfab_copy, amrex_sepunionmfab_copy_cgv ! This copies the data
     generic   :: parallel_copy => amrex_sepunionmfab_parallel_copy, amrex_sepunionmfab_parallel_copy_c, &
          amrex_sepunionmfab_parallel_copy_cg, amrex_sepunionmfab_parallel_copy_cgv
     generic   :: fill_boundary => amrex_sepunionmfab_fill_boundary, amrex_sepunionmfab_fill_boundary_c
     procedure, private :: amrex_sepunionmfab_copy
     procedure, private :: amrex_sepunionmfab_copy_cgv
     procedure, private :: amrex_sepunionmfab_fill_boundary
     procedure, private :: amrex_sepunionmfab_fill_boundary_c
     procedure, private :: amrex_sepunionmfab_parallel_copy
     procedure, private :: amrex_sepunionmfab_parallel_copy_c
     procedure, private :: amrex_sepunionmfab_parallel_copy_cg
     procedure, private :: amrex_sepunionmfab_parallel_copy_cgv
     procedure, private :: amrex_sepunionmfab_assign
     procedure, private :: amrex_sepunionmfab_install
     procedure, private :: amrex_sepunionmfab_dataptr_iter
     procedure, private :: amrex_sepunionmfab_dataptr_int
     final :: amrex_sepunionmfab_destroy
  end type amrex_sepunionmfab

  interface amrex_sepunionmfab_build
    module procedure amrex_sepunionmfab_build_s
    module procedure amrex_sepunionmfab_build_a
  end interface amrex_sepunionmfab_build

#ifdef __NVCOMPILER
  interface amrex_sepunionmfab_destroy
    module procedure amrex_sepunionmfab_destroy
  end interface amrex_sepunionmfab_destroy
#endif

  interface
     subroutine amrex_fi_new_sepunionmfab (mf,ba,dm,nc,ng,nodal) bind(c)
       import
       implicit none
       type(c_ptr) :: mf, ba, dm
       integer(c_int), value :: nc
       integer(c_int), intent(in) :: ng(3), nodal(3)
     end subroutine amrex_fi_new_sepunionmfab

     subroutine amrex_fi_new_sepunionmfab_alias (mf, srcmf, comp, ncomp) bind(c)
       import
       implicit none
       type(c_ptr) :: mf
       type(c_ptr), value :: srcmf
       integer(c_int), value :: comp, ncomp
     end subroutine amrex_fi_new_sepunionmfab_alias

     subroutine amrex_fi_delete_sepunionmfab (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end subroutine amrex_fi_delete_sepunionmfab

     integer(c_int) function amrex_fi_sepunionmfab_ncomp (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end function amrex_fi_sepunionmfab_ncomp

     subroutine amrex_fi_sepunionmfab_ngrow (mf, ngv) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
       integer(c_int) :: ngv(3)
     end subroutine amrex_fi_sepunionmfab_ngrow

     type(c_ptr) function amrex_fi_sepunionmfab_boxarray (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end function amrex_fi_sepunionmfab_boxarray

     type(c_ptr) function amrex_fi_sepunionmfab_distromap (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end function amrex_fi_sepunionmfab_distromap

     subroutine amrex_fi_sepunionmfab_dataptr_iter (mf, mfi, dp, lo, hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, mfi
       type(c_ptr) :: dp
       integer(c_int) :: lo(3), hi(3)
     end subroutine amrex_fi_sepunionmfab_dataptr_iter

     subroutine amrex_fi_sepunionmfab_dataptr_int (mf, igrd, dp, lo, hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
       integer, value :: igrd
       type(c_ptr) :: dp
       integer(c_int) :: lo(3), hi(3)
     end subroutine amrex_fi_sepunionmfab_dataptr_int

    !  subroutine amrex_fi_sepunionmfab_setval (mf, val, ic, nc, ng) bind(c)
    !    import
    !    implicit none
    !    type(c_ptr), value :: mf
    !    real(amrex_real), value :: val
    !    integer(c_int), value :: ic, nc
    !    integer(c_int), intent(in) :: ng(*)
    !  end subroutine amrex_fi_sepunionmfab_setval

     subroutine amrex_fi_sepunionmfab_copy (dstmf, srcmf, srccomp, dstcomp, nc, ng) bind(c)
       import
       implicit none
       type(c_ptr), value :: dstmf, srcmf
       integer(c_int), value :: srccomp, dstcomp, nc
       integer(c_int), intent(in) :: ng(*)
     end subroutine amrex_fi_sepunionmfab_copy

     subroutine amrex_fi_sepunionmfab_parallelcopy(dstmf, srcmf, srccomp, dstcomp, nc,&
          srcng, dstng, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: dstmf, srcmf, geom
       integer(c_int), value :: srccomp, dstcomp, nc, srcng, dstng
     end subroutine amrex_fi_sepunionmfab_parallelcopy

     subroutine amrex_fi_sepunionmfab_parallelcopy_gv(dstmf, srcmf, srccomp, dstcomp, nc,&
          srcng, dstng, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: dstmf, srcmf, geom
       integer(c_int), value :: srccomp, dstcomp, nc
       integer, intent(in) :: srcng(*), dstng(*)
     end subroutine amrex_fi_sepunionmfab_parallelcopy_gv

     subroutine amrex_fi_sepunionmfab_fill_boundary (mf, geom, c, nc, cross) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, geom
       integer(c_int), value :: c, nc, cross
     end subroutine amrex_fi_sepunionmfab_fill_boundary

     subroutine amrex_fi_write_sepunionmfab (mf, name) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
       character(kind=c_char), intent(in) :: name(*)
     end subroutine amrex_fi_write_sepunionmfab

     subroutine amrex_fi_read_sepunionmfab (mf, name) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
       character(kind=c_char), intent(in) :: name(*)
     end subroutine amrex_fi_read_sepunionmfab
  end interface

contains

  subroutine amrex_sepunionmfab_build_a (mf, ba, dm, nc, ng, nodal)
    type(amrex_sepunionmfab), intent(inout) :: mf
    type(amrex_boxarray), intent(in )   :: ba
    type(amrex_distromap),intent(in )   :: dm
    integer, intent(in) :: nc, ng(*)
    logical, intent(in), optional :: nodal(*)
    integer :: inodal(3), dir
    mf%owner = .true.
    mf%nc = nc
    mf%ng(1:ndims) = ng(1:ndims)
    inodal = 0
    if (present(nodal)) then
       do dir = 1, ndims
          if (nodal(dir)) inodal(dir) = 1
       end do
    end if
    mf%ba = ba
    mf%dm = dm
    call amrex_fi_new_sepunionmfab(mf%p, mf%ba%p, mf%dm%p, mf%nc, mf%ng, inodal)
  end subroutine amrex_sepunionmfab_build_a

  subroutine amrex_sepunionmfab_build_s (mf, ba, dm, nc, ng, nodal)
    type(amrex_sepunionmfab), intent(inout) :: mf
    type(amrex_boxarray), intent(in )   :: ba
    type(amrex_distromap),intent(in )   :: dm
    integer, intent(in) :: nc, ng
    logical, intent(in), optional :: nodal(*)
    call amrex_sepunionmfab_build_a(mf, ba, dm, nc, (/ng,ng,ng/), nodal)
  end subroutine amrex_sepunionmfab_build_s

  subroutine amrex_sepunionmfab_build_alias (mf, srcmf, comp, ncomp)
    type(amrex_sepunionmfab), intent(inout) :: mf
    type(amrex_sepunionmfab), intent(in   ) :: srcmf
    integer, intent(in) :: comp, ncomp
    call amrex_sepunionmfab_destroy(mf)
    mf%owner = .true.
    mf%nc    = ncomp
    mf%ng    = srcmf%ng
    mf%ba    = srcmf%ba
    mf%dm    = srcmf%dm
    call amrex_fi_new_sepunionmfab_alias(mf%p, srcmf%p, comp-1, ncomp)
  end subroutine amrex_sepunionmfab_build_alias

   !> Rebuild a sepunionmfab and set to zero
   subroutine amrex_sepunionmfab_rebuild (mf, ba, dm, nc, ng)
      implicit none
      type(amrex_sepunionmfab), intent(inout) :: mf
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      integer, intent(in) :: nc, ng
      call amrex_sepunionmfab_destroy(mf)
      call amrex_sepunionmfab_build(mf,ba,dm,nc,ng)
      ! call mf%setval(0.0_amrex_real) %%% TODO
   end subroutine amrex_sepunionmfab_rebuild

  impure elemental subroutine amrex_sepunionmfab_destroy (this)
    type(amrex_sepunionmfab), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_sepunionmfab(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
    call amrex_boxarray_destroy(this%ba)
    call amrex_distromap_destroy(this%dm)
  end subroutine amrex_sepunionmfab_destroy

  subroutine amrex_sepunionmfab_assign (dst, src)
    class(amrex_sepunionmfab), intent(inout) :: dst
    type (amrex_sepunionmfab), intent(in   ) :: src
    call amrex_sepunionmfab_destroy(dst)
    dst%owner = .false.
    dst%p     = src%p
    dst%nc    = src%nc
    dst%ng    = src%ng
    dst%ba    = src%ba
    dst%dm    = src%dm
  end subroutine amrex_sepunionmfab_assign

  subroutine amrex_sepunionmfab_install (this, p)
    class(amrex_sepunionmfab), intent(inout) :: this
    type(c_ptr), intent(in) :: p
    this%owner = .false.
    this%p     = p
    this%nc    = amrex_fi_sepunionmfab_ncomp(p)
    call amrex_fi_sepunionmfab_ngrow(p, this%ng)
    this%ba    = amrex_fi_sepunionmfab_boxarray(p)
    this%dm    = amrex_fi_sepunionmfab_distromap(p)
  end subroutine amrex_sepunionmfab_install

  subroutine amrex_sepunionmfab_move (dst, src)
    class(amrex_sepunionmfab), intent(inout) :: dst
    type (amrex_sepunionmfab), intent(inout) :: src
    call amrex_sepunionmfab_destroy(dst)
    dst%owner = src%owner
    dst%p     = src%p
    dst%nc    = src%nc
    dst%ng    = src%ng
    call dst%ba%move(src%ba)
    call dst%dm%move(src%dm)
    src%owner = .false.
    src%p     = c_null_ptr
  end subroutine amrex_sepunionmfab_move

  subroutine amrex_sepunionmfab_swap(mf1, mf2)
    type(amrex_sepunionmfab), intent(inout) :: mf1, mf2
    type(amrex_sepunionmfab) :: mftmp
    call mftmp%move(mf1)
    call mf1%move(mf2)
    call mf2%move(mftmp)
    call amrex_sepunionmfab_destroy(mftmp)
  end subroutine amrex_sepunionmfab_swap

  pure integer function amrex_sepunionmfab_ncomp (this)
    class(amrex_sepunionmfab), intent(in) :: this
    amrex_sepunionmfab_ncomp = this%nc
  end function amrex_sepunionmfab_ncomp

  pure integer function amrex_sepunionmfab_nghost (this)
    class(amrex_sepunionmfab), intent(in) :: this
    amrex_sepunionmfab_nghost = this%ng(1)
  end function amrex_sepunionmfab_nghost

  pure function amrex_sepunionmfab_nghost_vect (this) result(ngv)
    class(amrex_sepunionmfab), intent(in) :: this
    integer, dimension(3) :: ngv
    ngv = this%ng
  end function amrex_sepunionmfab_nghost_vect

  pure function amrex_sepunionmfab_nodal_type (this) result(nodal)
    class(amrex_sepunionmfab), intent(in) :: this
    logical, dimension(3) :: nodal
    nodal = this%ba%nodal_type()
  end function amrex_sepunionmfab_nodal_type

  function amrex_sepunionmfab_dataPtr_iter (this, mfi) result(dp)
    class(amrex_sepunionmfab), intent(in) :: this
    type(amrex_mfiter), intent(in) :: mfi
    type(SeparatorUnion_type_raw), contiguous, pointer, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    type(SeparatorUnion_type_raw), contiguous, pointer :: fp(:,:,:,:)
    integer(c_int) :: n(4)
    type(amrex_box) :: bx
    call amrex_fi_sepunionmfab_dataptr_iter(this%p, mfi%p, cp, bx%lo, bx%hi)
    n(1:3) = bx%hi - bx%lo + 1
    n(4)   = this%ncomp()
    call c_f_pointer(cp, fp, shape=n)
    dp(bx%lo(1):,bx%lo(2):,bx%lo(3):,1:) => fp
  end function amrex_sepunionmfab_dataPtr_iter

  function amrex_sepunionmfab_dataPtr_int (this, gid) result(dp)
    class(amrex_sepunionmfab), intent(in) :: this
    integer, intent(in) :: gid
    type(SeparatorUnion_type_raw), contiguous, pointer, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    type(SeparatorUnion_type_raw), contiguous, pointer :: fp(:,:,:,:)
    integer(c_int) :: n(4)
    type(amrex_box) :: bx
    call amrex_fi_sepunionmfab_dataptr_int(this%p, gid, cp, bx%lo, bx%hi)
    n(1:3) = bx%hi - bx%lo + 1
    n(4)   = this%ncomp()
    call c_f_pointer(cp, fp, shape=n)
    dp(bx%lo(1):,bx%lo(2):,bx%lo(3):,1:) => fp
  end function amrex_sepunionmfab_dataPtr_int

  subroutine amrex_sepunionmfab_copy (this, srcmf, srccomp, dstcomp, nc, ng)
    class(amrex_sepunionmfab) :: this
    type(amrex_sepunionmfab), intent(in) :: srcmf
    integer, intent(in) :: srccomp, dstcomp, nc, ng
    call amrex_fi_sepunionmfab_copy(this%p, srcmf%p, srccomp-1, dstcomp-1, nc, (/ng,ng,ng/))
  end subroutine amrex_sepunionmfab_copy

  subroutine amrex_sepunionmfab_copy_cgv (this, srcmf, srccomp, dstcomp, nc, ng)
    class(amrex_sepunionmfab) :: this
    type(amrex_sepunionmfab), intent(in) :: srcmf
    integer, intent(in) :: srccomp, dstcomp, nc, ng(*)
    call amrex_fi_sepunionmfab_copy(this%p, srcmf%p, srccomp-1, dstcomp-1, nc, ng)
  end subroutine amrex_sepunionmfab_copy_cgv

  subroutine amrex_sepunionmfab_parallel_copy (this, srcmf, geom)
    class(amrex_sepunionmfab) :: this
    type(amrex_sepunionmfab), intent(in) :: srcmf
    type(amrex_geometry), intent(in) :: geom
    call amrex_fi_sepunionmfab_parallelcopy(this%p, srcmf%p, 0, 0, this%nc, 0, 0, geom%p)
  end subroutine amrex_sepunionmfab_parallel_copy

  subroutine amrex_sepunionmfab_parallel_copy_c (this, srcmf, srccomp, dstcomp, nc, geom)
    class(amrex_sepunionmfab) :: this
    type(amrex_sepunionmfab), intent(in) :: srcmf
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: srccomp, dstcomp, nc
    call amrex_fi_sepunionmfab_parallelcopy(this%p, srcmf%p, srccomp-1, dstcomp-1, nc, 0, 0, geom%p)
  end subroutine amrex_sepunionmfab_parallel_copy_c

  subroutine amrex_sepunionmfab_parallel_copy_cg (this, srcmf, srccomp, dstcomp, nc, srcng, dstng, geom)
    class(amrex_sepunionmfab) :: this
    type(amrex_sepunionmfab), intent(in) :: srcmf
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: srccomp, dstcomp, nc, srcng, dstng
    call amrex_fi_sepunionmfab_parallelcopy(this%p, srcmf%p, srccomp-1, dstcomp-1, nc, srcng, dstng, geom%p)
  end subroutine amrex_sepunionmfab_parallel_copy_cg

  subroutine amrex_sepunionmfab_parallel_copy_cgv (this, srcmf, srccomp, dstcomp, nc, srcng, dstng, geom)
    class(amrex_sepunionmfab) :: this
    type(amrex_sepunionmfab), intent(in) :: srcmf
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: srccomp, dstcomp, nc
    integer, intent(in) :: srcng(*), dstng(*)
    call amrex_fi_sepunionmfab_parallelcopy_gv(this%p, srcmf%p, srccomp-1, dstcomp-1, nc, srcng, dstng, geom%p)
  end subroutine amrex_sepunionmfab_parallel_copy_cgv

  subroutine amrex_sepunionmfab_fill_boundary (this, geom, cross)
    class(amrex_sepunionmfab) :: this
    type(amrex_geometry), intent(in) :: geom
    logical, intent(in), optional :: cross
    call this%amrex_sepunionmfab_fill_boundary_c(geom, 1, this%nc, cross)
  end subroutine amrex_sepunionmfab_fill_boundary

  subroutine amrex_sepunionmfab_fill_boundary_c (this, geom, c, nc, cross)
    class(amrex_sepunionmfab) :: this
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: c, nc
    logical, intent(in), optional :: cross
    integer :: lcross
    lcross = 0
    if (present(cross)) then
       if (cross) then
          lcross = 1
       else
          lcross = 0
       end if
    end if
    call amrex_fi_sepunionmfab_fill_boundary(this%p, geom%p, c-1, nc, lcross)
  end subroutine amrex_sepunionmfab_fill_boundary_c

  subroutine amrex_sepunionmfab_write (mf, name)
    type(amrex_sepunionmfab), intent(in) :: mf
    character(*), intent(in) :: name
    call amrex_fi_write_sepunionmfab(mf%p, amrex_string_f_to_c(name))
  end subroutine amrex_sepunionmfab_write

  subroutine amrex_sepunionmfab_read (mf, name)
    type(amrex_sepunionmfab), intent(inout) :: mf
    character(*), intent(in) :: name
    call amrex_fi_read_sepunionmfab(mf%p, amrex_string_f_to_c(name))
    mf%owner = .true.
    mf%nc    = amrex_fi_sepunionmfab_ncomp(mf%p)
    call amrex_fi_sepunionmfab_ngrow(mf%p, mf%ng)
    mf%ba    = amrex_fi_sepunionmfab_boxarray(mf%p)
    mf%dm    = amrex_fi_sepunionmfab_distromap(mf%p)
  end subroutine amrex_sepunionmfab_read

end module amrex_sepunionmfab_module

