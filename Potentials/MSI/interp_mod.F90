
module interpolation
  implicit none
  !---------------------------------------------
  ! parameters which are read in from IN_INTERP
  !---------------------------------------------
  type interp_params
    integer :: ipart,neigh_update,nneigh,ngroup
    integer :: ipow,lowp, ipow2
    integer :: ndata,nforc
    real(kind=8) :: wtol,outer,vmax, vmin, sigma
    real(kind=8) :: avrads
    integer, dimension(:,:), pointer :: aperm, bperm
  end type interp_params
  !---------------------------------------
  ! an ab initio data point 
  !---------------------------------------
  type pot_data_point
    real(kind=8), pointer :: r(:)
    real(kind=8), pointer :: z(:)
    !Used when some of the system is held rigid. Precomputes the 
    !rigid components
    real(kind=8) :: UtR_rigid
    real(kind=8) :: Weight_rigid
    real(kind=8), pointer :: ut(:,:)
    real(kind=8) ::  v0
    real(kind=8), pointer :: v1(:),v2(:)
    real(kind=8) :: rad
    real(kind=8), pointer :: sumb(:)
    integer :: id,perm
  end type pot_data_point
  
  !---------------------------------------
  !  neighbour lists
  !---------------------------------------
  type neighbour_list
    real(kind=8), dimension(:), pointer :: r
    integer, dimension(:), pointer :: inner
    integer :: numInner
  end type
  !-------------------------------------
  ! constructor for PotDataPoint
  !-------------------------------------
  interface new
    module procedure PotDataPoint_init
  end interface
  interface destroy
    module procedure PotDataPoint_clean
  end interface

  !-------------------------------------
  ! constructor for interp_params
  !-------------------------------------
  interface perms
    module procedure InterpParams_perms
  end interface
contains
  
  subroutine InterpParams_perms(this,n,sys)
    use molecule_specs
    type(interp_params), intent(out) :: this
    integer, intent(in) :: n
    type(molsysdat), intent(in) :: sys
    integer :: ierr
    this%ngroup = n
    allocate(this%aperm(sys%natom,this%ngroup),stat=ierr)
    if (ierr.ne.0) stop ' interp_params allocation error: aperm '
    allocate(this%bperm(sys%nbond,this%ngroup),stat=ierr)
    if (ierr.ne.0) stop ' Interp_params allocation error: bperm '
  end subroutine InterpParams_perms

  !--------------------------------------
  ! create a pot_data_point
  !--------------------------------------
  subroutine PotDataPoint_init(this,n) 
    type (pot_data_point) :: this
    integer, intent(in) :: n
    integer :: nbond,nint,ierr
    nbond = n*(n-1)/2  ! n is number of atoms
    if (n > 2) then
       nint = 3*n-6
    else
       nint = 3*n-5
    endif
    allocate(this%r(nbond),stat=ierr)
    if (ierr.ne.0) stop ' PotDataPoint allocation error '
    allocate(this%z(nint),stat=ierr)
    if (ierr.ne.0) stop ' PotDataPoint allocation error '
    allocate(this%ut(nint,nbond),stat=ierr)
    if (ierr.ne.0) stop ' PotDataPoint allocation error '
    allocate(this%v1(nint),stat=ierr)
    if (ierr.ne.0) stop ' PotDataPoint allocation error '
    allocate(this%v2(nint),stat=ierr)
    if (ierr.ne.0) stop ' PotDataPoint allocation error '
    allocate(this%sumb(nbond),stat=ierr)
    if (ierr.ne.0) stop ' PotDataPoint allocation error '
    this%rad = 0.0
    return
  end subroutine PotDataPoint_init

  !--------------------------------------
  ! destroy a pot_data_point
  !--------------------------------------
  subroutine PotDataPoint_clean(this)
    type (pot_data_point):: this
    integer :: ierr
    deallocate(this%r,stat=ierr)
    if (ierr.ne.0) stop ' Error deallocating PotDataPoint '
    deallocate(this%z,stat=ierr)
    if (ierr.ne.0) stop ' Error deallocating PotDataPoint '
    deallocate(this%ut,stat=ierr)
    if (ierr.ne.0) stop ' Error deallocating PotDataPoint '
    deallocate(this%v1,stat=ierr)
    if (ierr.ne.0) stop ' Error deallocating PotDataPoint '
    deallocate(this%v2,stat=ierr)
    if (ierr.ne.0) stop ' Error deallocating PotDataPoint '
    deallocate(this%sumb,stat=ierr)
    if (ierr.ne.0) stop ' Error deallocating PotDataPoint '
    return
  end subroutine PotDataPoint_clean

end module interpolation


