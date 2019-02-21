!  parameters read in from IN_SYSTEM
module molecule_specs
   implicit none
   real(kind=8), parameter :: amu2au = 1822.888481
   real(kind=8), parameter :: bohr2m = 5.29177e-11
   real(kind=8), parameter :: twopi = 6.28318530717959
   real(kind=8), parameter :: kJ = 2625.4973
   real(kind=8), parameter :: kCal = kJ/4.12
   ! define the type
   type molsysdat
     integer natom,nbond,nint,iseed,ngroup
     character(len=70) :: title
     character(len=3), dimension(:), pointer ::  atom_label
     real(kind=8), dimension(:), pointer :: mass
     integer, dimension(:), pointer :: nb,mb
     integer :: dimen   !Stores the dimensionality of the system (3 for a molecule embedded in real space)
     real(kind=8), dimension(:,:), pointer :: EquilibriumGeom
   end type molsysdat
   ! declare the constructor
   interface new
     module procedure MolSysDat_init
   end interface

 contains
   ! the allocation routine
   subroutine MolSysDat_init(this,n,dimen,title)
     type (molsysdat), intent(out) :: this
     integer, intent(in) :: n
     integer, intent(in) :: dimen
     character(len=70) :: title
     integer :: ierr
     this%natom = n
     this%nbond = n*(n-1)/2
     if (n > 2) then
       this%nint  = 3*n-6
     else
       this%nint  = 3*n-5
     endif
     this%title = title
     this%dimen = dimen
     allocate(this%EquilibriumGeom(dimen,n),stat=ierr)
     if (ierr.ne.0) stop 'Error allocating pimc_par%EquilibriumGeom'
     allocate(this%atom_label(n),stat=ierr)
     if (ierr.ne.0) stop ' molsysdat allocation error: atom_label '
     allocate(this%mass(n),stat=ierr)
     if (ierr.ne.0) stop ' molsysdat allocation error: mass '
     allocate(this%nb(this%nbond),stat=ierr)
     if (ierr.ne.0) stop ' molsysdat allocation error: nb '
     allocate(this%mb(this%nbond),stat=ierr)
     if (ierr.ne.0) stop ' molsysdat allocation error: mb '
     return
   end subroutine MolSysDat_init

end module molecule_specs
