

program symplectic_int_test
  use molecule_specs
  use pimc_structures
  use potential_msi

  implicit none
  
  type(molsysdat) :: sys
  type(pimc_par) :: pimc
  !type(pimc_particle), dimension(:), pointer :: OldDiscret
  !type(pimc_particle), dimension(:), pointer :: OlderDiscret
  
  character(len=80) :: IN_SYSTEM, IN_FILE

  ! for MSI
  type(msi_params) :: pot
  character(len=80) :: POT_FILE, IN_INTERP, IN_ATOMPERM, OUT_COORD
  
  integer :: i 
  ! two parts, first I need Beads(i)%dVdx(k,j) to be parsed from pimc90
  ! second, I need to read a initial geometry from somewhere

  !two thoughts for generating the initial geometry
  ! 1: averaged from MC converged geometries
  ! 2: one of the converged geometries that gives the minimium energy 
  call getarg(1,IN_SYSTEM)
  call getarg(2,IN_FILE)
  call getarg(3,IN_INTERP)
  call getarg(4,POT_FILE)
  call getarg(5,IN_ATOMPERM)
  call getarg(6,OUT_COORD)


  !call system('mkdir -p '//trim(OUT_COORD)) 

  !read input file
  call read_system_data(sys,IN_SYSTEM) !reading initial geo is also implemented here
  !print *, sys%Beads1Geom
  !print *, ''
  !print *, sys%Bead2Geom
  call read_input(pimc,IN_FILE) ! this input should store the Trotter numnber size
  !print *, sys%dimen, sys%natom
  !print *, pimc%NumBeadsEff
   
  ! Intialise the MSI PES
  call MSI_INIT(pot,sys,IN_INTERP,POT_FILE,IN_ATOMPERM,pimc%NumDiscretisation)
  
  ! computing the potential derivatives w.r.t. beads need to be done in the integrator subroutine
  !do i=1,pimc%NumDiscretisation
      !MSI potential energy surfaces
  !     call potential(i,pot,Beads(i)%x,Beads(i)%r,Beads(i)%VCurr,Beads(i)%dVdx)
  !   print *, 'for bead ', i
  !   print *, 'cartesian coordinates of 5 atoms are '
  !   print *, Beads(i)%x
  !   print *, 'VCurr is '
  !   print *, Beads(i)%VCurr
  !   print *, 'derivative w.r.t. bead ', i, 'cartesian coordinate'
  !   print *, Beads(i)%dVdx
  !enddo
  !print *, 'exiting call potential module'
  
  !print *,  Beads(i)%VCurr
  ! this line enforces the periodicity
  !call copy(Beads(1),Beads(pimc%NumBeads+1)) 
  call numerical_integrator(sys,pimc,pot,OUT_COORD) 

  ! main routine - iterating the symplectic integrator
  !call symplectic_integrator()

end
