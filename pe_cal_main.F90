
program pe_calculator

  use molecule_specs
  use potential_msi

  implicit none

  character(len=80) :: IN_SYSTEM
  character(len=80) :: IN_INTERP
  character(len=80) :: POT_FILE
  character(len=80) :: IN_ATOMPERM

  integer :: i,j,k


  type (molsysdat) :: sys
  type (msi_params) :: pot

  call getarg(1,IN_SYSTEM)
  call getarg(2,IN_INTERP)
  call getarg(3,IN_ATOMPERM)
  call getarg(4,POT_FILE)

  call read_system_data(sys,IN_SYSTEM)
  !do i=1,sys%natom
  !   print *, (sys%EquilibriumGeom(k,i),k=1,sys%dimen)
  !enddo

  call MSI_INIT(pot,sys,IN_INTERP,POT_FILE,IN_ATOMPERM,1)

  call Potential(pot,sys%EquilibriumGeom)

end program
