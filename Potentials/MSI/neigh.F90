
! choose which data points to include in the neighbour list


subroutine neighbour(sys,interp,pot,RawWeight,r,neigh,RawWeightTemp,current_MCstep,update,pimc,neigh_copy,iatom,ind)
  use interpolation
  use molecule_specs
  use pimc_structures
  implicit none

  type (interp_params), intent(in) :: interp
  type (molsysdat), intent(in) :: sys
  type (pot_data_point), dimension(:), pointer :: pot
  real(kind=8), dimension(:), intent(out) :: RawWeight
  real(kind=8), dimension(:), intent(out) :: RawWeightTemp
  real(kind=8), dimension(:), intent(in) :: r
  type (neighbour_list), intent(out) :: neigh
  type(neighbour_list), dimension(:), pointer :: old_neigh
  integer, intent(in) :: current_MCstep, iatom, ind
  logical, intent(inout) :: update
  type(pimc_par), intent(in) :: pimc

  integer :: i,j
  real(kind=8) totsum,tol, tmpWeight

  ! neighlist copy 
  integer, dimension(pimc%atom_pass,pimc%NumBeadsEff,interp%ndata), intent(inout) :: neigh_copy
  integer :: counter

  neigh%numInner=0    ! number of neighbours
  neigh%inner = 0      ! list of (inner) neighbours
  counter = 0

  !----------------------------------------------------------
  ! calculate raw weights and totsum in two loops for speed
  !----------------------------------------------------------

  RawWeight = 0.d0
  !$OMP PARALLEL DO PRIVATE(i) SHARED(r,RawWeight,RawWeightTemp)
  do i=1,interp%ndata
     RawWeightTemp(i) = 1.0/(sum((r-pot(i)%r)**2))
     RawWeight(i) = RawWeightTemp(i)**interp%ipow
  enddo
  !$OMP END PARALLEL DO
  totsum = sum(RawWeight)

  !if (update) then

    !----------------------------------------------------------
    !  build the inner neighbour list
    !----------------------------------------------------------

    tol = interp%wtol*totsum
    do i=1,interp%ndata
       if (RawWeight(i) > tol) then
         neigh%numInner = neigh%numInner + 1
         neigh%inner(neigh%numInner) = i
       endif
    enddo
    print *, (neigh%inner(i), i=1,neigh%numInner)
    print *, '************'

    !do i=1,neigh%numInner
    !   neigh_copy(iatom,ind,i) = neigh%inner(i)
    !enddo

     !print *, (neigh_copy(iatom,ind,i),i=1,neigh%numInner)
     !print *, '------------------------------'
    
     update = .FALSE.

  !endif

  !if (update==.FALSE. .AND. current_MCstep .GT. 1) then 

  !   !=====================================
  !   ! Restore the neigh%inner array
  !   !===================================== 
  !   do i=1, interp%ndata
  !      if (neigh_copy(iatom,ind,i).GT.0) then
  !         counter = counter+1 
  !         neigh%inner(counter) = neigh_copy(iatom,ind,counter)
  !      endif     
  !   enddo
  !  do i=1, interp%ndata
  !     print *, neigh_copy(iatom,ind,i)
  !  enddo

  !endif
  
  return
end subroutine
