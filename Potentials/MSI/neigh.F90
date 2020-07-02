
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
  
    !----------------------------------------------------------
    !  build the inner neighbour list
    !----------------------------------------------------------

    !totsum = sum(RawWeight)
    !tol = interp%wtol*totsum
    !do i=1,interp%ndata
    !   if (RawWeight(i) > tol) then
    !     neigh%numInner = neigh%numInner + 1
    !     neigh%inner(neigh%numInner) = i
    !   endif
    !enddo
  if (interp%inneigh_update == 1) then
     call inner_neigh_update(interp,neigh,RawWeight) 
  else 
     if (update) then

         call inner_neigh_update(interp,neigh,RawWeight) 

         do i=1,neigh%numInner
            neigh_copy(iatom,ind,i) = neigh%inner(i)
         enddo

     elseif (update==.FALSE. .AND. current_MCstep.LE.1) then

         !print *, 'this should only appear before and during the first update'
         call inner_neigh_update(interp,neigh,RawWeight)

     elseif (update==.FALSE. .AND. current_MCstep.GT.1) then 

         !=====================================
         ! Restore the neigh%inner array
         !===================================== 
    
         !do i=1, interp%ndata
         !   if (neigh_copy(iatom,ind,i).GT.0) then
         !      neigh%numInner = neigh%numInner+1 
         !      neigh%inner(neigh%numInner) = neigh_copy(iatom,ind,neigh%numInner)
         !   endif     
         !enddo
         !print *, 'restoring neighbour list'

         neigh%numInner = count(neigh_copy(iatom,ind,:).ne.0)
         neigh%inner = neigh_copy(iatom,ind,:)

     else
         stop 'Unexpected exception in neighbour list update: neigh.F90'
         call exit(0)
     endif
     
     update =.FALSE.
  endif
  !print *, (neigh%inner(i),i=1,neigh%numInner)
  
  return
end subroutine

subroutine inner_neigh_update(interp,neigh,RawWeight)

    use interpolation  

    implicit none

    type (interp_params), intent(in) :: interp
    real(kind=8), dimension(interp%ndata), intent(in) :: RawWeight
    type (neighbour_list), intent(out) :: neigh
    
    real(kind=8) :: tol,totsum 
    integer :: i
    tol = 0.d0
    totsum = 0.d0

    !----------------------------------------------------------
    !  build the inner neighbour list
    !----------------------------------------------------------
    totsum = sum(RawWeight)
    tol = interp%wtol*totsum
    do i=1,interp%ndata
       if (RawWeight(i) > tol) then
         neigh%numInner = neigh%numInner + 1
         neigh%inner(neigh%numInner) = i
       endif
    enddo

    return
end subroutine
