
! choose which data points to include in the neighbour list


subroutine neighbour(sys,interp,pot,RawWeight,r,neigh,RawWeightTemp,current_MCstep)
  use interpolation
  use molecule_specs
  implicit none

  type (interp_params), intent(in) :: interp
  type (molsysdat), intent(in) :: sys
  type (pot_data_point), dimension(:), pointer :: pot
  real(kind=8), dimension(:), intent(out) :: RawWeight
  real(kind=8), dimension(:), intent(out) :: RawWeightTemp
  real(kind=8), dimension(:), intent(in) :: r
  type (neighbour_list), intent(out) :: neigh
  integer, intent(in) :: current_MCstep

  integer :: i,j
  real(kind=8) totsum,tol, tmpWeight

  neigh%numInner=0    ! number of neighbours
  neigh%inner = 0      ! list of (inner) neighbours

  !----------------------------------------------------------
  ! calculate raw weights and totsum in two loops for speed
  !----------------------------------------------------------
  
  ! the size of array r and pot(i)%r mismatch
  ! however it does not give the wrong answer but rather it's faster
  ! than the fixed version below
  RawWeight = 0.d0
  !$OMP PARALLEL DO PRIVATE(i) SHARED(r,RawWeight,RawWeightTemp)
  do i=1,interp%ndata
     RawWeightTemp(i) = 1.0/(sum((r-pot(i)%r)**2))
     RawWeight(i) = RawWeightTemp(i)**interp%ipow
  enddo
  !$OMP END PARALLEL DO
  totsum = sum(RawWeight)
  
  !RawWeight = 0.d0
  !do j = 1, interp%ndata
  !   tmpWeight = 0.d0
  !   do i = 1,sys%nbond
  !      tmpWeight = tmpWeight + (r(i) - pot(j)%r(i))**2
  !   enddo
  !      RawWeight(j) = tmpWeight
        
  !enddo
  !RawWeight = 1.0/ (RawWeight**interp%ipow)
        
  !totsum = sum(RawWeight)
  

  ! periodically update inner neighbour list
  if (mod(current_MCstep,interp%inneigh_update)==0) then
     
  endif

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

  !----------------------------------------------------------


  ! build the outer neighour list

  
  return
end subroutine

