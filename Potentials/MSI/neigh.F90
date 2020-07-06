
! choose which data points to include in the neighbour list


subroutine neighbour(sys,interp,pot,RawWeight,r,neigh,RawWeightTemp,inner_update,outer_update)
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

  integer :: i,j
  real(kind=8) totsum,tol, tmpWeight
  logical, intent(in) :: inner_update, outer_update

  if (outer_update) then 
     call update_outer_neighlist(interp,neigh,totsum,RawWeight)
  endif

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

  
  if (inner_update) then
     !----------------------------------------------------------
     !  build the inner neighbour list
     !----------------------------------------------------------
     neigh%numInner=0    ! number of neighbours
     neigh%inner = 0      ! list of (inner) neighbours

     tol = interp%wtol*totsum
     do i=1,interp%ndata
        if (RawWeight(i) > tol) then
          neigh%numInner = neigh%numInner + 1
          neigh%inner(neigh%numInner) = i
        endif
     enddo
  endif

  !----------------------------------------------------------

  return
end subroutine


subroutine update_outer_neighlist(interp,neigh,totsum,RawWeight)
     use interpolation
     
     type(interp_params), intent(in) :: interp
     type(neighbour_list), intent(out) :: neigh 
     real(kind=8), dimension(interp%ndata), intent(in) :: RawWeight
     real(kind=8), intent(in) :: totsum
     real(kind=8) :: outtol

     integer :: i

     !---------------------------------------------------------
     ! build the outer neighbour list
     !---------------------------------------------------------
     neigh%numOuter=0
     neigh%outer=0
     outtol=interp%outer*totsum
     do i=1,interp%ndata
        if (RawWeight(i) > outtol) then
            neigh%numOuter = neigh%numOuter + 1
            neigh%outer(neigh%numOuter) = i
        endif
     enddo

     return
end subroutine
