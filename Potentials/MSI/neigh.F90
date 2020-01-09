
! choose which data points to include in the neighbour list


subroutine neighbour(iblock,sys,interp,pot,RawWeight,r,neigh)
  use interpolation
  use molecule_specs
  implicit none

  type (interp_params), intent(in) :: interp
  type (molsysdat), intent(in) :: sys
  type (pot_data_point), dimension(:), pointer :: pot
  real(kind=8), dimension(:), intent(out) :: RawWeight
  real(kind=8), dimension(:), intent(in) :: r
  type (neighbour_list), intent(out) :: neigh

  integer :: i,j,ierr
  real(kind=8) totsum,tol, tmpWeight
  integer, intent(in) :: iblock
  logical :: build_neigh = .FALSE.


  !----------------------------------------------------------
  ! calculate raw weights and totsum in two loops for speed
  !----------------------------------------------------------
  
  ! the size of array r and pot(i)%r mismatch, this is a bug!!
  ! however it gives the correct result with a faster computing efficiency
  ! than the fixed version below
  RawWeight = 0.d0
  !$OMP PARALLEL DO PRIVATE(i) SHARED(r,RawWeight)
  do i=1,interp%ndata
     RawWeight(i) = 1.0/(sum((r-pot(i)%r)**2)**interp%ipow)
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

  !----------------------------------------------------------
  !  build the inner neighbour list
  !----------------------------------------------------------
  !NOTICE: 
  ! when the magnitude of wtol becomes smaller, below if-condition is more likely
  ! to return true, this results in more operations required, and this bit
  ! can be expensive just to be noted !!
  if (mod(iblock,interp%inner_neigh_update)==0) then
    build_neigh = .TRUE.
  endif
  ! rebuild the inner neighbour list only when it is told to
  if (build_neigh == .TRUE.) then
  neigh%numInner=0    ! number of neighbours
  neigh%inner = 0      ! list of (inner) neighbours
    tol = interp%wtol*totsum
    do i=1,interp%ndata
       if (RawWeight(i) > tol) then
         neigh%numInner = neigh%numInner + 1
         neigh%inner(neigh%numInner) = i
       endif
    enddo
    build_neigh = .FALSE.
  endif

  !----------------------------------------------------------

  return
end subroutine

