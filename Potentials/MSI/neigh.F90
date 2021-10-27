
! choose which data points to include in the neighbour list


subroutine neighbour(sys,interp,pot,RawWeight,r,neigh,RawWeightTemp)
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
  real(kind=8) :: tot

  integer :: i,j
  real(kind=8) totsum,tol, tmpWeight

  neigh%numInner=0    ! number of neighbours
  neigh%inner = 0      ! list of (inner) neighbours

  !----------------------------------------------------------
  ! calculate raw weights and totsum in two loops for speed
  !----------------------------------------------------------
  
  ! the size of array r and pot(i)%r mismatch, this is a bug!!
  ! however it gives the correct result with a faster computing efficiency
  ! than the fixed version below


!    RawWeightTemp = 0.d0
!    RawWeight = 0.d0
!    do i=1,interp%ndata
!       RawWeightTemp(i) = 1.0/(sum((r-pot(i)%r)**2))
!       RawWeight(i) = RawWeightTemp(i)**interp%ipow
!    enddo
!    totsum = sum(RawWeight)


  totsum = 0.d0
!  !$acc parallel loop reduction(+:totsum) 
  do j=1,interp%ndata
     RawWeightTemp(j) = 0.d0
     RawWeight(j) = 0.d0
     tot = 0.d0
!     !$acc loop 
     do i=1,size(r)
       tot = tot + (r(i) - pot(j)%r(i))**2 
     enddo
     RawWeightTemp(j) = 1.d0/tot
     RawWeight(j) = RawWeightTemp(j)**interp%ipow
     totsum = totsum + RawWeight(j)
  enddo
!  !$acc end parallel loop


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

  return
end subroutine

