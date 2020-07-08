
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
  real(kind=8) totsum,tol, outertotsum
  logical, intent(in) :: inner_update, outer_update
 
  RawWeightTemp=0.d0
  RawWeight=0.d0
  !print *, interp%outer_disable
  !call exit(0)

  if (interp%outer_disable) then

     !----------------------------------------------------------
     ! calculate raw weights and totsum in two loops for speed
     !----------------------------------------------------------
     !$OMP PARALLEL DO PRIVATE(i) SHARED(r,RawWeight,RawWeightTemp)
     do i=1,interp%ndata
        RawWeightTemp(i) = 1.0/(sum((r-pot(i)%r)**2))
        RawWeight(i) = RawWeightTemp(i)**interp%ipow
     enddo
     !$OMP END PARALLEL DO

     totsum = sum(RawWeight)
     
     if (outer_update) then
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

  else
     ! outer neighbour list defines the partial PES data points
     if (outer_update) then 
        !$OMP PARALLEL DO PRIVATE(i) SHARED(r,RawWeight,RawWeightTemp)
        do i=1,interp%ndata
           RawWeightTemp(i) = 1.0/(sum((r-pot(i)%r)**2))
           RawWeight(i) = RawWeightTemp(i)**interp%ipow
        enddo
        !$OMP END PARALLEL DO
        totsum = sum(RawWeight)
        call update_outer_neighlist(interp,neigh,totsum,RawWeight)
     endif
        !print *, 'outer:', neigh%numOuter

     !----------------------------------------------------------
     ! calculate raw weights and totsum in two loops for speed
     !----------------------------------------------------------
     !$OMP PARALLEL DO PRIVATE(i) SHARED(r,RawWeight,RawWeightTemp)
     do i=1,neigh%numOuter
        RawWeightTemp(neigh%outer(i)) = 1.0/(sum((r-pot(neigh%outer(i))%r)**2))
        RawWeight(neigh%outer(i)) = RawWeightTemp(neigh%outer(i))**interp%ipow
     enddo
     !$OMP END PARALLEL DO

     outertotsum = sum(RawWeight)
     
     if (inner_update) then
        !----------------------------------------------------------
        !  build the inner neighbour list
        !----------------------------------------------------------
        neigh%numInner=0    ! number of neighbours
        neigh%inner = 0      ! list of (inner) neighbours

        tol = interp%wtol*outertotsum
        do i=1,interp%ndata
           if (RawWeight(i) > tol) then
             neigh%numInner = neigh%numInner + 1
             neigh%inner(neigh%numInner) = i
           endif
        enddo
     endif

  endif

  return
end subroutine

subroutine append_array(neighlist,size_neigh,reset,vec)  
   use interpolation
   
   integer, intent(in) :: size_neigh
   integer, dimension(size_neigh), intent(in) ::neighlist
   integer, dimension(size_neigh+1), intent(inout) :: vec
   integer, dimension(:), allocatable :: temp_neigh
   integer :: i,j, first_ind, last_ind
   logical, intent(in) :: reset
   j = 0
   
   !----------------------------------------------------------- 
   !> append all non-duplicated MSI data points that contribute
   !> to the potential into an array 
   !----------------------------------------------------------- 
   if (reset) then
      vec=0
   endif

   allocate(temp_neigh(count(neighlist.ne.0)))
   temp_neigh=0

   do i=1, size_neigh
     if (ANY(vec .EQ. neighlist(i))==.FALSE.) then
        j = j + 1 
        temp_neigh(j) = neighlist(i) 
     endif
   enddo
   
   first_ind = count(vec.NE.0)+1
   last_ind = count(vec.NE.0)+count(temp_neigh.ne.0)
   do i=first_ind, last_ind 
      j=i-first_ind+1
      vec(i) = temp_neigh(j)
   enddo

   deallocate(temp_neigh)
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


!function(array,array_size) result(mean,std)
subroutine num_neigh(array,array_size)
    implicit none
    integer, intent(in) :: array(1) 
    integer, intent(in) :: array_size
    integer :: i
    real(kind=8) :: mean, std
    
    !do i=1, array_size
    !   print *, array(i)
    !enddo
    ! calculate the mean
    mean = 0.d0
    do i=1,array_size
       mean=mean+dble(array(i))
    enddo
    mean=mean/dble(array_size)
    
    std=0.d0
    ! calculate 2sigma
    do i=1,array_size
       std=std+dble(array(i)-mean)**2
    enddo
    std=2.0*Sqrt(std/(array_size-1))   
    
    print *, 'mean =', mean, '2sigma =', std 
    return
!end function
end subroutine
