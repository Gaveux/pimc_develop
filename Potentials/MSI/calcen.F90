

subroutine calcen(sys,interp,pot,neigh,Weight,r,V,dVdR,RawWeightTemp) 
    use molecule_specs
    use interpolation

    implicit none

    type (molsysdat), intent(in) :: sys
    type (interp_params), intent(in) :: interp
    type (pot_data_point), dimension(:), pointer :: pot
    type (neighbour_list), intent(in) :: neigh

    real(kind=8), dimension(:), intent(inout) :: Weight
    real(kind=8), dimension(neigh%numInner) :: testWeight
    real(kind=8), dimension(:), intent(in) :: RawWeightTemp
    real(kind=8), dimension(:), intent(in) :: r
    real(kind=8), intent(out) :: V
    real(kind=8), dimension(sys%nbond), intent(out) :: dVdR

    real(kind=8), dimension(neigh%numInner) :: Raw
    ! by defining the dimension of Raw to neigh%numInner, this restricts the
    ! size of array to the correct size w.r.t. each iteration within the dynamical array
    !real(kind=8), dimension(neigh%numInner) :: Raw
    real(kind=8), dimension(sys%nbond,neigh%numInner) :: DWeight
    real(kind=8), dimension(neigh%numInner,sys%nbond) :: dTaydR
    real(kind=8), dimension(sys%nint,neigh%numInner) :: DTay
    real(kind=8), dimension(sys%nbond) :: SumDWeight
    real(kind=8) :: totsum, temp, tempsum

    !Stores the value of the taylor series expansions
    real(kind=8), dimension(neigh%numInner) :: Tay
    !Stores the zeta coordinates of the system
    real(kind=8), dimension(sys%nint,neigh%numInner) :: z
    integer :: i,j,k
    
    !---------------------------------------------------
    !  Calculate the Weights 
    !---------------------------------------------------
    totsum = 0.d0
    do i=1,neigh%numInner
        Raw(i) = Weight(neigh%inner(i))
    enddo

    totsum = sum(Raw(1:neigh%numInner))

    do i=1,neigh%numInner
       testWeight(i) = Raw(i)/totsum
    enddo

    ! derivative of raw weights
    do k=1,neigh%numInner
    !print *, RawWeightTemp(k)
        temp = interp%ipow2*testWeight(k)*RawWeightTemp(k)
        do i=1,sys%nbond
            DWeight(i,k) = temp*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**2
        enddo
    enddo
    print *, DWeight

    SumDWeight = 0.0
    do k=1,neigh%numInner
        do i=1,sys%nbond
            SumDWeight(i) = SumDWeight(i) + DWeight(i,k)
        enddo
    enddo
    !print *, SumDWeight

    ! derivative of relative weights
    do k=1,neigh%numInner
        do i=1,sys%nbond
            DWeight(i,k) = DWeight(i,k) - testWeight(k)*SumDWeight(i)
        enddo
    enddo
    !print *, DWeight

    z = 0.d0
    do k=1,neigh%numInner
        do j = 1,sys%nbond
            do i = 1,sys%nint
                z(i,k) = z(i,k) + pot(neigh%inner(k))%ut(i,j)*r(j)
            enddo
        enddo
    enddo

    do k=1,neigh%numInner
        do i=1,sys%nint
            z(i,k) = z(i,k) - pot(neigh%inner(k))%z(i)
        enddo
!        print *, z(:,k)
    enddo

    !---------------------------------------------------
    !  Calculate the energy
    !---------------------------------------------------
    do k=1,neigh%numInner
        Tay(k) =  pot(neigh%inner(k))%v0
    enddo

    do k=1,neigh%numInner
        do i=1,sys%nint
            Tay(k) = Tay(k) + z(i,k)*(pot(neigh%inner(k))%v1(i) + 0.5*pot(neigh%inner(k))%v2(i)*z(i,k))
        enddo
    enddo
    !print *, Tay

    V =0.d0
    do k=1,neigh%numInner
        V = V + testWeight(k)*Tay(k)
    enddo
    !print *, V

    !----------------------------------------------------------
    !  Calculate the gradient of the energy (w.r.t internals)
    !----------------------------------------------------------
     
    ! derivative of Taylor polynomial in local internals
    do k=1,neigh%numInner
        do i=1,sys%nint
            DTay(i,k) = z(i,k)*pot(neigh%inner(k))%v2(i) + &
            pot(neigh%inner(k))%v1(i)
        enddo
    enddo
    !print *, DTay

    ! derivative of the Taylor polynomical w.r.t bondlengths (NOT inverses)
    do j=1,sys%nbond
        temp = -r(j)**2
        do k=1,neigh%numInner
            dTaydR(k,j) = 0.0
            do i=1,sys%nint
                dTaydR(k,j) = dTaydR(k,j) + DTay(i,k)*pot(neigh%inner(k))%ut(i,j)
            enddo
            dTaydR(k,j) = dTaydR(k,j)*temp
        enddo
    enddo
    !print *, dTaydR

    ! gradient of the energy
    dVdR = 0.0
    do k=1,neigh%numInner
        do i=1,sys%nbond
            dVdR(i) = dVdR(i) + Tay(k)*DWeight(i,k) + testWeight(k)*dTaydR(k,i)
        enddo
    enddo

    print *, ''

    totsum = 0.d0
    V = 0.d0
    DWeight = 0.d0

!$acc data copyin(neigh,sys,interp%ipow2,totsum,r(:),RawWeightTemp(:),Weight(:),pot(:)) create(SumDWeight(:),DWeight(1:sys%nbond,1:neigh%numInner),z(1:sys%nint,1:neigh%numInner),tempsum,Tay(1:neigh%numInner),DTay(1:sys%nint,1:neigh%numInner),dTaydR(1:neigh%numInner,1:sys%nbond)) copy(V) copyout(dVdR(1:sys%nbond))

    !$acc parallel loop reduction(+:totsum) async(1)
    do i=1,neigh%numInner
        !Raw(i) = Weight(neigh%inner(i))
        totsum = totsum + Weight(neigh%inner(i))
    enddo

    !totsum = sum(Raw(1:neigh%numInner))

    !$acc parallel loop async(1)
    do i=1,neigh%numInner
       Weight(i) = Weight(neigh%inner(i))/totsum
    enddo
   
    !---------------------------------------------------
    !  Calculate the derivatives of the weights
    !---------------------------------------------------
   
    ! derivative of raw weights
    !$acc parallel loop gang async(1)
    do k=1,neigh%numInner
        !$acc loop vector
        do i=1,sys%nbond
            DWeight(i,k) = interp%ipow2*Weight(k)*RawWeightTemp(k)*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**2
        enddo
    enddo

    !$acc loop
    do k=1,neigh%numInner
       do i=1,sys%nbond
          print *, i, k, DWeight(i,k)
       enddo
    enddo

    !$acc parallel loop gang async(1) 
    do k=1,sys%nbond
       !SumDWeight(k) = 0.d0
       temp = 0.d0
       !$acc loop reduction(+:temp) 
       do i=1,neigh%numInner
          !SumDWeight(k) = SumDWeight(k) + DWeight(k,i)
          temp = temp + DWeight(k,i)
       enddo
       SumDWeight(k) = temp
    enddo

    ! derivative of relative weights
    !$acc parallel loop gang async(1)
    do k=1,neigh%numInner
        !$acc loop vector
        do i=1,sys%nbond
            DWeight(i,k) = DWeight(i,k) - Weight(k)*SumDWeight(i)
        enddo
    enddo

!    do k=1,sys%nbond
!       do i=1,neigh%numInner
!          print *, i, k, DWeight(k,i)
!       enddo
!    enddo
    !---------------------------------------------------
    !  Evaluate (z-z0) the local internal coordinates
    !---------------------------------------------------
!    do k=1,neigh%numInner
!        do j = 1,sys%nbond
!            do i = 1,sys%nint
!                z(i,k) = z(i,k) + pot(neigh%inner(k))%ut(i,j)*r(j)
!            enddo
!        enddo
!    enddo

    !$acc parallel loop gang async(2) 
    do k=1,neigh%numInner
       !$acc loop worker
       do j=1,sys%nint
          tempsum = 0.d0
          !$acc loop reduction(+:tempsum)
          do i=1,sys%nbond
             !z(j,k) = z(j,k) + pot(neigh%inner(k))%ut(j,i)*r(i)
             tempsum = tempsum + pot(neigh%inner(k))%ut(j,i)*r(i)
          enddo
          z(j,k) = tempsum - pot(neigh%inner(k))%z(j)
       enddo
    enddo


!!    do k=1,neigh%numInner
!!          print *, z(:,k)
!!    enddo
    
    !---------------------------------------------------
    !  Calculate the energy
    !---------------------------------------------------
!    !$acc parallel loop
!    do k=1,neigh%numInner
!        Tay(k) =  pot(neigh%inner(k))%v0
!    enddo
    
    !$acc parallel loop gang async(2)
    do k=1,neigh%numInner
        Tay(k) =  pot(neigh%inner(k))%v0
        !$acc loop vector
        do i=1,sys%nint
           !$acc atomic
           Tay(k) = Tay(k) + z(i,k)*(pot(neigh%inner(k))%v1(i) + 0.5*pot(neigh%inner(k))%v2(i)*z(i,k))
        enddo
    enddo

   
    !$acc wait(1) async(2)
    !$acc parallel loop reduction(+:V) async(2)
    do k=1,neigh%numInner
        V = V + Weight(k)*Tay(k)
    enddo
     
    !----------------------------------------------------------
    !  Calculate the gradient of the energy (w.r.t internals)
    !----------------------------------------------------------
     
    ! derivative of Taylor polynomial in local internals
    !$acc parallel loop gang async(2) 
    do k=1,neigh%numInner
        !$acc loop vector
        do i=1,sys%nint
            DTay(i,k) = z(i,k)*pot(neigh%inner(k))%v2(i) + &
            pot(neigh%inner(k))%v1(i)
        enddo
    enddo
   
    ! derivative of the Taylor polynomical w.r.t bondlengths (NOT inverses)
    !$acc parallel loop gang worker async(2)
    do j=1,sys%nbond
        do k=1,neigh%numInner
            temp = 0.d0
            !$acc loop reduction(+:temp)
            do i=1,sys%nint
                !dTaydR(k,j) = dTaydR(k,j) + DTay(i,k)*pot(neigh%inner(k))%ut(i,j)
                temp = temp + DTay(i,k)*pot(neigh%inner(k))%ut(i,j)
            enddo
            !dTaydR(k,j) = dTaydR(k,j)*temp
            dTaydR(k,j) = -r(j)**2*temp
        enddo
    enddo

    ! gradient of the energy
    !$acc wait
    !$acc parallel loop gang async(2)
    do k=1,sys%nbond
        temp = 0.d0
        !$acc loop reduction(+:temp)
        do i=1,neigh%numInner
            temp = temp + Tay(i)*DWeight(k,i) + Weight(i)*dTaydR(i,k)
        enddo
        dVdR(k) = temp
    enddo

!$acc end data
    !!print *, dVdR
    !!print *, '============================'
    call exit(0)
    
  return
end subroutine
