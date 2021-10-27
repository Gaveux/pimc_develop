

subroutine calcen(sys,interp,pot,neigh,Weightmp,r,V,dVdR,RawWeightTemp) 
    use molecule_specs
    use interpolation

    implicit none

    type (molsysdat), intent(in) :: sys
    type (interp_params), intent(in) :: interp
    type (pot_data_point), dimension(:), pointer :: pot
    type (neighbour_list), intent(in) :: neigh

    real(kind=8), dimension(:), intent(inout) :: Weightmp
    real(kind=8), dimension(:), intent(in) :: RawWeightTemp
    real(kind=8), dimension(:), intent(in) :: r
    real(kind=8), intent(out) :: V
    real(kind=8), dimension(sys%nbond), intent(out) :: dVdR

    real(kind=8), dimension(neigh%numInner) :: Raw
    real(kind=8), dimension(neigh%numInner) :: Weight
    ! by defining the dimension of Raw to neigh%numInner, this restricts the
    ! size of array to the correct size w.r.t. each iteration within the dynamical array
    !real(kind=8), dimension(neigh%numInner) :: Raw
    real(kind=8), dimension(sys%nbond,size(Weight)) :: DWeight
    real(kind=8), dimension(size(Weight),sys%nbond) :: dTaydR
    real(kind=8), dimension(sys%nint,size(Weight)) :: DTay
    real(kind=8), dimension(sys%nbond) :: SumDWeight
    real(kind=8) :: totsum, energy, temp

    !Stores the value of the taylor series expansions
    real(kind=8), dimension(interp%ndata) :: Tay
    !Stores the zeta coordinates of the system
    real(kind=8), dimension(sys%nint,interp%ndata) :: z
    integer :: i,j,k
    
    !---------------------------------------------------
    !  Calculate the Weights 
    !---------------------------------------------------
    !$acc data copyin(Weightmp(:),neigh%inner(:),neigh%numInner) create(Raw(:)) copyout(Weight(:))
    !$acc parallel loop reduction(+:totsum)
    do i=1,neigh%numInner
       Raw(i) = Weightmp(neigh%inner(i))
       totsum = totsum + Raw(i)
    enddo

    !$acc parallel loop
    do i=1,neigh%numInner
       Weight(i) = Raw(i)/totsum
    enddo
    !$acc end data 

!    print *, Weight
!    print *, ''
!
!    Weight = 0.d0
!    do i=1,neigh%numInner
!        Raw(i) = Weightmp(neigh%inner(i))
!    enddo
!      
!    totsum = sum(Raw(1:neigh%numInner))
!    Weight = Raw/totsum
!    print *, Weight
!    call exit(0)
   
    !---------------------------------------------------
    !  Calculate the derivatives of the weights
    !---------------------------------------------------
   
    ! derivative of raw weights
    do k=1,neigh%numInner
        temp = interp%ipow2*Weight(k)*RawWeightTemp(neigh%inner(k))
        do i=1,sys%nbond
            DWeight(i,k) = temp*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**2
        enddo
    enddo
    
    SumDWeight = 0.0
    do k=1,neigh%numInner
        do i=1,sys%nbond
            SumDWeight(i) = SumDWeight(i) + DWeight(i,k)
        enddo
    enddo

    ! derivative of relative weights
    do k=1,neigh%numInner
        do i=1,sys%nbond
            DWeight(i,k) = DWeight(i,k) - Weight(k)*SumDWeight(i)
        enddo
    enddo
    !---------------------------------------------------
    !  Evaluate (z-z0) the local internal coordinates
    !---------------------------------------------------
    z = 0.0
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
   
    energy=0.0
    do k=1,neigh%numInner
        energy = energy + Weight(k)*Tay(k)
    enddo
     
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
   
    ! gradient of the energy
    dVdR = 0.0
    do k=1,neigh%numInner
        do i=1,sys%nbond
            dVdR(i) = dVdR(i) + Tay(k)*DWeight(i,k) + Weight(k)*dTaydR(k,i)
        enddo
    enddo
    
    V = energy 
  return
end subroutine





