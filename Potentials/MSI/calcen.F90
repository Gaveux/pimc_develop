

subroutine calcen(sys,interp,pot,neigh,Weight,r,V,dVdR,RawWeightTemp,ddr_Fsqr) 
    use molecule_specs
    use interpolation

    implicit none

    type (molsysdat), intent(in) :: sys
    type (interp_params), intent(in) :: interp
    type (pot_data_point), dimension(:), pointer :: pot
    type (neighbour_list), intent(in) :: neigh

    real(kind=8), dimension(:), intent(inout) :: Weight
    real(kind=8), dimension(:), intent(in) :: RawWeightTemp
    real(kind=8), dimension(:), intent(in) :: r
    real(kind=8), intent(out) :: V
    real(kind=8), dimension(sys%nbond), intent(out) :: dVdR
    real(kind=8), dimension(sys%nbond), intent(out) :: ddr_Fsqr

    real(kind=8), dimension(size(Weight)) :: Raw
    ! by defining the dimension of Raw to neigh%numInner, this restricts the
    ! size of array to the correct size w.r.t. each iteration within the dynamical array
    !real(kind=8), dimension(neigh%numInner) :: Raw
    real(kind=8), dimension(sys%nbond,size(Weight)) :: DWeight
    real(kind=8), dimension(sys%nbond,size(Weight)) :: DRawWeight
    real(kind=8), dimension(size(Weight),sys%nbond) :: dTaydR
    real(kind=8), dimension(size(Weight),sys%nbond) :: d2TaydR2
    real(kind=8), dimension(sys%nbond,size(Weight)) :: D2RawWeight
    real(kind=8), dimension(sys%nbond,size(Weight)) :: D2Weight
    real(kind=8), dimension(sys%nint,size(Weight)) :: DTay
    real(kind=8), dimension(sys%nbond) :: SumDRawWeight
    real(kind=8), dimension(sys%nbond) :: SumD2RawWeight
    real(kind=8) :: totsum, energy, temp, invTotsum, temp2

    !Stores the value of the taylor series expansions
    real(kind=8), dimension(interp%ndata) :: Tay
    !Stores the zeta coordinates of the system
    real(kind=8), dimension(sys%nint,interp%ndata) :: z
    real(kind=8), dimension(sys%nint,interp%ndata) :: rdzdr
    integer :: i,j,k
    
    !---------------------------------------------------
    !  Calculate the Weights 
    !---------------------------------------------------
    !!$OMP PARALLEL DO PRIVATE(i) SHARED(Raw, Weight)
    do i=1,neigh%numInner
        Raw(i) = Weight(neigh%inner(i))
    enddo
    !!$OMP END PARALLEL DO 
      
    totsum = sum(Raw(1:neigh%numInner))
    invTotsum = 1.0/totsum
    Weight = Raw/totsum
   
    !---------------------------------------------------
    !  Calculate the derivatives of the weights
    !---------------------------------------------------
   
    ! derivative of raw weights
    !!$OMP PARALLEL
    !!$OMP DO
    !!$OMP PARALLEL DO PRIVATE(i,k) SHARED(DWeight,temp,r)
    do k=1,neigh%numInner
        temp = interp%ipow2*Weight(k)*RawWeightTemp(k)
        do i=1,sys%nbond
            DRawWeight(i,k) = temp*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**2
            D2RawWeight(i,k) = temp*(2.0*pot(neigh%inner(k))%r(i)*r(i)**3 &
            - 3.0*r(i)**4) 
            D2RawWeight(i,k) = D2RawWeight(i,k) + temp*RawWeightTemp(k) &
            *2.0*(interp%ipow+1.0)*(r(i)**3 - pot(neigh%inner(k))%r(i)*r(i)**2)**2
        enddo
    enddo
    !!$OMP END PARALLEL DO
    
    SumDRawWeight = 0.0
    SumD2RawWeight = 0.0
    !!$OMP PARALLEL DO PRIVATE(i,k) SHARED(SumDWeight,DWeight)
    do k=1,neigh%numInner
        do i=1,sys%nbond
            SumDRawWeight(i) = SumDRawWeight(i) + DRawWeight(i,k)
            SumD2RawWeight(i) = SumD2RawWeight(i) + D2RawWeight(i,k)
        enddo
    enddo
    !!$OMP END PARALLE DO   

    ! derivative of relative weights
    !!$OMP PARALLEL DO PRIVATE(i,k) SHARED(DWeight,Weight,SumDWeight)
    do k=1,neigh%numInner
        do i=1,sys%nbond
            DWeight(i,k) = DRawWeight(i,k) - Weight(k)*SumDRawWeight(i)
            D2Weight(i,k) = 2.0*(Weight(k)*SumDRawWeight(i)**2 - DRawWeight(i,k)*SumDRawWeight(i)) &
                    + invTotsum*(D2RawWeight(i,k)-Weight(k)*SumD2RawWeight(i))
        enddo
    enddo

    !!$OMP END PARALLEL DO
    !!$OMP END PARALLEL
    !---------------------------------------------------
    !  Evaluate (z-z0) the local internal coordinates
    !---------------------------------------------------
    z = 0.0
    ! nbond =N(N-1)/2, nint=3N-6
    do k=1,neigh%numInner
        do j = 1,sys%nbond
            do i = 1,sys%nint
                z(i,k) = z(i,k) + pot(neigh%inner(k))%ut(i,j)*r(j)
                !rdzdr(i,k) = rdzdr(i,k) + pot(neigh%inner(k))%ut(i,j)*r(j)**2
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
    !  Calculate the gradient and normal of the energy (w.r.t internals)
    !----------------------------------------------------------
     
    ! derivative of Taylor polynomial in local internals
    do k=1,neigh%numInner
        do i=1,sys%nint
            DTay(i,k) = z(i,k)*pot(neigh%inner(k))%v2(i) + &
            pot(neigh%inner(k))%v1(i)
        enddo
    enddo

    ! ut(i,j) sum in nint
    do j=1,sys%nbond
        do k=1,neigh%numInner
        d2TaydR2(k,j) = 0.0
          do i=1,sys%nint
             d2TaydR2(k,j) = d2TaydR2(k,j) + pot(neigh%inner(k))%ut(i,j)
          enddo
        enddo
    enddo
   
    ! derivative of the Taylor polynomial w.r.t r^-1
    do j=1,sys%nbond
        temp = -r(j)**2
        temp2 = r(j)**3
        do k=1,neigh%numInner
            dTaydR(k,j) = 0.0
            do i=1,sys%nint
                dTaydR(k,j) = dTaydR(k,j) + DTay(i,k)*pot(neigh%inner(k))%ut(i,j)
                d2TaydR2(k,j) = d2TaydR2(k,j)*(pot(neigh%inner(k))%ut(i,j)*pot(neigh%inner(k))%v2(i)*r(j) + 2.0*DTay(i,k))
            enddo
            dTaydR(k,j) = dTaydR(k,j)*temp
            d2TaydR2(k,j) = d2TaydR2(k,j)*temp2
        enddo
    enddo
   
    ! gradient of the energy
    dVdR = 0.0
    do k=1,neigh%numInner
        do i=1,sys%nbond
            dVdR(i) = dVdR(i) + Tay(k)*DWeight(i,k) + Weight(k)*dTaydR(k,i)
        enddo
    enddo
    ddr_Fsqr = 0.0
    do k=1,neigh%numInner
        do i=1,sys%nbond
            ddr_Fsqr(i) = ddr_Fsqr(i) + Tay(k)*D2Weight(i,k)+Weight(k)*d2TaydR2(k,i)+2.0*dTaydR(k,i)*DWeight(i,k)
        enddo
    enddo

    do i=1,sys%nbond
           ddr_Fsqr(i) = 2.0*ddr_Fsqr(i)*dVdR(i)
    enddo
    
    V = energy 
  return
end subroutine





