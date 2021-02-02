

subroutine calcen(sys,interp,pot,neigh,Weight,r,V,dVdR,RawWeightTemp,dWTdr2,d2veightdr2tmp1,d2veightdr2tmp2,&
                TDWeightSumDWeight,WTSumDWeightDrSqr,WTaySumD2veightDr1,WTaySumD2veightDr2) 
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
    real(kind=8), dimension(sys%nbond), intent(out) :: dWTdr2

    real(kind=8), dimension(size(Weight)) :: Raw
    ! by defining the dimension of Raw to neigh%numInner, this restricts the
    ! size of array to the correct size w.r.t. each iteration within the dynamical array
    !real(kind=8), dimension(neigh%numInner) :: Raw
    real(kind=8), dimension(sys%nbond,size(Weight)) :: DWeight
    real(kind=8), dimension(size(Weight),sys%nbond) :: dTaydR
    real(kind=8), dimension(sys%nint,size(Weight)) :: DTay
    real(kind=8), dimension(sys%nbond) :: SumDWeight
    real(kind=8), dimension(sys%nbond) :: SumDWeightSqr
    real(kind=8), dimension(sys%nbond) :: Sumd2veightdr2tmp1
    real(kind=8), dimension(sys%nbond) :: Sumd2veightdr2tmp2
    real(kind=8) :: totsum, energy, temp, temp2, temp3

    ! for d2Vdx2
    real(kind=8), dimension(size(Weight),sys%nbond) :: d2TaydR2tmp1 ! intent(out)
    real(kind=8), dimension(size(Weight),sys%nbond) :: d2TaydR2tmp2 ! intent(out)

    ! for d2veightdx2
    real(kind=8), dimension(sys%nbond,size(Weight)), intent(out) :: d2veightdr2tmp1 ! intent(out)
    real(kind=8), dimension(sys%nbond,size(Weight)), intent(out) :: d2veightdr2tmp2 ! intent(out)

    ! 2 terms for Tayd2veightdr2
    real(kind=8), dimension(sys%nbond) :: Tayd2veightdr2tmp1 ! intent(out)
    real(kind=8), dimension(sys%nbond) :: Tayd2veightdr2tmp2 ! intent(out)

    ! 2 terms for WTaySumD2veightDr1
    real(kind=8), dimension(sys%nbond), intent(out) :: WTaySumD2veightDr1 ! intent(out)
    real(kind=8), dimension(sys%nbond), intent(out) :: WTaySumD2veightDr2 ! intent(out)

    ! Weight * Tay * (Sum dvdr)**2
    real(kind=8), dimension(sys%nbond), intent(out) :: WTSumDWeightDrSqr  ! intent(out)

    ! (Sum Tay*DWeight)*SumDWeight
    real(kind=8), dimension(sys%nbond), intent(out) :: TDWeightSumDWeight  ! intent(out)

    !Stores the value of the taylor series expansions
    real(kind=8), dimension(interp%ndata) :: Tay
    !Stores the zeta coordinates of the system
    real(kind=8), dimension(sys%nint,interp%ndata) :: z
    !Stores the d(eta)/dr
    real(kind=8), dimension(sys%nint,interp%ndata) :: detadr
    integer :: i,j,k
    
    !---------------------------------------------------
    !  Calculate the Weights 
    !---------------------------------------------------
    !$OMP PARALLEL DO PRIVATE(i) SHARED(Raw, Weight)
    do i=1,neigh%numInner
        Raw(i) = Weight(neigh%inner(i))
    enddo
    !$OMP END PARALLEL DO 
      
    totsum = sum(Raw(1:neigh%numInner))
    Weight = Raw/totsum
   
    !---------------------------------------------------
    !  Calculate the derivatives of the weights
    !---------------------------------------------------
   
    ! derivative of raw weights
    !!$OMP PARALLEL DO PRIVATE(i,k) SHARED(DWeight,temp,r)
    do k=1,neigh%numInner
        temp = interp%ipow2*Weight(k)*RawWeightTemp(k)
        temp2 = 2.0*temp*RawWeightTemp(k)*(interp%ipow+1.0)
        temp3 = -temp
        do i=1,sys%nbond
            DWeight(i,k) = temp*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**2
            d2veightdr2tmp1(i,k) = temp2*(r(i) - pot(neigh%inner(k))%r(i))**2*r(i)**4 
            d2veightdr2tmp2(i,k) = temp3*(r(i)**4 - (r(i) - pot(neigh%inner(k))%r(i))*r(i)**2)
        enddo
    enddo
    !!$OMP END PARALLEL DO
    
    SumDWeight = 0.0
    Sumd2veightdr2tmp1 = 0.0
    Sumd2veightdr2tmp2 = 0.0
    !!$OMP PARALLEL DO PRIVATE(i,k) SHARED(SumDWeight,DWeight)
    do k=1,neigh%numInner
        do i=1,sys%nbond
            SumDWeight(i) = SumDWeight(i) + DWeight(i,k)
            Sumd2veightdr2tmp1(i) = Sumd2veightdr2tmp1(i) + d2veightdr2tmp1(i,k)
            Sumd2veightdr2tmp2(i) = Sumd2veightdr2tmp2(i) + d2veightdr2tmp2(i,k)
        enddo
    enddo
    !!$OMP END PARALLE DO   

    SumDWeightSqr = 0.0
    do i=1,sys%nbond
        SumDWeightSqr(i) = SumDWeight(i)**2
    enddo

    ! derivative of relative weights
    !!$OMP PARALLEL DO PRIVATE(i,k) SHARED(DWeight,Weight,SumDWeight)
    do k=1,neigh%numInner
        do i=1,sys%nbond
            DWeight(i,k) = DWeight(i,k) - Weight(k)*SumDWeight(i)
        enddo
    enddo

    !!$OMP END PARALLEL DO
    !---------------------------------------------------
    !  Evaluate (z-z0) the local internal coordinates
    !---------------------------------------------------
    z = 0.0
    detadr = 0.0
    do k=1,neigh%numInner
        do j = 1,sys%nbond
            do i = 1,sys%nint
                z(i,k) = z(i,k) + pot(neigh%inner(k))%ut(i,j)*r(j)
                detadr(i,k) = detadr(i,k) + pot(neigh%inner(k))%ut(i,j)*r(j)**2
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

    do j=1,sys%nbond
       do i=1,neigh%numInner 
          d2TaydR2tmp1(i,j) = -2.0*dTaydR(i,j)*r(j)
          d2TaydR2tmp2(i,j) = dTaydR(i,j)
       enddo
    enddo

    ! 2 terms for Tay*d2veightdr2
    ! 2 terms for Tay*Weight*(Sum d2v/dr2)
    Tayd2veightdr2tmp1 = 0.0
    Tayd2veightdr2tmp2 = 0.0
    WTaySumD2veightDr1 = 0.0
    WTaySumD2veightDr2 = 0.0
    WTSumDWeightDrSqr = 0.0
    TDWeightSumDWeight = 0.0
    do k=1, neigh%numInner
       do i=1,sys%nbond
          Tayd2veightdr2tmp1(i) = Tayd2veightdr2tmp1(i) + Tay(k)*d2veightdr2tmp1(i,k)
          Tayd2veightdr2tmp2(i) = Tayd2veightdr2tmp2(i) + Tay(k)*d2veightdr2tmp2(i,k)

          WTaySumD2veightDr1(i) = WTaySumD2veightDr1(i) - Weight(k)*Tay(k)*Sumd2veightdr2tmp1(i)
          WTaySumD2veightDr2(i) = WTaySumD2veightDr2(i) - Weight(k)*Tay(k)*Sumd2veightdr2tmp2(i)

          ! weight*Tay*SumDWeight**2 remember to *2 later!
          WTSumDWeightDrSqr(i) = WTSumDWeightDrSqr(i) + Weight(k)*Tay(k)*SumDWeightSqr(i)

          ! (Sum Tay*DWeight)*SumDWeight
          TDWeightSumDWeight(i) =  TDWeightSumDWeight(i) - Tay(k)*DWeight(i,k)*SumDWeight(i)
       enddo
    enddo
    WTSumDWeightDrSqr = 2.0*WTSumDWeightDrSqr

    ! dWeight/dr * dTay/dr
    dWTdr2 = 0.0
    do k=1,neigh%numInner
       do i=1,sys%nbond
          dWTdr2(i) = dWTdr2(i) + 2.0*DWeight(i,k)*dTaydR(k,i)
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
