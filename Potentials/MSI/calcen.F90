

subroutine calcen(sys,interp,pot,neigh,Weight,r,V,dVdR,RawWeightTemp,dWTdr2,&
                d2TaydR2tmp1,&
                drdx,d2rdx2,d2rdxmdxn,d2TaydR2tmp2_mb,d2TaydR2tmp2_nb) 
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
    !real(kind=8), dimension(sys%nbond) :: Sumd2veightdr2tmp1
    !real(kind=8), dimension(sys%nbond) :: Sumd2veightdr2tmp2
    real(kind=8) :: totsum, energy, temp, temp2, temp3

    ! For d2veightdx2
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2vdrtmp1
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2vdrtmp2
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2vdrtmp3
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2vdrtmp4
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2vdr2tmp1

    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen) :: d2vdx2tmp1
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen) :: d2vdx2tmp2
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen) :: d2vdxmdxntmp1
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen) :: d2vdxmdxntmp2

    ! For sum_{j=1}^{Ndata} d2veightdx2
    real(kind=8), dimension(sys%dimen,sys%dimen) :: Sumd2veightdx2 
    real(kind=8), dimension(sys%dimen,sys%dimen) :: Sumd2veightdxmdxn 

    ! for d2Taydx2
    real(kind=8), dimension(sys%nbond), intent(out) :: d2TaydR2tmp1 ! intent(out)
    real(kind=8), dimension(sys%dimen,sys%nbond), intent(out) :: d2TaydR2tmp2_mb ! intent(out)
    real(kind=8), dimension(sys%dimen,sys%nbond), intent(out) :: d2TaydR2tmp2_nb ! intent(out)

    ! for d2veightdx2
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2veightdr2tmp1 ! intent(out)
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2veightdr2tmp2 ! intent(out)

    !derivative of bond lengths with respect to cartesians
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(in) :: drdx
    ! 2nd derivative of bondlengths w.r.t Cartesians
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%nbond), intent(in) :: d2rdx2
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%nbond), intent(in) :: d2rdxmdxn

    ! Sum_{j=1}^{Ndata} dveight/dx
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: Sumdveightdxm
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: Sumdveightdxn

    ! (Sum_{j=1}^{Ndata} dvdx)*(Sum_{j=1}^{Ndata} dvdx)
    real(kind=8), dimension(sys%dimen,sys%dimen) :: sqrSumdveightdx
    real(kind=8), dimension(sys%dimen,sys%dimen) :: sqrSumdveightdxmdxn

    ! Sum_{j=1}^{Ndata} Tay_j * d2veight/dx2
    real(kind=8), dimension(sys%dimen,sys%dimen) :: SumTayd2veightdx2
    real(kind=8), dimension(sys%dimen,sys%dimen) :: SumTayd2veightdxmdxn

    ! Sum_{j=1}^{Ndata} Tay*dvdx
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: SumTaydvdxm
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: SumTaydvdxn

    ! Sum_{j=1}^{Ndata} Tay*dvdx * Sum_{j=1}^{Ndata} dvdx
    real(kind=8), dimension(sys%dimen,sys%dimen) :: SumTaydvmSumdvm
    real(kind=8), dimension(sys%dimen,sys%dimen) :: SumTaydvmSumdvn

    ! d2Vdx2 
    real(kind=8), dimension(sys%dimen,sys%dimen) :: d2Vdx2
    real(kind=8), dimension(sys%dimen,sys%dimen) :: d2Vdxmdxn

    !Stores the value of r^2*drdx
    !real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: r2drdx
    real(kind=8), dimension(sys%dimen,sys%natom) :: r2drdx_mb
    real(kind=8), dimension(sys%dimen,sys%natom) :: r2drdx_nb

    ! ut*r^2*drdx
    real(kind=8), dimension(size(Weight),sys%nint,sys%dimen) :: utr2drdx_mb 
    real(kind=8), dimension(size(Weight),sys%nint,sys%dimen) :: utr2drdx_nb 
    ! second term in d2Taydx2, r^2 * (ut*r^2*drdx) * ut
    real(kind=8), dimension(size(Weight),sys%nbond,sys%dimen) :: v2detadxut_mb
    real(kind=8), dimension(size(Weight),sys%nbond,sys%dimen) :: v2detadxut_nb

    !Stores the value of the taylor series expansions
    real(kind=8), dimension(interp%ndata) :: Tay
    !Stores the zeta coordinates of the system
    real(kind=8), dimension(sys%nint,interp%ndata) :: z
    integer :: i,j,k,l
    
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
        temp2 = (2.0+interp%ipow2)*temp*RawWeightTemp(k)
        do i=1,sys%nbond
            DWeight(i,k) = temp*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**2

            ! For evaluating d2vdx2
            d2vdrtmp1(i,k) = temp2*((r(i) - pot(neigh%inner(k))%r(i))*r(i)**2)**2
            d2vdrtmp2(i,k) = -temp*r(i)**4
            d2vdrtmp3(i,k) = -2.0*temp*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**3
            d2vdrtmp4(i,k) = DWeight(i,k)
            
            d2vdr2tmp1(i,k) = d2vdrtmp1(i,k) + d2vdrtmp2(i,k) + d2vdrtmp3(i,k)
        enddo
    enddo
    !!$OMP END PARALLEL DO

    ! d2veightdx2
    d2vdx2tmp1 = 0.0
    d2vdxmdxntmp1 = 0.0
    d2vdx2tmp2 = 0.0
    d2vdxmdxntmp2 = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do i=1,sys%nbond
                d2vdx2tmp1(k,l,j) = d2vdx2tmp1(k,l,j) + d2vdr2tmp1(i,k)*drdx(j,sys%mb(i),i)*drdx(l,sys%mb(i),i)
                d2vdxmdxntmp1(k,l,j) = d2vdxmdxntmp1(k,l,j) + d2vdr2tmp1(i,k)*drdx(j,sys%mb(i),i)*drdx(l,sys%nb(i),i)

                d2vdx2tmp2(k,l,j) = d2vdx2tmp2(k,l,j) + d2vdrtmp4(i,k)*d2rdx2(l,j,i)
                d2vdxmdxntmp2(k,l,j) = d2vdxmdxntmp2(k,l,j) + d2vdrtmp4(i,k)*d2rdxmdxn(l,j,i)
             enddo
          enddo
       enddo
    enddo

    ! Sum_{j=1}^{Ndata} d2veightdx2
    Sumd2veightdx2 = 0.0
    !Sumd2veightdxmdxn = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%dimen
             Sumd2veightdx2(j,i) = Sumd2veightdx2(j,i) + d2vdx2tmp1(k,j,i) + d2vdx2tmp2(k,j,i)
             !Sumd2veightdxmdxn(j,i) = Sumd2veightdxmdxn(j,i) + d2vdxmdxntmp1(k,j,i) + d2vdxmdxntmp2(k,j,i)
          enddo
       enddo
    enddo
    Sumd2veightdxmdxn = - Sumd2veightdx2
    
    SumDWeight = 0.0
    !!$OMP PARALLEL DO PRIVATE(i,k) SHARED(SumDWeight,DWeight)
    do k=1,neigh%numInner
        do i=1,sys%nbond
            SumDWeight(i) = SumDWeight(i) + DWeight(i,k)
        enddo
    enddo
    !!$OMP END PARALLE DO   

    ! sum_{j=1}^{Ndata} dveight_j/dx
    Sumdveightdxm = 0.0
    !Sumdveightdxn = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
             Sumdveightdxm(j,sys%mb(i),i) = Sumdveightdxm(j,sys%mb(i),i) + DWeight(i,k)*drdx(j,sys%mb(i),i)
             !Sumdveightdxn(j,sys%nb(i),i) = Sumdveightdxn(j,sys%nb(i),i) + DWeight(i,k)*drdx(j,sys%nb(i),i)
          enddo
       enddo
    enddo
    do j=1,sys%dimen
       do i=1,sys%nbond
          Sumdveightdxn(j,sys%nb(i),i) = -Sumdveightdxm(j,sys%mb(i),i)
       enddo
    enddo

    ! sum_{j=1}^{Ndata} dveight_j/dx * sum_{j=1}^{Ndata} dveight_j/dx = 3*3 matrix
    sqrSumdveightdx = 0.0
    do k=1,sys%dimen
       do j=1,sys%dimen
          do i=1,sys%nbond
             sqrSumdveightdx(k,j) = sqrSumdveightdx(k,j) + Sumdveightdxm(k,sys%mb(i),i)*Sumdveightdxm(j,sys%mb(i),i)
          enddo
       enddo
    enddo
    sqrSumdveightdxmdxn = -sqrSumdveightdx

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
    
    !--------------------------------------------------
    !  Evaluate v2*ut*r^2*drdx
    !--------------------------------------------------
    do j=1,sys%dimen
       do i=1,sys%nbond
          r2drdx_mb(j,sys%mb(i)) = drdx(j,sys%mb(i),i)*r(i)**2
          r2drdx_nb(j,sys%nb(i)) = drdx(j,sys%nb(i),i)*r(i)**2
       enddo
    enddo
    !do i=1,sys%nbond
    !   print *, r2drdx(:,sys%mb(i),i)
    !enddo
    !call exit(0)
    !print *, sys%mb
    !print *, sys%nb
    !print *, '_________________________________________________________________'

    ! ut*r2drdx matrix multiplication 
    utr2drdx_mb = 0.0
    utr2drdx_nb = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%nint
             do i=1,sys%nbond
                utr2drdx_mb(k,j,l) = utr2drdx_mb(k,j,l) + pot(neigh%inner(k))%ut(j,i)*r2drdx_mb(l,sys%mb(i)) 
                utr2drdx_nb(k,j,l) = utr2drdx_nb(k,j,l) + pot(neigh%inner(k))%ut(j,i)*r2drdx_nb(l,sys%nb(i)) 
             enddo
             !print *, utr2drdx_mb(k,j,l), utr2drdx_nb(k,j,l)
          enddo
          !print *, ''
       enddo
       !call exit(0)
    enddo

    ! v2* (ut*r^2*drdx) * ut
    v2detadxut_mb = 0.0
    v2detadxut_nb = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do l=1,sys%nbond
             do i=1,sys%nint
                v2detadxut_mb(k,l,j) = v2detadxut_mb(k,l,j) + Weight(k)*pot(neigh%inner(k))%v2(i)*&
                        utr2drdx_mb(k,i,j)*pot(neigh%inner(k))%ut(i,l)
                v2detadxut_nb(k,l,j) = v2detadxut_nb(k,l,j) + Weight(k)*pot(neigh%inner(k))%v2(i)*&
                        utr2drdx_nb(k,i,j)*pot(neigh%inner(k))%ut(i,l)
             enddo
             v2detadxut_mb(k,l,j) = v2detadxut_mb(k,l,j)*r(l)**2
             v2detadxut_nb(k,l,j) = v2detadxut_nb(k,l,j)*r(l)**2
          enddo
       enddo
    enddo

    d2TaydR2tmp2_mb = 0.0
    d2TaydR2tmp2_nb = 0.0
    do k=1,sys%nbond
       do j=1, sys%dimen
          do i=1,neigh%numInner
             d2TaydR2tmp2_mb(j,k) = d2TaydR2tmp2_mb(j,k) + v2detadxut_mb(i,k,j)
             d2TaydR2tmp2_nb(j,k) = d2TaydR2tmp2_nb(j,k) + v2detadxut_nb(i,k,j)
          enddo
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
          d2TaydR2tmp1(j) = d2TaydR2tmp1(j) - Weight(i)*dTaydR(i,j)
          !d2TaydR2tmp1(j) = d2TaydR2tmp1(j) - Weight(i)*dTaydR(i,j)
          !d2TaydR2tmp2(j) = d2TaydR2tmp2(j) + Weight(i)*dTaydR(i,j)
       enddo
       !d2TaydR2tmp1(j) = 2.0*d2TaydR2tmp1(j)*r(j)
       d2TaydR2tmp1(j) = (2.0*r(j)-1.0)*d2TaydR2tmp1(j)
    enddo

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

    !----------------------------------------------------------
    !  Calculate the Sum_{j=1}^{Ndata} (Tay_j * d2weight_j/dx2)
    !----------------------------------------------------------
    SumTayd2veightdx2 = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%dimen
           SumTayd2veightdx2(j,i) = SumTayd2veightdx2(j,i) + (d2vdx2tmp1(k,j,i) + d2vdx2tmp2(k,j,i))*Tay(k)
           !SumTayd2veightdxmdxn(j,i) = SumTayd2veightdxmdxn(j,i) + (d2vdxmdxntmp1(k,j,i) + d2vdxmdxntmp2(k,j,i))*Tay(k)
          enddo
       enddo
    enddo
    SumTayd2veightdxmdxn = -SumTayd2veightdx2 
    
    !---------------------------------------------------------------------------
    !  Calculate the (Sum_j w_j*T_j) *Sum_{j=1}^{Ndata} (Tay_j * d2weight_j/dx2)
    !---------------------------------------------------------------------------
    Sumd2veightdx2 = -Sumd2veightdx2*V
    Sumd2veightdxmdxn = -Sumd2veightdxmdxn*V

    !------------------------------------------------------------------------------------
    !  Calculate the (Sum_{j=1}^{Ndata} -Tay_j * dveightdx)*(Sum_{j=1}^{Ndata} dveightdx)
    !------------------------------------------------------------------------------------
    SumTaydvdxm = 0.0
    !SumTaydvdxn = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
              SumTaydvdxm(j,sys%mb(i),i) = SumTaydvdxm(j,sys%mb(i),i) - Tay(k)*DWeight(i,k)*drdx(j,sys%mb(i),i)
              !SumTaydvdxn(j,i) = SumTaydvdxn(j,i) - Tay(k)*DWeight(i,k)*drdx(j,sys%nb(i),i)
          enddo
       enddo
    enddo
    do j=1,sys%dimen
       do i=1,sys%nbond
          SumTaydvdxn(j,sys%nb(i),i) = -SumTaydvdxm(j,sys%mb(i),i)
       enddo
    enddo

    SumTaydvmSumdvm = 0.0
    !SumTaydvmSumdvn = 0.0
    do k=1,sys%dimen
       do j=1,sys%dimen
          do i=1,sys%nbond
             SumTaydvmSumdvm(k,j) = SumTaydvmSumdvm(k,j)  - SumTaydvdxm(k,sys%mb(i),i)*Sumdveightdxm(j,sys%mb(i),i)
             !SumTaydvmSumdvn(k,j) = SumTaydvmSumdvn(k,j)  - SumTaydvdxm(k,sys%mb(i),i)*Sumdveightdxn(j,sys%nb(i),i)
          enddo
       enddo
    enddo
    SumTaydvmSumdvn = -SumTaydvmSumdvm

    !---------------------------------------------------------------------------
    !  Calculate (Sum_{j=1}^{Ndata} 2.0 * w_j * T_j) * (Sum dvdx)*(Sum dvdx)
    !---------------------------------------------------------------------------
    sqrSumdveightdx = 2.0*V*sqrSumdveightdx    !FIXME 
    sqrSumdveightdxmdxn = 2.0*V*sqrSumdveightdxmdxn   !FIXME
    
    !---------------------------------------------------------------------------
    !  Calculate (Sum_{j=1}^{Ndata} 2.0 * w_j * T_j) * (Sum dvdx)*(Sum dvdx)
    !---------------------------------------------------------------------------
    d2Vdx2 = SumTayd2veightdx2 + Sumd2veightdx2 + 2.0*SumTaydvmSumdvm + sqrSumdveightdx
    d2Vdxmdxn = SumTayd2veightdxmdxn + Sumd2veightdxmdxn + 2.0*SumTaydvmSumdvn + sqrSumdveightdxmdxn 

  return
end subroutine
