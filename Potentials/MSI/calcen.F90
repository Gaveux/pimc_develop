

subroutine calcen(sys,interp,pot,neigh,Weight,r,V,dVdR,RawWeightTemp,&
                drdx,d2rdx2,d2Vdx2) 
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
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2vdr2tmp2

    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2vdx2tmp1
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2vdx2tmp2
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2vdx2tmp3

    real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: d2veightdxtmp1
    real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: d2veightdxtmp2

    ! For sum_{j=1}^{Ndata} d2veightdx2
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: Sumd2veightdx2 

    ! For d2Taydx2
    real(kind=8), dimension(size(Weight),sys%nbond) :: d2Taydx2tmp1
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Taydx2_tmp1
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Taydx2_tmp2
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Taydx2_tmp3
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: SumWd2Taydx2

    ! For d2veightdx2
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2veightdr2tmp1 ! intent(out)
    real(kind=8), dimension(sys%nbond,size(Weight)) :: d2veightdr2tmp2 ! intent(out)

    ! Derivative of bond lengths with respect to cartesians
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(in) :: drdx
    ! 2nd derivative of bondlengths w.r.t Cartesians
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom,sys%nbond), intent(in) :: d2rdx2

    ! Sum_{j=1}^{Ndata} dveight/dx
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: Sumdveightdxm
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: Sumdveightdxn
    real(kind=8), dimension(sys%dimen,sys%natom) :: Sumdveightdx

    ! (Sum_{j=1}^{Ndata} dvdx)*(Sum_{j=1}^{Ndata} dvdx)
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: sqrSumdveightdx

    ! Sum_{j=1}^{Ndata} Tay_j * d2veight/dx2
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: SumTayd2veightdx2

    ! Sum_{j=1}^{Ndata} Tay*dvdx
    real(kind=8), dimension(sys%dimen,sys%natom) :: SumTaydvdx

    ! Sum_{j=1}^{Ndata} Tay*dvdx * Sum_{j=1}^{Ndata} dvdx
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: SumTaydvSumdv
    !------------------------------------------------------------------------
    ! Sum_{j=1}^{Ndata} Tay*dvdx_m * Sum_{j=1}^{Ndata} dvdx_n -
    ! Sum_{j=1}^{Ndata} Tay*dvdx_n * Sum_{j=1}^{Ndata} dvdx_m 
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: SumTaydvSumdv2
    !------------------------------------------------------------------------

    ! d2Weightdx2 
    real(kind=8), dimension(sys%dimen,sys%dimen) :: d2Weightdx2    !intent(out)
    real(kind=8), dimension(sys%dimen,sys%dimen) :: d2Weightdxmdxn !intent(out)

    ! d2Vdx2
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom), intent(out) :: d2Vdx2

    !DWeightdxm
    real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: DWeightdx
    !dTaydxm and dTaydxn
    real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: dTaydx
    ! dWTdx2
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: dWTdx2

    ! Sum_{l=1}^{N choose 2} UT*r^2 * drdx
    real(kind=8), dimension(size(Weight),sys%nint,sys%dimen,sys%natom) :: SumUTrsqrdrdx
    ! v2 * ut * detadx
    real(kind=8), dimension(size(Weight),sys%nbond,sys%dimen,sys%natom) :: v2utdetadx

    !Stores the value of the taylor series expansions
    real(kind=8), dimension(interp%ndata) :: Tay
    !Stores the zeta coordinates of the system
    real(kind=8), dimension(sys%nint,interp%ndata) :: z
    integer :: i,j,k,l,m
    
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
            d2vdrtmp1(i,k) = temp2*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**2
            d2vdrtmp2(i,k) = -temp*r(i)**4
            d2vdrtmp3(i,k) = -2.0*temp*(r(i) - pot(neigh%inner(k))%r(i))*r(i)**3
            d2vdrtmp4(i,k) = DWeight(i,k)
            
            d2vdr2tmp2(i,k) = d2vdrtmp2(i,k) + d2vdrtmp3(i,k)
        enddo
    enddo
    !!$OMP END PARALLEL DO

    ! For d2vdx2tmp1
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
             d2veightdxtmp1(k,j,sys%mb(i)) = d2veightdxtmp1(k,j,sys%mb(i)) + &
                   d2vdrtmp1(i,k)*drdx(j,sys%mb(i),i)
             d2veightdxtmp1(k,j,sys%nb(i)) = d2veightdxtmp1(k,j,sys%nb(i)) + &
                   d2vdrtmp1(i,k)*drdx(j,sys%nb(i),i)

             ! Calculates (z_k-z_j)*dz_k/dx
             d2veightdxtmp2(k,j,sys%mb(i)) = d2veightdxtmp2(k,j,sys%mb(i)) + &
                   (r(i)-pot(neigh%inner(k))%r(i))*r(i)**2*drdx(j,sys%mb(i),i)
             d2veightdxtmp2(k,j,sys%nb(i)) = d2veightdxtmp2(k,j,sys%nb(i)) + &
                   (r(i)-pot(neigh%inner(k))%r(i))*r(i)**2*drdx(j,sys%nb(i),i)
          enddo
       enddo
    enddo

    ! d2veightdx2 a (3*3*5*5*numInner) matrix
    d2vdx2tmp1 = 0.0 
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do m=1,sys%natom
                do i=1,sys%natom
                   d2vdx2tmp1(k,l,j,m,i) = d2veightdxtmp1(k,l,m)*d2veightdxtmp2(k,j,i) 
                enddo
             enddo
          enddo
       enddo
    enddo

    d2vdx2tmp2 = 0.0
    d2vdx2tmp3 = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do i=1,sys%nbond
                !same dr_i/dxm dr_i/dxn same r_i
                d2vdx2tmp2(k,l,j,sys%mb(i),sys%mb(i)) = d2vdx2tmp2(k,l,j,sys%mb(i),sys%mb(i)) + &
                             d2vdr2tmp2(i,k)*drdx(l,sys%mb(i),i)*drdx(j,sys%mb(i),i)
                d2vdx2tmp2(k,l,j,sys%nb(i),sys%mb(i)) = d2vdx2tmp2(k,l,j,sys%nb(i),sys%mb(i)) + &
                             d2vdr2tmp2(i,k)*drdx(l,sys%nb(i),i)*drdx(j,sys%mb(i),i)
                d2vdx2tmp2(k,l,j,sys%mb(i),sys%nb(i)) = d2vdx2tmp2(k,l,j,sys%mb(i),sys%nb(i)) + &
                             d2vdr2tmp2(i,k)*drdx(l,sys%mb(i),i)*drdx(j,sys%nb(i),i)
                d2vdx2tmp2(k,l,j,sys%nb(i),sys%nb(i)) = d2vdx2tmp2(k,l,j,sys%nb(i),sys%nb(i)) + &
                             d2vdr2tmp2(i,k)*drdx(l,sys%nb(i),i)*drdx(j,sys%nb(i),i)

                d2vdx2tmp3(k,l,j,sys%mb(i),sys%mb(i)) = d2vdx2tmp3(k,l,j,sys%mb(i),sys%mb(i)) + &
                             d2vdrtmp4(i,k)*d2rdx2(l,j,sys%mb(i),sys%mb(i),i)
                d2vdx2tmp3(k,l,j,sys%nb(i),sys%mb(i)) = d2vdx2tmp3(k,l,j,sys%nb(i),sys%mb(i)) + &
                             d2vdrtmp4(i,k)*d2rdx2(l,j,sys%nb(i),sys%mb(i),i)
                d2vdx2tmp3(k,l,j,sys%mb(i),sys%nb(i)) = d2vdx2tmp3(k,l,j,sys%mb(i),sys%nb(i)) + &
                             d2vdrtmp4(i,k)*d2rdx2(l,j,sys%mb(i),sys%nb(i),i)
                d2vdx2tmp3(k,l,j,sys%nb(i),sys%nb(i)) = d2vdx2tmp3(k,l,j,sys%nb(i),sys%nb(i)) + &
                             d2vdrtmp4(i,k)*d2rdx2(l,j,sys%nb(i),sys%nb(i),i)
                !d2vdx2tmp2(k,l,j) = d2vdx2tmp2(k,l,j) + d2vdrtmp4(i,k)*d2rdx2(l,j,i)
                !d2vdxmdxntmp2(k,l,j) = d2vdxmdxntmp2(k,l,j) + d2vdrtmp4(i,k)*d2rdxmdxn(l,j,i)
             enddo
          enddo
       enddo
    enddo

    !------------------------------------------
    ! Sum_{j=1}^{Ndata} d2veightdx2
    !------------------------------------------
    Sumd2veightdx2 = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do m=1,sys%natom
                do i=1,sys%natom
                   Sumd2veightdx2(l,j,m,i) = Sumd2veightdx2(l,j,m,i) + d2vdx2tmp1(k,l,j,m,i) &
                        + d2vdx2tmp2(k,l,j,m,i) + d2vdx2tmp3(k,l,j,m,i)
                enddo
             enddo
          enddo
       enddo
    enddo
             !do l=1,sys%natom
             !   do i=1,sys%natom
             !      print *, Sumd2veightdx2(:,:,l,i)
             !      print *, ''
             !   enddo
             !   print *, '----------'
             !enddo
             !call exit(0)
    
    SumDWeight = 0.0
    !!$OMP PARALLEL DO PRIVATE(i,k) SHARED(SumDWeight,DWeight)
    do k=1,neigh%numInner
        do i=1,sys%nbond
            SumDWeight(i) = SumDWeight(i) + DWeight(i,k)
        enddo
    enddo
    !!$OMP END PARALLE DO   

    !------------------------------------------
    ! sum_{j=1}^{Ndata} dveight_j/dx
    !------------------------------------------
    Sumdveightdx = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
             Sumdveightdx(j,sys%mb(i)) = Sumdveightdx(j,sys%mb(i)) + DWeight(i,k)*drdx(j,sys%mb(i),i)
             Sumdveightdx(j,sys%nb(i)) = Sumdveightdx(j,sys%nb(i)) + DWeight(i,k)*drdx(j,sys%nb(i),i)
          enddo
       enddo
    enddo
             !do i=1,sys%dimen
             !   do j=1,sys%dimen
             !      print *, Sumdveightdx(i,5)*Sumdveightdx(j,5)
             !   enddo
             !   print *, ''
             !enddo
             !print *, '--------------------------------'

    !---------------------------------------------------------------------
    ! (sum_{j=1}^{Ndata} dveight_j/dx) * (sum_{j=1}^{Ndata} dveight_j/dx) 
    !---------------------------------------------------------------------
    sqrSumdveightdx = 0.0
    do l=1,sys%dimen
       do k=1,sys%dimen
          do j=1,sys%natom
             do i=1,sys%natom
                sqrSumdveightdx(l,k,j,i) = Sumdveightdx(l,j)*Sumdveightdx(k,i)
             enddo
          enddo
       enddo
    enddo
            !do i=1,sys%natom
            !    print *, sqrSumdveightdx(:,:,5,5)
            !    print *, ''
            !enddo
            !call exit(0)

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
    
    !==================================================================
    !===========================d2Taydx2===============================
    !==================================================================

    !------------------------------------------------------------
    !  Evaluate Sum_{l=1}^{N choose 2} ut*r^2*drdx  For d2Taydx2
    !------------------------------------------------------------
    SumUTrsqrdrdx = 0.0
    do l=1,neigh%numInner
       do k=1,sys%nint
          do j=1,sys%dimen
             do i=1,sys%nbond
                SumUTrsqrdrdx(l,k,j,sys%mb(i)) = SumUTrsqrdrdx(l,k,j,sys%mb(i)) + &
                      pot(neigh%inner(l))%ut(k,i)*r(i)**2*drdx(j,sys%mb(i),i)
                SumUTrsqrdrdx(l,k,j,sys%nb(i)) = SumUTrsqrdrdx(l,k,j,sys%nb(i)) + &
                      pot(neigh%inner(l))%ut(k,i)*r(i)**2*drdx(j,sys%nb(i),i)
             enddo
          enddo
       enddo
    enddo


    !-------------------------------------------------------
    ! Sum_{i=1}^{3N-6} SumUTrsqrdrdx * v2_i * ut_il
    !-------------------------------------------------------
    v2utdetadx = 0.0
    do j=1,sys%dimen
       do i=1,sys%natom
          do l=1,neigh%numInner
             do m=1,sys%nbond
                do k=1,sys%nint
                   v2utdetadx(l,m,j,i) = v2utdetadx(l,m,j,i) + &
                      SumUTrsqrdrdx(l,k,j,i)*pot(neigh%inner(l))%v2(k)*pot(neigh%inner(l))%ut(k,m)
                enddo
             enddo
          enddo
       enddo
    enddo

    !-------------------------------------------------------
    ! sum_{l=1}^{N choose 2} v2utdetadx * r_l dr_ldx
    !-------------------------------------------------------
    d2Taydx2_tmp2 = 0.0
    do l=1,neigh%numInner
       do k=1,sys%dimen
          do j=1,sys%dimen
             do m=1,sys%natom
                do i=1,sys%nbond
                   d2Taydx2_tmp2(l,k,j,m,sys%mb(i)) = d2Taydx2_tmp2(l,k,j,m,sys%mb(i)) + &
                       v2utdetadx(l,i,k,m)*r(i)**2*drdx(j,sys%mb(i),i)
                   d2Taydx2_tmp2(l,k,j,m,sys%nb(i)) = d2Taydx2_tmp2(l,k,j,m,sys%nb(i)) + &
                       v2utdetadx(l,i,k,m)*r(i)**2*drdx(j,sys%nb(i),i)
                enddo
             enddo
          enddo
       enddo
    enddo

    !=================================================================
    !=================================================================
    !=================================================================


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

    !---------------------------------------------------------
    ! Calculate the 1st term of d2Taydx2
    !---------------------------------------------------------

    do j=1,neigh%numInner
       do i=1,sys%nbond
          d2Taydx2tmp1(j,i) = -2.0*dTaydR(j,i)*r(i)
       enddo
    enddo

    ! dr_ldx*dr_ldx same l
    d2Taydx2_tmp1 = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do i=1,sys%nbond
                 d2Taydx2_tmp1(k,l,j,sys%mb(i),sys%mb(i)) = d2Taydx2_tmp1(k,l,j,sys%mb(i),sys%mb(i)) + &
                        d2Taydx2tmp1(k,i)*drdx(l,sys%mb(i),i)*drdx(j,sys%mb(i),i)
                 d2Taydx2_tmp1(k,l,j,sys%mb(i),sys%nb(i)) = d2Taydx2_tmp1(k,l,j,sys%mb(i),sys%nb(i)) + &
                        d2Taydx2tmp1(k,i)*drdx(l,sys%mb(i),i)*drdx(j,sys%nb(i),i)
                 d2Taydx2_tmp1(k,l,j,sys%nb(i),sys%nb(i)) = d2Taydx2_tmp1(k,l,j,sys%nb(i),sys%nb(i)) + &
                        d2Taydx2tmp1(k,i)*drdx(l,sys%nb(i),i)*drdx(j,sys%nb(i),i)
                 d2Taydx2_tmp1(k,l,j,sys%nb(i),sys%mb(i)) = d2Taydx2_tmp1(k,l,j,sys%nb(i),sys%mb(i)) + &
                        d2Taydx2tmp1(k,i)*drdx(l,sys%nb(i),i)*drdx(j,sys%mb(i),i)
             enddo
          enddo
       enddo
    enddo
               !do j=1,sys%dimen
               !   do i=1,sys%dimen
               !      print *,  d2Taydx2_tmp1(1,j,i,1,1)
               !   enddo
               !enddo
               !call exit(0)
    
    !---------------------------------------------------------
    ! Calculate the 3rd term of d2Taydx2
    !---------------------------------------------------------
    d2Taydx2_tmp3 = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do i=1,sys%nbond
                d2Taydx2_tmp3(k,l,j,sys%mb(i),sys%mb(i)) = d2Taydx2_tmp3(k,l,j,sys%mb(i),sys%mb(i)) + &
                       dTaydR(k,i)*d2rdx2(l,j,sys%mb(i),sys%mb(i),i)
                d2Taydx2_tmp3(k,l,j,sys%mb(i),sys%nb(i)) = d2Taydx2_tmp3(k,l,j,sys%mb(i),sys%nb(i)) + &
                       dTaydR(k,i)*d2rdx2(l,j,sys%mb(i),sys%nb(i),i)
                d2Taydx2_tmp3(k,l,j,sys%nb(i),sys%mb(i)) = d2Taydx2_tmp3(k,l,j,sys%nb(i),sys%mb(i)) + &
                       dTaydR(k,i)*d2rdx2(l,j,sys%nb(i),sys%mb(i),i)
                d2Taydx2_tmp3(k,l,j,sys%nb(i),sys%nb(i)) = d2Taydx2_tmp3(k,l,j,sys%nb(i),sys%nb(i)) + &
                       dTaydR(k,i)*d2rdx2(l,j,sys%nb(i),sys%nb(i),i)
             enddo
          enddo
       enddo
    enddo
    
    !---------------------------------------------------------
    ! Calculate Sum_{j=1}^{Ndata} weight_j * d2Taydx2
    !---------------------------------------------------------
    ! Sum_{j=1}^{Ndata} Weight_j * dTay_j/dx
    SumWd2Taydx2 = 0.0
    do m=1,sys%dimen
       do l=1,sys%dimen
          do k=1,sys%natom
             do j=1,sys%natom
                do i=1,neigh%numInner
                   SumWd2Taydx2(m,l,k,j) = SumWd2Taydx2(m,l,k,j) + &
                      Weight(i)*(d2Taydx2_tmp1(i,m,l,k,j) + d2Taydx2_tmp2(i,m,l,k,j) + &
                      d2Taydx2_tmp3(i,m,l,k,j))
                enddo
             enddo
          enddo
       enddo
    enddo

    !----------------------------------------------------------
    !  Calculate sum_{j}^{Ndata} 2.0 * dWeightdx * dTaydx
    !----------------------------------------------------------
    ! dWeightdx
    DWeightdx = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
            DWeightdx(k,j,sys%mb(i)) = DWeightdx(k,j,sys%mb(i)) + DWeight(i,k)*drdx(j,sys%mb(i),i)
            DWeightdx(k,j,sys%nb(i)) = DWeightdx(k,j,sys%nb(i)) + DWeight(i,k)*drdx(j,sys%nb(i),i)
          enddo
       enddo
    enddo

    ! dTaydx
    dTaydx = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
             dTaydx(k,j,sys%mb(i)) = dTaydx(k,j,sys%mb(i)) + dTaydR(k,i)*drdx(j,sys%mb(i),i)
             dTaydx(k,j,sys%nb(i)) = dTaydx(k,j,sys%nb(i)) + dTaydR(k,i)*drdx(j,sys%nb(i),i)
          enddo
       enddo
    enddo

    ! DWeightdx * dTaydx
    do m=1,neigh%numInner
       do l=1,sys%dimen
          do k=1,sys%dimen
             do j=1,sys%natom
                do i=1,sys%natom 
                   dWTdx2(l,k,j,i) = dWTdx2(l,k,j,i) + DWeightdx(m,l,j)*dTaydx(m,k,i)
                enddo
             enddo
          enddo
       enddo
    enddo
    dWTdx2 = 2.0*dWTdx2
                  !print *, dWTdx2
                  !call exit(0)

    !------------------------------------------
    ! Sum_{j=1}^{Ndata} Tay_j*d2veight_jdx2
    !------------------------------------------
    SumTayd2veightdx2 = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do m=1,sys%natom
                do i=1,sys%natom
                   SumTayd2veightdx2(l,j,m,i) = SumTayd2veightdx2(l,j,m,i) + Tay(k)*(d2vdx2tmp1(k,l,j,m,i) &
                        + d2vdx2tmp2(k,l,j,m,i) + d2vdx2tmp3(k,l,j,m,i))
                enddo
             enddo
          enddo
       enddo
    enddo
    
    !------------------------------------------------------------------------------------------------
    !  Calculate the (Sum_j w_j*T_j) * (Sum_{j=1}^{Ndata} d2weight_j/dx2) For Tay*d2Weightdx2
    !------------------------------------------------------------------------------------------------
    Sumd2veightdx2 = -Sumd2veightdx2*V

    !---------------------------------------------------------------------------------------------------------
    !  Calculate the (Sum_{j=1}^{Ndata} -Tay_j * dveightdx)*(Sum_{j=1}^{Ndata} dveightdx) For Tay*d2Weightdx2
    !---------------------------------------------------------------------------------------------------------
    SumTaydvdx = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
              SumTaydvdx(j,sys%mb(i)) = SumTaydvdx(j,sys%mb(i)) - Tay(k)*DWeight(i,k)*drdx(j,sys%mb(i),i)
              SumTaydvdx(j,sys%nb(i)) = SumTaydvdx(j,sys%nb(i)) - Tay(k)*DWeight(i,k)*drdx(j,sys%nb(i),i)
          enddo
       enddo
    enddo

    do l=1,sys%dimen
       do k=1,sys%dimen
          do j=1,sys%natom
             do i=1,sys%natom
                SumTaydvSumdv(l,k,j,i) = - SumTaydvdx(l,j)*Sumdveightdx(k,i)
             enddo
          enddo
       enddo
    enddo
    
    ! (Sum_{j=1}^{Ndata} Tay_j*dvdx_m)*(Sum_{j=1}^{Ndata} dvdx_n) + 
    ! (Sum_{j=1}^{Ndata} Tay_j*dvdx_n)*(Sum_{j=1}^{Ndata} dvdx_m)
    do l=1,sys%dimen
       do k=1,sys%dimen
          do j=1,sys%natom
             do i=1,sys%natom
                SumTaydvSumdv2(l,k,j,i) = SumTaydvSumdv(l,k,j,i) + SumTaydvSumdv(l,k,i,j)
             enddo
          enddo
       enddo
    enddo

    !------------------------------------------------------------------------------------------------------
    !  Calculate (Sum_{j=1}^{Ndata} 2.0 * w_j * T_j) * (Sum dvdx)*(Sum dvdx)  For Tay*d2Weightdx2
    !-------------------------------------------------------------------------------------------------------
    sqrSumdveightdx = 2.0*V*sqrSumdveightdx 
    
    !---------------------------------------------------------------------------
    ! Calculate d2Vdx2
    !---------------------------------------------------------------------------
    d2Vdx2 = 0.0
    do l=1,sys%dimen
       do k=1,sys%dimen
          do j=1,sys%natom
             do i=1,sys%natom
                d2Vdx2(l,k,j,i) = Sumd2veightdx2(l,k,j,i) +  SumTaydvSumdv2(l,k,j,i) + &
                        sqrSumdveightdx(l,k,j,i) + SumTayd2veightdx2(l,k,j,i) + &
                        SumWd2Taydx2(l,k,j,i) + dWTdx2(l,k,j,i) 
             enddo
          enddo
       enddo
    enddo
    !print *, d2Vdx2
    !call exit(0)

  return
end subroutine
