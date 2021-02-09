

subroutine calcen(sys,interp,pot,neigh,Weight,r,V,dVdR,RawWeightTemp,&
                drdx,d2rdx2) 
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

    ! For Sum_{l=1}^{N choose 2} drdx * r^2 d2veightdx2
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen) :: v2r4ut2d2rdx2
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen) :: v2r4ut2d2rdxmdxn

    ! For d2Taydx2
    real(kind=8), dimension(sys%dimen,sys%dimen) :: d2Taydx2_tmp2 ! intent(out)
    real(kind=8), dimension(sys%dimen,sys%dimen) :: d2Taydxmdxn_tmp2 ! intent(out)
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen) :: d2Taydx2_tmp1
    real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen) :: d2Taydxmdxn_tmp1
    real(kind=8), dimension(sys%dimen,sys%dimen) :: SumWd2Taydx2_tmp1
    real(kind=8), dimension(sys%dimen,sys%dimen) :: SumWd2Taydxmdxn_tmp1

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
!    real(kind=8), dimension(sys%dimen,sys%dimen), intent(out) :: d2Vdx2
!    real(kind=8), dimension(sys%dimen,sys%dimen), intent(out) :: d2Vdxmdxn

    !DWeightdxm
    real(kind=8), dimension(size(Weight),sys%dimen) :: DWeightdxm
    real(kind=8), dimension(size(Weight),sys%dimen) :: DWeightdxn
    !dTaydxm and dTaydxn
    real(kind=8), dimension(size(Weight),sys%dimen) :: dTaydxm
    real(kind=8), dimension(size(Weight),sys%dimen) :: dTaydxn
    ! dWTdx2
    real(kind=8), dimension(sys%dimen,sys%dimen) :: dWTdx2
    real(kind=8), dimension(sys%dimen,sys%dimen) :: dWTdxmdxn

    !Stores the value of r^2*drdx
    !real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: r2drdx
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: r2drdx_mb
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: r2drdx_nb

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
    
    !----------------------------------------------------
    !  Evaluate v2*ut*r^2*drdx i.e. 2nd term in d2Taydx2
    !----------------------------------------------------
    do j=1,sys%dimen
       do i=1,sys%nbond
          r2drdx_mb(j,sys%mb(i),i) = drdx(j,sys%mb(i),i)*r(i)**2
          r2drdx_nb(j,sys%nb(i),i) = -r2drdx_mb(j,sys%mb(i),i)
       enddo
    enddo

    ! ut*r2drdx matrix multiplication 
    utr2drdx_mb = 0.0
    !utr2drdx_nb = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%nint
             do i=1,sys%nbond
                utr2drdx_mb(k,j,l) = utr2drdx_mb(k,j,l) + pot(neigh%inner(k))%ut(j,i)*r2drdx_mb(l,sys%mb(i),i) 
                !utr2drdx_nb(k,j,l) = utr2drdx_nb(k,j,l) + pot(neigh%inner(k))%ut(j,i)*r2drdx_nb(l,sys%nb(i),i) 
             enddo
          enddo
       enddo
    enddo
    utr2drdx_nb = -utr2drdx_mb

    ! v2* (ut*r^2*drdx) * ut
    v2detadxut_mb = 0.0
    !v2detadxut_nb = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do l=1,sys%nbond
             do i=1,sys%nint
                v2detadxut_mb(k,l,j) = v2detadxut_mb(k,l,j) + pot(neigh%inner(k))%v2(i)*&
                        utr2drdx_mb(k,i,j)*pot(neigh%inner(k))%ut(i,l)
                !v2detadxut_nb(k,l,j) = v2detadxut_nb(k,l,j) + pot(neigh%inner(k))%v2(i)*&
                !        utr2drdx_nb(k,i,j)*pot(neigh%inner(k))%ut(i,l)
             enddo
             v2detadxut_mb(k,l,j) = v2detadxut_mb(k,l,j)*r(l)**2
             !v2detadxut_nb(k,l,j) = v2detadxut_nb(k,l,j)*r(l)**2
          enddo
       enddo
    enddo
    v2detadxut_nb = -v2detadxut_mb

    ! Sum_{l=1}^{N choose 2} drdx * r^-2 * v2detadxut
    v2r4ut2d2rdx2 = 0.0
    !v2r4ut2d2rdxmdxn = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do i=1,sys%nbond
                v2r4ut2d2rdx2(k,l,j) = v2r4ut2d2rdx2(k,l,j) + v2detadxut_mb(k,i,j)*r(i)**2*drdx(l,sys%mb(i),i)
                !v2r4ut2d2rdxmdxn(k,l,j) = v2r4ut2d2rdxmdxn(k,l,j) + v2detadxut_nb(k,i,j)*r(i)**2*drdx(l,sys%mb(i),i)
             enddo
          enddo
       enddo
    enddo
    v2r4ut2d2rdxmdxn = -v2r4ut2d2rdx2

    ! Sum_{j=1}^{Ndata} w_j * v2r4ut2d2rdx2
    d2Taydx2_tmp2 = 0.0
    !d2Taydxmdxn_tmp2 = 0.0
    do k=1,sys%dimen
       do j=1,sys%dimen
          do i=1,neigh%numInner
             d2Taydx2_tmp2(k,j) = d2Taydx2_tmp2(k,j) + Weight(i)*v2r4ut2d2rdx2(i,k,j)
             !d2TaydR2tmp2_nb(k,j) = d2TaydR2tmp2_nb(k,j) + Weight(i)*v2r4ut2d2rdxmdxn(i,k,j)
          enddo
       enddo
    enddo
    d2Taydxmdxn_tmp2 = -d2Taydx2_tmp2


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
    ! Calculate the 1st and 3rd term of d2Taydx2
    !---------------------------------------------------------
    ! Sum_{l=1}^{N choose 2} d2Tay_j/dr_l * dr_l/dx
    d2Taydx2_tmp1 = 0.0
    !d2Taydxmdxn_tmp1 = 0.0
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do j=1,sys%dimen
             do i=1,sys%nbond
!                d2Taydx2_tmp1(k,l,j) = d2Taydx2_tmp1(k,l,j) + (2.0*r(i)-1.0)*dTaydR(k,i)*d2rdx2(l,j,i)
                !d2Taydxmdxn_tmp1(k,l,j) = d2Taydxmdxn_tmp1(k,l,j) + (2.0*r(i)-1.0)*dTaydR(k,i)*d2rdxmdxn(l,j,i)
             enddo
          enddo
       enddo
    enddo
    d2Taydxmdxn_tmp1 = -d2Taydx2_tmp1

    ! Sum_{j=1}^{Ndata} Weight_j * dTay_j/dx
    SumWd2Taydx2_tmp1 = 0.0
    !SumWd2Taydxmdxn_tmp1 = 0.0
    do k=1,sys%dimen
       do j=1,sys%dimen
          do i=1,neigh%numInner
             SumWd2Taydx2_tmp1(k,j) = SumWd2Taydx2_tmp1(k,j) + d2Taydx2_tmp1(i,k,j)*Weight(i)
             !SumWd2Taydxmdxn_tmp1(k,j) = SumWd2Taydxmdxn_tmp1(k,j) + d2Taydxmdxn_tmp1(i,k,j)*Weight(i)
          enddo
       enddo
    enddo
    SumWd2Taydxmdxn_tmp1 = -SumWd2Taydx2_tmp1

    !----------------------------------------------------------
    !  Calculate sum_{j}^{Ndata} 2.0 * dWeightdx * dTaydx
    !----------------------------------------------------------
    ! dWeight/dx_m and dWeight/dx_n
    DWeightdxm = 0.0
    !DWeightdxn = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
            DWeightdxm(k,j) = DWeightdxm(k,j) + DWeight(i,k)*drdx(j,sys%mb(i),i)
            !DWeightdxn(k,j) = DWeightdxn(k,j) + DWeight(i,k)*drdx(j,sys%nb(i),i)
          enddo
       enddo
    enddo
    DWeightdxn = -DWeightdxm

    ! dTaydxm dTaydxn
    dTaydxm = 0.0
    do k=1,neigh%numInner
       do j=1,sys%dimen
          do i=1,sys%nbond
             dTaydxm(k,j) = dTaydxm(k,j) + dTaydR(k,i)*drdx(j,sys%mb(i),i)
          enddo
       enddo
    enddo
    dTaydxn = -dTaydxm

    ! DWeightdx * dTaydx
    do k=1,neigh%numInner
       do l=1,sys%dimen
          do i=1,sys%dimen
             dWTdx2(l,i) = dWTdx2(l,i) + DWeightdxm(k,l)*dTaydxm(k,i)
             dWTdxmdxn(l,i) = dWTdxmdxn(l,i) + DWeightdxm(k,l)*dTaydxn(k,i) 
          enddo
       enddo
    enddo
    dWTdx2 = 2.0*dWTdx2
    dWTdxmdxn = 2.0*dWTdxmdxn

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
    !  Calculate (Sum_{j=1}^{Ndata} 2.0 * w_j * T_j) * (Sum dvdx)*(Sum dvdx)
    !---------------------------------------------------------------------------
!    d2Weightdx2 = SumTayd2veightdx2 + Sumd2veightdx2 + 2.0*SumTaydvmSumdvm + sqrSumdveightdx
!    d2Weightdxmdxn = SumTayd2veightdxmdxn + Sumd2veightdxmdxn + 2.0*SumTaydvmSumdvn + sqrSumdveightdxmdxn 

    !---------------------------------------------------------------------------
    ! Calculate d2Vdx2
    !---------------------------------------------------------------------------
!    do j=1,sys%dimen
!       do i=1,sys%dimen
!          d2Vdx2(j,i) = d2Weightdx2(j,i) + sqrSumdveightdx(j,i) + SumTaydvmSumdvm(j,i) + Sumd2veightdx2(j,i) +&
!                SumTayd2veightdx2(j,i) + SumWd2Taydx2_tmp1(j,i) + d2Taydx2_tmp2(j,i) 
!          !d2Vdxmdxn(j,i) = d2Weightdxmdxn(j,i) + sqrSumdveightdxmdxn(j,i) + SumTaydvmSumdvn(j,i) + &
!          !     Sumd2veightdxmdxn(j,i) +  SumTayd2veightdxmdxn(j,i) + SumWd2Taydxmdxn_tmp1(j,i) + &
!          !     + d2Taydxmdxn_tmp2(j,i) 
!       enddo
!    enddo
!    d2Vdxmdxn = -d2Vdx2

  return
end subroutine
