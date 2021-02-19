

  subroutine cal_d2Vdx2(sys,interp,pot,neigh,r,drdx,d2rdx2,dr2dxdx,Weight,RawWeightTemp&
  ,DWeight,Tay,V,dTaydR,d2Vdx2)
     use molecule_specs
     use interpolation
     
     implicit none

     type (molsysdat), intent(in) :: sys
     type (interp_params), intent(in) :: interp
     type (neighbour_list), intent(in) :: neigh
     type (pot_data_point), dimension(:), pointer :: pot

     real(kind=8), dimension(interp%ndata), intent(in) :: Tay
     real(kind=8), dimension(:), intent(in) :: Weight
     real(kind=8), dimension(size(Weight),sys%nbond), intent(in) :: dTaydR
     real(kind=8), intent(in) :: V
     real(kind=8), dimension(sys%nbond,size(Weight)), intent(in) :: DWeight
     real(kind=8), dimension(:), intent(in) :: RawWeightTemp
     real(kind=8), dimension(:), intent(in) :: r ! inverse bondlength
     real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(in) :: drdx
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom,sys%nbond), intent(in):: d2rdx2
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom,sys%nbond), intent(in):: dr2dxdx

     ! Temporate terms for d2v_jdx2
     real(kind=8), dimension(sys%nbond,size(Weight)) :: D2veightdr_tmp1
     real(kind=8), dimension(sys%nbond,size(Weight)) :: D2veightdr_tmp2
     real(kind=8), dimension(sys%nbond,size(Weight)) :: D2veightdr_tmp3
     real(kind=8), dimension(sys%nbond,size(Weight)) :: D2veightdr_tmp4
     real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: dveightdx_tmp1
     real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: zdzdx_tmp1
     real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2veightdx2_tmp1
     real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2veightdx2_tmp2
     real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2veightdx2_tmp3
     real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2veightdx2
     real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: dveightdx
     real(kind=8), dimension(sys%dimen,sys%natom) :: Sumdveightdx
     real(kind=8), dimension(sys%dimen,sys%natom) :: SumTaydveightdx
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: Sumd2veightdx2 
     real(kind=8) :: SumWeight

     ! Temp terms for Tay*d2Weightdx2
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Weightdx2_tmp1
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Weightdx2_tmp2
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Weightdx2_tmp3
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Weightdx2_tmp4
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Weightdx2

     ! Temp terms for dTaydx * dWeightdx and dWeightdx * dTaydx
     real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: DWeightdx
     real(kind=8), dimension(size(Weight),sys%dimen,sys%natom) :: DTaydx 
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: SumdWeightdxdTaydx
     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: SumdTaydxdWeightdx

     ! Temp terms for d2Taydx2
     real(kind=8), dimension(size(Weight),sys%nint,sys%dimen,sys%natom) :: SumUTr2drdx
     real(kind=8), dimension(size(Weight),sys%nbond,sys%dimen,sys%natom) :: Sumv2utdetadx

     real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond) :: drdx_tmp
     real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Taydx2_tmp1
     real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Taydx2_tmp2
     real(kind=8), dimension(size(Weight),sys%dimen,sys%dimen,sys%natom,sys%natom) :: d2Taydx2_tmp3

     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom) :: Weightd2Taydx2

     real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom), intent(out) :: d2Vdx2


     integer :: i,j,k,l,m,n
     real(kind=8) :: temp,temp2,temp3

     !-------------------------------------------------
     ! Construct d2v_jdx2 * (sum_{j=1}^{Ndata} v_j)
     !-------------------------------------------------
     do k=1,neigh%numInner
        temp = interp%ipow2*Weight(k)*RawWeightTemp(k)
        temp2 = (interp%ipow2+2.0)*temp*RawWeightTemp(k)
        temp3 = -2.0*temp
        do i=1,sys%nbond
           D2veightdr_tmp1(i,k) = temp2*(r(i)-pot(neigh%inner(k))%r(i))*r(i)**2
           D2veightdr_tmp2(i,k) = -temp*r(i)**4
           D2veightdr_tmp3(i,k) = temp3*(r(i)-pot(neigh%inner(k))%r(i))*r(i)**3
           D2veightdr_tmp4(i,k) = temp*(r(i)-pot(neigh%inner(k))%r(i))*r(i)**2
        enddo
     enddo

     dveightdx_tmp1 = 0.0
     zdzdx_tmp1 = 0.0
     do k=1,neigh%numInner
        do j=1,sys%dimen
           do i=1,sys%nbond
              dveightdx_tmp1(k,j,sys%mb(i)) = dveightdx_tmp1(k,j,sys%mb(i)) + &
                   D2veightdr_tmp1(i,k)*drdx(j,sys%mb(i),i)
              dveightdx_tmp1(k,j,sys%nb(i)) = dveightdx_tmp1(k,j,sys%nb(i)) + &
                   D2veightdr_tmp1(i,k)*drdx(j,sys%nb(i),i)
              zdzdx_tmp1(k,j,sys%mb(i)) = zdzdx_tmp1(k,j,sys%mb(i)) + &
                   (r(i)-pot(neigh%inner(k))%r(i))*drdx(j,sys%mb(i),i)*r(i)**2
              zdzdx_tmp1(k,j,sys%nb(i)) = zdzdx_tmp1(k,j,sys%nb(i)) + &
                   (r(i)-pot(neigh%inner(k))%r(i))*drdx(j,sys%nb(i),i)*r(i)**2
           enddo
        enddo
     enddo

!     zdzdx_tmp1 = 0.0
!     do k=1,neigh%numInner
!        do j=1,sys%dimen
!           do i=1,sys%nbond
!              zdzdx_tmp1(k,j,sys%mb(i)) = zdzdx_tmp1(k,j,sys%mb(i)) + &
!                   (r(i)-pot(neigh%inner(k))%r(i))*drdx(j,sys%mb(i),i)*r(i)**2
!              zdzdx_tmp1(k,j,sys%nb(i)) = zdzdx_tmp1(k,j,sys%nb(i)) + &
!                   (r(i)-pot(neigh%inner(k))%r(i))*drdx(j,sys%nb(i),i)*r(i)**2
!           enddo
!        enddo
!     enddo

!     do k=1,neigh%numInner
!        do m=1,sys%dimen
!           do l=1,sys%dimen
!              do j=1,sys%natom
!                 do i=1,sys%natom
!                    d2veightdx2_tmp1(k,m,l,j,i) = zdzdx_tmp1(k,l,i)*dveightdx_tmp1(k,m,j)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
         !do m=1,sys%dimen
         !   do l=1,sys%dimen
         !      !print *, zdzdx_tmp1(1,m,1)*dveightdx_tmp1(1,l,5)
         !      print *, d2veightdx2_tmp1(1,m,l,1,1)
         !      !print *, ''
         !   enddo
         !enddo
         !call exit(0)

     d2veightdx2_tmp2 = 0.0
     d2veightdx2_tmp3 = 0.0
     do k=1,neigh%numInner
        do m=1,sys%dimen
           do n=1,sys%dimen
              do l=1,sys%natom
                 do j=1,sys%natom
                    d2veightdx2_tmp1(k,m,n,l,j) = zdzdx_tmp1(k,n,j)*dveightdx_tmp1(k,m,l)
                    do i=1,sys%nbond
                       d2veightdx2_tmp2(k,m,n,l,j) = d2veightdx2_tmp2(k,m,n,l,j) + &
                           dr2dxdx(m,n,l,j,i)*(D2veightdr_tmp2(i,k)+D2veightdr_tmp3(i,k))
                       d2veightdx2_tmp3(k,m,n,l,j) = d2veightdx2_tmp3(k,m,n,l,j) + &
                           d2rdx2(m,n,l,j,i)*D2veightdr_tmp4(i,k)
                    enddo
                    ! For single threading this is slower somehow
                    !d2veightdx2(k,m,n,l,j) = d2veightdx2_tmp1(k,m,n,l,j) + &
                    !     d2veightdx2_tmp2(k,m,n,l,j) + d2veightdx2_tmp3(k,m,n,l,j)
                 enddo
              enddo
           enddo
        enddo
     enddo

     do k=1,neigh%numInner
        do m=1,sys%dimen
           do l=1,sys%dimen
              do j=1,sys%natom
                 do i=1,sys%natom
                    d2veightdx2(k,m,l,j,i) = d2veightdx2_tmp1(k,m,l,j,i) + &
                         d2veightdx2_tmp2(k,m,l,j,i) + d2veightdx2_tmp3(k,m,l,j,i)
                 enddo
              enddo
           enddo
       enddo
     enddo
         !do j=1,sys%dimen
         !   do i=1,sys%dimen
         !      print *, d2veightdx2(1,j,i,4,3)
         !   enddo
         !enddo
         !call exit(0)

     !---------------------------------------------------------
     ! Evaluate Sum_{j=1}^{Ndata} Tay_j * d2veight_jdx2 
     !---------------------------------------------------------
     d2Weightdx2_tmp1 = 0.0
     do k=1,neigh%numInner
        do m=1,sys%dimen
           do l=1,sys%dimen
              do j=1,sys%natom
                 do i=1,sys%natom
                    d2Weightdx2_tmp1(m,l,j,i) = d2Weightdx2_tmp1(m,l,j,i) + &
                          Tay(k)*d2veightdx2(k,m,l,j,i)
                 enddo
              enddo
           enddo
        enddo
     enddo
          !print *, d2Weightdx2_tmp1(:,:,1,1)

     !---------------------------------------------------
     ! Evaluate Sum_{j=1}^{Ndata} d2v_jdx2
     !---------------------------------------------------
     Sumd2veightdx2 = 0.0
     do m=1,sys%dimen
        do l=1,sys%dimen
           do k=1,sys%natom
              do j=1,sys%natom
                 do i=1,neigh%numInner
                    Sumd2veightdx2(m,l,k,j) = Sumd2veightdx2(m,l,k,j) + d2veightdx2(i,m,l,k,j)
                 enddo
              enddo
           enddo
        enddo
     enddo
         !do j=1,sys%dimen
         !   do i=1,sys%dimen
         !      print *, Sumd2veightdx2(j,i,1,1)
         !   enddo
         !enddo
         !call exit(0)
         !print *, Sumd2veightdx2(:,:,1,1)
         
     !-------------------------------------------------
     ! Evaluate dv_jdx*( Sum_j^Ndata v_j)^-1
     !-------------------------------------------------
     dveightdx = 0.0
     do k=1,neigh%numInner
        do j=1,sys%dimen
           do i=1,sys%nbond
              dveightdx(k,j,sys%mb(i)) = dveightdx(k,j,sys%mb(i)) + &
                     DWeight(i,k)*drdx(j,sys%mb(i),i)
              dveightdx(k,j,sys%nb(i)) = dveightdx(k,j,sys%nb(i)) + &
                     DWeight(i,k)*drdx(j,sys%nb(i),i)
           enddo
        enddo
     enddo

     !----------------------------------------------------------
     ! Evaluate Sum_{j=1}^{Ndata} dv_jdx *( Sum_j^Ndata v_j)^-1
     !----------------------------------------------------------
     Sumdveightdx = 0.0
     do k=1,neigh%numInner
        do j=1,sys%dimen
           do i=1,sys%natom
              Sumdveightdx(j,i) =  Sumdveightdx(j,i) + dveightdx(k,j,i) 
           enddo
        enddo
     enddo

     !---------------------------------------------------------------
     ! Evaluate (Sum_{j=1}^{Ndata} Tay_j*dveightdx) * Sumdveightdx
     !---------------------------------------------------------------
     SumTaydveightdx = 0.0
     do k=1,neigh%numInner
        do j=1,sys%dimen
           do i=1,sys%nbond
             SumTaydveightdx(j,sys%mb(i)) = SumTaydveightdx(j,sys%mb(i)) - &
                 dveightdx(k,j,sys%mb(i))*Tay(k)
             SumTaydveightdx(j,sys%nb(i)) = SumTaydveightdx(j,sys%nb(i)) - &
                 dveightdx(k,j,sys%nb(i))*Tay(k)
           enddo
        enddo
     enddo

     !--------------------------------------------------------------
     ! Evaluate SumTaydveightdx * Sumdveightdx 
     !                         &
     ! Evaluate (Sum_{j=1}^{Ndata} 2.0*Tay_j*Weight_j) Sumdveightdx  
     !                         &
     ! Evaluate (Sum_{j=1}^{Ndata} Weight_j) * Sumd2veightdx2  
     !--------------------------------------------------------------
     do l=1,sys%dimen
        do k=1,sys%dimen
           do j=1,sys%natom
              do i=1,sys%natom
                 d2Weightdx2_tmp2(l,k,j,i) = SumTaydveightdx(l,j)*Sumdveightdx(k,i) &
                    + SumTaydveightdx(k,i)*Sumdveightdx(l,j)
                 d2Weightdx2_tmp3(l,k,j,i) = 2.0*V*Sumdveightdx(l,j)*Sumdveightdx(k,i)
                 d2Weightdx2_tmp4(l,k,j,i) = -V*Sumd2veightdx2(l,k,j,i)
              enddo
           enddo
        enddo
     enddo
            !do l=1,sys%dimen
            !   do k=1,sys%dimen
            !      print *, 2.0*SumTaydveightdx(l,5)*Sumdveightdx(k,4)
            !   enddo
            !enddo
            !do i=1, sys%dimen
            !   print *, d2Weightdx2_tmp2(:,i,5,4)
            !enddo
            !print *, d2Weightdx2_tmp2(2,1,5,4)
            !call exit(0)

     !------------------------------------------------------------
     ! Evaluate (Sum_{j=1}^{Ndata} 2.0*Tay_j*Weight_j) Sumdveightdx  
     !------------------------------------------------------------
!     do l=1,sys%dimen
!        do k=1,sys%dimen
!           do j=1,sys%natom
!              do i=1,sys%natom
!                 d2Weightdx2_tmp3(l,k,j,i) = 2.0*V*Sumdveightdx(l,j)*Sumdveightdx(k,i)
!              enddo
!           enddo
!        enddo
!     enddo
            !print *, d2Weightdx2_tmp3(:,:,1,1)

     !-----------------------------------------------------------
     ! Evaluate (Sum_{j=1}^{Ndata} Weight_j) * Sumd2veightdx2  
     !-----------------------------------------------------------
!     do l=1,sys%dimen
!        do k=1,sys%dimen
!           do j=1,sys%natom
!              do i=1,sys%natom
!                 d2Weightdx2_tmp4(l,k,j,i) = -V*Sumd2veightdx2(l,k,j,i)
!              enddo
!           enddo
!        enddo
!     enddo
          !print *, d2Weightdx2_tmp4(:,:,1,1)
          !print *, ''

     !----------------------------------------------------------
     ! Evaluate Sum_{j=1}^{Ndata} d2Weight_jdx2
     !----------------------------------------------------------
     do l=1,sys%dimen
        do k=1,sys%dimen
           do j=1,sys%natom
              do i=1,sys%natom
                 d2Weightdx2(l,k,j,i) = d2Weightdx2_tmp1(l,k,j,i) + d2Weightdx2_tmp2(l,k,j,i) &
                       + d2Weightdx2_tmp3(l,k,j,i) + d2Weightdx2_tmp4(l,k,j,i)
              enddo
           enddo
        enddo
     enddo
          !print *, d2Weightdx2(:,:,1,1)
          !print *, ''
          !call exit(0)

     !----------------------------------------------------------
     ! Evaluate Sum_{j=1}^{Ndata} dWeight_jdx * dTay_jdx
     !----------------------------------------------------------
     DWeightdx = 0.0
     DTaydx = 0.0
     do k=1,neigh%numInner
        do j=1,sys%dimen
           do i=1,sys%nbond
              DWeightdx(k,j,sys%mb(i)) = DWeightdx(k,j,sys%mb(i)) + DWeight(i,k)*drdx(j,sys%mb(i),i)
              DWeightdx(k,j,sys%nb(i)) = DWeightdx(k,j,sys%nb(i)) + DWeight(i,k)*drdx(j,sys%nb(i),i)
              DTaydx(k,j,sys%mb(i)) = DTaydx(k,j,sys%mb(i)) + DTaydR(k,i)*drdx(j,sys%mb(i),i)
              DTaydx(k,j,sys%nb(i)) = DTaydx(k,j,sys%nb(i)) + DTaydR(k,i)*drdx(j,sys%nb(i),i)
           enddo
        enddo
     enddo

     SumdWeightdxdTaydx = 0.0
     SumdTaydxdWeightdx = 0.0
     do m=1,sys%dimen
        do l=1,sys%dimen
           do k=1,sys%natom
              do j=1,sys%natom
                 do i=1,neigh%numInner
                    SumdWeightdxdTaydx(m,l,k,j) = SumdWeightdxdTaydx(m,l,k,j) + &
                        DWeightdx(i,m,k)*DTaydx(i,l,j)
                    SumdTaydxdWeightdx(m,l,k,j) = SumdTaydxdWeightdx(m,l,k,j) + &
                        DTaydx(i,m,k)*DWeightdx(i,l,j)
                 enddo
              enddo
           enddo
        enddo
     enddo
            !print *, SumdWeightdxdTaydx(:,:,1,1)
            !print *, SumdTaydxdWeightdx(:,:,1,1)
            !call exit(0)
     
     !=========================================================
     ! Below session is for d2Taydx2 
     !=========================================================
     !-----------------------------------------------
     ! Sum_{l=1}^{N choose 2} detadz * r^{-2} * drdx
     !-----------------------------------------------
     SumUTr2drdx = 0.0
     do k=1,neigh%numInner
        do j=1,sys%nint
           do l=1,sys%dimen
              do i=1,sys%nbond
                 SumUTr2drdx(k,j,l,sys%mb(i)) = SumUTr2drdx(k,j,l,sys%mb(i)) + &
                     pot(neigh%inner(k))%ut(j,i)*r(i)**2*drdx(l,sys%mb(i),i)
                 SumUTr2drdx(k,j,l,sys%nb(i)) = SumUTr2drdx(k,j,l,sys%nb(i)) + &
                     pot(neigh%inner(k))%ut(j,i)*r(i)**2*drdx(l,sys%nb(i),i)
              enddo
           enddo
        enddo
     enddo

     !------------------------------------------------
     ! Sum_{i=1}^{3N-6} v2_i * ut_{il} * SumUTr2drdx
     !------------------------------------------------
     Sumv2utdetadx = 0.0
     do k=1,neigh%numInner
        do m=1,sys%nbond
           do j=1,sys%nint
              do l=1,sys%dimen
                 do i=1,sys%natom
                    Sumv2utdetadx(k,m,l,i) = Sumv2utdetadx(k,m,l,i) + & 
                        SumUTr2drdx(k,j,l,i)*pot(neigh%inner(k))%ut(j,m) &
                        * pot(neigh%inner(k))%v2(j)
                 enddo
              enddo
           enddo
        enddo
     enddo

     !----------------------------------------------------
     ! Sum_{l=1}^{N choose 2} r_l * drdx * Sumv2utdetadx
     !----------------------------------------------------
     d2Taydx2_tmp1 = 0.0
     do k=1,neigh%numInner
        do m=1,sys%dimen
           do l=1,sys%dimen
              do n=1,sys%natom
                 do i=1,sys%nbond
                    d2Taydx2_tmp1(k,m,l,n,sys%mb(i)) = d2Taydx2_tmp1(k,m,l,n,sys%mb(i)) & 
                      + Sumv2utdetadx(k,i,m,n)*drdx(l,sys%mb(i),i)*r(i)**2
                    d2Taydx2_tmp1(k,m,l,n,sys%nb(i)) = d2Taydx2_tmp1(k,m,l,n,sys%nb(i)) & 
                      + Sumv2utdetadx(k,i,m,n)*drdx(l,sys%nb(i),i)*r(i)**2
                 enddo
              enddo
           enddo
        enddo
     enddo
           !do j=1,sys%dimen
           !   do i=1,sys%dimen
           !      print *, d2Taydx2_tmp1(1,j,i,2,1), d2Taydx2_tmp1(1,j,i,1,2)
           !   enddo
           !enddo
           !call exit(0)

     !--------------------------------------------------
     ! Second term of d2Taydx2
     !          &
     ! Third term of d2Taydx2
     !--------------------------------------------------
     d2Taydx2_tmp2 = 0.0
     d2Taydx2_tmp3 = 0.0
     do k=1,neigh%numInner
        do m=1,sys%dimen
           do n=1,sys%dimen
              do l=1,sys%natom
                 do j=1,sys%natom
                    do i=1,sys%nbond
                       d2Taydx2_tmp2(k,m,n,l,j) = d2Taydx2_tmp2(k,m,n,l,j) &
                          - 2.0*dr2dxdx(m,n,l,j,i)*dTaydR(k,i)*r(i)
                       d2Taydx2_tmp3(k,m,n,l,j) = d2Taydx2_tmp3(k,m,n,l,j) &
                          + d2rdx2(m,n,l,j,i)*dTaydR(k,i)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
         !do j=1,sys%dimen
         !   do i=1,sys%dimen
         !      print *, d2Taydx2_tmp1(1,j,i,1,1), d2Taydx2_tmp2(1,j,i,1,1)
         !   enddo
         !enddo
         !call exit(0)
     
     !-------------------------------------------------
     ! Third term of d2Taydx2
     !-------------------------------------------------
!     d2Taydx2_tmp3 = 0.0
!     do k=1,neigh%numInner
!        do m=1,sys%dimen
!           do n=1,sys%dimen
!              do l=1,sys%natom
!                 do j=1,sys%natom
!                    do i=1,sys%nbond
!                       d2Taydx2_tmp3(k,m,n,l,j) = d2Taydx2_tmp3(k,m,n,l,j) &
!                          + d2rdx2(m,n,l,j,i)*dTaydR(k,i)
!                    enddo
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo

          !do j=1,sys%dimen
          !   do i=1,sys%dimen
          !      print *, d2Taydx2_tmp1(1,j,i,1,1), d2Taydx2_tmp2(1,j,i,1,1), d2Taydx2_tmp3(1,j,i,1,1)
          !   enddo
          !enddo
          !do j=1,sys%dimen
          !   do i=1,sys%dimen
          !      print *, d2Taydx2_tmp3(1,j,i,1,2)
          !   enddo
          !enddo
          !call exit(0)

     !---------------------------------------------
     ! Evaluate Sum_{j=1}^{Ndata} d2Tay_jdx2 
     !---------------------------------------------
     Weightd2Taydx2 = 0.0
     do k=1,neigh%numInner
        do m=1,sys%dimen
           do n=1,sys%dimen
              do l=1,sys%natom
                 do i=1,sys%natom
                    Weightd2Taydx2(m,n,l,i) = Weightd2Taydx2(m,n,l,i) + &
                       (d2Taydx2_tmp3(k,m,n,l,i) + d2Taydx2_tmp1(k,m,n,l,i) + &
                       d2Taydx2_tmp2(k,m,n,l,i))*Weight(k)
                 enddo
              enddo
           enddo
        enddo
     enddo
            !print *, Weightd2Taydx2(:,:,1,1)
            !call exit(0)

     !-------------------------------------------
     ! Evaluate d2Vdx2
     !-------------------------------------------
     do l=1,sys%dimen
        do k=1,sys%dimen
           do j=1,sys%natom
              do i=1,sys%natom
                 d2Vdx2(l,k,j,i) = Weightd2Taydx2(l,k,j,i) + SumdWeightdxdTaydx(l,k,j,i) & 
                 !d2Vdx2(l,k,j,i) = SumdWeightdxdTaydx(l,k,j,i) & 
                    + SumdTaydxdWeightdx(l,k,j,i) + d2Weightdx2(l,k,j,i)
              enddo
           enddo
        enddo
     enddo

            !print *, d2Vdx2(:,:,1,1)
            !print *, d2Vdx2(:,:,5,1)
            !print *, ''
            !call exit(0)

     return
  end subroutine
