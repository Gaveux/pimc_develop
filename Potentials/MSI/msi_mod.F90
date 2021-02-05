      !> @file
      !! The potential_msi module is used to calculate the value of a Modified Shepard Interpolated 
      !> Potential Energy Surface
      !
      ! DESCRIPTION:
      !!@details 
      
      module potential_msi

        use interpolation
        use molecule_specs

      
        implicit none


        type msi_params
            type(interp_params) :: interp
            type(pot_data_point), dimension(:), pointer :: pot
            type (molsysdat) :: sys
            type(neighbour_list), dimension(:), pointer :: neighlist
        end type msi_params
        

        contains

        !Evaluate the potential energy of the system and optionally the cartesian first and second derivatives
        !of the potential.  
        subroutine potential(ind,param,x,r,V,dV)
            integer, intent(in) :: ind
            !declaration of the variables passed to the subroutine
            type (msi_params) :: param
            ! Cartesian positions of atoms in system
            real(kind=8), dimension(param%sys%dimen,param%sys%natom), intent(in) :: x
            !bond lengths
            real(kind=8), dimension(param%sys%nbond), intent(inout) :: r
            !Calculated potential energy
            real(kind=8), intent(out) :: V          
            !Calculated derivatives of the potential energy with respect to cartesian coordinates - optional
            real(kind=8), dimension(param%sys%dimen,param%sys%natom), intent(out) :: dV 
            real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%natom) :: d2Taydx2 ! intent(out)
            real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%natom) :: d2Taydxmdxn ! intent(out)
            real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%natom) :: d2wdx2
            real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%natom) :: d2wdxmdxn
            real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%natom) :: dWTdx2
            real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%natom) :: dWTdxmdxn
            !real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%natom) :: d2Vdx2
            !real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%natom) :: d2Vdxmdxn
            

            !Variables used internally by the modified shepard code
              
            !derivative of bond lengths with respect to cartesians
            real(kind=8), dimension(param%sys%dimen,param%sys%natom&
            ,param%sys%nbond) :: dr
            !second derivative matrix of bondlengths w.r.t. Cartesians
            real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%nbond) :: d2rdx2
            real(kind=8), dimension(param%sys%dimen,param%sys%dimen,param%sys%nbond) :: d2rdxmdxn
            !Derivatives of the potential with respect to internal coordinates
            real(kind=8), dimension(param%sys%nbond) :: dVdr
            !dWeightdr*dTaydr
            real(kind=8), dimension(param%sys%nbond) :: dWTdr2
            !Stores the value of the weight function
            real(kind=8), dimension(param%interp%ndata) :: Weight
            real(kind=8), dimension(param%interp%ndata) :: RawWeightTemp

            !Variables used for d2Vdr2
            real(kind=8), dimension(param%sys%nbond) :: Tayd2veightdr2tmp1
            real(kind=8), dimension(param%sys%nbond) :: Tayd2veightdr2tmp2

            ! (Sum Tay*DWeight)*SumDWeight
            real(kind=8), dimension(param%sys%nbond) :: TDWeightSumDWeight

            ! Weight * Tay * (Sum dvdr)**2
            real(kind=8), dimension(param%sys%nbond) :: WTSumDWeightDrSqr

            real(kind=8), dimension(param%sys%nbond) :: d2TaydR2tmp1
            real(kind=8), dimension(param%sys%dimen,param%sys%nbond) :: d2TaydR2tmp2_mb
            real(kind=8), dimension(param%sys%dimen,param%sys%nbond) :: d2TaydR2tmp2_nb

            ! second term in d2Taydx2, r^2 * (ut*r^2*drdx) * ut
            real(kind=8), dimension(size(Weight),param%sys%nbond,param%sys%dimen) :: v2detadxut_mb
            real(kind=8), dimension(size(Weight),param%sys%nbond,param%sys%dimen) :: v2detadxut_nb            
            
            integer :: i,j,k
            include 'intern.int'
            include 'neigh.int'
            include 'calcen.int'

            call intern(param%sys,x,r,dr,d2rdx2,d2rdxmdxn)

            !Update the inner neighbour list each potential evaluation
            call neighbour(param%sys,param%interp,param%pot,Weight,r,&
            param%neighlist(ind),RawWeightTemp)

            !Interpolate the surface - currently only the one part weight function is implemented
            !the two part weight function can easily be included, however there will need to be 
            !a slight shuffling of the variables in the calcen2w.f90 file to make them compatible
            !with the rest of the code

            !if (param%interp%ipart == 1) then
            call calcen(param%sys,param%interp,param%pot,&
            param%neighlist(ind),Weight,r,V,dVdr,RawWeightTemp,&
            dWTdr2,Tayd2veightdr2tmp1,Tayd2veightdr2tmp2,TDWeightSumDWeight,&
            WTSumDWeightDrSqr,&
            d2TaydR2tmp1,dr,d2rdx2,d2rdxmdxn,d2TaydR2tmp2_mb,d2TaydR2tmp2_nb)
            !endif

            !print *, dVdr
            !print *, dWTdr2
            !call exit(0)

            V = V - param%interp%vmin
            do j=1,param%sys%natom
                do k=1,param%sys%dimen
                    dV(k,j) = 0.d0
                enddo
            enddo

            do j=1,param%sys%nbond
                do k=1,param%sys%dimen
                    dV(k,param%sys%mb(j))=dV(k,param%sys%mb(j))+dVdR(j)*dr(k,param%sys%mb(j),j)
                    dV(k,param%sys%nb(j))=dV(k,param%sys%nb(j))+dVdR(j)*dr(k,param%sys%nb(j),j)
                enddo
            enddo

            !do i=1, param%sys%nbond
            !   print *, dV(:,param%sys%mb(i))
            !   print *, dV(:,param%sys%nb(i))
            !   print *, ''
            !enddo
            !call exit(0)

            !==========================================
            ! Evaluate Sum weight*d2Taydx2
            !==========================================

            ! note: d2r(:,sys%mb) is diagonal elements, d2r(:,sys%mb)
            ! is off-diagonal elements
!            do j=1,param%sys%nbond
!               do k=1,param%sys%dimen
!                  do i=1,param%sys%dimen
!                     d2Taydx2(k,i,param%sys%mb(j)) = d2Taydx2(k,i,param%sys%mb(j)) &
!                          + d2TaydR2tmp1(j)*d2r(k,i,param%sys%mb(j),j)
!                     d2Taydxmdxn(k,i,param%sys%nb(j)) = d2Taydxmdxn(k,i,param%sys%nb(j)) &
!                          + d2TaydR2tmp1(j)*d2r(k,i,param%sys%nb(j),j)
!                  enddo
!               enddo
!            enddo

            !==========================================
            ! Evaluate Sum Tay*d2Weightdx2
            !==========================================
!            d2wdx2 = 0.0
!            d2wdxmdxn = 0.0
!            do k=1, param%sys%dimen
!               do j=1,param%sys%dimen
!                  do i=1,param%sys%nbond
!                     d2wdx2(k,j,param%sys%mb(j)) = d2wdx2(k,j,param%sys%mb(j)) & 
!                        + (Tayd2veightdr2tmp1(i) + 2.0*TDWeightSumDWeight(i) &
!                        + WTSumDWeightDrSqr(i) + WTaySumD2veightDr1(i))*dr(j,param%sys%mb(i),i)&
!                        * dr(k,param%sys%mb(i),i) + (Tayd2veightdr2tmp2(i) &
!                        + WTaySumD2veightDr2(i))*d2r(k,j,param%sys%mb(i),i)
!
!
!                     d2wdxmdxn(k,j,param%sys%nb(j)) = d2wdxmdxn(k,j,param%sys%nb(j)) & 
!                        + (Tayd2veightdr2tmp1(i) + 2.0*TDWeightSumDWeight(i) &
!                        + WTSumDWeightDrSqr(i) + WTaySumD2veightDr1(i))*dr(j,param%sys%mb(i),i)&
!                        * dr(k,param%sys%nb(i),i) + (Tayd2veightdr2tmp2(i) &
!                        + WTaySumD2veightDr2(i))*d2r(k,j,param%sys%nb(i),i)
!                  enddo
!               enddo
!            enddo


               !do j=1,param%sys%dimen
               !   do i=1, param%sys%dimen
               !      print *, d2wdx2(j,i,param%sys%mb(1))
               !   enddo
               !enddo
               !print *, ''

               !do j=1,param%sys%dimen
               !   do i=1, param%sys%dimen
               !      print *, d2wdxmdxn(j,i,param%sys%nb(1))
               !   enddo
               !enddo

               !call exit(0)


            !==========================================
            ! Evaluate Sum dWeight * dTay
            !==========================================
!            dWTdx2 = 0.0
!            dWTdxmdxn = 0.0
!            do k=1,param%sys%dimen
!               do j=1,param%sys%dimen
!                  do i=1,param%sys%nbond
!                     dWTdx2(k,j,param%sys%mb(i)) = dWTdx2(k,j,param%sys%mb(i)) + &
!                        dWTdr2(i)*dr(k,param%sys%mb(i),i)*dr(j,param%sys%mb(i),i)
!
!                     dWTdxmdxn(k,j,param%sys%mb(i)) = dWTdxmdxn(k,j,param%sys%mb(i)) + &
!                        dWTdr2(i)*dr(k,param%sys%mb(i),i)*dr(j,param%sys%nb(i),i)
!                  enddo
!               enddo
!            enddo


            !do j=1, param%sys%dimen
            !   do i=1, param%sys%dimen
            !      print *, dWTdx2(j,i,param%sys%mb(1)), dWTdxmdxn(j,i,param%sys%mb(1))
            !   enddo
            !   print *, ''
            !enddo
            !call exit(0)

            !==========================================
            ! Evaluate Sum d2Vdx2
            !==========================================
            !d2Vdx2 = 0.0
            !d2Vdxmdxn = 0.0
            !do k=1,param%sys%dimen
            !   do j=1,param%%sys%dimen
            !      do i=1, param%sys%nbond
            !         d2Vdx2(k,j,param%sys%mb(i)) = d2Vdx2(k,j,param%sys%mb(i)) + &
            !            d2Taydx2(k,j,param%sys%mb(i)) + d2wdx2(k,j,param%sys%mb(j)) &
            !            + dWTdx2(k,j,param%sys%mb(i))
            !      enddo
            !   enddo
            !enddo
            

            r=1/r

        end subroutine potential

        !Initialise variables used in the modified shepard interpolation 
        !requires the molsysdat object to be initialised before hand and passed
        !to this initialisation subroutine
        subroutine MSI_INIT(this, sys, interp_file, pot_file, atom_perm_file, num_copies)
            type (msi_params), intent(out) :: this
            type(molsysdat), intent(in) ::sys
            character(len=80), intent(in) :: interp_file
            character(len=80), intent(in) :: pot_file
            character(len=80), intent(in) :: atom_perm_file
            integer, intent(in) :: num_copies
            integer :: ierr,i

            include 'bondperms.int'
            include 'read_pot.int'

            this%sys = sys
            call read_interp(this%interp,interp_file)
            call bondperm(this%interp,this%sys,atom_perm_file)
            call read_pot(this%sys,this%interp,this%pot,pot_file)
            
            allocate(this%neighlist(num_copies),stat=ierr)
            if(ierr.ne.0) stop 'Error allocating the neighbour list array in MSI_INIT'
            do i=1,num_copies
                allocate(this%neighlist(i)%inner(this%interp%ndata*this%interp%ngroup),stat=ierr)
                if(ierr.ne.0) stop 'Error allocating inner neighbour lists'
                this%neighlist(i)%numInner=0
            enddo
            
        end subroutine MSI_INIT









      end module potential_msi
