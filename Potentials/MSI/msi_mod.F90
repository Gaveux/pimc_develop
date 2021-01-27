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
            real(kind=8), dimension(param%sys%dimen,param%sys%natom)  :: dFsqr 

            !Variables used internally by the modified shepard code
              
            !derivative of bond lengths with respect to cartesians
            real(kind=8), dimension(param%sys%dimen,param%sys%natom,param%sys%nbond) :: dr
            !second derivative of bondlengths w.r.t. Cartesians
            real(kind=8), dimension(param%sys%dimen,param%sys%natom,param%sys%nbond) :: d2r
            !derivative of bondlengths w.r.t. Cartesians squared
            real(kind=8), dimension(param%sys%dimen,param%sys%natom,param%sys%nbond) :: dr2
            !Derivatives of the potential with respect to internal coordinates
            real(kind=8), dimension(param%sys%nbond) :: dVdr
            !Dummy variables for Weight, RawWeightTemp is Weight^-p
            real(kind=8), dimension(param%interp%ndata) :: Weight
            real(kind=8), dimension(param%interp%ndata) :: RawWeightTemp
            !derivative of |F|^2 w.r.t. inverse bonde lengths
            !real(kind=8), dimension(param%sys%nbond) :: ddr_Fsqr

            ! sum d2wdx2Tay
            !real(kind=8), dimension(param%sys%dimen,param%sys%natom),intent(out) :: d2wdx2Tay
            real(kind=8), dimension(param%sys%dimen,param%sys%natom) :: d2wdx2Tay
            
            ! 6 terms in d2wdx2tay
            real(kind=8), dimension(param%sys%nbond) :: SumD2Weight1
            real(kind=8), dimension(param%sys%nbond) :: SumD2Weight2
            real(kind=8), dimension(param%sys%nbond) :: SumD2Weight3
            real(kind=8), dimension(param%sys%nbond) :: SumD2Weight4
            real(kind=8), dimension(param%sys%nbond) :: SumD2Weight5
            real(kind=8), dimension(param%sys%nbond) :: SumD2Weight6

            ! 3 terms in w*d2Taydx2
            real(kind=8), dimension(param%sys%nbond) :: SumWeightd2Taydx2tmp1
            real(kind=8), dimension(param%sys%nbond) :: SumWeightd2Taydx2tmp2
            real(kind=8), dimension(param%sys%nbond) :: SumWeightd2Taydx2tmp3

            ! w*d2Taydx2
            real(kind=8), dimension(param%sys%dimen,param%sys%natom) :: wd2Tay

            ! 2 terms in dwdx*dTaydx
            real(kind=8), dimension(param%sys%nbond) :: dwdrdTaydrtmp1
            real(kind=8), dimension(param%sys%nbond) :: dwdrdTaydrtmp2
            ! dwdx*dTaydx
            real(kind=8), dimension(param%sys%dimen,param%sys%natom) :: dwdxdTaydx

            ! d2Vdx2
            real(kind=8), dimension(param%sys%dimen,param%sys%natom) :: d2Vdx2

            ! dF^2/dx
            real(kind=8), dimension(param%sys%dimen,param%sys%natom) :: dFsqrdx

            
            integer :: i,j,k
            include 'intern.int'
            include 'neigh.int'
            include 'calcen.int'

            call intern(param%sys,x,r,dr,d2r,dr2)

            !Update the inner neighbour list each potential evaluation
            call neighbour(param%sys,param%interp,param%pot,Weight,r,param%neighlist(ind),RawWeightTemp)

            !Interpolate the surface - currently only the one part weight function is implemented
            !the two part weight function can easily be included, however there will need to be 
            !a slight shuffling of the variables in the calcen2w.f90 file to make them compatible
            !with the rest of the code
            !if (param%interp%ipart == 1) then
            call calcen(param%sys,param%interp,param%pot,param%neighlist(ind), &
            Weight,r,V,dVdr,RawWeightTemp,SumD2Weight1,SumD2Weight2,SumD2Weight3, &
            SumD2Weight4,SumD2Weight5,SumD2Weight6,SumWeightd2Taydx2tmp1,&
            SumWeightd2Taydx2tmp2,SumWeightd2Taydx2tmp3,dwdrdTaydrtmp1,&
            dwdrdTaydrtmp2)
            !endif

            V = V - param%interp%vmin
            !do j=1,param%sys%natom
            !    do k=1,param%sys%dimen
            !        dV(k,j) = 0.d0
            !    enddo
            !enddo

            dV = 0.0
            d2wdx2Tay = 0.0
            wd2Tay = 0.0
            dwdxdTaydx = 0.0
            d2Vdx2 = 0.0
            do j=1,param%sys%nbond
                do k=1,param%sys%dimen
                    dV(k,param%sys%mb(j)) = dV(k,param%sys%mb(j))+dVdR(j)* & 
                         dr(k,param%sys%mb(j),j)
                    dV(k,param%sys%nb(j)) = dV(k,param%sys%nb(j))+ & 
                         dVdR(j)*dr(k,param%sys%nb(j),j)

                    ! d2wdx2*Tay
                    d2wdx2Tay(k,param%sys%mb(j)) = d2wdx2Tay(k,param%sys%mb(j)) &
                    + SumD2Weight1(j)*dr2(k,param%sys%mb(j),j) &
                    + d2r(k,param%sys%mb(j),j)*SumD2Weight2(j) &
                    + 2.d0*dr2(k,param%sys%mb(j),j)*SumD2Weight3(j) &
                    - 2.d0*dr2(k,param%sys%mb(j),j)*SumD2Weight4(j) &
                    + dr2(k,param%sys%mb(j),j)*SumD2Weight5(j) &
                    - d2r(k,param%sys%mb(j),j)*SumD2Weight6(j)

                    d2wdx2Tay(k,param%sys%nb(j)) = d2wdx2Tay(k,param%sys%nb(j)) &
                    + SumD2Weight1(j)*dr2(k,param%sys%nb(j),j) &
                    + d2r(k,param%sys%nb(j),j)*SumD2Weight2(j) &
                    + 2.d0*dr2(k,param%sys%nb(j),j)*SumD2Weight3(j) &
                    - 2.d0*dr2(k,param%sys%nb(j),j)*SumD2Weight4(j) &
                    + dr2(k,param%sys%nb(j),j)*SumD2Weight5(j) &
                    - d2r(k,param%sys%nb(j),j)*SumD2Weight6(j)

                    ! w*d2Taydx2
                    wd2Tay(k,param%sys%mb(j))=wd2Tay(k,param%sys%mb(j))&
                    +dr2(k,param%sys%mb(j),j)*SumWeightd2Taydx2tmp1(j) &
                    +dr2(k,param%sys%mb(j),j)*SumWeightd2Taydx2tmp2(j) &
                    +d2r(k,param%sys%mb(j),j)*SumWeightd2Taydx2tmp3(j)

                    wd2Tay(k,param%sys%nb(j))=wd2Tay(k,param%sys%nb(j))&
                    +dr2(k,param%sys%nb(j),j)*SumWeightd2Taydx2tmp1(j) &
                    +dr2(k,param%sys%nb(j),j)*SumWeightd2Taydx2tmp2(j) &
                    +d2r(k,param%sys%nb(j),j)*SumWeightd2Taydx2tmp3(j)

                    ! dwdx*dTaydx
                    dwdxdTaydx(k,param%sys%mb(j))=dwdxdTaydx(k,param%sys%mb(j))&
                    + dr2(k,param%sys%mb(j),j)*dwdrdTaydrtmp1(j) &
                    - dr2(k,param%sys%mb(j),j)*dwdrdTaydrtmp2(j)

                    dwdxdTaydx(k,param%sys%nb(j))=dwdxdTaydx(k,param%sys%nb(j))&
                    + dr2(k,param%sys%nb(j),j)*dwdrdTaydrtmp1(j) &
                    - dr2(k,param%sys%nb(j),j)*dwdrdTaydrtmp2(j)

                    ! d2Vdx2
                    d2Vdx2(k,param%sys%mb(j)) = d2Vdx2(k,param%sys%mb(j)) &
                    + d2wdx2Tay(k,param%sys%mb(j)) + wd2Tay(k,param%sys%mb(j)) &
                    + dwdxdTaydx(k,param%sys%mb(j))

                    d2Vdx2(k,param%sys%nb(j)) = d2Vdx2(k,param%sys%nb(j)) &
                    + d2wdx2Tay(k,param%sys%nb(j)) + wd2Tay(k,param%sys%nb(j)) &
                    + dwdxdTaydx(k,param%sys%nb(j))
                enddo
            enddo
            !call exit(0)

            ! d/dx (dV/dx)^2 = 2(dV/dx) d2V/dx2
            do j=1,param%sys%natom
               do i=1,param%sys%dimen
                  dFsqrdx(i,j) = 2.0*dV(i,j)*d2Vdx2(i,j)
               enddo
            enddo
            print *, dFsqrdx
            print *, '------------------'
            print *, dV

            r=1/r
           ! print *, dV
           ! print *, '-------------------------'
           ! print *, d2wdx2Tay
           ! print *, '-------------------------'
           ! print *, wd2Tay
           !print *, d2Vdx2
           call exit(0)
           

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
