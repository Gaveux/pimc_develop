       !a module to perform binning of internal coordinates. Can be used
       !to bin an arbitrary set of bond lengths, bond angles and
       !dihedral angles
       module binning

         implicit none
         type bin_type
           integer, dimension(:,:), allocatable :: r_bins
           integer, dimension(:,:), allocatable :: angle_bins
           integer, dimension(:,:), allocatable :: dihedral_bins

           !2d binning (I could make this general for all systems but I shalln't and so it will only work with SSSH
           integer, dimension(:,:,:), allocatable :: dihedral_r_bins

           !stores the number of times that the update Binning subroutine has been called
           integer(kind=8) :: normalise

         end type bin_type

         type binning_params
            integer, dimension(:,:), allocatable :: dihedral_indices !stores a list four atoms to be considered for each dihedral angle calculation
            integer, dimension(:,:), allocatable :: angle_indices !stores a list of the 3 atoms to be considered in each angle calculation
            integer, dimension(:), allocatable :: bondlength_indices !stores a list of the bondlengths to be binned

            integer :: num_dihedrals
            integer :: num_angles
            integer :: num_bondlengths

            integer :: num_dihedral_bins
            integer :: num_angle_bins
            integer :: num_bondlength_bins

            !stores the smallest length to be binned the upper bound of
            !the lengths to be binned, the bin step size, and 1/bin_step_size
            real(kind=8), dimension(4) :: dihedral_bin_size
            real(kind=8), dimension(4) :: angle_bin_size
            real(kind=8), dimension(4) :: length_bin_size
         end type binning_params

         interface new
            module procedure init_bin_type
            module procedure init_binning_params
         end interface new

         contains
            !initialises the binning parameter data type
            subroutine init_binning_params(this,n_bl,n_ba,n_d)
                type(binning_params), intent(out) :: this
                integer, intent(in) :: n_bl
                integer, intent(in) :: n_ba
                integer, intent(in) :: n_d
                integer :: ierr

                this%num_dihedrals = n_d
                this%num_angles = n_ba
                this%num_bondlengths = n_bl

                allocate(this%dihedral_indices(4,n_d),stat=ierr)
                if(ierr.ne.0) stop 'Error allocating the dihedral indices array'
                allocate(this%angle_indices(3,n_ba),stat=ierr)
                if(ierr.ne.0) stop 'Error allocating the bond angle indices array'
                allocate(this%bondlength_indices(n_bl),stat=ierr)
                if(ierr.ne.0) stop 'Error allocating the bond length indices array'

            end subroutine init_binning_params

            !initialise the bins data type
            subroutine init_bin_type(this,params)
                type(bin_type), intent(out) :: this
                type(binning_params), intent(in) :: params
                integer :: ierr,i,j,k
                this%normalise = 0

                allocate(this%r_bins(params%num_bondlengths,params%num_bondlength_bins), stat=ierr)
                if(ierr.ne.0) stop 'Error allocating the bond length bins'
                allocate(this%angle_bins(params%num_angles,params%num_angle_bins), stat=ierr)
                if(ierr.ne.0) stop 'Error allocating the bond angle bins'
                allocate(this%dihedral_bins(params%num_dihedrals,params%num_dihedral_bins), stat=ierr)
                if(ierr.ne.0) stop 'Error allocating the dihedral angle bins'

                !can be generalised if more variables are added to the data structures above to allow for the generalisation
                allocate(this%dihedral_r_bins(1,params%num_dihedral_bins,params%num_bondlength_bins),stat=ierr)
                if(ierr.ne.0) stop 'Error allocating dihedral r'

                do i=1,params%num_bondlength_bins
                    do j=1,params%num_bondlengths
                        this%r_bins(j,i)=0
                    enddo
                enddo
                do i=1,params%num_dihedral_bins
                    do j=1,params%num_dihedrals
                        this%dihedral_bins(j,i)=0
                    enddo
                enddo
                do i=1,params%num_angle_bins
                    do j=1,params%num_angles
                        this%angle_bins(j,i)=0
                    enddo
                enddo
                !can be generalised if more data structures are added to the system
                do i=1,params%num_bondlength_bins
                    do j=1,params%num_dihedral_bins
                        do k=1,1
                            this%dihedral_r_bins(k,j,i)=0
                        enddo
                    enddo
                enddo

            end subroutine init_bin_type

            subroutine read_number_hist(input,n_bond,n_bl,n_ba,n_d)
                integer, intent(in) :: n_bond
                character(len=80), intent(in) :: input
                integer, intent(out):: n_bl, n_ba, n_d
                character(len=80) :: icomm

281      format(a80)
                open (unit=245,file=trim(input),status='old')

                !read header line
                read(245,281) icomm

                !read in the number of bondlengths, bond angles and dihedral angles and initialise the binning parameters
                read(245,281) icomm
                read(245,*) n_bl, n_ba, n_d

                if(n_bl.eq.-1) then
                    n_bl = n_bond
                endif

                close(unit=245)
                return

            end subroutine read_number_hist

            !reads the binning params from an input file
            subroutine read_binning_params(this,input,n_bond)
                type(binning_params), intent(inout) :: this
                character(len=80), intent(in) :: input
                integer, intent(in) :: n_bond
                character(len=80) :: icomm
                real(kind=8) :: valmin,valmax
                integer :: i
                integer :: a1,a2,a3,a4
                integer :: n_bl, n_ba, n_d

                open (unit=245,file=trim(input),status='old')
280      format(a80)

                !read header line
                read(245,280) icomm
                !read in the number of bondlengths, bond angles and dihedral angles and initialise the binning parameters
                read(245,280) icomm
                read(245,*) n_bl, n_ba, n_d

                n_bl = this%num_bondlengths
                n_ba = this%num_angles
                n_d = this%num_dihedrals

                !read in the bond length bin sizes
                read(245,280) icomm
                read(245,*) valmin,valmax,this%num_bondlength_bins
                call bin_shape(this,valmin,valmax,this%num_bondlength_bins,1)

                !read in the bond angle bin sizes
                read(245,280) icomm
                read(245,*) valmin,valmax,this%num_angle_bins
                call bin_shape(this,valmin,valmax,this%num_angle_bins,2)

                !read in the dihedral angle bin sizes
                read(245,280) icomm
                read(245,*) valmin,valmax,this%num_dihedral_bins
                call bin_shape(this,valmin,valmax,this%num_dihedral_bins,3)

                !read newline
                read(245,280) icomm

                !read in the specific bond lengths to be binned
                read(245,280) icomm
                read(245,280) icomm
                if(n_bl.eq.n_bond) then
                    do i=1,n_bond
                        call add_bondlength(this,i,i)
                    enddo
                else
                    do i=1,n_bl
                        read(245,*) a1
                        call add_bondlength(this,i,a1)
                    enddo
                endif
                read(245,280) icomm

                !read in the specific bond angles to be binned
                read(245,280) icomm
                do i=1,n_ba
                    read(245,*) a1,a2,a3
                    call add_angle(this,a1,a2,a3,i)
                enddo
                read(245,280) icomm

                !read in the specific dihedral angles to be binned
                read(245,280) icomm
                do i=1,n_d
                    read(245,*) a1,a2,a3,a4
                    call add_dihedral(this,a1,a2,a3,a4,i)
                enddo

                !now the bond length parameters should all have been read in and the params object initialised
                !after returning from this the actual bin_type object can be initialised

                close (unit=245)
                return

            end subroutine read_binning_params


            !a couple of setters for the params type

            !adds a dihedral angle to the list of dihedrals to be looked at
            !assumes a1,a2,a3,a4 have valid atom indices
            subroutine add_dihedral(this,a1,a2,a3,a4,ind)
                type(binning_params)  :: this
                integer, intent(in) :: a1,a2,a3,a4,ind
                if(ind>this%num_dihedrals) stop 'Error ind out of bound'

                this%dihedral_indices(1,ind)=a1
                this%dihedral_indices(2,ind)=a2
                this%dihedral_indices(3,ind)=a3
                this%dihedral_indices(4,ind)=a4
                return
            end subroutine add_dihedral

            !adds a bond angle to the list of dihedrals to be looked at
            subroutine add_angle(this,a1,a2,a3,ind)
                type(binning_params)  :: this
                integer, intent(in) :: a1,a2,a3,ind
                if(ind>this%num_angles) stop 'Error ind out of bound'

                this%angle_indices(1,ind)=a1
                this%angle_indices(2,ind)=a2
                this%angle_indices(3,ind)=a3
                return
            end subroutine add_angle

            !adds a bond length to the list of angles to be looked at
            subroutine add_bondlength(this,a1,ind)
                type(binning_params)  :: this
                integer, intent(in) :: a1,ind
              if(ind>this%num_bondlengths) stop 'Error ind out of bound'

                this%bondlength_indices(ind)=a1
                return
            end subroutine add_bondlength

            !ind specifies whether this is for bondlength (1), bondangle (2), or dihedral(3)
            subroutine bin_shape(this,valmin,valmax,num, ind)
                type(binning_params)  :: this
                real(kind=8), intent(in) :: valmin,valmax
                integer, intent(in) :: num,ind
                real(kind=8) :: dx
                dx = (valmax-valmin)/dble(num)
                if(ind==1) then
                    this%length_bin_size(1)=valmin
                    this%length_bin_size(2)=valmax
                    this%length_bin_size(3)=dx
                    this%length_bin_size(4)=1.0d0/dx
                elseif (ind==2) then
                    this%angle_bin_size(1)=valmin
                    this%angle_bin_size(2)=valmax
                    this%angle_bin_size(3)=dx
                    this%angle_bin_size(4)=1.0d0/dx
                elseif(ind ==3 ) then
                    this%dihedral_bin_size(1)=valmin
                    this%dihedral_bin_size(2)=valmax
                    this%dihedral_bin_size(3)=dx
                    this%dihedral_bin_size(4)=1.0d0/dx
                else
                    stop 'bin_shape ind out of bounds'
                endif

            end subroutine bin_shape


            !need to add a subroutine for initialising theses and
            !reading the parameters in from file

            !performs the binning of the bondlengths, bond angles and
            !dihedral angles
            subroutine updateBinning(bins, params, x, r)
                type(bin_type) :: bins
                type(binning_params), intent(in) :: params
                real(kind=8), dimension(:,:), intent(in) :: x
                real(kind=8), dimension(:), intent(in) :: r
                integer :: i, inc, inc2,ind, NumBins, NumBins2
                integer, dimension(3) :: aA
                integer, dimension(4) :: dA
                real(kind=8) :: minx,invdx, angle, dihedral, minx2, invdx2

                minx = params%length_bin_size(1)
                invdx = params%length_bin_size(4)

                NumBins = params%num_bondlength_bins
                do i=1,params%num_bondlengths

                    ind = params%bondlength_indices(i)
                    inc = int((r(ind)-minx)*invdx) + 1
                    if(inc > NumBins .or. inc < 1) then
                    else
                        bins%r_bins(i,inc) = bins%r_bins(i,inc) + 1
                    endif

                enddo

                minx = params%angle_bin_size(1)
                invdx = params%angle_bin_size(4)
                NumBins = params%num_angle_bins

                do i=1,params%num_angles
                    aA = params%angle_indices(:,i)
                    call computeAngle(x(:,aA(1)),x(:,aA(2)),x(:,aA(3)),angle)
                    angle = angle*57.2957795131
                    inc = int((angle-minx)*invdx) + 1

                    if(inc > NumBins .or. inc < 1) then
                    else
                        bins%angle_bins(i,inc) = bins%angle_bins(i,inc)+ 1
                    endif

                enddo

                minx = params%dihedral_bin_size(1)
                invdx = params%dihedral_bin_size(4)
                NumBins = params%num_dihedral_bins

                !this code here is specific to HSSS remove for other things
                minx2 = params%length_bin_size(1)
                invdx2 = params%length_bin_size(4)

                NumBins2 = params%num_bondlength_bins

                !HSSS specific dihedrals (shalln't need)
                do i=1,params%num_dihedrals
                    dA = params%dihedral_indices(:,i)
                    !here is some super specific HSSS code to deal with the probability that the H can migrate to the
                    if(r(4) > r(6)) then
                        call computeDihedral(x(:,dA(1)),x(:,dA(2)),x(:,dA(3)),x(:,dA(4)),dihedral)
                        inc2 = int((r(4)-minx2)*invdx2) + 1
                    else
                        call computeDihedral(x(:,dA(3)),x(:,dA(2)),x(:,dA(1)),x(:,dA(4)),dihedral)
                        inc2 = int((r(1)-minx2)*invdx2) + 1
                    endif


                    dihedral = dihedral*57.2957795131
                    inc = int((dihedral-minx)*invdx) + 1

                    if(inc > NumBins .or. inc < 1) then
                    else
                        bins%dihedral_bins(i,inc) = bins%dihedral_bins(i,inc) + 1
                        if(inc2 > NumBins2 .or. inc2 < 1) then

                        else
                            bins%dihedral_r_bins(i,inc,inc2) = bins%dihedral_r_bins(i,inc,inc2)+1
                        endif
                    endif
                enddo

                !update the variable counting the number of bins
                bins%normalise = bins%normalise + 1

                return
            end subroutine updateBinning

            !returns the dihedral formed between the atoms at position
            !x1,x2,x3,x4 using the alternative definition given
            !assumes that the vectors formed from x1,x2,x3,x4 aren't
            !collinear as when they are the dihedral angle is ill
            !defined
            subroutine computeDihedral(x1,x2,x3,x4,dihedral)
              real(kind=8), dimension(3), intent(in) :: x1,x2,x3,x4
              real(kind=8), intent(out) :: dihedral
              real(kind=8), dimension(3) :: r1,r2,r3
              real(kind=8), dimension(3) :: b12, b23
              real(kind=8) :: norm_r2
              r1 = x2-x1
              r2 = x3-x2
              r3 = x4-x3
              b12 = cross(r1,r2)
              b23 = cross(r2,r3)
              norm_r2 = sqrt(dot_product(r2,r2))
              dihedral = atan2(dot_product(cross(b12,b23),r2)/norm_r2,dot_product(b12,b23))
              !if(dihedral < 0) dihedral = 6.28318530718+dihedral
              return
            end subroutine computeDihedral

            !computes the angle formed by 3 atoms at positions x1,x2,x3
            subroutine computeAngle(x1,x2,x3,angle)
                real(kind=8), dimension(3), intent(in) :: x1,x2,x3
                real(kind=8), intent(out) :: angle
                real(kind=8), dimension(3) :: r1,r2
                r1 = x1-x2
                r2 = x3-x2
                angle = acos(dot_product(r1,r2)/sqrt(dot_product(r1,r1)*dot_product(r2,r2)))
                return
            end subroutine computeAngle

            function cross(x1,x2) result(res)
              real(kind=8), dimension(3), intent(in) :: x1,x2
              real(kind=8), dimension(3) :: res

              res(1) = x1(2)*x2(3) - x1(3)*x2(2)
              res(2) = x1(3)*x2(1) - x1(1)*x2(3)
              res(3) = x1(1)*x2(2) - x1(2)*x2(1)
              return
            end function cross


            !doesn't print a header to the binning as the bins are
            !printed in the order specified in the IN_BINNING input file
            subroutine writeHistogram(bins,params,OUT_DIR)
              type(bin_type), intent(in) :: bins
              type(binning_params), intent(in) :: params
              character(len=80), intent(in) :: OUT_DIR
              integer :: i,j,k
              real(kind=8) :: minx, dx, minx2,dx2

                open(unit=3,file=trim(OUT_DIR)//'/hist_out.bond',status='unknown')
                minx = params%length_bin_size(1)
                dx = params%length_bin_size(3)
                do j=1,params%num_bondlength_bins
                    write(3,*) minx+(j-0.5)*dx,(dble(bins%r_bins(i,j))/dble(bins%normalise*dx),i=1,params%num_bondlengths)
                enddo
                close(unit=3)

                open(unit=3,file=trim(OUT_DIR)//'/hist_out.angle',status='unknown')
                minx = params%angle_bin_size(1)
                dx = params%angle_bin_size(3)
                do j=1,params%num_angle_bins
                    write(3,*) minx+(j-0.5)*dx,(dble(bins%angle_bins(i,j))/dble(bins%normalise*dx),i=1,params%num_angles)
                enddo
                close(unit=3)

                open(unit=3,file=trim(OUT_DIR)//'/hist_out.dihedral',status='unknown')
                minx = params%dihedral_bin_size(1)
                dx = params%dihedral_bin_size(3)
                do j=1,params%num_dihedral_bins
                    write(3,*) minx+(j-0.5)*dx,(dble(bins%dihedral_bins(i,j))/dble(bins%normalise*dx),i=1,params%num_dihedrals)
                enddo

                close(unit=3)


                open(unit=3,file=trim(OUT_DIR)//'/hist_out.dihedral_r',status='unknown')
                minx=params%dihedral_bin_size(1)
                minx2=params%length_bin_size(1)
                dx = params%dihedral_bin_size(3)
                dx2 = params%length_bin_size(3)
                do k=1,params%num_bondlength_bins
                    do j=1,params%num_dihedral_bins
                        write(3,*) minx+(j-0.5)*dx, minx2+(k-0.5)*dx2,(dble(bins%dihedral_r_bins(1,j,k))/dble(bins%normalise*dx))
                    enddo
                    write(3,*) ""
                enddo
                close(unit=3)

            end subroutine writeHistogram

            !called at the end of each block. If this is the first time
            !the subroutine is called it overwrites the old TOUT otherwise
            !it appends to TOUT
            subroutine writeTOUT(x, V, OUT_DIR, close_file,natom,dimen)
              real(kind=8), dimension(:,:), intent(in) :: x
              real(kind=8), intent(in) :: V
              character(len=80), intent(in) :: OUT_DIR
              logical, intent(in) :: close_file
              integer, intent(in) :: natom,dimen
              logical, save :: first_time = .TRUE.
              integer :: i,j

              !if we have not been told to close the file write to tout
              if(.not.close_file) then
              !if the first time this is called overwrite the old tout file
              if(first_time) then
                first_time=.FALSE.
              open(unit=999,file=trim(OUT_DIR)//'/TOUT',status='unknown',action='write')
              endif

              do i=1,natom
                write(999,*) (x(j,i),j=1,dimen)
              enddo
              write(999,*) V

              !if it is time to close the file close the file
              else
                close(unit=999)
              endif

            end subroutine writeTOUT

            ! Write the beads geometry configuration into a checkpoint file
            ! at the end of each block, this is called in pimc_monte.F90	
            subroutine writeCheckpoint(x,CHECKPOINT_DIR,start,WritingCheckpoint,close_file,natom,dimen)
               real(kind=8), dimension(:,:), intent(in) :: x
               character(len=80), intent(in) :: CHECKPOINT_DIR, start
               logical, intent(in) :: close_file
               integer, intent(in) :: natom, dimen
               integer :: i,j
               character(len=1), intent(in) :: WritingCheckpoint

               if(WritingCheckpoint == 'y') then
                  ! if it is told not to close the file write to checkpint
                  if(.not.close_file) then
                    !open(unit=599,file=trim(CHECKPOINT_DIR)//'/checkpoint',status='unknown',action='write')
                    open(unit=599,file=trim(checkpoint_dir)//trim(start),status='unknown',action='write')
                    ! print *, 'writing bead configuration into', trim(checkpoint_dir)//trim(start)
                    ! this is Cartesian coordinates of atom in a single bead form
                    ! see pimc_monte.F90 for writing all beads configurations
                    do i=1,natom
                       !do j=1,dimen
                          write(599,*) (x(j,i),j=1,dimen)
                       !enddo
                    enddo
            
                  else 
                     close(unit=599)
                  endif

               elseif(WritingCheckpoint == 'n') then
                  !print *, 'print rubbish'
               endif
            end subroutine
       end module binning
