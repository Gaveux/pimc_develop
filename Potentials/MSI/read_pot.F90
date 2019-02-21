
      subroutine read_pot(sys,interp,allpot,filename)

      use molecule_specs
      use interpolation

      implicit none
      type (molsysdat), intent(in) :: sys
      type (interp_params), intent(inout) :: interp
      type (pot_data_point), dimension(:), pointer :: pot
      type (pot_data_point), dimension(:), pointer :: allpot
      character(len=80), intent(in) :: filename

      integer  :: ierr,i,j,k,counter,n, reject
      real(kind=8) :: vmin_pot
      character(len=80) :: comment_line

80    format(a80)

!----------------------------------------------------------------
! read in the derivative data from which the PES is interpolated
!----------------------------------------------------------------

      write(11,*) '----------------------------------------------------'
      write(11,*) ' '
      write(11,*) '  read potential energy data from ', trim(filename)
      write(11,*) ' '

      open(unit=3,file=trim(filename),status='old')
      
!  read through the file to find out how many data points it contains

      interp%ndata = 0
      read(3,80) comment_line
      do 
        read(3,80,iostat=ierr,end=11) comment_line
        if (index(comment_line,'data') /= 0) then
           interp%ndata = interp%ndata + 1
        endif
      enddo

11    continue

!  allocate space

      allocate(pot(interp%ndata),stat=ierr)
      if (ierr.ne.0) stop '  Error allocating pot '

      do i=1,interp%ndata
        !print *, ' allocating pot(',i,')'
         call new(pot(i),sys%nbond)
      enddo

!  rewind and do a formatted read

      rewind(unit=3,iostat=ierr)
      if (ierr /= 0) stop ' Error rewinding POT file '
      
      read(3,80) comment_line
      reject = 0
      i=0
      do while (i <= interp%ndata)
        i=i+1
        ! read in the potential energy data, including the 
        ! definitions of the normal local coordinates
        read(3,80,err=20,end=21) comment_line
        read(3,*)  (pot(i)%r(j),j=1,sys%nbond)
        do j=1,sys%nint
           read(3,*) (pot(i)%ut(j,k),k=1,sys%nbond)
        enddo
        read(3,*)   pot(i)%v0
        read(3,*)  (pot(i)%v1(j),j=1,sys%nint)
        read(3,*)  (pot(i)%v2(j),j=1,sys%nint)
        ! if the energy is too high skip this point
        if (pot(i)%v0 > interp%vmax) then
           i = i - 1
           interp%ndata = interp%ndata - 1
           reject = reject + 1
        endif
      enddo

20    if (ierr >= 0) stop ' Error reading from POT file '
21    continue
      
      close(unit=3)

      write(11,*) ' Read in',interp%ndata,' data points from ', trim(filename)
      write(11,*) ' rejected ',reject,' data points as higher than vmax'

! set vmin to be the minimum of the value read in from IN_INTERP and the
! minimum found in the POT file 

      vmin_pot = 10.0
      do i=1,interp%ndata
         if (pot(i)%v0 < vmin_pot) then
            vmin_pot = pot(i)%v0
         endif
      enddo

      write(11,*)' lowest potential energy in the POT file          = ',vmin_pot
      write(11,*)' lowest potential energy expected from IN_INTERP  = ',interp%vmin


!     if (vmin_pot < interp%vmin) then
         interp%vmin = vmin_pot
!     endif

      write(11,*)' lowest potential energy used (ie value of vmin)  = ',interp%vmin

!  invert data bonds once and for all

      do i=1,interp%ndata
         do j=1,sys%nbond
            pot(i)%r(j) = 1.0/pot(i)%r(j)
         enddo
      enddo

      !  construct the local coords for each data point
      !  remembering that rda are read in as bond lengths NOT inverses 
      !  and that we have inverted them 

      do i=1,interp%ndata
         do j=1,sys%nint
         pot(i)%z(j) = 0.0
         do k=1,sys%nbond
            pot(i)%z(j) =  pot(i)%z(j) + pot(i)%ut(j,k)*pot(i)%r(k)
         enddo
         enddo
      enddo

! expand permutations once and for all

      ! make space for new array
      allocate(allpot(interp%ngroup*interp%ndata),stat=ierr)
      if (ierr.ne.0) stop ' Error allocating allpot '

      do i=1,interp%ndata*interp%ngroup
         call new(allpot(i),sys%nbond)
      enddo

      ! copy permuted data into new array, checking for and 
      ! removing identical geometries
      counter = 0
      do i=1,interp%ndata
      do n=1,interp%ngroup
!        diff = 0.0
!        do j=1,sys%nbond
!           diff = diff + abs(pot(i)%r(j) -  pot(i)%r(sys%bperm(j,n)))
!        enddo
!        if ((n > 1).and.(diff < 1.d-6)) cycle   ! skip perms with identical bonds
         counter = counter + 1
         do j=1,sys%nbond
            allpot(counter)%r(j) = pot(i)%r(interp%bperm(j,n))
         enddo
         do j=1,sys%nbond
         do k=1,sys%nint
            allpot(counter)%ut(k,j) = pot(i)%ut(k,interp%bperm(j,n))
         enddo
         enddo
         allpot(counter)%v0 = pot(i)%v0
         allpot(counter)%v1 = pot(i)%v1
         allpot(counter)%v2 = pot(i)%v2
         allpot(counter)%z  = pot(i)%z 
         allpot(counter)%id = i
         allpot(counter)%perm = n
      enddo
      enddo

      ! update global counter 
      interp%ndata = interp%ndata*interp%ngroup

      print *, '# Number of data points used for interpolation:',interp%ndata
         
      ! deallocate the pot array now that we're finished with it
      do i=1,interp%ndata/interp%ngroup
         call destroy(pot(i))
      enddo

      deallocate(pot,stat=ierr)
      if (ierr.ne.0) stop ' Error deallocating pot '


      return
      end

