
      subroutine bondperm(interp,sys,filename)

      use interpolation
      use molecule_specs

      type (interp_params), intent(inout) :: interp
      type (molsysdat), intent(in) :: sys
      character(len=80), intent(in) :: filename

      character(len=80) icomm
      integer idummy
      integer i,j,k,n,ic

!  read atomic permutations into the aperm array

      open(unit=8, file=trim(filename),status='old')

      read(8,80)icomm
80    format(a80)

!  read in a comment line and the order of the group

      read(8,80)icomm
      read(8,*) n

!  allocate sys%aperm and sys%bperm 

      call perms(interp,n,sys)

!  read in the atomic perm

      do n=1,interp%ngroup
        read(8,80)icomm
        do k=1,sys%natom
          read(8,*)idummy,interp%aperm(k,n)   ! idummy will equal k
        enddo
      enddo

      close(unit=8)

!  now set up the bperm array of bond permutations by looking at aperm

      do n=1,interp%ngroup

        ic=0
        do i=1,sys%natom-1
        do j=i+1,sys%natom

          ic=ic+1
          do k=1,sys%nbond
            if (sys%mb(k).eq.interp%aperm(i,n).and.sys%nb(k).eq.interp%aperm(j,n)) then
               interp%bperm(k,n)=ic
            endif
            if (sys%nb(k).eq.interp%aperm(i,n).and.sys%mb(k).eq.interp%aperm(j,n)) then
               interp%bperm(k,n)=ic
            endif
          enddo

        enddo
        enddo

      enddo

      return
      end
