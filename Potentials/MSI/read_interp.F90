

subroutine read_interp(interp, filename)

   use interpolation

   implicit none
   type (interp_params), intent(out) :: interp
   character(len=80), intent(in) :: filename

   character(len=80) comment_line

80    format(a80)
81    format(1x,a79)

!----------------------------------------------------------------
! read in the parameters for the interpolation
!----------------------------------------------------------------
   print *, trim(filename)
   open(unit=7,file=trim(filename),status='old')
   write(11,*) '----------------------------------------------------'
   write(11,*) ''
   write(11,*) ' read interpolation parameters from ', trim(filename)
   write(11,*) ''
   
!  read title line
   
   read(7,80)   comment_line
   write(11,81) comment_line

!  determine one or twp part weight function is used, ipart=1,2

   read(7,80)   comment_line
   write(11,81) comment_line
   read(7,*)    interp%ipart
   write(11,91) interp%ipart
91    format(4x,'we will use the ',i1,' part weight function')

   if (interp%ipart.ne.1.and.interp%ipart.ne.2) then
      write(6,*)' ERROR: you must choose 1 or 2 part wt function'
      call exit(1)
   endif

!  read in the parameters defining the weight function, the neighbour
!  list and the number of timesteps between calls to neighbour

   read(7,80)   comment_line
   write(11,81) comment_line
   read(7,*)    interp%lowp, interp%ipow
   write(11,*)  interp%lowp, interp%ipow

   read(7,80)   comment_line
   write(11,81) comment_line
   read(7,*)    interp%wtol, interp%outer
   write(11,*)  interp%wtol, interp%outer

   read(7,80)   comment_line
   write(11,81) comment_line
   read(7,*)    interp%neigh_update
   write(11,*)  interp%neigh_update

!  read in energy maximum and minimum with which to screen data

   read(7,80)   comment_line
   write(11,81) comment_line
   read(7,*)    interp%vmax, interp%vmin
   write(11,*)  interp%vmax, interp%vmin

! read in the multipicitive factor used in radius.f, confidrad.f

   read(7,80)   comment_line
   write(11,81) comment_line
   read(7,*)    interp%nneigh
   write(11,*)  interp%nneigh

!  Enter the energy error tolerance for confidrad.f

   read(7,80)   comment_line
   write(11,81) comment_line
   read(7,*)    interp%sigma
   write(11,*)  interp%sigma

   interp%ipow2 = interp%ipow*2

! close IN_INTERP

   close(unit=7)

   return
   end
