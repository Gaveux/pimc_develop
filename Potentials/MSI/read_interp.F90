

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
   read(7,*)    interp%wtol, interp%wouter
   write(11,*)  interp%wtol, interp%wouter

   if (interp%wtol==1.d0) then
      print *, 'inner and outer neighbour list share the same weight cutoff'
      print *, 'hence outer neighbour list is disabled for better performance'
      interp%outer_disable = .TRUE.
   else
      interp%outer_disable = .FALSE.
   endif

   read(7,80)   comment_line
   write(11,81) comment_line
   read(7,*)    interp%inneigh_update, interp%outneigh_update
   write(11,*)  interp%inneigh_update, interp%outneigh_update
   
   !if (interp%inneigh_update.GT.interp%outneigh_update) then
   !   print *, 'ERROR: inner neighbour list must update more often than outer neighbour list' 
   !   call exit(1)
   if (interp%inneigh_update.EQ.interp%outneigh_update) then
      !if (interp%inneigh_update.EQ.1) then
         print *, 'update frequency for inner and outer neighbour list cannot be equal'
         call exit(1)
         !print *, 'inner and outer neighbour lists are updated every step'
      !else
        
        !call exit(0)
         !print *, 'outer neighbour list is disabled, neighbour list is updated every',interp%inneigh_update, 'step' 
      !endif
   endif

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
