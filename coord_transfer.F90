


   module coordinate_transformation
     use molecule_specs
     use pimc_structures
     use quaternion_mod
     use debug

     implicit none

     interface trans
        module procedure coord_translation
     end interface

     contains
     

     subroutine coord_translation(Beads,pimc,sys)
      implicit none
      
      type(molsysdat),intent(in) :: sys
      type(pimc_par) :: pimc
      type(pimc_particle), dimension(:), pointer :: Beads
      
      real(kind=8),dimension(sys%dimen) :: diff
 
      integer :: i,j,k
      
       
      ! calculate and store the translational variable into an array
      diff = Beads(1)%x(:,1)
      !print *, diff

      do k=1,sys%natom
           do i=1,pimc%NumBeadsEff
               Beads(i)%x(:,k) = Beads(i)%x(:,k) - diff
           enddo
      enddo
      !do i=1,pimc%NumBeadsEff
         !do j=1,sys%natom
         !   call print_Beads_distance(i,j,Beads,sys,pimc)
         !enddo
      !enddo

      !call print_beads_coordinates(Beads,sys,pimc)
      !print *, ''


     end subroutine coord_translation
    
     subroutine coord_rotation(Beads,pimc,sys)
     !subroutine coord_rotation(Beads,pimc,sys,seedval)
       implicit none

       type(molsysdat), intent(in) :: sys
       type(pimc_par) :: pimc
       type(pimc_particle),dimension(:), pointer :: Beads
       type(quaternion) :: rot
       !type(mod_seed), intent(inout) :: seedval
       type(quaternion) :: quad_x, quad_xy
       
       real(kind=8),dimension(sys%dimen) :: x_axis, y_axis
       
       integer :: i,j,k, ind
       real(kind=8) :: theta, pi
       pi=acos(-1.d0)
       x_axis= [1.d0,0.d0,0.d0]
       
       !call print_beads_coordinates(Beads,sys,pimc)
        
       ! map bead 1 of atom 2 onto x-axis
       quad_x = get_rotation_between(Beads(1)%x(:,2), x_axis)

       ! rotate the entire system 

       do k =1, pimc%NumBeadsEff
          do i = 1, sys%natom
           Beads(k)%x(:,i) = quaternion_rotate(quad_x, Beads(k)%x(:,i))
          enddo
       enddo
       !call print_beads_coordinates(Beads,sys,pimc)
       !do i = 1, pimc%NumBeadsEff
       !   call print_Beads_distance(i,Beads,sys,pimc)
       !enddo

       ! map bead 1 of atom 2 onto the xy-plane that rotates around x-axis
       ! calculates the rotation angle around x-axis
       if (Beads(1)%x(3,3) .gt. 0.d0) then
         if (Beads(1)%x(2,3) .gt. 0.d0) then
            theta = -atan(Beads(1)%x(3,3)/Beads(1)%x(2,3))
         elseif (Beads(1)%x(2,3) .gt. 0.d0) then
            theta = -acos(pi*0.5)
         else
            theta = -pi - atan(Beads(1)%x(3,3)/Beads(1)%x(2,3))
         endif
       elseif (Beads(1)%x(3,3) .lt. 0.d0) then
         if (Beads(1)%x(2,3) .gt. 0.d0) then
            theta = -2.0*pi - atan(Beads(1)%x(3,3)/Beads(1)%x(2,3))
         elseif (Beads(1)%x(3,3) .eq. 0.d0) then
            theta = -acos(pi*3.0/2.0)
         else
            theta = -pi + atan(Beads(1)%x(3,3)/Beads(1)%x(2,3))
         endif
       else 
         if (Beads(1)%x(2,3) .ge. 0.d0) then
            theta = 0.d0
         elseif (Beads(1)%x(2,3) .le. 0.d0) then
            theta = -pi
         endif
       endif
       
       ! construct the quaternion using the theta
       call quaternion_val(quad_xy, cos(theta*0.5), sin(theta*0.5), 0.d0, 0.d0)
       
       ! rotate the entire system
       do k =1, pimc%NumBeadsEff
          do i = 1, sys%natom
             Beads(k)%x(:,i) = quaternion_rotate(quad_xy, Beads(k)%x(:,i))
          enddo
       enddo
      !call print_beads_coordinates(Beads,sys,pimc)
      ! do i = 1, pimc%NumBeadsEff
      !    call print_Beads_distance(i,Beads,sys,pimc)
      ! enddo


     end subroutine coord_rotation

     subroutine reset_reference(Beads)
       type(pimc_particle), dimension(:), pointer :: Beads
      
       Beads(1)%x(:,1) = 0.d0

     end subroutine reset_reference

   end
