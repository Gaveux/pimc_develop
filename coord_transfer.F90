


   module coordinate_transformation
     use molecule_specs
     use pimc_structures
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
      real(kind=8),dimension(2) :: x_axis_diff
 
      integer :: i,j,k
       
      ! calculate and store the translational variable into an array
      diff = Beads(1)%x(:,1)
      !print *, diff

      do k=1,sys%natom
           do i=1,pimc%NumBeadsEff
               Beads(i)%x(:,k) = Beads(i)%x(:,k) - diff
           enddo
      enddo
      call print_beads_coordinates(Beads,sys,pimc)


     end subroutine coord_trans

     subroutine reset_reference(Beads)
       !type(molsysdat), intent(in) :: sys
       !type(pimc_par) :: pimc
       type(pimc_particle), dimension(:), pointer :: Beads
      
       Beads(1)%x(:,1) = 0.d0

     end subroutine reset_reference

   end
