


   module coordinate_transformation
     use molecule_specs
     use pimc_structures

     implicit none

     interface trans
        module procedure coord_trans
     end interface

     contains
     

     ! take the first bead of the first atom as the origin
     subroutine coord_trans(Beads,pimc,sys)
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
       
     end subroutine coord_trans

   end
