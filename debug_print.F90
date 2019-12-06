
   module debug

     use molecule_specs
     use pimc_structures

     implicit none

     interface beads_coordinate
           module procedure print_beads_coordinates
     end interface
    
     contains

     subroutine print_beads_coordinates(Beads,sys,pimc)
        implicit none
        
        type(molsysdat), intent(in) :: sys
        type(pimc_par) :: pimc
        type(pimc_particle), dimension(:), pointer :: Beads
        integer :: i,j,k

          do k =1, sys%natom
             print *, 'atom', k
             do i=1, pimc%NumBeadsEff
                print *, Beads(i)%x(:,k)
             enddo
          enddo  
          print *, '----------------------------------' 
         !call exit(1) 
     end subroutine print_beads_coordinates


   end
