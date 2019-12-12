
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
     
     subroutine print_Beads_distance(ind,Beads,sys,pimc)

         implicit none
         
         type(molsysdat), intent(in) :: sys
         type(pimc_par) :: pimc
         type(pimc_particle), dimension(:), pointer :: Beads
         !real(kind=8), dimension(sys%dimen) :: norm
         integer :: i,ind,natom,j,k
         real(kind=8) :: norm
         norm = 0.d0
         
         do k=1,sys%natom-1
           do i=1,sys%dimen
             norm = norm + (Beads(ind)%x(i,k+1) - Beads(ind)%x(i,k))**2
           enddo
         norm = sqrt(norm)
         print *, norm
         enddo
         !print *, 'distance between bead',ind,'=', norm, 'for atom ', natom
     end subroutine print_beads_distance

   end
