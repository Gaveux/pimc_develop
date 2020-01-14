 

   program coord_conv
     use quaternion_mod
     use pimc_structures
     use coordinate_transformation

     implicit none

     type(pimc_particle),dimension(:), pointer :: Beads
     
  
     integer :: i,j,k,line,ierr,skip, lines
     integer :: move, natom, dimen, NumBeadsEff
     real(kind=8), dimension(3) :: y    
     character(len=100) :: filename
     character(len=3) :: arg2
     integer :: BeadNumber

     call getarg(1,filename)
     call getarg(2,arg2)
     read (arg2,*) BeadNumber
     allocate(Beads(BeadNumber),stat=ierr)
     if (ierr.ne.0) stop 'Error allocating Beads array'

     do i=1,size(Beads)
        call new(Beads(i),5,3)
     enddo
    
     open(unit=7,file=trim(filename),status='unknown',action='read',position='rewind')
       do move = 1, 1000000

         do NumBeadsEff = 1,BeadNumber
            do natom = 1,5
               !read(7,*) Beads(1)%x(:,1) 
               read(7,*) Beads(NumBeadsEff)%x(:,natom)
            enddo
         enddo
         ! do transformation operations here
         call coord_translation(Beads,BeadNumber)
         call coord_rotation(Beads,BeadNumber)
        do NumBeadsEff = 1,BeadNumber
           do natom = 1,5
              print *, Beads(NumBeadsEff)%x(:,natom)
           enddo
        enddo 
       enddo

     
     close(unit=7)
     deallocate(Beads)

   end program

