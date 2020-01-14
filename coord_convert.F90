 

   program coord_conv
     use quaternion_mod
     use pimc_structures
     use coordinate_transformation

     implicit none

     type(pimc_particle),dimension(:), pointer :: Beads, BeadSum
     
  
     integer :: i,j,k,line,ierr,skip, lines
     integer :: move, natom, dimen, NumBeadsEff
     real(kind=8), dimension(3) :: y    
     character(len=100) :: filename, atom1, atom2, atom3, atom4, atom5
     character(len=3) :: arg2
     integer :: BeadNumber
     integer :: a1,a2,a3,a4,a5

     call getarg(1,filename)
     call getarg(2,arg2)
     call getarg(3,atom1)
     call getarg(4,atom2)
     call getarg(5,atom3)
     call getarg(6,atom4)
     call getarg(7,atom5)
     read (arg2,*) BeadNumber
     allocate(Beads(BeadNumber),stat=ierr)
     if (ierr.ne.0) stop 'Error allocating Beads array'
     allocate(BeadSum(BeadNumber),stat=ierr)
     if (ierr.ne.0) stop 'Error allocating BeadSum array'

     do i=1,size(Beads)
        call new(Beads(i),5,3)
        call new(BeadSum(i),5,3)
     enddo
    
     open(unit=7,file=trim(filename),status='unknown',action='read',position='rewind')
     open(newunit=a1,file=trim(atom1),status='unknown',action='write')
     open(newunit=a2,file=trim(atom2),status='unknown',action='write')
     open(newunit=a3,file=trim(atom3),status='unknown',action='write')
     open(newunit=a4,file=trim(atom4),status='unknown',action='write')
     open(newunit=a5,file=trim(atom5),status='unknown',action='write')
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
        !do NumBeadsEff = 1,BeadNumber
        !   do natom = 1,5
        !      print *, Beads(NumBeadsEff)%x(:,natom)
        !   enddo
        !enddo 

        ! do moving average of Beads coordinates
        do NumBeadsEff =1,BeadNumber
           do natom = 1,5
             BeadSum(NumBeadsEff)%x(:,natom) = BeadSum(NumBeadsEff)%x(:,natom) + Beads(NumBeadsEff)%x(:,natom) 
           enddo
           write(a1,*) BeadSum(NumBeadsEff)%x(:,1)/dble(move)
           write(a2,*) BeadSum(NumBeadsEff)%x(:,2)/dble(move)
           write(a3,*) BeadSum(NumBeadsEff)%x(:,3)/dble(move)
           write(a4,*) BeadSum(NumBeadsEff)%x(:,4)/dble(move)
           write(a5,*) BeadSum(NumBeadsEff)%x(:,5)/dble(move)
        enddo
       enddo

     
     close(unit=7)
     deallocate(Beads)
     deallocate(BeadSum)

   end program

