
! prepare the initial collection of beads

subroutine pimc_setup(seedval,sys,pimc,Beads,OldBeads)
    !Description:
    !! @detail
    !> Subroutine pimc_setup generates the initial configuration of the
    !> ring polynmer beads used in the path integral simulations from
    !> the molsysdat and pimc_par objects
    !
    !\param[in] sys The system data object
    !\param[in] pimc The pimc simulation parameters
    !\param[out] Beads The array storing the current value of the beads of the ring polynmer
    !>has size pimc%NumBeads+1 with the final element being a copy of the first element
    !\param[out] OldBeads The array storing the previous value of the beads of the ring polynmer
    !>has size pimc%NumBeads+1 with the final element being a copy of the first element

    use molecule_specs
    use pimc_structures
    !use mt19937
    use prng
    use seed
    implicit none
    
    type(molsysdat), intent(in) :: sys
    type(mod_seed), intent(inout) :: seedval
    type(pimc_par), intent(inout) :: pimc
    type(pimc_particle), dimension(:), pointer :: Beads
    type(pimc_particle), dimension(:), pointer :: OldBeads
 
    integer :: ierr, i,j,k
    real(kind=8) :: rand, total_mass
    real(kind=8), dimension(sys%dimen) :: centre
 
    !------------------------------------------------------------
    !  declare Beads array
    !------------------------------------------------------------
 
    allocate(Beads(pimc%NumBeadsEff+1),stat=ierr)
    if (ierr.ne.0) stop ' Error allocating Beads array '
 
    do i=1,size(Beads)
        call new(Beads(i),sys%natom,sys%dimen)
    enddo
 
    !------------------------------------------------------------
    !  declare OldBeads array
    !------------------------------------------------------------
 
    allocate(OldBeads(pimc%NumBeadsEff+1),stat=ierr)
    if (ierr.ne.0) stop ' Error allocating OldBeads array '
 
    do i=1,size(OldBeads)
        call new(OldBeads(i),sys%natom,sys%dimen)
    enddo
 
    !------------------------------------------------------------
    ! initialise coordinates for each Bead
    !------------------------------------------------------------
 
    if (pimc%Restart == 'n') then
        ! not restarting so generate as random displacements from equil
        do i=1,pimc%NumBeadsEff
            do j=1,sys%natom
                do k=1,sys%dimen
                    rand=genrand_real(seedval%seedvalue)*2.0 - 1.0
                    !rand=genrand_real3()*2.0 - 1.0
                    Beads(i)%x(k,j) = sys%EquilibriumGeom(k,j) + rand*pimc%IniDisp
                enddo
            enddo
        enddo
 
    else if (pimc%Restart == 'y') then
        
        ! flag the error if the file does not already exist
        print *,'reading checkpoint file: ', adjustl(trim(pimc%resume))
        open(unit=599,file=adjustl(trim(pimc%resume)),status='old',action='read')
        ! read in the latest beads configuration in 3D Cartesian coordinate
        read(599,*) 
        do i=1,pimc%NumBeadsEff
           do j=1,sys%natom
              !do k=1,sys%dimen
                 read(599,*) (Beads(i)%x(k,j),k=1,sys%dimen)
                  !Beads(i)%x(k,j) = sys%EquilibriumGeom(k,j) + rand*pimc%IniDisp ! enable for debugging
              !enddo
           enddo
        enddo
        close(unit=599)
        
        !checking dimensionality - enable for debugging only
        !do i=1,pimc%NumBeadsEff
         !  do j=1,sys%natom
          !       print *, (Beads(i)%x(k,j),k=1,sys%dimen)
           !enddo
       ! enddo
       ! call exit(1)
       
    else
        print *,' Error: restart parameter must be "y" or "n" '
        print *,'      : restart is ',pimc%Restart
        call exit(1)
    endif
 
    !------------------------------------------------------------
    ! translate Equilibrium geoms to centre of mass frame
    ! for later use in pimc_output
    !------------------------------------------------------------
 
    centre = 0.0
 
    do j=1,sys%natom
        do k=1,sys%dimen
            centre(k) = centre(k) +  sys%EquilibriumGeom(k,j)*sys%mass(j)
        enddo
    enddo
 
    total_mass = sum(sys%mass)
    centre = centre/total_mass
 
    do j=1,sys%natom
        do k=1,sys%dimen
            sys%EquilibriumGeom(k,j) = sys%EquilibriumGeom(k,j) - centre(k)
        enddo
    enddo
 
    return
end subroutine



