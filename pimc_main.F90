
program pimc90  

    !> @file
    !Description:
    !!@details
    !>program pimc90: The core path integral monte carlo program. Sets up the core data structures
    !>                from input files specified as command line arguments and runs the core monte
    !>                carlo subroutine
    use seed
    use molecule_specs
    use pimc_structures
    use path_integral_monte_carlo
    use prng
    !use file_names

#if POT == 0
    use potential_msi
#elif POT == 1
    use potential_h2o
#elif POT == 2
    use potential_nh3
#elif POT == 3
    use potential_hcn
#endif

    implicit none
    
    type (mod_seed) :: seedval
    type (molsysdat) :: sys
    type (pimc_par) :: pimc
    type (pimc_particle), dimension(:), pointer :: Beads
    type (pimc_particle), dimension(:), pointer :: OldBeads

    real :: t1,t2,hours,minutes,seconds
    character(len=80) :: sfilename, rfilename, bfilename
    integer :: i,iargs
 
    character(len=80) :: OUT_DIRECTORY, IN_PIMC, IN_SYSTEM, IN_ISEED, IN_BINNING
    character(len=80) :: CHECKPOINT_DIRECTORY
#if POT == 0
    type (msi_params) :: pot
    character(len=80) :: POT_FILE, IN_INTERP, IN_ATOMPERM
    integer :: ini_MCstep = 0
    integer :: iatom, imove
    logical :: inner_update = .FALSE.
    integer, dimension(:,:,:), allocatable :: old_neigh
    integer :: ierr
#endif

    include 'pimc_setup.int'

    !Get the number of command line arguments
    iargs = IARGC()
  
    !If the number of command line arguments specified is not the number required
    !quit execution and print an error message

#if POT == 0
    if(iargs.ne.12) then
        write(*,*) 'Error: incorrect number of command line arguments'
        write(*,*) iargs
        write(*,*) 'Correct Usage:'
        write(*,*) '    pimc90 pimc.in system.in iseed.in binning.in interp.in pot.in atomperms.in output_dir checkpoint_dir checkpoint_output_filename checkpoint_input_filename blocking_filename'
        stop
    endif
#else

    if(iargs.ne.9) then
        write(*,*) 'Error: incorrect number of command line arguments'
        write(*,*) iargs
        write(*,*) 'Correct Usage:'
        write(*,*) '    pimc90 pimc.in system.in iseed.in binning.in output_dir checkpoint_dir checkpoint_output_filename checkpoint_input_filename blocking_filename'
        stop
    endif
#endif

    !read the command line arguments in

    call getarg(1,IN_PIMC)
    call getarg(2,IN_SYSTEM)
    call getarg(3,IN_ISEED)
    call getarg(4,IN_BINNING)

#if POT == 0
    call getarg(5,IN_INTERP)
    call getarg(6,POT_FILE)
    call getarg(7,IN_ATOMPERM)
    call getarg(8,OUT_DIRECTORY)
    call getarg(9,CHECKPOINT_DIRECTORY)
    call getarg(10,SFILENAME)
    call getarg(11,RFILENAME)
    call getarg(12,BFILENAME)
#else
    call getarg(5,OUT_DIRECTORY)
    call getarg(6,CHECKPOINT_DIRECTORY)
    call getarg(7,SFILENAME)
    call getarg(8,RFILENAME)
    call getarg(9,BFILENAME)
#endif
    !print which type of potential energy surface is being used

#ifdef FREE_ENERGY
    write(*,*) "Free Energy Simulation Mode is On"
#endif

#if POT == 0
    write(*,*) "Modified Shepard Potential" 
#elif POT == 1
    write(*,*) "PJT2 Water Potential" 
#elif POT == 2
    write(*,*) "AMMPOT4 Ammonia Potential" 
#elif POT == 3
    write(*,*) "Murrell, Carter, Halonen Hydrogen Cyanide Potential" 
#endif

    !---------------------------------------------------------------
    ! Start timing CPU time used
    !---------------------------------------------------------------
  
    call cpu_time(t1)
    
    !make the output directory
    call system('mkdir -p '//trim(OUT_DIRECTORY))
    call system('mkdir -p '//trim(CHECKPOINT_DIRECTORY))
    
    !read the input files
    call read_system_data(sys,IN_SYSTEM)
    call read_pimc(sys,pimc,IN_PIMC)
    pimc%resume=rfilename
    print *, 'Write a checkpoint file for the current job? ', pimc%WritingCheckpoint
    call read_iseed(seedval,pimc,IN_ISEED)
    call pimc_setup(seedval,sys,pimc,Beads,OldBeads)
    pimc%blk=bfilename
 
#if POT == 0
    !Initialise the modified shepard potential energy surface
    call MSI_INIT(pot,sys,IN_INTERP,POT_FILE,IN_ATOMPERM,pimc%numBeadsEff)
    !print *, pot%interp%ndata
    !call neigh_copy(pot,pimc,old_neigh)
    allocate(old_neigh(pimc%atom_pass,pimc%NumBeadsEff,pot%interp%ndata),stat=ierr)
    if (ierr.ne.0) stop 'Error allocating old_neigh array in pimc_main.F90'
    old_neigh = 0
#endif

    ! calculate potential energies for the initial geometry
    iatom = 0
    imove = 0
    do i=1,pimc%NumBeadsEff
#if POT == 0
      !MSI potential energy surfaces
      call potential(i,pot,Beads(i)%x,Beads(i)%r,Beads(i)%VCurr,Beads(i)%dVdx, ini_MCstep,pimc,iatom,imove,inner_update,old_neigh)
#else
      !Analytic potential energy surfaces
      call potential(sys,Beads(i)%x,Beads(i)%r,Beads(i)%VCurr,Beads(i)%dVdx)
#endif
    enddo
    call exit(0)
 
    ! copy the initial geometry to the end of the array (NumGeoms + 1)
    ! to make summing the action over the whole closed system easier
 
    call copy(Beads(1),Beads(pimc%NumBeads+1))
   
    
    ! main routine - move beads, evaluate action and monte carlo step
    !Create file target for writing variables at end of simulation
    pimc%start=sfilename
#if POT == 0
    call pimc_monte(seedval,sys,pimc,Beads,OldBeads,pot,OUT_DIRECTORY,CHECKPOINT_DIRECTORY,IN_BINNING,old_neigh)
    deallocate(old_neigh,stat=ierr)
    if (ierr.ne.0) stop 'Error deallocating old_neigh in pimc_main.F90'
   
#else 
    call pimc_monte(seedval,sys,pimc,Beads,OldBeads,OUT_DIRECTORY,CHECKPOINT_DIRECTORY,IN_BINNING)
#endif
    ! stop timer and print CPU time used
    call cpu_time(t2)
    print *,' '
    hours = (t2-t1)/3600.0
    minutes = (hours - int(hours))*60.0
    seconds =  (minutes - int(minutes))*60.0
    print *,'CPU time: ',int(hours),&
    &       'hour(s)',int(minutes),'minute(s)',int(seconds),'second(s)'

end
