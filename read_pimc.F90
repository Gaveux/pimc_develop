
subroutine read_pimc(sys,pimc,in_file)
    ! Example input file: pimc.in
    !
    ! #Specify the core simulation parameters
    ! Number of beads
    ! 6
    ! Temperature
    ! 300
    ! Number of blocks, Blocks to equilibrium, steps per block
    ! 500,10,50000
    ! Run Flyvbjerg blocking algorithm (y/n)
    ! n
    ! Write out geometries at the end of each block (0=no, 1=yes)
    ! 0
    ! Do we use a restart file
    ! n
    ! Enable/Disable writing checkpoint file (y/n)
    ! n
    ! Maximum step size for initial displacements
    ! 0.01
    !
    ! #speciy the trial move parameters
    ! Trial Move Type (0=bead + atom move, 1=staging+atom move)
    ! 1
    ! Atomic displacement parameter
    ! 6.0
    ! Displacement factor for bead moves
    ! 0.5
    ! Number of beads to displace per step
    ! 9
    !
    ! #Specify the action used
    ! Action type (0=Primitive Action, 1=Takahashi-Imada Action, 2=Chin Action, 3=Suzuki-Chin Action)
    ! 2
    ! t0 Parameter for chin action
    ! 0.12
    ! a1 Parameter for chin action
    ! 0.33
    ! alpha Parameter for Suzuki-Chin action
    ! 0.5
    !
    ! #Read in the Free Energy Calculation Parameters - only used if built with FREE_ENERGY set
    ! Compute Free Energy (0 = no, 1=yes)
    ! 0
    ! Compute Normal (0=no,1=yes)
    ! 0
    ! Free Energy Type (TI_Class_Quan = 0, TI_Ref_Pot=1, RS = 2)
    ! 2
    ! lambda (Only used for the TI techniques) (takes values on [0,1]
    ! 0.0
    ! Initial Temperature (used for the reversible scaling simulations)
    ! 100
    ! Final Temperature(used for the reversible scaling simulations)
    ! 300
    !
    !\params[in] sys The molsysdat object storing the system parameters
    !\params[out] pimc The pimc_par object storing the PIMC simulation parameters
    !\params[in] in_file The input file containing the PIMC simulation parameters
    !

    use molecule_specs
    use pimc_structures
    implicit none

    type (molsysdat), intent(in) :: sys
    type (pimc_par), intent(out) :: pimc
    character(len=80), intent(in) :: in_file

    integer :: i,k
    character(len=80) :: icomm

    write(11,*) '----------------------------------------------------'
    write(11,*) ''
    write(11,*) ' read path integral MC parameters from IN_PIMC '
    write(11,*) ''

    open(unit=7,file=trim(in_file),status='old')

    ! read comment; a header for the IN_PIMC file

    read(7,80) icomm
    !write (11,80) icomm
80  format(a80)

    !##########################################################################
    !read in the core simulation parameters
    !##########################################################################
    !  read in the Number of beads

    read(7,80)icomm
    write(11,80) icomm
    read(7,*)  pimc%NumBeads
    write(11,*)  pimc%NumBeads

    read(7,80)icomm
    write(11,80) icomm
    read(7,*) pimc%Temperature
    write(11,*) pimc%Temperature

    !  read in the number of blocks, and number until equilibration
    read(7,80)icomm
    write(11,80) icomm
    read(7,*) pimc%NumBlocks, pimc%BlocksToEquil, pimc%StepsPerBlock
    write(11,*) pimc%NumBlocks, pimc%BlocksToEquil, pimc%StepsPerBlock

    ! run Flyvbjerg blocking algorithm (y/n)
    read(7,80)icomm
    write(11,80)icomm
    read(7,12) pimc%blocking
    write(11,12) pimc%blocking
12  format(a1)

    !  read in the probability of writing a walker to TOUT
    read(7,80)icomm
    write(11,80) icomm
    read(7,*) pimc%Sample
    write(11,*) pimc%Sample

    !  read in flag for starting from qdmc geometries (read in either y or n)
    read(7,80)icomm
    write(11,80) icomm
    read(7,8) pimc%Restart
    write(11,8) pimc%Restart
8   format(a1)
    
    ! Whether or not wirting checkpoint
    read(7,80) icomm
    write(11,80) icomm
    read(7,8) pimc%WritingCheckpoint
    write(11,8) pimc%WritingCheckpoint    
 
    !  read in the maximum step size for initial displacements
    read(7,80)icomm
    write(11,80) icomm
    read(7,*) pimc%IniDisp
    write(11,*) pimc%IniDisp

    !##########################################################################

    !reads blank line
    read(7,80)icomm
    write(11,80) icomm

    !read move parameters header
    read(7,80)icomm
    write(11,80) icomm

    !##########################################################################
    !parameters for the trial moves
    !##########################################################################
    !read in the move parameters
    read(7,80) icomm
    write(11,80) icomm
    read(7,*) pimc%move%move_type
    write(11,*) pimc%move%move_type

    !  read in the atomic displacement parameter (moving whole polymer chains)
    read(7,80)icomm
    write(11,80) icomm
    read(7,*) pimc%move%AtomDisp
    write(11,*) pimc%move%AtomDisp

    !  read in the displacement factor for bead moves
    read(7,80)icomm
    write(11,80) icomm
    read(7,*) pimc%move%BeadDisp
    write(11,*) pimc%move%BeadDisp

    !  read in number of beads to displace per step
    read(7,80)icomm
    write(11,80) icomm
    read(7,*) pimc%move%MovesPerStep
    write(11,*) pimc%move%MovesPerStep

    !##########################################################################

    !reads blank line
    read(7,80)icomm
    write(11,80) icomm

    !read action parameters header
    read(7,80)icomm
    write(11,80) icomm

    !##########################################################################
    !parameters for the action
    !##########################################################################
    !read in the action parameters
    read(7,80) icomm
    write(11,80) icomm
    read(7,*) pimc%act%act_type
    write(11,*) pimc%act%act_type

    read(7,80) icomm
    write(11,80) icomm
    read(7,*) pimc%act%t0
    write(11,*) pimc%act%t0

    read(7,80) icomm
    write(11,80) icomm
    read(7,*) pimc%act%a1
    write(11,*) pimc%act%a1

    read(7,80) icomm
    write(11,80) icomm
    read(7,*) pimc%act%alpha
    write(11,*) pimc%act%alpha

    !##########################################################################

#ifdef FREE_ENERGY
    !reads blank line
    read(7,80)icomm
    write(11,80) icomm

    !read free parameters header
    read(7,80)icomm
    write(11,80) icomm

    !##########################################################################
    !parameters for the free energy calculations
    !##########################################################################
    !read whether to calculate the free energy
    read(7,80)icomm
    write(11,80) icomm

    read(7,*) pimc%doFree
    write(11,*) pimc%doFree

    !read whether to calculate the probability distributions and energies
    read(7,80)icomm
    write(11,80) icomm

    read(7,*) pimc%doSample
    write(11,*) pimc%doSample

    !read in the free energy calculation type
    read(7,80)icomm
    write(11,80) icomm

    read(7,*) pimc%free%free_type
    write(11,*) pimc%free%free_type

    !read in the free energy lambda parameters
    read(7,80)icomm
    write(11,80) icomm

    read(7,*) pimc%free%lambda
    write(11,*) pimc%free%lambda

    !read in the initial temperature for the reversible scaling
    read(7,80)icomm
    write(11,80) icomm

    read(7,*) pimc%free%temp_init
    write(11,*) pimc%free%temp_init

    pimc%Temperature = pimc%free%temp_init

    !read in the final temperature for the reversible scaling
    read(7,80)icomm
    write(11,80) icomm

    read(7,*) pimc%free%temp_fin
    write(11,*) pimc%free%temp_fin

    !compute the amount beta must change each step
    if(pimc%doFree==1 .and. pimc%free%free_type == 2) then
      pimc%Temperature = pimc%free%temp_init
    endif

    pimc%free%beta_step = ((1.0/pimc%free%temp_fin) - (1.0/pimc%free%temp_init) )
    pimc%free%beta_step = pimc%free%beta_step/3.16681520371153d-6
    pimc%free%beta_step = pimc%free%beta_step/dble(((pimc%NumBlocks-pimc%BlocksToEquil)*pimc%StepsPerBlock)-1)
    write(*,*) pimc%free%beta_step
    write(*,*) 1.0/(pimc%free%temp_init*3.16681520371153d-6), 1.0/(pimc%free%temp_fin*3.16681520371153d-6)
    write(*,*) pimc%free%temp_init, pimc%free%temp_fin
    pimc%OldBeta = 1.0/(pimc%Temperature*3.16681520371153d-6)
    !##########################################################################
#endif

    !close the in_pimc file
    close (unit=7)

    !  setting the number of beads to equal the number of seed geometries
    !  times the number of beads per geometry
    if(pimc%act%act_type == 2) then
        if(pimc%act%t0==0) then
            pimc%NumBeadsEff=2*pimc%NumBeads
            stop 'currently broken when t0 = 0'
        else
            pimc%NumBeadsEff=3*pimc%NumBeads
        endif

        !illegal values of the t0 and a1 parameter in the chin action
        if(pimc%act%t0<0 .or. pimc%act%t0 > (0.5*(1.0-1.0/sqrt(3.0)))) then
            stop 't0 value outside of range'
        endif
        if(pimc%act%a1<0 .or. pimc%act%a1 > 1) then
            stop 'a1 value outside of range'
        endif

        !calculate the rest of the chin action parameters from the two input values
        pimc%act%u0=(1.0/12.0)*(1.0-(1.0/(1.0-2.0*pimc%act%t0))+(1.0/(6.0*(1.0-2.0*pimc%act%t0)**3)))
        pimc%act%v1=1.0/(6.0*(1.0-2.0*pimc%act%t0)**2)
        pimc%act%v2=1.0-2.0*pimc%act%v1
        pimc%act%t1=0.5-pimc%act%t0
        pimc%act%t1inv = 1.0/pimc%act%t1
        pimc%act%t0inv = 1.0/pimc%act%t0

! suzuki-chin action:

    elseif(pimc%act%act_type == 3) then
        !illegal values of the alpha parameter in the suzuki-chin action
        if(pimc%act%alpha<0 .or. pimc%act%alpha > 1) then
            stop 'alpha value outside of range'
        endif

! for suzuki-chin action need number of beads to be even

        if(mod(pimc%NumBeads,2)> 0)then
          write(11,*)' ****************************************************'
          write(11,*)' suzuki-chin action, resetting number of beads as even'
          write(11,*)' ****************************************************'
          print *, ' ****************************************************'
          print *, ' suzuki-chin action, resetting number of beads as even'
          print *, ' ****************************************************'
          pimc%NumBeads=pimc%NumBeads+1
          write(11,*)' NumBeads:',pimc%NumBeads
          write(11,*)' ****************************************************'
          print *, ' NumBeads:',pimc%NumBeads
          print *, ' ****************************************************'
        endif
         
        pimc%NumBeadsEff=pimc%NumBeads

! primitive action and takahashi-imada action

    else
        pimc%NumBeadsEff=pimc%NumBeads
    endif

    if(pimc%move%move_type==1 .and. pimc%NumBeadsEff .le. pimc%move%MovesPerStep  ) then
        stop 'Staging more beads than there are beads to move'
    endif

!determine the number of moves that need to be made per monte carlo pass
        pimc%atom_pass = 0
        if(pimc%move%move_type.eq.0) then
            pimc%num_moves=1
            pimc%atom_pass=1
        else if(pimc%move%move_type.eq.1) then
            pimc%num_moves=ceiling(dble(pimc%NumBeadsEff)/dble(pimc%move%MovesPerStep))+1
            pimc%atom_pass=sys%natom
        endif


    pimc%invBeta=3.16681520371153d-6*pimc%Temperature
    pimc%Beta = 1.0/pimc%invBeta
    pimc%invNumBeads = 1.0d0/dble(pimc%NumBeads)

    return
end
