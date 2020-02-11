module path_integral_monte_carlo
    use Estimator_class
    use molecule_specs
    use pimc_structures
    !use mt19937
    use actions
    use moves
    use binning
    use prng
    use seed
    use blocking
    use annealing_schedule

#ifdef FREE_ENERGY
    use free_energy
#endif

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

    contains

    subroutine print_pimc_header(sys,pimc,num_moves)
        type(molsysdat), intent(in) :: sys
        type(pimc_par), intent(inout) :: pimc
        integer, intent(in) :: num_moves
        integer :: i,j,k, NumBlocksLeft, BlocksToEquilLeft
        type (estimator) :: est
        print *, '==============================================================================='
        print *, '                      PIMC  the temperature effect  - '
        print *, '==============================================================================='
        print *, ' System Dimensionality                       ',sys%dimen
        print *, ' Total number of beads                       ',pimc%NumBeads 
        print *, ' Number of Effective Beads                   ',pimc%NumBeadsEff
        print *, ' Step size for initial displacement          ',pimc%IniDisp

        if (pimc%Restart == 'n') then
           print * ,' Number of blocks                            ',pimc%NumBlocks
           print *, ' Number of blocks to equilibrium             ',pimc%BlocksToEquil
        else if (pimc%Restart == 'y') then
             open(unit=599,file=adjustl(trim(pimc%resume)),status='old',action='read')
               read(599,*) ! skip the seedvalue line
               do j=1,pimc%NumBeadsEff
                  do k=1,sys%natom
                   read(599,*) !skip beads configurations
                  enddo
               enddo
               read(599,*) est
               read(599,*) ! skip the block number line
               !read(599,*) pimc%NumBlocks, pimc%BlocksToEquil
               read(599,*) NumBlocksLeft, BlocksToEquilLeft
            close(unit=599)
            
            !if NumBlocks = 0 then take the value from pimc.in
            if (NumBlocksLeft .NE. 0) then
               pimc%NumBlocks = NumBlocksLeft
               pimc%BlocksToEquil = BlocksToEquilLeft
            endif
            print *, ' Number of blocks                            ',pimc%NumBlocks 
            print *, ' Number of blocks to equilibrium             ',pimc%BlocksToEquil
        endif

        print *, ' Number of Monte Carlo steps per block       ',pimc%StepsPerBlock
        print *, ' Temperature in Kelvin                       ',pimc%Temperature
        print *, ' Trial Moves type                            ',pimc%move%move_type
        print *, ' Atomic displacement parameter (MC move)     ',pimc%move%AtomDisp
        print *, ' Displacement factor for moving beads        ',pimc%move%BeadDisp
        print *, ' Number of beads displaced per step (staging)',pimc%move%MovesPerStep
        print *, ' Moves Per Monte Carlo Pass =                ',num_moves
        print *, ' Action type                                 ',pimc%act%act_type
        print *, ' t0 =                                        ',pimc%act%t0
        print *, ' a1 =                                        ',pimc%act%a1
        print *, ' alpha =                                     ',pimc%act%alpha

#ifdef FREE_ENERGY
        print *, ' Compute Free Energy =                   ',pimc%doFree
        if(pimc%doFree==1) then
            if(pimc%free%free_type == 2) then
                print *, ' Free Energy Type =                      ','Reversible Scaling'
                print *, ' Initial Temperature =                   ',pimc%free%temp_init
                print *, ' Final Temperature =                     ',pimc%free%temp_fin
            elseif(pimc%free%free_type == 1) then
                print *, ' Free Energy Type =                      ','Thermodynamic Integration : ref Potential'
                print *, ' Currently not implemented'
                stop ' Exiting'
            elseif(pimc%free%free_type == 0) then
                print *, ' Free Energy Type =                      ','Thermodynamic Integration : Classical to Quantum'
                print *, ' Currently not implemented'
                stop ' Exiting'
            endif
        endif
#endif
        
            print *, '--------------------------------------------------------------------'
        if (pimc%Restart == 'y'.and.pimc%blocking == 'y') then
            print *, '     Resume numerical convergence test using blocking algorithm     ' 
        elseif(pimc%Restart =='y'.and.pimc%blocking == 'n' ) then
            print *, '              Resuming path integral MC calculation                 '
        elseif(pimc%Restart =='n'.and.pimc%blocking == 'n')then
            print *, '              Starting path integral MC calculation                 '
        elseif(pimc%Restart == 'n'.and. pimc%blocking == 'y') then 
            print *, '     Start numerical convergence test using blocking algorithm      '
        endif
            print *, '--------------------------------------------------------------------'
    end subroutine print_pimc_header


! main loop for path integral Monte Carlo
#if POT == 0
    subroutine pimc_monte(seedval,sys,pimc,Beads,OldBeads, pot,out_dir, checkpoint_dir, binning_in)
#else
    subroutine pimc_monte(seedval,sys,pimc,Beads,OldBeads, out_dir, checkpoint_dir, binning_in)
#endif
        type (molsysdat), intent(in)  :: sys
        type (mod_seed), intent(inout) :: seedval
        type (pimc_par) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        type (pimc_particle), dimension(:), pointer :: OldBeads
#if POT == 0
        type(msi_params) :: pot
#endif
        character(len=80), intent(in) :: out_dir
        character(len=80), intent(in) :: checkpoint_dir
        character(len=80), intent(in) :: binning_in
        character(len=80) :: IN_ISEED
#ifdef FREE_ENERGY
        type (free_energy_type) :: free
#endif
        type (estimator) :: est
        type (action_vals) :: act
        type(bin_type) :: bins
        type(binning_params) :: params
    
        integer :: iblock, iter, ibead
        real(kind=8) :: accept, acctot, moveacc, moveacctot
        real(kind=8) deltae, rand, B_init
        integer :: i, ind
        integer :: imove,iatom,num_moves,atom_pass
        integer :: n_bl, n_ba, n_d
        logical :: equil = .TRUE.
        logical :: first_time = .TRUE.
        logical :: atom_move
        integer :: first_moved,last_moved
        integer :: j,k, NumBlocksLeft, BlocksToEquilLeft
        real(kind=8) ::AtomDisp

        atom_pass=0
        !determine the number of moves that need to be made per monte carlo pass
        if(pimc%move%move_type.eq.0) then
            num_moves=1
            atom_pass=1
        else if(pimc%move%move_type.eq.1) then
            num_moves=ceiling(dble(pimc%NumBeadsEff)/dble(pimc%move%MovesPerStep))+1
            atom_pass=sys%natom
        endif
    
        !print the pimc parameters out
        call print_pimc_header(sys,pimc,num_moves)
        

        if (pimc%Restart == 'n') then
        
             call new(est)

        else if (pimc%Restart == 'y') then
          ! need a flag for just in case the previous job did not write a checkpoint file 
             open(unit=599,file=adjustl(trim(pimc%resume)),status='old',action='read')
             read(599,*) ! skip the seedvalue line
             do j=1,pimc%NumBeadsEff
                do k=1,sys%natom
                   read(599,*) !skip beads configurations
                enddo 
             enddo
             read(599,*) est 
             read(599,*) ! skip the block number line
             read(599,*) NumBlocksLeft, BlocksToEquilLeft
             close(unit=599)
             
             if (NumBlocksLeft .NE. 0) then
                pimc%NumBlocks = NumBlocksLeft
                pimc%BlocksToEquil = BlocksToEquilLeft
             endif
             !print *, pimc%NumBlocks, pimc%BlocksToEquil
        endif
         
        !initialise the action type
        act%act_old=0.0
        act%act=0.0

#ifdef FREE_ENERGY
        if(pimc%doFree==1) then
            !initialise the free energy type
            call new(free,pimc,sys%natom,sys%dimen)
            if(pimc%free%free_type==2) then
                call openFileOut(out_dir)
            endif
        endif
#endif
    
        !initialise the types for binning parameters of the system
        call read_number_hist(binning_in, sys%nbond, n_bl, n_ba, n_d)
        call new(params,n_bl,n_ba,n_d)
        call read_binning_params(params,binning_in,sys%nbond)
        call new(bins,params)
        
        !evaluate the action to start
        call eval_action(sys,pimc,Beads,act)

        acctot=0.0
        moveacctot=0.0
 
        !Start the main monte carlo loop
        do iblock=1,pimc%NumBlocks
            accept = 0.0
            moveacc = 0.0
            if (iblock.eq.pimc%BlocksToEquil+1) then
                equil=.FALSE.
            endif


            B_init = pimc%Beta
            do iter=1,pimc%StepsPerBlock
           
                do iatom=1,atom_pass
                    do imove=1,num_moves
                        atom_move=.false.
                        if(imove.eq.num_moves) then
                            atom_move=.true.
                        endif

              ! introduce a larger than bondlength atomic displacement every
              ! given number of blocks till equilibration complete
              if (atom_move .and. equil) then
                  call annealing_condition(pimc,AtomDisp,iblock,iter)
              endif
              !pimc%move%AtomDisp = AtomDisp
              !print *,  pimc%move%AtomDisp



                        do ibead=1,pimc%NumBeadsEff+1
                            call copy(Beads(ibead),OldBeads(ibead))
                        enddo
                        act%act_old = act%act
#ifdef FREE_ENERGY
                        !scale the coordinates for thermodynamic integration between the classical andquantum result
                        if (pimc%free%free_type == 0 .and. pimc%doFree==1) then
                            call scaleCoords(free,Beads,pimc,sys)
                        endif
#endif                  
                        call trial_move(seedval,sys,pimc,Beads,atom_move,first_moved,last_moved)
#ifdef FREE_ENERGY
                        !If using classical to quantum scaling all the beads have to be updated
                        if (pimc%free%free_type == 0.and. pimc%doFree==1) then
                            first_moved=1
                            last_moved = pimc%NumBeadsEff
                        endif
#endif
                        do i=first_moved,last_moved ! see the definition in staging_move subroutine 
                            ind = mod(i-1,pimc%NumBeadsEff)+1
#if POT == 0
                            !MSI potential energy surfaces
                            call potential(ind,pot,Beads(ind)%x,Beads(ind)%r,Beads(ind)%VCurr,Beads(ind)%dVdx)
#else
                            !Analytic potential energy surfaces
                            call potential(sys,Beads(ind)%x,Beads(ind)%r,Beads(ind)%VCurr,Beads(ind)%dVdx)
#endif
                        enddo
#ifdef FREE_ENERGY
                        !Unscale the coordinates so the rest proceeds as normal
                        if (pimc%free%free_type == 0.and. pimc%doFree==1) then
                            call unscaleCoords(free,Beads,pimc,sys)
                        endif
#endif
                
                        call copy (Beads(1),Beads(pimc%NumBeadsEff+1))

                        ! calculate the path integral action
                        call eval_action(sys,pimc,Beads,act)

                        ! now the actual Monte Carlo step!
#ifdef FREE_ENERGY
                        if(pimc%free%free_type==2.and. pimc%doFree==1) then
                            deltae = pimc%Beta*act%act - pimc%OldBeta*act%act_old
                        else
                            deltae = pimc%Beta*(act%act - act%act_old)
                        endif
#else
                        deltae = pimc%Beta*(act%act - act%act_old)
#endif
                        !rand = genrand_real3()
                        rand = genrand_real(seedval%seedvalue)

                        if (deltae.lt.0.0) then
                            if(atom_move) then
                                accept = accept + 1.0
                            else
                                moveacc = moveacc + 1.0 
                            endif
                        else if ((deltae.gt.0.0).and.(exp(-deltae).lt.rand)) then
                            do ibead=1,pimc%NumBeadsEff+1
                                call copy(OldBeads(ibead),Beads(ibead))
                            enddo
                            act%act = act%act_old
                        else
                            if(atom_move) then
                                accept = accept + 1.0
                            else
                                moveacc = moveacc + 1.0 
                            endif
                        endif

                        ! do some calculating of energy estimators and binning
                        if(equil) then
                        else
#ifdef FREE_ENERGY
                            if(pimc%doSample==1) then
#endif
                                do i=1,pimc%NumBeads
                                    if(pimc%NumBeads*3==pimc%NumBeadsEff) then
                                        call updateBinning(bins,params,Beads(3*(i-1)+1)%x,Beads(i)%r)
                                    else
                                        call updateBinning(bins,params,Beads(i)%x,Beads(i)%r)
                                    endif
                                enddo
                              
                               ! call update energy(sys,pimc,Beads,results) see "module estimator_class" 
                                call update(sys,pimc,Beads,est)
#ifdef FREE_ENERGY
                            endif

                            if(pimc%free%free_type==2.and. pimc%doFree==1) then
                                call reversibleScaling(free,Beads,pimc,sys)
                            endif
                            
#endif
                        endif
                    enddo
                enddo
#ifdef FREE_ENERGY
                if(equil) then
                else
                    !If reversible scaling update the temperature at the end of each step
                    if(pimc%free%free_type==2.and. pimc%doFree==1) then
                        call reversibleScalingStep(free,pimc,iter,dble(num_moves*atom_pass))
                        call tempStep(pimc)
                    endif
                endif
#endif      
            enddo

            !End of the step
            if (pimc%WritingCheckpoint =='y'.and.iblock.lt.pimc%BlocksToEquil+1) then 
               open(unit=599,file=trim(checkpoint_dir)//trim(pimc%start),status='unknown',action='write',position='rewind')
               ! save the seed value at the end of each block
               write(599,*) seedval%seedvalue
               ! Save the beads configuration at the end of each block
               do i=1,pimc%NumBeadsEff
                  call writeCheckpoint(Beads(i)%x, checkpoint_dir, pimc%start, pimc%WritingCheckpoint, .False.,sys%natom,sys%dimen)
               enddo
               !close(unit=599)

            endif

            !Process the results at the end of the block
            accept = accept / (pimc%StepsPerBlock*atom_pass)
            acctot = acctot + accept
     
            if(pimc%move%move_type.eq.1) then
                moveacc = moveacc / (pimc%StepsPerBlock*(num_moves-1)*sys%natom)
                moveacctot = moveacctot + moveacc
            endif
            write(*,*) 'Block: ', iblock, 'Acceptance Ratio: ', accept
         
            if(pimc%move%move_type.eq.1) then
                write(*,*) 'Block: ', iblock, 'Staging Acceptance Ratio: ', moveacc
            endif

            if(pimc%Sample==1) then
                !at the end of each block write some geometries to tout
                do i=1,pimc%NumBeadsEff
                    call writeTOUT(Beads(i)%x,Beads(i)%VCurr, out_dir,.False.,sys%natom,sys%dimen)
                enddo
            endif
            
            pimc%NumBlocksLeft = pimc%NumBlocks - iblock
            pimc%BlocksToEquilLeft = pimc%BlocksToEquil - iblock
 
            if(equil) then
                open(unit=599,file=trim(checkpoint_dir)//trim(pimc%start),status='unknown',action='write',position='append')
                write(599,*) est
                write(599,*) 'block number: ', iblock
                write(599,*) pimc%NumBlocksLeft, pimc%BlocksToEquilLeft 
                !close(unit=599)

            else
#ifdef FREE_ENERGY
                if(pimc%free%free_type==2.and. pimc%doFree==1) then
                    call writeIntInternal(free,pimc,B_init)
                    write(*,*) "Block: ", iblock, 'Final Temperature', pimc%Temperature
                endif
                if(pimc%doSample==1) then
#endif
                call update_block(pimc,est)
                
                if (pimc%WritingCheckpoint =='y'.and.iblock.ge.pimc%BlocksToEquil+1) then
                   ! writing checkpoint for energy estimator and errors
                   open(unit=599,file=trim(checkpoint_dir)//trim(pimc%start),status='unknown',action='write',position='rewind')
                   write(599,*) seedval%seedvalue
                   do i=1, pimc%NumBeadsEff
                      call writeCheckpoint(Beads(i)%x, checkpoint_dir, pimc%start, pimc%WritingCheckpoint, .False.,sys%natom,sys%dimen)
                   enddo
                   
                   write(599,*) est
                   write(599,*) 'block: ', iblock
                   write(599,*) pimc%NumBlocksLeft, '0'   
                endif

#ifdef FREE_ENERGY
                endif
#endif
            endif
            
            !write(*,*)  'End of the block ', iblock, ' seed value', seedval%seedvalue 
        enddo
#ifdef FREE_ENERGY
        if(pimc%doSample==1) then
#endif
        call writeHistogram(bins,params,out_dir) 
#ifdef FREE_ENERGY
        endif
#endif
        !close tout
        if(pimc%Sample==1) then
            call writeTOUT(Beads(1)%x,Beads(1)%VCurr,out_dir,.True.,sys%natom,sys%dimen)
        endif
        ! close checkpoint
        call close_file(599)

        if(pimc%blocking=='y') then
          call blk_count(pimc%blk)
        endif
        ! calculate total average energy and acceptance probability
        acctot=acctot/dble(pimc%NumBlocks)
        moveacctot=moveacctot/dble(pimc%NumBlocks)
#ifdef FREE_ENERGY
        if(pimc%doSample==1) then
#endif
        call end_sim(est)
#ifdef FREE_ENERGY
        endif

        if(pimc%free%free_type==2.and. pimc%doFree==1) then
            call closeFileOut()
        endif
#endif
        return
    end subroutine

end module path_integral_monte_carlo
