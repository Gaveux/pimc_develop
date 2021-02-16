module moves
    use molecule_specs
    use pimc_structures
    !use mt19937
    use prng
    use seed
  
    implicit none
    

    interface move
        module procedure trial_move
    end interface move

    contains 


    !-------------------------------------------
    ! gasdev from numerical recipes (Page 280)
    ! returns a random number from the 
    ! distribution
    !
    !     p(y)dy = exp{-0.5*y^2}/sqrt(2*pi)
    !
    !-------------------------------------------
    function gauss_dev(seedval) result(harvest)
         
        implicit none
        type (mod_seed) :: seedval
        real(kind=8) :: harvest
        real(kind=8) :: rsq,v1,v2
        real(kind=8), save :: g
        logical, save :: gauss_stored = .false.
 
        if (gauss_stored) then
            harvest = g
            gauss_stored = .false.
        else
            do
                !v1 = genrand_real1()
                v1 = genrand_real(seedval%seedvalue)
                !v2 = genrand_real1()
                v2 = genrand_real(seedval%seedvalue)
                v1 = 2.0*v1-1.0
                v2 = 2.0*v2-1.0
                rsq = v1**2 + v2**2
               if (rsq > 0.0 .and. rsq < 1.0) exit
            enddo

            rsq = sqrt(-2.0*log(rsq)/rsq)
            harvest = v1*rsq
            g = v2*rsq
            gauss_stored = .true.

        endif
    end function gauss_dev

    !general trial move function selects which type of move is to be
    !made given the pimc parameters

    subroutine trial_move(seedval,sys,pimc,Beads,atom_move, start, last)
    !subroutine trial_move(sys,pimc,Beads,atom_move, start, last)
        type (molsysdat), intent(in)  :: sys
        type (mod_seed), intent(inout) :: seedval
        type (pimc_par) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        logical :: atom_move
        integer, intent(out) :: start,last
        
        ! Trial move type 0 = beads + atom move
        ! type 1 = staging + atom move
 
        if(pimc%move%move_type.eq.0) then
            call move_beads(seedval,sys,pimc,Beads)
            !call move_beads(sys,pimc,Beads)
            start=1
            last=pimc%NumBeadsEff
 
        else if (pimc%move%move_type.eq.1) then
            if(pimc%act%act_type.eq.0 .or. pimc%act%act_type.eq.1 .or. pimc%act%act_type.eq.3) then
                call staging_move(seedval,sys,pimc,Beads,atom_move,start,last)
                !call staging_move(sys,pimc,Beads,atom_move,start,last)
 
            !If the chin action call the chin action staging move
            else if(pimc%act%act_type.eq.2) then
                call staging_move_ca(seedval,sys,pimc,Beads,atom_move,start,last)
                !call staging_move_ca(sys,pimc,Beads,atom_move,start,last)
            endif
        endif

    end subroutine trial_move

    !  randomly choose which bead/beads to move, then move the
    !  bead, using the Box-Muller method to generate random numbers
    !  with a spherical Gaussian distribution (as per walking in qdmc)
    !  only recalculate the potential energy for the beads that have
    !  been moved each step
    ! doi: 10.1119/1.18168 for gaussian distribution 
    !subroutine move_beads(sys,pimc,Beads)
    subroutine move_beads(seedval,sys,pimc,Beads)
        
        type (mod_seed), intent(inout) :: seedval
        type (molsysdat), intent(in)  :: sys
        type (pimc_par) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads

        integer :: i,n,j,mvs
        real(kind=8) :: sigma,rho
        real(kind=8) :: rand
       
        !----------------------------------------------------------------
        ! start by moving all beads a large amount
        !----------------------------------------------------------------
       
        do n=1,sys%natom
            do j=1,sys%dimen
                rho = gauss_dev(seedval)
                sigma = sqrt(pimc%move%AtomDisp/sys%mass(n))
                do i=1,pimc%NumBeadsEff
                    Beads(i)%x(j,n) = Beads(i)%x(j,n) + sigma*rho
                enddo
            enddo
        enddo

        !----------------------------------------------------------------
        ! start small displacements by choosing which bead to move
        !----------------------------------------------------------------

        do mvs=1,pimc%move%MovesPerStep

            !rand=genrand_real1()
            rand = genrand_real(seedval%seedvalue)
            i = 1 + int(rand*pimc%NumBeadsEff)
 
            !--------------------------------------------------------------
            ! then move it...
            !--------------------------------------------------------------

            do n=1,sys%natom
                do j=1,sys%dimen

                    rho   = gauss_dev(seedval)
                    sigma = sqrt(pimc%move%BeadDisp/sys%mass(n))
                    Beads(i)%x(j,n) = Beads(i)%x(j,n) + sigma*rho

                enddo
            enddo

        enddo

        return
    end subroutine


    ! Staging move algorithm for N particles.  Randomly selects
    ! an atom from the list and applies a monte carlo pass to
    ! it.  This pass consists of P/j staging moves where P are
    ! the number of beads and j are the number of beads moved
    ! followed by a move of the entire ring of the atom
    !subroutine staging_move(sys,pimc,Beads,atom_move,start,last)
    subroutine staging_move(seedval,sys,pimc,Beads,atom_move,start,last)
        type (molsysdat), intent(in)  :: sys
        type (mod_seed), intent(inout) :: seedval
        type (pimc_par) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        logical, intent(in) :: atom_move
        integer, intent(out) :: start,last
        
        integer :: l, j,n,k, atom, dimen, i
        real(kind=8) :: sigma,rho
        real(kind=8) :: rand
        real(kind=8) :: nextPos, firstPos!,x,y,z
       
        !Run the staging move if told to by the monte carlo subroutine
        if(atom_move.eqv..false.) then
            
            !Get index of pseudorandom atom to move
            !rand=genrand_real1()
            rand=genrand_real(seedval%seedvalue)
            atom = 1 + int(rand*sys%natom)

            
            !get index of first fixed end point bead, the other is l+pimc%NumBeads+1
            !rand=genrand_real1()
            rand=genrand_real(seedval%seedvalue)
            l = 1 + int(rand*pimc%NumBeads)
       
            !Need to sample the gaussian distribution for each bead here start 
            !from largest index bead then go down to smallest index bead so that
            !The recursive formula is accounted for calculates the new bead position
            !by transforming from the staging coordinates randomly sampled from the  gaussian with
            !variance 
       
            !mod(l-1,pimc%NumBeads)+1 should return l provided it is less than or equal 
            !to pimc%NumBeads and should wrap around correctly if l
            !is greater than pimc%NumBeads
       
            do i=0,pimc%move%MovesPerStep-1
                do dimen=1,sys%dimen
                    k=pimc%move%MovesPerStep-i
                    rho=gauss_dev(seedval)
                    sigma=sqrt((dble(k)*pimc%Beta)/(dble(k+1)*pimc%NumBeads*sys%mass(atom)))
                    nextPos = Beads(mod(l+k,pimc%NumBeads)+1)%x(dimen,atom)
                    firstPos = Beads(l)%x(dimen,atom)
                    Beads(mod(l-1+k,pimc%NumBeads)+1)%x(dimen,atom)=sigma*rho+dble(k)*nextPos/dble(k+1)+firstPos/dble(k+1)
               enddo
            enddo 
        
            start = l
            last=pimc%move%MovesPerStep+l

        !Otherwise move the entire atom (and hence all beads) by an amount
        else

            do n=1,sys%natom
                do j=1,sys%dimen
                    rho = gauss_dev(seedval)
                    sigma = sqrt(pimc%move%AtomDisp/sys%mass(n))
 
                    do i=1,pimc%NumBeads
                        Beads(i)%x(j,n) = Beads(i)%x(j,n) + sigma*rho
                    enddo
                enddo
            enddo

            start=1
            last=pimc%NumBeads

        endif

        return
    end subroutine



    !==============================================================================
    !The following subroutines are the chin action version of the above trial moves
    ! Staging move algorithm for N particles.  Randomly selects
    ! an atom from the list and applies a monte carlo pass to
    ! it.  This pass consists of P/j staging moves where P are
    ! the number of beads and j are the number of beads moved
    ! followed by a move of the entire ring of the atom
   
    !subroutine staging_move_ca(sys,pimc,Beads,atom_move,start,last)
    subroutine staging_move_ca(seedval,sys,pimc,Beads,atom_move,start,last)
        type (molsysdat), intent(in)  :: sys
        type (mod_seed), intent(inout) :: seedval
        type (pimc_par) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        logical, intent(in) :: atom_move
        real(kind=8), dimension(pimc%NumBeadsEff) :: q_vals, b_vals, a_vals
        integer, intent(out) :: start,last

        integer :: i,n,j, start_bead, atom,bead_ind
        real(kind=8) :: sigma,rho
        real(kind=8) :: rand
  
        !A1 corresponds to the alpha+1 value and A to the alpha value
        real(kind=8) :: qPrev, aPrev, bPrev, cPrev, qCurr, bCurr, aCurr, cCurr, cStart
        real(kind=8) :: rStart, rNext, omega
       
        !If not atom move stage the beads
        if(atom_move.eqv..false.)then
            !----------------------------------------------------------------
            ! start by calculating the staging coefficients for each bead
            ! The recursive formula for these staging coefficients can be
            ! found in my Workbook page for april 29th - in addition a
            ! derivation of this expression is included
            !
            ! Rather than recursively computing the values of the coefficients
            ! when they are needed they are precomputed in a single pass and 
            ! stored in the arrays q_vals, b_vals, and a_vals 
            !----------------------------------------------------------------
        
            !Initialise the endpoint q,b,a values to zero
            qPrev=0.0
            bPrev=0.0
            aPrev=0.0
            
            !Choose the random starting bead
            !rand=genrand_real1()
            rand=genrand_real(seedval%seedvalue)
            
            !The first endpoint (not moved)
            start_bead=1+int(rand*pimc%NumBeadsEff)
        
            !Choose the atom to be staged
            !rand=genrand_real1()
            rand=genrand_real(seedval%seedvalue)
            atom=1+int(rand*sys%natom)
           
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !This chunk is most likely to have an error
            omega = sys%mass(atom)*dble(pimc%NumBeads)*(0.5*pimc%invBeta)
           
            !The beads at position which are non-zero modulo 3 have the
            !same weighting where as the other has a different weighting
            if(mod(start_bead,3).eq.0) then
                cPrev=omega/(2.0d0*pimc%act%t0)
            else
                cPrev=omega/(0.5d0-pimc%act%t0)
            endif
            
            cStart=cPrev
            
            !Compute the staging coefficients for each value of k
            do i=1,pimc%move%MovesPerStep
                bead_ind=mod(start_bead+i-1,pimc%NumBeadsEff)+1 !+1,-1 to ensure correct looping when edge of array reached

                !The beads at position which are non-zero modulo 3 have the same 
                !weighting where as the other has a different weighting
                if(mod(bead_ind,3).eq.0) then
                    cCurr=omega/dble(2.0*pimc%act%t0)
                else
                    cCurr=omega/dble(0.5-pimc%act%t0)
                endif

           
                !Calculate and store the current q value
                qCurr=cCurr+cPrev-qPrev*aPrev**2
                q_vals(bead_ind)=qCurr
              
                !Calculate and store the current a value
                aCurr=cCurr/qCurr
                a_vals(bead_ind)=aCurr
            
                !Calculate and store the current b value
                bCurr=qPrev*aPrev*bPrev
                
                !If we are at the first staged point then the kronecker delta is 1 and 
                !we need to include the c term otherwise we do not use this term
                if(i.eq.1) then
                    bCurr=bCurr+cStart
                endif
                bCurr=bCurr/qCurr
                b_vals(bead_ind)=bCurr
        
                !Update the previous values of these so that the next step has the required previous values
                qPrev=qCurr
                bPrev=bCurr
                aPrev=aCurr
                cPrev=cCurr
           
            enddo
            
            !Do the staging Move - here we are iterating from end_bead - 1 to start_bead + 1
            do i=pimc%move%MovesPerStep,1,-1
                bead_ind=mod(start_bead+i-1,pimc%NumBeadsEff)+1 !correctly wraps 0 to pimc%NumBeads*3
                do n=1,sys%dimen
                    !rho random value from gaussian of variance 1 - sigma standard deviation of the gaussian  
                    rho=gauss_dev(seedval)
                    sigma=sqrt(1.0/(2.0*q_vals(bead_ind)))
           
                    rStart=Beads(start_bead)%x(n,atom)
           
                    rNext=Beads(mod(bead_ind,pimc%NumBeadsEff)+1)%x(n,atom)
                    Beads(bead_ind)%x(n,atom)=sigma*rho + a_vals(bead_ind)*rNext + b_vals(bead_ind)*rStart
               
                enddo
            enddo
           
            start = start_bead+1
            last = start_bead+pimc%move%MovesPerStep
           
        !If atom_move move atoms
        else
            do n=1,sys%natom
                do j=1,sys%dimen
                    rho = gauss_dev(seedval)
                    sigma = sqrt(pimc%move%AtomDisp/sys%mass(n))
       
                    do i=1,pimc%NumBeadsEff
                        Beads(i)%x(j,n) = Beads(i)%x(j,n) + sigma*rho
                    enddo
       
                enddo
            enddo
       
            start = 1
            last=pimc%NumBeadsEff
        endif

        return
    end subroutine

end module moves

