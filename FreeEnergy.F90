!Helper module used in free energy calcutions 

#ifdef FREE_ENERGY

module free_energy
    use molecule_specs
    use pimc_structures

    implicit none
    
    type free_energy_type
        real(kind=8), dimension(:,:,:), allocatable :: coords
        real(kind=8), dimension(:), allocatable :: int_internal_1
        real(kind=8), dimension(:), allocatable :: int_internal_2
        real(kind=8), dimension(:), allocatable :: int_internal_3
        real(kind=8), dimension(3) :: curr,prev, cum_sum
    end type free_energy_type

    interface new
        module procedure init_free_energy
    end interface new

    contains
    !###############################################################################
    !initialise the free energy object differently depending on the free energy type
    subroutine init_free_energy(this,pimc,atoms,dimen)
        type (free_energy_type), intent(out)  :: this
        type (pimc_par), intent(in) :: pimc
        integer, intent(in)  ::  atoms, dimen
        integer :: ierr,i
        if(pimc%free%free_type==0) then
            allocate(this%coords(dimen,atoms,pimc%NumBeadsEff), stat=ierr)
            if(ierr.ne.0) stop 'Error allocating freeEnergyCoords'

        elseif(pimc%free%free_type==2) then
            allocate(this%int_internal_1(pimc%StepsPerBlock), stat=ierr)
            if(ierr.ne.0) stop 'Error allocating int_internal'

            if(pimc%act%act_type==0) then
                allocate(this%int_internal_2(pimc%StepsPerBlock), stat=ierr)
                if(ierr.ne.0) stop 'Error allocating int_internal'
                allocate(this%int_internal_3(pimc%StepsPerBlock), stat=ierr)
                if(ierr.ne.0) stop 'Error allocating int_internal'
            endif
        endif
        
        do i=1,3
            this%prev(i)=0.0
            this%curr(i)=0.0
            this%cum_sum(i)=0.0
        enddo

        return
    end subroutine init_free_energy 
    !###############################################################################




    !###############################################################################
    !open the output file for reversible scaling
    subroutine openFileOut(OUT_DIR)
        character(len=80), intent(in) :: OUT_DIR
        open(unit=302,file=trim(OUT_DIR)//'/reversible_scaling.out',status='unknown')
    end subroutine openFileOut

    !close the output file for reversible scaling
    subroutine closeFileOut()
        close(unit=302)
    end subroutine closeFileOut

    !updates the temperature changing beta by the amount beta_step and 
    !inverting to find the temperature

    subroutine tempStep(pimc)
        type(pimc_par) :: pimc
        pimc%OldBeta = pimc%Beta
        pimc%Beta = pimc%Beta + pimc%free%beta_step
        pimc%invBeta =1.0d0/pimc%Beta
        pimc%Temperature = pimc%invBeta/3.16681520371153d-6
    end subroutine tempStep


    subroutine writeIntInternal(this,pimc,B_init)
        type (free_energy_type), intent(in) :: this
        type(pimc_par), intent(in) :: pimc
        real(kind=8), intent(in) :: B_init
        integer :: i
        
        if(pimc%act%act_type==0) then
            do i=1,pimc%StepsPerBlock
                write(302,*) B_init + (i-1)*pimc%free%beta_step, this%int_internal_1(i)*pimc%free%beta_step, &
                &            this%int_internal_2(i)*pimc%free%beta_step, this%int_internal_3(i)*pimc%free%beta_step
            enddo
        else
            do i=1,pimc%StepsPerBlock
                write(302,*) B_init + (i-1)*pimc%free%beta_step, this%int_internal_1(i)*pimc%free%beta_step
            enddo
        endif

    end subroutine writeIntInternal
    
    subroutine reversibleScaling(this,Beads,pimc,sys)
        use Estimator_class

        type (free_energy_type) :: this
        type (pimc_particle), dimension(:), pointer :: Beads
        type (pimc_par), intent(in) :: pimc
        type (molsysdat), intent(in)  :: sys

        real(kind=8), dimension(3) :: results
        if(pimc%act%act_type==0) then
            call energy_thermo_pa(sys,pimc,Beads,results)
            this%curr(1) = this%curr(1) + results(1)
            results = 0

            call energy_virial_pa(sys,pimc,Beads,results)
            this%curr(2) = this%curr(2) + results(1)
            results = 0

            call energy_centvirial_pa(sys,pimc,Beads,results)
            this%curr(3) = this%curr(3) + results(1)
            results = 0

        else if (pimc%act%act_type==1) then
            call energy_thermo_ta(sys,pimc,Beads,results)
            this%curr(1) = this%curr(1) + results(1)
            results = 0

        else if(pimc%act%act_type==2) then
            call energy_thermo_ca(sys,pimc,Beads,results)
            this%curr(1) = this%curr(1) + results(1)
            results = 0
        else
            stop 'Invalid action type reached'
        endif

    end subroutine reversibleScaling
    
    !update the values of the value of the integral at the end of each monte carlo pass
    subroutine reversibleScalingStep(this,pimc,indx,normalise)
        type (free_energy_type) :: this
        type (pimc_par), intent(in) :: pimc
        integer, intent(in) :: indx
        real(kind=8), intent(in) :: normalise

        this%curr = this%curr/normalise
        this%cum_sum = this%cum_sum+(this%curr+this%prev)/2.0

        if(pimc%act%act_type==0) then
            this%int_internal_1(indx) = this%cum_sum(1)

            this%int_internal_2(indx) = this%cum_sum(2)

            this%int_internal_3(indx) = this%cum_sum(3)
            
        else if (pimc%act%act_type==1) then
            this%int_internal_1(indx) = this%cum_sum(1)
        else if(pimc%act%act_type==2) then
            this%int_internal_1(indx) = this%cum_sum(1)
        else
            stop 'In valid action type reached'
        endif

        this%prev=this%curr
        this%curr=0

    end subroutine reversibleScalingStep

    !###############################################################################








    !###############################################################################
    !Subroutines used for thermodynamic integration from classical to quantum mechanics
    !###############################################################################
    subroutine scaleCoords(this,Beads,pimc,sys)
        type (free_energy_type), intent(inout) :: this
        type (pimc_particle), dimension(:), pointer :: Beads
        type (pimc_par), intent(in) :: pimc
        type (molsysdat), intent(in)  :: sys
        integer :: i,j,k
        real(kind=8), dimension(sys%dimen, sys%natom) :: cv

        do i=1,pimc%NumBeadsEff
            do j=1,sys%natom
                do k=1,sys%dimen
                    this%coords(k,j,i) = Beads(i)%x(k,j) 
                    cv(k,j) = cv(k,j) + Beads(i)%x(k,j)
                enddo
            enddo
        enddo
        cv=cv/dble(pimc%NumBeadsEff)

        do i=1,pimc%NumBeadsEff
            do j=1,sys%natom
                do k=1,sys%dimen
                    Beads(i)%x(k,j) = Beads(i)%x(k,j)*pimc%free%lambda + &
                    &                 (1-pimc%free%lambda)*cv(k,j)
                enddo
            enddo
        enddo
        return
    end subroutine scaleCoords

    subroutine unscaleCoords(this,Beads,pimc,sys)
        type (free_energy_type), intent(inout) :: this
        type (pimc_particle), dimension(:), pointer :: Beads
        type (pimc_par), intent(in) :: pimc
        type (molsysdat), intent(in)  :: sys
        integer :: i,j,k

        do i=1,pimc%NumBeadsEff
            do j=1,sys%natom
                do k=1,sys%dimen
                    Beads(i)%x(k,j) = this%coords(k,j,i) 
                enddo
            enddo
        enddo
        return
    end subroutine unscaleCoords
    !###############################################################################


end module free_energy

#endif

