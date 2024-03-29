module pimc_structures
    implicit none
    !-----------------------
    !Action parameters
    !-----------------------
    type act_par
        integer :: act_type
        !parameters for the
        real(kind=8) :: t0,a1,u0,v1,v2,t1, t1inv, t0inv, alpha
    end type act_par

    !-----------------------
    !Move parameters
    !-----------------------
    type move_par
        integer :: move_type
        !parameters
        real(kind=8) :: AtomDisp
        integer :: MovesPerStep
        real(kind=8) :: BeadDisp
    end type move_par

#ifdef FREE_ENERGY
    !-----------------------
    !Free Energy Parameters
    !-----------------------
    type free_eng
        integer :: free_type

        real(kind=8) :: lambda
        real(kind=8) :: temp_init
        real(kind=8) :: temp_fin
        real(kind=8) :: beta_step
    end type free_eng
#endif
    !-----------------------
    ! pimc parameters
    !-----------------------
    type pimc_par
        !Number of initial geometries, Beads for each geometry, Number of Beads, Number of effective beads (3*NumBeads for chin action), Effective Beads per geom = 3*Beads per Geom for chin action
        integer :: NumBeads, NumBeadsEff
        integer :: NumBlocks, BlocksToEquil, StepsPerBlock, NumBlocksLeft,BlocksToEquilLeft
        integer :: Sample
        character(len=80) :: start, resume, blk

        !stores the action parameters
        type (act_par) :: act

        !stores the move parameters
        type (move_par) :: move
#ifdef FREE_ENERGY
        !stores the free energy parameters
        integer :: doFree, doSample
        type (free_eng) :: free
        real(kind=8) :: OldBeta
#endif
        real(kind=8) :: IniDisp, Temperature
        real(kind=8) :: Beta
        real(kind=8) :: invBeta !1/beta = kb*T
        real(kind=8) :: invNumBeads !1/P
        character(len=1) :: Restart 
        character(len=1) :: blocking, WritingCheckpoint  

    end type pimc_par
    !-----------------------
    ! the bead data type
    !-----------------------
    type pimc_particle
        real(kind=8) :: Vcurr
        real(kind=8), dimension(:,:), pointer :: x    ! position
        real(kind=8), dimension(:), pointer :: r      ! bond lengths
        real(kind=8), dimension(:,:), pointer :: dVdx     ! derivs w.r.t. cartesians
        !real(kind=8), dimension(:,:,:,:), pointer :: d2Vdx2   ! derivs w.r.t. cartesians
    end type pimc_particle
    !-------------------------
    ! overload constructor
    !-------------------------
    interface new
        module procedure pimc_particle_init
    end interface
    interface copy
        module procedure pimc_particle_copy
    end interface

    contains
    !----------------------------
    !  pimc_particle constructor
    !----------------------------
    subroutine pimc_particle_init(this,natom,ndimen)
        type (pimc_particle) :: this
        integer, intent(in) :: natom
        integer, intent(in) :: ndimen
        integer :: ierr, nbond
        nbond = natom*(natom-1)/2
        this%Vcurr = 0.0
        allocate(this%x(ndimen,natom),stat=ierr)
        if (ierr.ne.0) stop ' Error allocating pimc_particle%x '
        allocate(this%r(nbond),stat=ierr)
        if (ierr.ne.0) stop ' Error allocating pimc_particle%r '
        allocate(this%dVdx(ndimen,natom),stat=ierr)
        if (ierr.ne.0) stop ' Error allocating pimc_particle%dVdx'
        !allocate(this%d2Vdx2(ndimen,ndimen,natom,natom),stat=ierr)
        !if (ierr.ne.0) stop ' Error allocating pimc_particle%d2Vdx2'
        return
    end subroutine pimc_particle_init
    !----------------------------
    !  pimc_particle copy
    !----------------------------
    subroutine pimc_particle_copy(Win,Wout)
        type (pimc_particle), intent(in) :: Win
        type (pimc_particle), intent(out) :: Wout
        integer :: i,j
        Wout%Vcurr = Win%Vcurr
        do i=1,size(Win%x,dim=1)
            do j=1,size(Win%x,dim=2)
                Wout%x(i,j) = Win%x(i,j)
            enddo
        enddo
        do i=1,size(Win%r)
            Wout%r(i) = Win%r(i)
        enddo
        do i=1,size(Win%dVdx,dim=1)
            do j=1,size(win%dVdx,dim=2)
                Wout%dVdx(i,j) = Win%dVdx(i,j)
            enddo
        enddo
        !do i=1,size(Win%d2Vdx2,dim=1)
            !do j=1,size(win%d2Vdx2,dim=2)
                !do k=1,size(win%d2Vdx2,dim=3)
                    !do l=1,size(win%d2Vdx2,dim=4)
                        !Wout%d2Vdx2(i,j,k,l) = Win%d2Vdx2(i,j,k,l)
                    !enddo
                !enddo
            !enddo
        !enddo

        return
    end subroutine pimc_particle_copy

end module pimc_structures

