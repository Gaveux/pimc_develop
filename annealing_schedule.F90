module annealing_schedule
  use pimc_structures 
  implicit none

  contains

  subroutine annealing_condition(pimc,AtomDisp,iblock,first_time)
    implicit none
  
    type (pimc_par) :: pimc
    real(kind=8), intent(out) :: AtomDisp
    logical :: first_time

    integer :: iblock
    real(kind=8) :: totalIteration, currentIteration, move
    real(kind=8) :: min_AtomDisp
    
     currentIteration=dble(iblock*pimc%StepsPerBlock)
     totalIteration=dble(pimc%NumBlocks*pimc%StepsPerBlock)

     if (first_time.eq..TRUE.) then
     AtomDisp = pimc%move%AtomDisp
     else
     !determine at what iteration k updates
        min_AtomDisp = 0.3d0
        move = 0.1d0*totalIteration
        

        if (mod(currentIteration,move).eq.0) then
           AtomDisp = pimc%move%AtomDisp*(1.0 - currentIteration/totalIteration)**4.0

           ! ideally we don't want AtomDisp to be too small
           if (AtomDisp .le. min_AtomDisp) then
              AtomDisp = min_AtomDisp
           endif
        endif
     endif


    return
  end subroutine

end module
