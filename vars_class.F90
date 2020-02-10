!A module holding all variables and containing all functions required for computing the average and
!variance of a quantity using the online variance algotirhm
module vars_class
    use pimc_structures
    implicit none

    !Stores the variables required for averaging a quantity over a simulation run
    type vars
        real(kind=8) :: mean_block, var_block, mean_tot, var_tot, curr
        real(kind=8) :: delta, diffsqr, equil_sum, equil_move_avg
        integer(kind=4) :: n_block, n_tot

    end type vars

    !Define functions to be used by the vars class
    interface new
        module procedure init_vars
    end interface 

    interface update
        module procedure update_step_var
    end interface

    interface update_block
        module procedure update_block_var
    end interface

    interface print_block
        module procedure print_var_block
    end interface

    interface print_end
        module procedure print_var_end
    end interface

    interface equilibrated_sum
        module procedure update_equilibrated_sum
    end interface

    interface reset
        module procedure reset_block_var
    end interface 

    contains
        subroutine init_vars(this)!,var_name
            type (vars), intent(out) :: this
            this%mean_block=0.0
            this%var_block=0.0
            this%mean_tot=0.0
            this%var_tot=0.0
            this%delta=0.0
            this%diffsqr=0.0
            this%curr=0.0
            this%n_block=0
            this%n_tot=0
            this%equil_sum=0.0
            this%equil_move_avg = 0.0

            return
        end subroutine init_vars

        subroutine update_step_var(this,val)
            type (vars), intent(inout) :: this
            real(kind=8) :: val
            this%curr=val
            !print *, 'this%curr = ', val
            this%n_block=this%n_block+1
!            print *, 'curr', this%curr, 'mean_block', this%mean_block
            this%delta=this%curr-this%mean_block
            this%mean_block=this%mean_block+this%delta/dble(this%n_block)
            
            this%diffsqr=this%diffsqr+(this%curr-this%mean_block)**2

            !print *, 'n_block = ', this%n_block

            return
        end subroutine update_step_var

        subroutine reset_block_var(this)
            type (vars), intent(inout) ::this
            this%n_block=0
            this%var_block=0.0
            this%mean_block=0.0
            this%diffsqr=0.0
            this%equil_move_avg = 0.0

        end subroutine reset_block_var

        subroutine update_block_var(this)
            type (vars), intent(inout) :: this
            !print *, 'this%n_block = ', this%n_block
            !If there was more than one step in the block
            if(this%n_block.gt.1) then
               ! !Calculate the block variance given it is well defined
                this%var_block=this%diffsqr/dble(this%n_block-1)

                !If we have previously had a block of length greater than one
                if(this%n_tot.gt.1) then
            
                    !Calculate the combined variance of the two data sets
                    this%var_tot=(dble(this%n_block-1)*this%var_block+dble(this%n_tot-1)*this%var_tot&
                    &           +((this%mean_block-this%mean_tot)**2)*dble(this%n_block*this%n_tot)&
                    &           /dble(this%n_block+this%n_tot))/(dble(this%n_block+this%n_tot-1))
               !     print *, 'calculated via n_tot > 1'
                    !Calculate the weighted average of the means
                    this%mean_tot=(this%mean_tot*dble(this%n_tot)+this%mean_block*dble(this%n_block))&
                    &               /dble(this%n_block+this%n_tot)

                else
                    !Set the total variance to the block variance
                    this%var_tot=this%var_block
                !    print *, 'calculated via else'
                    !Update the mean for the total simulation
                    this%mean_tot=(this%mean_tot*dble(this%n_tot)+this%mean_block*dble(this%n_block))&
                    &               /dble(this%n_block+this%n_tot)
                 endif
                 !   print *, 'this%var_tot = ', this%var_tot
            endif
            !If we have only one step in the block then the variance of the block is not well defined
            if(this%n_block.eq.1) then
                !set block variance to zero
                this%var_block=0.0
                !update the total mean
                this%mean_tot=(this%mean_tot*dble(this%n_tot)+this%mean_block*dble(this%n_block))&
                &               /dble(this%n_block+this%n_tot)
                
                !update the total variance
                this%var_tot=(dble(this%n_tot-1)*this%var_tot+((this%mean_block-this%mean_tot)**2)&
                &           *dble(this%n_block*this%n_tot)/dble(this%n_block+this%n_tot))&
                &           /(dble(this%n_block+this%n_tot-1))

            endif   
            !Update the total number of times that this has been called
            this%n_tot=this%n_tot+this%n_block
        end subroutine update_block_var




        subroutine print_var_block(pimc,this)
            type (vars), intent(inout) :: this
            type (pimc_par), intent(in) :: pimc  
               
              if(this%n_block.ne.0) then
                if (pimc%blocking == 'n') then
                   write(*,*) 2625.5*this%mean_block, '+/-', 2625.5*sqrt(this%var_block/this%n_block)!, &
                   !&           'Block Size: ', this%n_block
                else
                   open(unit=199,file=trim(pimc%blk),status='unknown',action='write',position='append') 
                   ! Do not make a dynamical array for this, because you cannot
                   ! checkpoint easily
                   write(199,*) this%mean_block
                   close(unit=199)
                endif
              endif
            !this%n_block=0
            !this%var_block=0.0
            !this%mean_block=0.0
            !this%diffsqr=0.0
        end subroutine print_var_block

        subroutine update_equilibrated_sum(this)
            type (vars), intent(inout) :: this
            real(kind=8) :: nblk = 0.0
            
            nblk=nblk + 1.0
           
            this%equil_sum = this%equil_sum + this%mean_block 
            this%equil_move_avg = this%equil_sum/nblk
           
            !print *, this%equil_move_avg

        end subroutine update_equilibrated_sum



        subroutine print_var_end(this)
            type (vars), intent(in) :: this
            if(this%n_tot.ne.0) then
            write(*,*) 2625.5*this%mean_tot, '+/-', 2625.5*sqrt(this%var_tot/this%n_tot), &
            &           'Averages: ', this%n_tot
            endif
        end subroutine print_var_end


end module vars_class
