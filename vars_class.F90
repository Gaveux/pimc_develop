!A module holding all variables and containing all functions required for computing the average and
!variance of a quantity using the online variance algotirhm
module vars_class
    implicit none

    !Stores the variables required for averaging a quantity over a simulation run
    type vars
        real(kind=8) :: mean_block, var_block, mean_tot, var_tot, curr
        real(kind=8) :: delta, diffsqr
        integer(kind=4) :: n_block, n_tot, n_block_save

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

    !interface reset
     !   module procedure reset_block_var
    !end interface 

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

            return
        end subroutine init_vars

        subroutine update_step_var(this,val)
            type (vars), intent(inout) :: this
            real(kind=8) :: val
            

            ! we use sample standard deviation for error propagation within each 
            ! block since each evaluation here only estimates the energy instead 
            ! of being exact

            this%curr=val
            !this%n_block=this%n_block+1
            ! use recursive calculation to calculate mean
            !this%mean_block=((this%mean_block*dble(this%n_block-1))+this%curr)/dble(this%n_block)
            !print *, this%mean_block, this%curr
            !this%delta=this%curr-this%mean_block
            !print *,this%delta
            !this%diffsqr=this%diffsqr+this%delta*this%delta
            !print *,this%diffsqr
           
            !=====original code===================
            this%n_block=this%n_block+1
            ! this is cumulative mean average
            this%delta=this%curr-this%mean_block
            this%mean_block=this%mean_block+this%delta/dble(this%n_block)
            this%diffsqr=this%diffsqr+(this%curr-this%mean_block)*(this%curr-this%mean_block)
            !==================================
            !print *, 'n_block = ', this%n_block
            return
        end subroutine update_step_var

        !subroutine reset_block_var(this) ! this is never called in estimator_class
         !   type (vars), intent(inout) ::this
          !  this%n_block=0
           ! this%var_block=0.0
           ! this%mean_block=0.0
           ! this%diffsqr=0.0

        !end subroutine reset_block_var

        subroutine update_block_var(this)
            type (vars), intent(inout) :: this
            !this%var_block=this%diffsqr/dble(this%n_block-1)

            !print *,'n_tot', this%n_tot
            !this%n_tot = this%n_tot + 1
            !print *, 'n_tot = ', this%n_tot
            !this%mean_tot=((this%mean_tot*dble(this%n_tot-1))+this%mean_block)/dble(this%n_tot)
            !this%var_tot= this%var_tot+(this%mean_block-this%mean_tot)*(this%mean_block-this%mean_tot)
            !If there was more than one step in the block
            !if(this%n_block.gt.1) then
                !Calculate the sample standard deviation for the block
                !print *, this%n_block, this%diffsqr
                this%var_block=this%diffsqr/dble(this%n_block-1)


                ! using the definition of the pooled variance
                this%var_tot = this%var_tot + this%diffsqr
                !print *, 'var_tot = ', this%var_tot, 'diffsqr = ', this%diffsqr
                !print *, 'this%n_tot', this%n_tot
 
                !If we have previously had a block of length greater than one
                !if(this%n_tot.gt.1) then
                    
                    !Calculate the average of the means
            !         print *, 'mean_block', this%mean_block
                     
                    !======= original code=================
                    !Calculate the combined variance of the two data sets
                !    print *, 'this%n_block = ', this%n_block, 'this%n_tot = ', this%n_tot 
                !    this%var_tot=(dble(this%n_block-1)*this%var_block+dble(this%n_tot-1)*this%var_tot&
                !    &           +((this%mean_block-this%mean_tot)**2)*dble(this%n_block*this%n_tot)&
                !    &           /dble(this%n_block+this%n_tot))/(dble(this%n_block+this%n_tot-1))
                !    print *, this%var_tot
                    !Calculate the weighted average of the means
                    this%mean_tot=(this%mean_tot*dble(this%n_tot)+this%mean_block*dble(this%n_block))&
                    &               /dble(this%n_block+this%n_tot)
                    !print *, 'mean_tot = ', this%mean_tot
                    !======================================

                !else
                    !Set the total variance to the block variance
                !    this%var_tot=this%var_block
                    
                    !Update the mean for the total simulation
                !    this%mean_tot=(this%mean_tot*dble(this%n_tot)+this%mean_block*dble(this%n_block))&
                !    &               /dble(this%n_block+this%n_tot)
                ! endif
            !endif
            !If we have only one step in the block then the variance of the block is not well defined
            !if(this%n_block.eq.1) then
                !set block variance to zero
            !    this%var_block=0.0
                !update the total mean
            !    this%mean_tot=(this%mean_tot*dble(this%n_tot)+this%mean_block*dble(this%n_block))&
            !    &               /dble(this%n_block+this%n_tot)
                
                !update the total variance
            !    this%var_tot=(dble(this%n_tot-1)*this%var_tot+((this%mean_block-this%mean_tot)**2)&
            !    &           *dble(this%n_block*this%n_tot)/dble(this%n_block+this%n_tot))&
            !    &           /(dble(this%n_block+this%n_tot-1))

            !endif   
            !Update the total number of times that this has been called
            this%n_tot=this%n_tot+this%n_block
            !print *, 'end', this%n_tot
        end subroutine update_block_var




        subroutine print_var_block(this)
            type (vars), intent(inout) :: this
            if(this%n_block.ne.0) then
            ! returns the mean of block and its standard deviation of the mean
            !write(*,*) this%mean_block, '+/-', sqrt(this%var_block/(this%n_block*this%n_block)), &
            !&           'Block Size: ', this%n_block
            !print *, this%var_block
            endif
            this%n_block_save = this%n_block
            this%n_block=0
            this%var_block=0.0
            this%mean_block=0.0
            this%diffsqr=0.0
        end subroutine print_var_block




        subroutine print_var_end(this)
            type (vars), intent(in) :: this
            if(this%n_tot.ne.0) then
            !print *, "n_block_save = ", this%n_block_save, 'MC steps = ', this%n_tot/this%n_block_save
            ! this calculates the mean of MC-block and its population standard deviation of the mean
            write(*,*) this%mean_tot, '+/-', sqrt(this%var_tot/(dble(this%n_tot-(this%n_tot/this%n_block_save))*dble(this%n_tot)))
            !&           'Average+)
            !print *, 'MC block = ', this%n_tot/this%n_block_save
            !print *, this%n_tot - (this%n_tot/this%n_block_save)
            !print *, dble(this%n_tot)*dble(this%n_tot-(this%n_tot/this%n_block_save))
            endif
        end subroutine print_var_end

end module vars_class
