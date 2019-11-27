

  subroutine numerical_integrator(sys,pimc,pot)
   use molecule_specs
   use pimc_structures
   use potential_msi
 
   implicit none
   
   type(molsysdat) :: sys
   type(pimc_par) :: pimc
   type (msi_params) :: pot

   real(kind=8), dimension(sys%dimen,sys%natom) :: x_old, x_older, x_temp, dVdx
   real(kind=8), dimension(sys%nbond) :: r
   real(kind=8) :: VCurr
    
   
   
   integer :: i,j,k

   ! initial conditions
   x_old = sys%Bead2Geom
   x_older = sys%Bead1Geom
   !OlderDiscret%x = sys%Bead1Geom
   ! numerical integrator
   do k=1,pimc%NumDiscretisation
       ! call potential

       !call potential(i-1,pot,OldDiscret(i-1)%x,OldDiscret(i-1)%r,OldDiscret(i-1)%VCurr,OldDiscret(i-1)%dVdx)
       call potential(k,pot,x_old,r,VCurr,dVdx)
       !print *, 'dVdx' 
       !print *, dVdx 
       ! integrator
       !do j=1,sys%natom
       !   do k=1,sys%dimen
       !       OldDiscret(i)%x(k,j) = 2.0*OldDiscret(i)%x(k,j) - OldDiscret(i-1)%x(k,j) + 2.0*pimc%invNumDiscretisation**2*sys%invMass(j)*pimc%Beta**2*OldDiscret(i-1)%dVdx(k,j)
       !   enddo
       !enddo
       do j=1,sys%natom
          do i=1,sys%dimen
             x_temp(i,j) = 2.0*x_old(i,j) - x_older(i,j) + 2.0*pimc%invNumDiscretisation**2*sys%invMass(j)*pimc%Beta**2*dVdx(i,j)
          enddo
       enddo 
       !x_temp = 2.0*x_old - x_older
       x_older = x_old
       x_old = x_temp 
   enddo
   print *, x_temp
    

  end
