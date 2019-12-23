

  subroutine numerical_integrator(sys,pimc,pot,OUT_COORD)
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
   character(len=80), intent(in) :: OUT_COORD
   CHARACTER(*), PARAMETER :: fileplace = "/home/gaveuxifort/PIMC_DEV/inputs/OUT_COORD/"
    
   
   
   integer :: i,j,k
   print *, 'beta**2 = ', pimc%Beta**2
   ! initial conditions
   x_old = sys%Bead2Geom
   x_older = sys%Bead1Geom
   !OlderDiscret%x = sys%Bead1Geom
   ! numerical integrator
   !do k=1,pimc%NumDiscretisation
   open(unit=333,file=fileplace//trim(OUT_COORD),status='unknown',action='write',position='append')
   write(333,*) x_older(:,1)
   write(333,*) x_old(:,1)
   do k=1, pimc%NumDiscretisation
       ! call potential

       call potential(k,pot,x_old,r,VCurr,dVdx)
       !print *, dVdx
       !do i=1,sys%natom
       !  print *, sys%invMass(i)*(pimc%Beta*pimc%invNumDiscretisation)**2
       !enddo
       !print *, sys%invMass(1)
       !print *, dble(pimc%invNumDiscretisation*pimc%Beta)**2
       ! integrator
       do j=1,sys%natom
          do i=1,sys%dimen
             x_temp(i,j) = 2.0*x_old(i,j) - x_older(i,j) - sys%invMass(j)*dVdx(i,j)*(pimc%invNumDiscretisation*pimc%Beta)**2
          enddo
          !print *, dVdx
       enddo
       !do i=1, sys%dimen
          !print *, x_temp(i,1), x_old(i,1)
          !print *, dVdx(i,1)
          !print *, (2.0*(pimc%invNumDiscretisation**2)*sys%invMass(1)*(pimc%Beta**2)*dVdx(i,1),i=1,sys%dimen)
       !enddo
       write(333,*) (x_temp(i,1), i = 1, sys%dimen)
       !print *, (x_temp(i,1), i = 1, sys%dimen)
      ! print *, x_old 
      ! print *, x_temp
       x_older = x_old
       x_old = x_temp 
   enddo
       close(333)
   
   print *, x_temp

  end
