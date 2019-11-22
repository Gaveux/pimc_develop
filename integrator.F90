

  subroutine numerical_integrator(sys,pimc,pot)
   use molecule_specs
   use pimc_structures
   use potential_msi
 
   implicit none
   
   type(molsysdat) :: sys
   type(pimc_par) :: pimc
   type (msi_params) :: pot
   type(pimc_particle), dimension(:), pointer :: OldDiscret
   type(pimc_particle), dimension(:), pointer :: OlderDiscret
   type(pimc_particle), dimension(:), pointer :: OldDiscret_temp

   !type(pimc_particle), dimension(:), pointer :: Beads
   !type(pimc_particle), dimension(:), pointer :: OldBeads

   
   
   integer :: i,j,k,ierr
   !print *, 'iterating', pimc%NumDiscretisation, 'steps'
   !do i=1,sys%natom 
   !print *, 'atom', i
   !print *, (sys%Bead1Geom(k,i),k=1,sys%dimen)
   !enddo
   !print * ,'bead1Geom'
   !print *, sys%Bead1Geom
   !print *, 'bead2geom'
   !print *, sys%Bead2Geom
   !print *, '2*bead2geom - bead1geom'
   !print *, 2.0*sys%Bead2Geom - sys%Bead1Geom
   
   ! declare dynamical arrays
   allocate(OldDiscret(pimc%NumDiscretisation),stat=ierr)
   if (ierr.ne.0) stop 'Error allocating OldDiscret array'
   allocate(OlderDiscret(pimc%NumDiscretisation),stat=ierr)
   if (ierr.ne.0) stop 'Error allocating OlderDiscret array'
   allocate(OldDiscret_temp(pimc%NumDiscretisation),stat=ierr)
   if (ierr.ne.0) stop 'Error allocating temporal OlderDiscret array'
   do i=0,pimc%NumDiscretisation
      call new(OldDiscret(i),sys%natom,sys%dimen)
      call new(OlderDiscret(i),sys%natom,sys%dimen)
      call new(OldDiscret_temp(i),sys%natom,sys%dimen)
   enddo
   
   ! create the  geometries as the initial conditions
   
   !OldDiscret(0)%x = sys%Bead2Geom
   !print *, 'oldDiscret(1)' 
   !print *, OldDiscret(1)%x
   !OlderDiscret(1)%x = sys%Bead1Geom
   !print *, 'OldDiscret(0)'
   !print *, OldDiscret(0)%x
   print *, 'Bead1Geom'
   print *, sys%Bead1Geom
   print *, 'Bead2Geom'
   print *, sys%Bead2Geom


   ! compute dV/dr_1 for the initial conditions w.r.t. potential
   !call potential(1,pot,OldDiscret(1)%x,OldDiscret(1)%r,OldDiscret(1)%VCurr,OldDiscret(1)%dVdx)
   !print *,'dVdx for x_1'
   !print *, OldDiscret(1)%dVdx
   !print *, 'invNumDiscretisation = ', pimc%invNumDiscretisation
   !print *, 'beta = ', pimc%beta
   !print *, 'temperature = ', pimc%Temperature
   
   ! recursion of integrator
   !do i=2,pimc%NumDiscretisation
       ! initial conditions
       !OldDiscret(1)%x = sys%Bead2Geom
       !OldDiscret(0)%x = sys%Bead1Geom
       
       ! call potential
   !    call potential(i-1,pot,OldDiscret(i-1)%x,OldDiscret(i-1)%r,OldDiscret(i-1)%VCurr,OldDiscret(i-1)%dVdx)
       !print *, 'dVdx for r_{k}'
       !print *, OldDiscret(i-1)%dVdx       

       ! integrator
       !do j=1,sys%natom
       !   do k=1,sys%dimen
       !       OldDiscret(i)%x(k,j) = 2.0*OldDiscret(i)%x(k,j) - OldDiscret(i-1)%x(k,j) + 2.0*pimc%invNumDiscretisation**2*sys%invMass(j)*pimc%Beta**2*OldDiscret(i-1)%dVdx(k,j)
       !   enddo
       !enddo
       
       ! a cheaper method of doing integration
       !OldDiscret_temp(i)%x = OldDiscret(i)%x 
       !OldDiscret(i)%x = 2.0*OldDiscret(i)%x - OlderDiscret(i)%x
       !OlderDiscret(i)%x = OldDiscret_temp(i)%x
       
   !enddo
    
   OldDiscret(1)%x=sys%Bead2Geom
   OlderDiscret(1)%x=sys%Bead1Geom
   do i = 2, pimc%NumDiscretisation
     print *, OldDiscret(1)%x
     OldDiscret(i)%x=2.0*OldDiscret(i-1)%x-OlderDiscret(i-1)%x
     OlderDiscret(i)%x=OldDiscret(i-1)%x
   enddo 

  end
