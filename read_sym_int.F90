   subroutine read_input(pimc,in_file)

       use pimc_structures
     
       implicit none

       type(pimc_par), intent(out) :: pimc
       character(len=80), intent(in) :: in_file
       character(len=80) :: icomm
  
       open(unit=17,file=trim(in_file),status='old')

80     format(a80)

       ! read in Trotter number
       read(17,80) icomm
       write(21,80) icomm
       read(17,*) pimc%NumDiscretisation
       write(21,*) pimc%NumDiscretisation
       
       ! read in Temperature
       read(17,80) icomm
       write(21,80) icomm
       read(17,*) pimc%Temperature
       write(21,*) pimc%Temperature
 
       close(unit=17)
       
       pimc%invBeta=3.16681520371153d-6*dble(pimc%Temperature)
       pimc%Beta=1.0/pimc%invBeta
       pimc%invNumDiscretisation = 1.0d0/dble(pimc%NumDiscretisation)
      
       return
   end
