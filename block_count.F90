!-------------------------------------------------------------------------
!  Countint data size and parse data into blocking transformation
!-------------------------------------------------------------------------

  module block_count
      use pimc_structures
      implicit none
      
      contains

      subroutine blk_count(b_filename)
          type(pimc_par), intent(in) :: pimc
          real(kind=8), dimension(:), allocatable :: blockdata
          integer i
          character(len=80), intent(in) :: b_filename
          

          ! read and store all data into blockdata
          open (unit=13, file='b_file',status='old', action='read')
          rows=0 !Count the number of lines in the file
          
          do
            read(1,*,iostat=io)
            if (io/=0) exit
               rows=rows+1
          enddo

          rewind(13)

          print *, 'number of rows = ', rows

          allocate(blockdata(rows))

          do i=1,rows,1
            read(1,*) blockdata(i)
          enddo
          close(13)
          !print *, blockdata
          

          ! do blocking transformation 
          call blocking(blockdata)
      end subroutine blk_count


  end module block_count
