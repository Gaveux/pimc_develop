!-------------------------------------------------------------------------
!  Countint data size and parse data into blocking transformation
!-------------------------------------------------------------------------

  module blocking
      use pimc_structures
      implicit none
      
      
      !type blocking_related
      !     type(block_data) :: blockdata
      !end type blocking_related
        
      contains

      subroutine blk_count(blk)
          !type(pimc_par), intent(in) :: pimc
          real(kind=8), dimension(:), allocatable :: blockdata
          integer i, rows, io
          character(len=80), intent(in) :: blk
          
          include 'blocking_transformation.int'
          

          ! read and store all data into blockdata
          open (unit=13, file=trim(blk),status='old', action='read')
          rows=0 !Count the number of lines in the file
          
          do
            read(13,*,iostat=io)
            if (io/=0) exit
               rows=rows+1
          enddo

          rewind(13)

          !print *, 'number of rows = ', rows

          allocate(blockdata(rows))

          do i=1,rows
            read(13,*) blockdata(i)
          enddo
          close(13)
          !print *, blockdata
          !print *, 'size of blockdata = ', size(blockdata) 

          ! do blocking transformation 
          call blocking_transformation(blockdata)
      end subroutine blk_count


  end module blocking
