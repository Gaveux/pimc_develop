!------------------------------------------------------------------
! calculating blocking algorithm
!------------------------------------------------------------------

subroutine blocking(blockdata,b_file)
!J. Chem. Phys. 91, 461 (1989); http://dx.doi.org/10.1063/1.457480

  real(kind=8), dimension(:), allocatable :: blockdata
  real(kind=8) :: c0, xbar, sigma, sigma2, c0divn, sigmadev, sigma2dev
  integer i,j, ndata, nblock, rows
  character(len=80), intent(in) :: b_file

! note that number of data must be even !
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
!-----------------------------------------
! do stats for untransformed data first
!-----------------------------------------

  ndata = size(blockdata)

  xbar = 0.0
  do i=1,ndata
    xbar = xbar + blockdata(i)
  enddo
  xbar = xbar/ndata

  c0 = 0.0
  do i=1,ndata
    c0 = c0 + (blockdata(i) - xbar)**2
  enddo
  c0 = c0/ndata

  c0divn = c0/(ndata-1)

  sigma2 = c0divn
  sigma2dev = sqrt(2.0/(dble(ndata-1)))*sigma2

  sigma = sqrt(c0divn)
  sigmadev = sigma/(sqrt(2*dble(ndata-1)))

  write (*,*) ''
  write (*,*) ' statistics for untransformed data'
! write (*,*) 'sigma^2 =', sigma2, '+/-', sigma2dev
  write (*,*) ' xbar = ',xbar
  write (*,*) ' sigma =', sigma, '+/-', sigmadev
  write (*,*) ''

!------------------------------------------------------------------
! calculate the number of block transformations we'll be able to do
!------------------------------------------------------------------

  nblock = int(log(dble(ndata))/log(2.0)) - 1

!-------------------------------------------------------
! then dive right in to the blocking algorithm
!-------------------------------------------------------

  print *, ' blocking transformation '
  print *, ' j,   sigma,   sigma dev '

  do j=1,nblock

    ndata = ndata/2
    do i=1,ndata
      blockdata(i) = (blockdata(2*i-1) + blockdata(2*i))/2.0
    enddo

   ! note here that xbar is constant under blocking transformation
   ! therefore not recalculated

    c0 = 0.0
    do i=1,ndata
      c0 = c0 + (blockdata(i) - xbar)**2
    enddo
    c0 = c0/ndata

    c0divn = c0/(ndata-1)

    sigma2 = c0divn
    sigma2dev = (sqrt(2.0/dble(ndata-1)))*c0divn

    sigma = sqrt(c0divn)
    sigmadev = (1.0/(sqrt(2*dble(ndata-1))))*sigma

!   write (*,*) 'blocking transformation number', j
!   write (*,*) 'sigma^2 =', sigma2, '+/-', sigma2dev
!   write (*,*) 'sigma =', sigma, '+/-', sigmadev

    print 10, j, sigma,sigmadev

  enddo

10  format(2x,i4,2x,2(2x,g12.5))

return

end subroutine
