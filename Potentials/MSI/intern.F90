!Computes the bond length coordinates and the first and second derivatives of these coordinates with respect
!to cartesians if these values are specified
subroutine intern(sys,x,r,dr,d2rdxdx,dr2)
    use molecule_specs
    type (molsysdat), intent(in) :: sys
    real(kind=8), dimension(sys%dimen,sys%natom), intent(in) :: x
    real(kind=8), dimension(sys%nbond), intent(out) :: r
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(out) :: dr  ! dr/dx
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(out) :: dr2 !(dr/dx)^2
    real(kind=8), dimension(sys%dimen,sys%dimen) :: dxkdxj ! a diagonal matrix for calculating d2rdx2 matrix
    !real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%nbond), intent(out) :: d2r !d2r/dxdx
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%nbond), intent(out) :: d2rdxdx !d2r/dxdx

    integer :: i,j,k

    do i=1,sys%nbond
        r(i)=0.0d0
        do j=1,sys%dimen
            r(i) = r(i) + (x(j,sys%mb(i))-x(j,sys%nb(i)))**2
        enddo
        r(i) = 1.0/sqrt(r(i))
    enddo

    dr=0.0d0
    dr2=0.0d0
    
    do j=1,sys%dimen
        do i=1,sys%nbond
            dr(j,sys%mb(i),i)=(x(j,sys%mb(i))-x(j,sys%nb(i)))*r(i)
            dr2(j,sys%mb(i),i)=dr(j,sys%mb(i),i)**2
            dr(j,sys%nb(i),i)=-dr(j,sys%mb(i),i)
            dr2(j,sys%nb(i),i)=dr(j,sys%nb(i),i)**2
        enddo
    enddo

    ! d2r/dx2
    dxkdxj=0.0 ! identity matrix
    do i=1,sys%dimen
       dxkdxj(i,i)=1.0
    enddo
    ! for d2r/dx2 matrix:
    ! sys%mb(i) is d2r/dx_m*dx_m = d2r/dx_n*dx_n
    ! sys%nb(i) is d2r/dx_m*dx_n = d2r/dx_n*dx_m
    d2rdxdx=0.0
    do k=1,sys%dimen
       do j=1,sys%dimen
          do i=1,sys%nbond 
             d2rdxdx(k,j,sys%mb(i),i)=r(i)*dxkdxj(k,j) - (x(j,sys%mb(i))-x(j,sys%nb(i)))*(x(k,sys%mb(i))-x(k,sys%nb(i)))*r(i)**3 
             d2rdxdx(k,j,sys%nb(i),i)=-d2rdxdx(k,j,sys%mb(i),i)
          enddo
       enddo
    enddo

    r=1/r
end subroutine intern
