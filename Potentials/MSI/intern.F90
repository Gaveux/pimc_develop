!Computes the bond length coordinates and the first and second derivatives of these coordinates with respect
!to cartesians if these values are specified
subroutine intern(sys,x,r,dr,dr2,d2rdxdx)
    use molecule_specs
    type (molsysdat), intent(in) :: sys
    real(kind=8), dimension(sys%dimen,sys%natom), intent(in) :: x
    real(kind=8), dimension(sys%nbond), intent(out) :: r
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(out) :: dr
    real(kind=8), dimension(sys%dimen,sys%dimen) :: dxkdxj ! a diagonal matrix serves for evaluating d2rdx2 matrix
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(out) :: dr2
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%nbond), intent(out) :: d2rdxdx

    integer :: i,j,k

    do i=1,sys%nbond
        r(i)=0.0d0
        do j=1,sys%dimen
            r(i) = r(i) + (x(j,sys%mb(i))-x(j,sys%nb(i)))**2
        enddo
        r(i) = sqrt(r(i))
    enddo

    dr=0.0d0
    
    do j=1,sys%dimen
        do i=1,sys%nbond
            dr(j,sys%mb(i),i)=(x(j,sys%mb(i))-x(j,sys%nb(i)))/r(i)
            dr(j,sys%nb(i),i)=-dr(j,sys%mb(i),i)

            !dr2(j,sys%mb(i),i) is (dr/dx_k)^2 and (dr/dx_j)^2
            dr2(j,sys%mb(i),i)=dr(j,sys%mb(i),i)**2
            !dr2(j,sys%nb(i),i) is (dr/dx_k)*(dr/dx_j) and (dr/dx_j)*(dr/dx_k)
            dr2(j,sys%nb(i),i)=-dr2(j,sys%mb(i),i)
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
