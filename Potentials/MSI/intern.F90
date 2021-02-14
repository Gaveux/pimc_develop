!Computes the bond length coordinates and the first and second derivatives of these coordinates with respect
!to cartesians if these values are specified
subroutine intern(sys,x,r,dr,d2r,dr2dxdx)
    use molecule_specs
    type (molsysdat), intent(in) :: sys
    real(kind=8), dimension(sys%dimen,sys%natom), intent(in) :: x
    real(kind=8), dimension(sys%nbond), intent(out) :: r
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(out) :: dr
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom,sys%nbond), intent(out) :: d2r
    real(kind=8), dimension(sys%dimen,sys%dimen,sys%natom,sys%natom,sys%nbond), intent(out) :: dr2dxdx
    real(kind=8), dimension(sys%dimen,sys%dimen) :: dxdx

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
        enddo
    enddo

    dr2dxdx = 0.0
    do m=1,sys%nbond
       do k=1,sys%natom
          do l=1,sys%natom
             do j=1,sys%dimen
                do i=1,sys%dimen
                   dr2dxdx(j,i,k,l,m) = dr(j,k,m)*dr(i,l,m)
                enddo
             enddo
          enddo
       enddo
    enddo

    dxdx = 0.0
    do i=1,sys%dimen
       dxdx(i,i) = 1.0
    enddo

    d2r=0.0
    do k=1,sys%nbond
       do j=1,sys%dimen
          do i=1,sys%dimen
             d2r(j,i,sys%mb(k),sys%mb(k),k) = dxdx(j,i)/r(k) - &
                (x(j,sys%mb(k))-x(j,sys%nb(k)))*(x(i,sys%mb(k))-x(i,sys%nb(k)))/r(k)**3
             d2r(j,i,sys%nb(k),sys%mb(k),k) = -d2r(j,i,sys%mb(k),sys%mb(k),k)
             d2r(j,i,sys%mb(k),sys%nb(k),k) = -d2r(j,i,sys%mb(k),sys%mb(k),k)
             d2r(j,i,sys%nb(k),sys%nb(k),k) = d2r(j,i,sys%mb(k),sys%mb(k),k)
          enddo
       enddo
    enddo

    r=1/r
end subroutine intern
