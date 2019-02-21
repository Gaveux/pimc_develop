!Computes the bond length coordinates and the first and second derivatives of these coordinates with respect
!to cartesians if these values are specified
subroutine intern(sys,x,r,dr)
    use molecule_specs
    type (molsysdat), intent(in) :: sys
    real(kind=8), dimension(sys%dimen,sys%natom), intent(in) :: x
    real(kind=8), dimension(sys%nbond), intent(out) :: r
    real(kind=8), dimension(sys%dimen,sys%natom,sys%nbond), intent(out) :: dr

    integer :: i,j

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
    r=1/r
end subroutine intern
