module vec_mod

implicit none

contains

FUNCTION cross(a, b)
  real(kind=8), DIMENSION(3) :: cross
  real(kind=8), DIMENSION(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)

END FUNCTION cross

function dot(a, b)
  real(kind=8) :: dot
  real(kind=8), dimension(:), intent(in) :: a, b

  dot = sum(a(:)*b(:))

end function dot

function normalise(v)
  real(kind=8), dimension(3) :: normalise
  real(kind=8), dimension(:), intent(in)  :: v
  real(kind=8) :: dist

  dist = sqrt( sum(v**2) )
  normalise(:) = v(:)/dist

end function normalise

real(kind=8) function length(vec)
      implicit none
      real(kind=8), dimension(:), intent(in) :: vec
      integer :: i

      length = 0.d0
      do i=1,3
        length = length + vec(i)**2
      enddo
      length = sqrt(length)
      return
end function length


end module vec_mod



