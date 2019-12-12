!A quaternion rotation module. The euler angle used have the ZYZ convention 11 in http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
!This convention is the same one that is used in zare (defined on pages 77-81).
!All quaternions are normalised on creation

module quaternion_mod
    use vec_mod
    !use mt19937
    use prng
    use seed
    use moves

    implicit none

    real(kind=8),parameter :: twopi_quat = 8.D0*ATAN(1.D0)
    real(kind=8),parameter :: pi_quat = 4.D0*ATAN(1.D0)
    type quaternion
        real(kind=8), dimension(4) :: q
    end type quaternion

    !All of these normalise the quaternion
    interface new
        module procedure quaternion_empty
        module procedure quaternion_val
        module procedure quaternion_axisangle
        module procedure quaternion_euler
    end interface new

    interface operator(*)
        module procedure qq_mul
        module procedure vq_mul
    end interface

    interface assignment(=)
        module procedure qq_copy
    end interface

    interface rotate
        module procedure quaternion_rotate
    end interface rotate
    contains

    subroutine quaternion_empty(this)
        type(quaternion), intent(out) :: this
        this%q(1) = 1.d0
        this%q(2) = 0.d0
        this%q(3) = 0.d0
        this%q(4) = 0.d0
    end subroutine quaternion_empty

    subroutine quaternion_val(this,w,x,y,z)
        type(quaternion), intent(out) :: this
        real(kind=8), intent(in) :: w,x,y,z
        this%q(1) = w
        this%q(2) = x
        this%q(3) = y
        this%q(4) = z
        call normalise_quat(this)
    end subroutine quaternion_val

    subroutine quaternion_axisangle(this,v,omega)
        type(quaternion), intent(out) :: this
        real(kind=8), dimension(3), intent(in) :: v
        real(kind=8), intent(in) :: omega
        real(kind=8) :: cosw, sinw
        cosw = cos(omega*0.5)
        sinw = sin(omega*0.5)
        this%q(1) = cosw
        this%q(2) = sinw*v(1)
        this%q(3) = sinw*v(2)
        this%q(4) = sinw*v(3)
        call normalise_quat(this)
    end subroutine quaternion_axisangle

    !Creates a quaternion that represents the ZYZ euler angles (Zare pg 80).
    !The form of the quaternion is taken from http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
    subroutine quaternion_euler(this,a,b,c)
        type(quaternion), intent(out) :: this
        real(kind=8), intent(in) :: a,b,c
        real(kind=8) :: c1, c2, c3, s1, s2, s3
        c1 = cos(0.5*(a+c))
        c2 = cos(0.5*b)
        c3 = cos(0.5*(a-c))
        s1 = sin(0.5*(a+c))
        s2 = sin(0.5*b)
        s3 = sin(0.5*(a-c))
        this%q(1) =  c2*c1
        this%q(2) = -s2*s3
        this%q(3) =  s2*c3
        this%q(4) =  c2*s1
        call normalise_quat(this)
    end subroutine quaternion_euler

    !copies quaternion b to quaternion a
    subroutine qq_copy(a,b)
        type(quaternion), intent(out) :: a
        type(quaternion), intent(in) :: b
        a%q(1) = b%q(1)
        a%q(2) = b%q(2)
        a%q(3) = b%q(3)
        a%q(4) = b%q(4)
    end subroutine qq_copy


    !Returns the euler angle (ZYZ) representation of the quaternion.  This code deals with the singularities in the euler angles by
    !setting the chi variable to zero at these angles and setting phi to the entire angle spanned by chi+phi
    subroutine getEuler(t, phi,theta,chi)
        type(quaternion), intent(in) :: t
        real(kind=8), intent(out) :: phi,theta,chi
        real(kind=8) :: norm, test
        real(kind=8) :: tanphix, tanphiy, costheta, tanchix, tanchiy
        norm = length2_quat(t)
        tanphiy = 2*(t%q(3)*t%q(4)-t%q(1)*t%q(2))
        tanphix = 2*(t%q(2)*t%q(4) + t%q(1)*t%q(3))
        costheta = t%q(1)**2 - t%q(2)**2 - t%q(3)**2 + t%q(4)**2
        tanchiy = 2*(t%q(3)*t%q(4)+t%q(1)*t%q(2))
        tanchix = -2*(t%q(2)*t%q(4) - t%q(1)*t%q(3))
        !deals with singularities
        if(costheta >= 1.0) then
            phi =  2*atan2(t%q(4),t%q(1))
            theta = 0.0
            chi = 0.0
        elseif(costheta <= -1.0) then
            phi =  2*atan2(-t%q(3),t%q(2))
            theta = pi_quat
            chi = 0.0
        else
            phi = atan2(tanphiy,tanphix)
            theta = acos(costheta)
            chi = atan2(tanchiy,tanchix)
        endif
    end subroutine getEuler

    !normalises a quaternion
    subroutine normalise_quat(this)
        type(quaternion), intent(inout) :: this
        real(kind=8) :: norm
        integer :: i
        norm = sqrt(sum(this%q(:)**2))
        this%q = this%q/norm
    end subroutine normalise_quat

    function conj(this) result(co)
        type(quaternion), intent(in) :: this
        type(quaternion) :: co
        co%q(1) = this%q(1)
        co%q(2) = -this%q(2)
        co%q(3) = -this%q(3)
        co%q(4) = -this%q(4)
    end function conj

    !multiply two quaternions
    function qq_mul(a,b) result(res)
        type(quaternion), intent(in) :: a,b
        type(quaternion) :: res
        res%q(1) = a%q(1)*b%q(1) - a%q(2)*b%q(2) - a%q(3)*b%q(3) - a%q(4)*b%q(4)
        res%q(2) = a%q(1)*b%q(2) + a%q(2)*b%q(1) + a%q(3)*b%q(4) - a%q(4)*b%q(3)
        res%q(3) = a%q(1)*b%q(3) - a%q(2)*b%q(4) + a%q(3)*b%q(1) + a%q(4)*b%q(2)
        res%q(4) = a%q(1)*b%q(4) + a%q(2)*b%q(3) - a%q(3)*b%q(2) + a%q(4)*b%q(1)
    end function qq_mul

    !multiplies a vector * quaternion and returns the resultant quaternion
    function vq_mul(a,b) result(res)
        real(kind=8), dimension(3), intent(in) :: a
        type(quaternion), intent(in) :: b
        type(quaternion) :: res
        res%q(1) =  - a(1)*b%q(2) - a(2)*b%q(3) - a(3)*b%q(4)
        res%q(2) =  + a(1)*b%q(1) + a(2)*b%q(4) - a(3)*b%q(3)
        res%q(3) =  - a(1)*b%q(4) + a(2)*b%q(1) + a(3)*b%q(2)
        res%q(4) =  + a(1)*b%q(3) - a(2)*b%q(2) + a(3)*b%q(1)
    end function vq_mul

    !returns the vector that is obtained when a is rotated by b.  Assumes that all of the quaternions are normalised, a is the quaternion, b is the vector
    function quaternion_rotate(a, b) result(arr)
        type(quaternion), intent(inout) :: a
        real(kind=8), dimension(3), intent(in) :: b
        type(quaternion) :: co
        real(kind=8), dimension(3) :: arr
        type(quaternion) :: temp
        real(kind=8) :: leng
        leng = length2_quat(a)
        if(leng /= 1.0) then
            leng = sqrt(leng)
            a%q = a%q/leng
        endif
        co = conj(a)
        temp = qq_mul(a,vq_mul(b,co))
        arr(1) = temp%q(2)
        arr(2) = temp%q(3)
        arr(3) = temp%q(4)
    end function quaternion_rotate

    function length2_quat(this) result(norm)
        type(quaternion), intent(in) :: this
        real(kind=8) :: norm
        integer :: i
        norm = 0.0
        do i=1,4
            norm = norm + this%q(i)**2
        enddo 
    end function length2_quat

    !gets the quaternion that will rotate a onto b
    function get_rotation_between(a,b) result(quat)
        real(kind=8), dimension(3), intent(in) :: a,b
        real(kind=8), dimension(3) :: other
        type(quaternion) :: quat
        real(kind=8) :: cos_theta,k, temp

        cos_theta = dot_product(a,b)
        !print *, cos_theta
        !call exit(1)
        k = sqrt(dot_product(a,a)*dot_product(b,b))
        !print *, 'cos(theta) = ', cos_theta/k
        if(cos_theta/k == -1) then
            temp = abs(dot_product(a,[1.0,0.0,0.0]))
            if(temp < 1.0) then
                other=[1.0,0.0,0.0]
            else
                other = [0.0,1.0,0.0]
            endif
            call quaternion_axisangle(quat, normalise(cross(a,other)),pi_quat)
        endif
        other = cross(a,b)
        call quaternion_val(quat,cos_theta+k,other(1), other(2), other(3)) 
        call normalise_quat(quat)
    end function get_rotation_between

    !generates a random quaternion that uniformly samples SO(3)
    !https://www-preview.ri.cmu.edu/pub_files/pub4/kuffner_james_2004_1/kuffner_james_2004_1.pdf
    function random_quaternion(dang,seedval) result(quat)
        real(kind=8), intent(in) :: dang
        type(mod_seed), intent(inout) :: seedval
        type(quaternion) :: quat
        real(kind=8) :: x1, x2, x3, xnorm, ang
        real(kind=8), dimension(3) :: axis
        x1 = gauss_dev(seedval)
        x2 = gauss_dev(seedval)
        x3 = gauss_dev(seedval)
        do while (abs(x1) < 1e-14 .and. abs(x2) < 1e-14 .and. abs(x3) < 1e-14) 
            x1 = gauss_dev(seedval)
            x2 = gauss_dev(seedval)
            x3 = gauss_dev(seedval)
        enddo
        xnorm = 1.0/sqrt(x1**2 + x2**2 + x3**2)
        x1 = x1*xnorm
        x2 = x2*xnorm
        x3 = x3*xnorm
        axis = [x1,x2,x3]
        !ang = genrand_real3()*dang
        ang = genrand_real(seedval%seedvalue)*dang
        call new(quat,axis,ang)
    end function random_quaternion
end module quaternion_mod
