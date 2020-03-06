
      subroutine normod(r,nbond)
         implicit none
         
         integer :: nbond, i
         real(kind=8), dimension(nbond) :: r
         real(kind=8) :: s, as1, as2, as3, as4

         !print *, nbond
         !print *, r

         ! calculate symmetric vibrational normal mode
         ! r1+r2+r3+r4
         s = 0.d0
         as1 = 0.d0
         as2 = 0.d0
         as3 = 0.d0
         as4 = 0.d0

         do i=1,4
           s = s + r(i)
         enddo
         
         ! non-symmetrical vibrational normal mode
         as1 = 3.0*r(1) - r(2) - r(3) - r(4)
         as2 = 3.0*r(2) - r(1) - r(3) - r(4)
         as3 = 3.0*r(3) - r(1) - r(2) - r(4)
         as4 = 3.0*r(4) - r(1) - r(2) - r(3)

         write(149,*) s
         write(159,*) as1
         write(169,*) as2
         write(179,*) as3
         write(189,*) as4
                  
      end subroutine
