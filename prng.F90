 module prng
 
  implicit none
  ! see Numerical Recipes in Fortran, Ch7- Quick and Dirty Generators 
  contains 
   real*8 function genrand_real( ir )
     implicit real*8 (a-h,o-z)
 
     integer, intent(inout) :: ir
     parameter(da=16807.d0,db=2147483647.d0,dc=2147483648.d0)
     ir = abs(mod(da*ir,db)+0.5d0)
     genrand_real = dfloat(ir)/dc
     
     return 
   end function genrand_real
 
 end module prng
