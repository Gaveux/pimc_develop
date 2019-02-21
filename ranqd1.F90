!A quick and dirty pseudo random number generator
!see numerical recipes in Fortran 77, ed. 2. Page 275~276
 module ranqd1
   implicit none
   
   contains

   integer idum, itemp, jflone, jflmsk

   real ftemp
   equivalence (itemp,ftemp)
   data jflone /Z'3F800000'/, jflmsk /Z'007FFFFF'/

   idum=1664525*idum+1013904223
   itemp=ior(jflone,iand(jflmsk,idum))
   ran=ftemp-1.0

 end module ranqd1
