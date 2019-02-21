MODULE POTENTIAL_NH3
#if POT==2
    use molecule_specs
      
    CONTAINS
       

    !Subroutine for calculating the value of the potential of methane with atoms at cartesian positions
    !given in X with order N,H,H,H
    SUBROUTINE Potential(sys,x,r,V,dV)
      implicit none
      type (molsysdat), intent(in) :: sys
      REAL(KIND=8), DIMENSION(sys%dimen,sys%natom), INTENT(IN) :: x
      REAL(KIND=8), INTENT(OUT) :: V
      REAL(KIND=8), DIMENSION(1,sys%nbond), intent(out) :: r
      real(kind=8), dimension(sys%dimen,sys%natom), intent(out) :: dV

      REAL(KIND=8), DIMENSION(1) :: VTEMP
      real(kind=8) :: step
      step = 1.0d-5
      
      CALL INTERN(sys,x,r)
      CALL VNH3(r,VTEMP)
      CALL FINITEDIFFERENCEDERIVS(x,sys,step,dV)
      V=VTEMP(1)

    END SUBROUTINE Potential

    !Converts from cartesian coordinates to internal coordinates used by VNH3
    SUBROUTINE INTERN(sys,x,r)
      implicit none
      type (molsysdat), intent(in) :: sys
      REAL(KIND=8), DIMENSION(sys%dimen,sys%natom), INTENT(IN) :: x
      REAL(KIND=8), DIMENSION(1,sys%nbond), INTENT(OUT) :: r
      REAL(KIND=8), DIMENSION(sys%nbond) :: bl
      INTEGER :: i,j,k
      !Numbering for the bonds to make everything easier
      do i=1,sys%nbond
        bl(i)=0
        do j=1,sys%dimen
          bl(i) = bl(i) + (x(j,sys%mb(i))-x(j,sys%nb(i)))**2
        enddo
        bl(i) = sqrt(bl(i))
      enddo

      !N-H distances (in angstrom)
      r(1,1:3)=bl(1:3)*0.529177249

      !H-N-H angles (in degrees) - number out the front is 180/pi
      r(1,4) = 57.2957795131d0*acos((bl(1)**2+bl(2)**2-bl(4)**2)/(2.0d0*bl(1)*bl(2))) !H1, H2
      r(1,5) = 57.2957795131d0*acos((bl(1)**2+bl(3)**2-bl(5)**2)/(2.0d0*bl(1)*bl(3))) !H1, H3
      r(1,6) = 57.2957795131d0*acos((bl(2)**2+bl(3)**2-bl(6)**2)/(2.0d0*bl(2)*bl(3))) !H2, H3

      RETURN
    END SUBROUTINE INTERN

    !Compute the finite difference derivatives of the ammonia potential
    SUBROUTINE FINITEDIFFERENCEDERIVS(x,sys,step,grad)
      implicit none
      type (molsysdat), intent(in) :: sys
      REAL(KIND=8), DIMENSION(sys%dimen,sys%natom), INTENT(IN) :: x
      REAL(KIND=8), DIMENSION(sys%dimen,sys%natom), INTENT(OUT) :: grad
      REAL(KIND=8), INTENT(IN) :: step
      REAL(KIND=8), DIMENSION(sys%dimen,sys%natom) :: temp_x
      REAL(KIND=8), DIMENSION(1,sys%nbond) :: RINT
      REAL(KIND=8), DIMENSION(1) :: V1, V2
      INTEGER :: i,j
      do i=1,sys%natom
        do j=1,sys%dimen
          temp_x(j,i)=x(j,i)
        enddo
      enddo

      do i=1,sys%natom
        do j=1,sys%dimen
          temp_x(j,i) = x(j,i) + step
          CALL INTERN(sys,temp_x,RINT)
          CALL VNH3(RINT,V1)
          temp_x(j,i) = x(j,i) - step
          CALL INTERN(sys,temp_x,RINT)
          CALL VNH3(RINT,V2)
          temp_x(j,i) = x(j,i)
          grad(j,i) = (V2(1)-V1(1))/(2.0d0*step)

        enddo
      enddo
      RETURN
    END SUBROUTINE FINITEDIFFERENCEDERIVS


! AMMPOT4
!
! Library with routines for the definition of parameters and sections of the
! potential surface of ammonia and products.
!
! This program was adapted from a program written by
!          R. Marquardt
!          Laboratorium fuer Physikalische Chemie
!          ETH Zuerich (Zentrum)
!          8092 Zuerich - Switzerland
!
! Coordinates: R1, R2, R3, A12, A13, A23 (the 3 bond lengths and angles)
!
! Units are: angstroem for ri, and degrees for aij.
! 
! Evaluation of the PES at a given set of these coordinates is performed with
!
!    CALL VNH3(RINT,V)
!
! where RINT is a real*8 vector of length 6 containing the
! aforementioned coordinates in the given units. 
!
! The potential energy V is returned in the corresponding 
! units of a wavenumber (hc cm-1).
!
! Before any call to subroutine VNH3, subroutine VINI must be called once,
! to initiallize parameters.
!
!-------------------------------------------------------------------------------
      SUBROUTINE VINI
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      PARAMETER(NINT=400,NEXT=200)
      CHARACTER(len=6) LPINT(NINT)
      COMMON /CP/ P(NEXT)
      COMMON /CPINT/ PINT(NINT)
      COMMON /CLINT/ LPINT
! version
      LPINT(NINT)='VERS. '
      PINT(NINT)=AMMPOT4
! assignment of model parameters
      CALL POTPAR(PINT,NINT)
      RETURN
      END
!-----------------------------------------------------------------------

      SUBROUTINE VNH3(RINT,V)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      PARAMETER(N=1,NI=3)         
      PARAMETER(MI=NI-1,NNI=NI*(NI-1)/2,NAT=NI+1,NDIM=3*NAT-6)
      PARAMETER(NORD=4,NABORD=10)
! conversion parameters to reciprocal centimeters
      PARAMETER(AJOULE=50341.12503,QJMOL=83.59346114)
! nominal parameters
! units are aJ=10**(-18)J for energies and A=100pm for lengths
      COMMON /PAR0/   V0(1:NI)
      COMMON /PARS/   AR(1:NI),RE(1:NI),FS(1:NI),AS(1:NI),BS(1:NI),     &
     &                E6(1:NI),E8(1:NI),RS(1:NI)
      COMMON /PARB/   AC(1:MI),CE(1:MI),FB1(1:MI),FB2(1:MI),            &
     &                AB1(1:NORD,1:NABORD,1:MI),                        &    
     &                AB2(1:NORD,1:NABORD,1:MI),                        & 
     &                N1(1:NORD),N2(1:NORD)
      COMMON /PARD/   DL1C(1:NORD,1:MI),DL1O(1:NORD,1:MI),              &
     &                DL2C(1:NORD,1:MI),DL2O(1:NORD,1:MI),              &
     &                DX1C(1:NORD,1:MI),DX1O(1:NORD,1:MI),              &
     &                DX2C(1:NORD,1:MI),DX2O(1:NORD,1:MI),              &
     &                DLPC(1:MI),DLPO(1:MI),                            &
     &                DXPC(1:MI),DXPO(1:MI),                            &
     &                ND1,ND2,NDCOS
      COMMON /PARP/   DII(1:MI),AII(1:MI),                              &
     &                RIII,AIII,DIII        
      COMMON /PARSW/  RSWSB,RSWDB,NSWSB,NSWDB
      COMMON /DECOMP/ IFDEC
! local parameters
      DIMENSION V(N),VS(N),VB(N),VP(N)
      DIMENSION RINT(N,NDIM) 
      DIMENSION R(N,NI),DR(N,NI),R2(N,NI)
      DIMENSION C(N,NNI),DC(N,NNI)
      DIMENSION D(N,NNI),DD(N,NNI)
      DIMENSION SPSR(N,NI),SPDR(N,NI),SPDD(N,NNI)
      DIMENSION SQSR(N,NI),SQDR(N,NI),SQDD(N,NNI)
      DIMENSION V0J(N),ARJ(N),REJ(N),ACJ(N),CEJ(N),FB1J(N),FB2J(N)
      DIMENSION FSJ(N,NI),ASJ(N,NI),BSJ(N,NI),E6J(N,NI),E8J(N,NI),      &
     &          RSJ(N,NI)
      DIMENSION DIIJ(N,NNI),AIIJ(N,NNI)           
      DIMENSION DLPCJ(N,NNI),DLPOJ(N,NNI),DLPJ(N,NNI)
      DIMENSION DXPCJ(N,NNI),DXPOJ(N,NNI),DXPJ(N,NNI)
      DIMENSION DL1CJ(N,NNI,NORD),DL1OJ(N,NNI,NORD),DL1J(N,NNI,NORD)
      DIMENSION DX1CJ(N,NNI,NORD),DX1OJ(N,NNI,NORD),DX1J(N,NNI,NORD)
      DIMENSION DL2CJ(N,NNI,NORD),DL2OJ(N,NNI,NORD),DL2J(N,NNI,NORD)
      DIMENSION DX2CJ(N,NNI,NORD),DX2OJ(N,NNI,NORD),DX2J(N,NNI,NORD)
      DIMENSION AB1J(N,NORD,NABORD)
      DIMENSION AB2J(N,NORD,NABORD)
      DIMENSION SP6(N,NI),SP8(N,NI),SPRE(N,NI),SQRE(N,NI),              &
     &          SPCE(N,NI),SQCE(N,NI)
      DIMENSION REQJ(N,NI),CEQJ(N,NNI)
      DIMENSION Y(N,NI)
      DIMENSION YD1(N,NNI,NORD),X1(N,NNI,NORD)
      DIMENSION YD2(N,NNI,NORD),X2(N,NNI,NORD)
      DIMENSION YDP(N,NNI)
      DIMENSION SB1(N,NORD),SB2A(N,NORD),SB2B(N,NORD)
      DIMENSION SBA1(N,NORD,NABORD),SBEa(N,NORD,NABORD),                &
     &          SBEb(N,NORD,NABORD)
      DIMENSION SBgA1(N),SBgEa(N),SBgEb(N)
      DIMENSION RJJ(N,NNI),AJJ(N,NNI),DJJ(N,NNI)
      DIMENSION NA(NORD),NE(NORD)
      DIMENSION D3F(N,NNI)
! the number pi
      PI=ACOS(-1.)
      URF=PI/180.
! map internal coordinates
      DO J=1,N
         R(J,1)=RINT(J,1)
         R(J,2)=RINT(J,2)
         R(J,3)=RINT(J,3)
         C(J,1)=COS(RINT(J,4)*URF)
         C(J,2)=COS(RINT(J,5)*URF)
         C(J,3)=COS(RINT(J,6)*URF)
      ENDDO
      DO K=1,NI
         DO J=1,N
            R2(J,K)=R(J,K)**2
         ENDDO
      ENDDO
      KK=0
      DO K=1,NI
         DO L=(K+1),NI
            KK=KK+1
            DO J=1,N
               D(J,KK)=SQRT(R2(J,K)+R2(J,L)-2.*R(J,K)*R(J,L)*C(J,KK))
            ENDDO
         ENDDO
      ENDDO
! evaluate generic one dimensional switching functions
      DO K=1,NI
         DO J=1,N
            SPSR(J,K)=SP(NSWSB,RSWSB,R(J,K))
            SPDR(J,K)=SP(NSWDB,RSWDB,R(J,K))
            SQSR(J,K)=1.-SPSR(J,K)
            SQDR(J,K)=1.-SPDR(J,K)
         ENDDO
      ENDDO
      DO KK=1,NNI
         DO J=1,N
            SPDD(J,KK)=SP(NSWDB,RSWDB,D(J,KK))
            SQDD(J,KK)=1.-SPDD(J,KK)
         ENDDO
      ENDDO
! multidimensional switching
      CALL SW1G3(N,NI,NI,SPSR,SQSR,V0,V0J)
      CALL SW1G3(N,NI,NI,SPSR,SQSR,AR,ARJ)
      CALL SW1G3(N,NI,NI,SPSR,SQSR,RE,REJ)
      CALL SW2G3(N,NI,MI,SPSR,SQSR,AC,ACJ)
      CALL SW2G3(N,NI,MI,SPSR,SQSR,CE,CEJ)
      CALL SW2G3(N,NI,MI,SPSR,SQSR,FB1,FB1J)
      CALL SW2G3(N,NI,MI,SPSR,SQSR,FB2,FB2J)
      CALL SW1L3(N,NI,NI,NI,SPSR,SQSR,FS,FSJ)
      CALL SW1L3(N,NI,NI,NI,SPSR,SQSR,AS,ASJ)
      CALL SW1L3(N,NI,NI,NI,SPSR,SQSR,BS,BSJ)
      CALL SW1L3(N,NI,NI,NI,SPSR,SQSR,E6,E6J)
      CALL SW1L3(N,NI,NI,NI,SPSR,SQSR,E8,E8J)
      CALL SW1L3(N,NI,NI,NI,SPSR,SQSR,RS,RSJ)
      CALL SW2L3(N,NI,MI,NNI,SPSR,SQSR,DII,DIIJ)
      CALL SW2L3(N,NI,MI,NNI,SPSR,SQSR,AII,AIIJ)
      CALL SW2L3(N,NI,MI,NNI,SPSR,SQSR,DLPC,DLPCJ)
      CALL SW2L3(N,NI,MI,NNI,SPSR,SQSR,DLPO,DLPOJ)
      CALL SW2L3(N,NI,MI,NNI,SPSR,SQSR,DXPC,DXPCJ)
      CALL SW2L3(N,NI,MI,NNI,SPSR,SQSR,DXPO,DXPOJ)
      CALL SP2L3(N,NI,MI,NNI,NORD,SPSR,SQSR,DL1C,DL1CJ)
      CALL SP2L3(N,NI,MI,NNI,NORD,SPSR,SQSR,DL1O,DL1OJ)
      CALL SP2L3(N,NI,MI,NNI,NORD,SPSR,SQSR,DL2C,DL2CJ)
      CALL SP2L3(N,NI,MI,NNI,NORD,SPSR,SQSR,DL2O,DL2OJ)
      CALL SP2L3(N,NI,MI,NNI,NORD,SPSR,SQSR,DX1C,DX1CJ)
      CALL SP2L3(N,NI,MI,NNI,NORD,SPSR,SQSR,DX1O,DX1OJ)
      CALL SP2L3(N,NI,MI,NNI,NORD,SPSR,SQSR,DX2C,DX2CJ)
      CALL SP2L3(N,NI,MI,NNI,NORD,SPSR,SQSR,DX2O,DX2OJ)
      CALL SP2G3(N,NI,MI,NORD,NABORD,SPSR,SQSR,AB1,AB1J)
      CALL SP2G3(N,NI,MI,NORD,NABORD,SPSR,SQSR,AB2,AB2J)
! angular dependency of damping parameters
      CALL DANG(N,NNI,NORD,NDCOS,C,DL1CJ,DL1OJ,DL1J)
      CALL DANG(N,NNI,NORD,NDCOS,C,DX1CJ,DX1OJ,DX1J)
      CALL DANG(N,NNI,NORD,NDCOS,C,DL2CJ,DL2OJ,DL2J)
      CALL DANG(N,NNI,NORD,NDCOS,C,DX2CJ,DX2OJ,DX2J)
      CALL DANG(N,NNI,1,NDCOS,C,DLPCJ,DLPOJ,DLPJ)
      CALL DANG(N,NNI,1,NDCOS,C,DXPCJ,DXPOJ,DXPJ)
!
      DO K=1,NI
         DO J=1,N
            SP6(J,K)=SP(6,RSJ(J,K),R(J,K))
            SP8(J,K)=SP(8,RSJ(J,K),R(J,K))
            SPRE(J,K)=MTANH(ARJ(J)*(R(J,K)-REJ(J)))
            SQRE(J,K)=1.-SPRE(J,K)
            SPCE(J,K)=MTANH(ACJ(J)*(R(J,K)-REJ(J)))
            SQCE(J,K)=1.-SPCE(J,K)
         ENDDO
      ENDDO
! equilibrium structure
      CALL SW1L3(N,NI,NI,NI,SPRE,SQRE,RE,REQJ)
      CALL SW2L3(N,NI,MI,NNI,SPCE,SQCE,CE,CEQJ)
! stretching coordinates
      DO K=1,NI
         DO J=1,N
            DR(J,K)=R(J,K)-REQJ(J,K)
            Y(J,K)=(1.+EXP(-BSJ(J,K)*DR(J,K)))/2.
            Y(J,K)=(1.-EXP(-ASJ(J,K)*Y(J,K)*DR(J,K)))/ASJ(J,K)          &       
     &             *(1.+E6J(J,K)*SP6(J,K)+E8J(J,K)*SP8(J,K))
         ENDDO
      ENDDO
! damping functions
      CALL DAMPG(N,NI,NNI,NORD,R,REQJ,C,DX1J,DL1J,ND1,ND2,YD1)
      CALL DAMPG(N,NI,NNI,NORD,R,REQJ,C,DX2J,DL2J,ND1,ND2,YD2)
      CALL DAMPG(N,NI,NNI,1   ,R,REQJ,C,DXPJ,DLPJ,ND1,ND2,YDP)
! bending coordinates
      DO KK=1,NNI
         DO J=1,N
            DC(J,KK)= C(J,KK)-CEQJ(J,KK)
         ENDDO
      ENDDO
      DO NO=1,NORD
         DO KK=1,NNI
            DO J=1,N
               X1(J,KK,NO) = DC(J,KK)*YD1(J,KK,NO)
               X2(J,KK,NO) = DC(J,KK)*YD2(J,KK,NO)
            ENDDO
         ENDDO
      ENDDO
! symmetrized anharmonic bending coordinates xy3 type molecules
      DO NO=1,NORD
         DO J=1,N
            sb1(J,NO)=(X1(J,1,NO)+X1(J,2,NO)+X1(J,3,NO))/sqrt(3.)
            sb2a(J,NO)=(2*X2(J,1,NO)-(X2(J,2,NO)+X2(J,3,NO)))/sqrt(6.)
            sb2b(J,NO)=(X2(J,2,NO)-X2(J,3,NO))/sqrt(2.)
         ENDDO
      ENDDO
      NA(1)=1
      NA(2)=2
      NA(3)=3
      NA(4)=4
      NE(1)=1
      NE(2)=2
      NE(3)=3
      NE(4)=5
      DO J=1,N
         SBA1(J,4,4) = (sb2a(J,4)**2+sb2b(J,4)**2)**2
         SBA1(J,4,3) = (sb2a(J,4)**2-3*sb2b(J,4)**2)*sb2a(J,4)*sb1(J,4)
         SBA1(J,4,2) = (sb2a(J,4)**2+sb2b(J,4)**2)*sb1(J,4)**2
         SBA1(J,4,1) = sb1(J,4)**4
         SBA1(J,3,3) = sb2a(J,3)**3-3*sb2a(J,3)*sb2b(J,3)**2
         SBA1(J,3,2) = (sb2a(J,3)**2+sb2b(J,3)**2)*sb1(J,3)
         SBA1(J,3,1) = sb1(J,3)**3
         SBA1(J,2,2) = sb2a(J,2)**2+sb2b(J,2)**2
         SBA1(J,2,1) = sb1(J,2)**2
         SBA1(J,1,1) = sb1(J,1)
         SBEa(J,4,5) = sb2a(J,4)**4-6*sb2a(J,4)**2*sb2b(J,4)**2                 &
     &        +sb2b(J,4)**4
         SBEa(J,4,4) = sb2a(J,4)**4-sb2b(J,4)**4
         SBEa(J,4,3) = (sb2a(J,4)**2+sb2b(J,4)**2)*sb2a(J,4)*sb1(J,4)
         SBEa(J,4,2) = (sb2a(J,4)**2-sb2b(J,4)**2)*sb1(J,4)**2/2
         SBEa(J,4,1) = sb2a(J,4)*sb1(J,4)**3
         SBEa(J,3,3) = (sb2a(J,3)**2+sb2b(J,3)**2)*sb2a(J,3)
         SBEa(J,3,2) = (sb2a(J,3)**2-sb2b(J,3)**2)*sb1(J,3)/2
         SBEa(J,3,1) = sb2a(J,3)*sb1(J,3)**2
         SBEa(J,2,2) = sb2a(J,2)**2/2-sb2b(J,2)**2/2
         SBEa(J,2,1) = sb2a(J,2)*sb1(J,2)
         SBEa(J,1,1) = sb2a(J,1)
         SBEb(J,4,5) = 4*sb2a(J,4)**3*sb2b(J,4)-4*sb2a(J,4)*sb2b(J,4)**3
         SBEb(J,4,4) =-2*sb2a(J,4)**3*sb2b(J,4)-2*sb2a(J,4)*sb2b(J,4)**3
         SBEb(J,4,3) = (sb2a(J,4)**2+sb2b(J,4)**2)*sb2b(J,4)*sb1(J,4)
         SBEb(J,4,2) = -sb2a(J,4)*sb2b(J,4)*sb1(J,4)**2
         SBEb(J,4,1) = sb2b(J,4)*sb1(J,4)**3
         SBEb(J,3,3) = (sb2a(J,3)**2+sb2b(J,3)**2)*sb2b(J,3)
         SBEb(J,3,2) = -sb2a(J,3)*sb2b(J,3)*sb1(J,3)
         SBEb(J,3,1) = sb2b(J,3)*sb1(J,3)**2
         SBEb(J,2,2) = -sb2a(J,2)*sb2b(J,2)
         SBEb(J,2,1) = sb2b(J,2)*sb1(J,2)
         SBEb(J,1,1) = sb2b(J,1)
         SBgA1(J)=AB1J(J,1,1)*SBA1(J,1,1)
         SBgEa(J)=AB2J(J,1,1)*SBEa(J,1,1)
         SBgEb(J)=AB2J(J,1,1)*SBEb(J,1,1)
      ENDDO
      DO NO=2,NORD
         DO L=1,NA(NO)
            DO J=1,N
               SBgA1(J)=SBgA1(J)+AB1J(J,NO,L)*SBA1(J,NO,L)
            ENDDO
         ENDDO
         DO M=1,NE(NO)
            DO J=1,N
               SBgEa(J)=SBgEa(J)+AB2J(J,NO,M)*SBEa(J,NO,M)
               SBgEb(J)=SBgEb(J)+AB2J(J,NO,M)*SBEb(J,NO,M)
            ENDDO
         ENDDO
      ENDDO
! calculate potential: initialization
      DO J=1,N
         VS(J) = 0.0
         VB(J) = 0.0   
         VP(J) = 0.0   
      ENDDO
! stretching potential
      DO K=1,NI
         DO J=1,N
            VS(J) = VS(J)+0.5*FSJ(J,K)*Y(J,K)**2
         ENDDO
      ENDDO
! bending potential
      DO J=1,N
         VB(J) = VB(J) + FB1J(J)*SBgA1(J)**2
         VB(J) = VB(J) + FB2J(J)*(SBgEa(J)**2+SBgEb(J)**2)
      ENDDO
! asymptotic values for the pair potential
      KK=0
      DO K=1,NI
         DO L=(K+1),NI
            KK=KK+1
            DO J=1,N
               RJJ(J,KK)=SQRT(REQJ(J,K)**2+REQJ(J,L)**2                     &
     &            -2.*REQJ(J,K)*REQJ(J,L)*CEQJ(J,KK))                       &
     &                  *(1.-SPDR(J,K)*SPDR(J,L))                           &
     &          +RIII*      (SPDR(J,K)*SPDR(J,L))
               AJJ(J,KK)=AIIJ(J,KK)*(1.-SPDR(J,K)*SPDR(J,L))                &
     &          +AIII*      (SPDR(J,K)*SPDR(J,L))
               DJJ(J,KK)=DIIJ(J,KK)*YDP(J,KK)                            
               DD(J,KK)=D(J,KK)-RJJ(J,KK)
            ENDDO
         ENDDO
      ENDDO
! pair potential
      DO KK=1,NNI
         DO J=1,N
            VP(J)=VP(J)+DJJ(J,KK)*(1.-EXP(-AJJ(J,KK)*DD(J,KK)))**2
         ENDDO
      ENDDO
! substract DIII binding energy at the rate of increase of DJJ
! to avoid 3-fold clustering, substraction is inhibited whenever
! the distances of more than two YY bonds become too small
      KK=0
      DO K=1,NI
         DO L=(K+1),NI
            KK=KK+1
            DO J=1,N
               D3F(J,KK)=SPDR(J,K)*SPDR(J,L)
            ENDDO
            DO MM=1,NNI
               IF (MM.NE.KK) THEN 
                  DO J=1,N
                     D3F(J,KK)=D3F(J,KK)*SPDD(J,MM)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO               
      DO KK=1,NNI
         DO J=1,N
            VP(J)=VP(J)+DIII*((1.-EXP(-AJJ(J,KK)*DD(J,KK)))**2-1)           &
     &           *D3F(J,KK)
         ENDDO
      ENDDO
! energy conversion and sum
      DO J=1,N
         VS(J)=VS(J)*AJOULE
         VB(J)=VB(J)*AJOULE
         VP(J)=VP(J)*AJOULE
         V(J) =V0J(J)*AJOULE
      ENDDO                  
      IF     (IFDEC.EQ.1) THEN
         DO J=1,N
            V(J)=V(J)+VS(J)
         ENDDO
      ELSEIF (IFDEC.EQ.2) THEN
         DO J=1,N
            V(J)=V(J)+VB(J)
         ENDDO
      ELSEIF (IFDEC.EQ.3) THEN
         DO J=1,N
            V(J)=V(J)+VP(J)
         ENDDO
      ELSE
         DO J=1,N
            V(J)=V(J)+VS(J)+VB(J)+VP(J)
         ENDDO
      ENDIF
      RETURN
      END
!--------------------------------------------------------------------------
      REAL(kind=8) FUNCTION SP(N,S,R)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      PARAMETER(ONE=1.D0)
      PARAMETER(STHRE=200.D0,SMALL=1.0D-06)
      SP=0.
      IF(N.LE.0) STOP ' Nsw.LE.0'
      IF(R.LE.0.) RETURN
      IF((S/R)**N.GT.STHRE) THEN
       SP=0.
      ELSEIF((S/R)**N.LT.SMALL) THEN
       SP=1.
      ELSE
       SP=EXP(-(S/R)**N)
      ENDIF
      RETURN
      END
!--------------------------------------------------------------------------
      SUBROUTINE SW1G3(N,NB,NL,SP,SQ,P,PJ)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      DIMENSION P(NL),PJ(N),SP(N,NB),SQ(N,NB)
      DO J=1,N
         PJ(J)=P(1)*(SQ(J,1)*SQ(J,2)*SQ(J,3))                               &
     &        +P(2)*(SP(J,1)*SQ(J,2)*SQ(J,3)+                               &
     &               SQ(J,1)*SP(J,2)*SQ(J,3)+                               &
     &               SQ(J,1)*SQ(J,2)*SP(J,3))                               &
     &        +P(3)*(SP(J,1)*SP(J,2)*SQ(J,3)+                               &
     &               SP(J,1)*SQ(J,2)*SP(J,3)+                               &
     &               SQ(J,1)*SP(J,2)*SP(J,3))
      ENDDO
      RETURN
      END
!--------------------------------------------------------------------------
      SUBROUTINE SW2G3(N,NB,NL,SP,SQ,P,PJ)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      DIMENSION P(NL),PJ(N),SP(N,NB),SQ(N,NB)
      DO J=1,N
         PJ(J)=P(1)*(SQ(J,1)*SQ(J,2)*SQ(J,3))                               &
     &        +P(2)*(SP(J,1)*SQ(J,2)*SQ(J,3)+                               &
     &               SQ(J,1)*SP(J,2)*SQ(J,3)+                               &
     &               SQ(J,1)*SQ(J,2)*SP(J,3))                               &
     &        +P(2)*(SP(J,1)*SP(J,2)*SQ(J,3)+                               &
     &               SP(J,1)*SQ(J,2)*SP(J,3)+                               &
     &               SQ(J,1)*SP(J,2)*SP(J,3))
      ENDDO
      RETURN
      END
!--------------------------------------------------------------------------------------------------
      SUBROUTINE SW1L3(N,NB,NL,NPJ,SP,SQ,P,PJ)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      DIMENSION P(NL),PJ(N,NB),SP(N,NB),SQ(N,NB)
      DO J=1,N
         PJ(J,1)=P(1)*(SQ(J,2)*SQ(J,3))                                     &
     &          +P(2)*(SP(J,2)*SQ(J,3)+                                     &
     &                 SQ(J,2)*SP(J,3))                                     &
     &          +P(3)*(SP(J,2)*SP(J,3))
         PJ(J,2)=P(1)*(SQ(J,1)*SQ(J,3))                                     &
     &          +P(2)*(SP(J,1)*SQ(J,3)+                                     &
     &                 SQ(J,1)*SP(J,3))                                     &
     &          +P(3)*(SP(J,1)*SP(J,3))
         PJ(J,3)=P(1)*(SQ(J,2)*SQ(J,1))                                     &
     &          +P(2)*(SP(J,2)*SQ(J,1)+                                     &
     &                 SQ(J,2)*SP(J,1))                                     &
     &          +P(3)*(SP(J,2)*SP(J,1))
      ENDDO
      RETURN
      END
!--------------------------------------------------------------------------------------------------
      SUBROUTINE SW2L3(N,NB,NL,NPJ,SP,SQ,P,PJ)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      DIMENSION P(NL),PJ(N,NPJ),SP(N,NB),SQ(N,NB)
      DO J=1,N
         PJ(J,1)=P(1)*SQ(J,3)                                               &
     &          +P(2)*SP(J,3)
         PJ(J,2)=P(1)*SQ(J,2)                                               &
     &          +P(2)*SP(J,2)
         PJ(J,3)=P(1)*SQ(J,1)                                               &
     &          +P(2)*SP(J,1)
      ENDDO
      RETURN
      END
!--------------------------------------------------------------------------------------------------
      SUBROUTINE SP2L3(N,NB,NL,NPJ,NO,SP,SQ,P,PJ)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      DIMENSION P(NO,NL),PJ(N,NPJ,NO),SP(N,NB),SQ(N,NB)
      DO L=1,NO
         DO J=1,N
            PJ(J,1,L)=P(L,1)*SQ(J,3)                                        &
     &               +P(L,2)*SP(J,3)
            PJ(J,2,L)=P(L,1)*SQ(J,2)                                        &
     &               +P(L,2)*SP(J,2)
            PJ(J,3,L)=P(L,1)*SQ(J,1)                                        &
     &               +P(L,2)*SP(J,1)
         ENDDO
      ENDDO
      RETURN
      END
!--------------------------------------------------------------------------
      SUBROUTINE SP2G3(N,NB,NL,N1,N2,SP,SQ,P,PJ)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      DIMENSION P(N1,N2,NL),PJ(N,N1,N2),SP(N,NB),SQ(N,NB)
      DO L1=1,N1
         DO L2=1,N2
            DO J=1,N
               PJ(J,L1,L2)=P(L1,L2,1)*(SQ(J,1)*SQ(J,2)*SQ(J,3))             &
     &                    +P(L1,L2,2)*(SP(J,1)*SQ(J,2)*SQ(J,3)+             &
     &                                 SQ(J,1)*SP(J,2)*SQ(J,3)+             &
     &                                 SQ(J,1)*SQ(J,2)*SP(J,3))             &
     &                    +P(L1,L2,2)*(SP(J,1)*SP(J,2)*SQ(J,3)+             &
     &                                 SP(J,1)*SQ(J,2)*SP(J,3)+             &
     &                                 SQ(J,1)*SP(J,2)*SP(J,3))
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
!--------------------------------------------------------------------------------------------------
      REAL(kind=8) FUNCTION MTANH(X)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      PARAMETER(TOL=100.)
      IF(X.GE.0.) THEN
       MTANH=1.
       IF(X.LT.TOL) THEN
        EP=EXP(-X)
        MTANH=(1.-EP**2)/(1.+EP**2)
       ENDIF
      ELSE
       MTANH=-1.
       IF(X.GT.-TOL) THEN
        EP=EXP(X)
        MTANH=-(1.-EP**2)/(1.+EP**2)
       ENDIF
      ENDIF
      RETURN
      END
!-----------------------------------------------------------------------------
      SUBROUTINE POTPAR(P,NP)
! definition of molecule specific potential parameters
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      PARAMETER(NI=3,MI=NI-1,NORD=4,NABORD=10)
      DIMENSION P(NP)
! nominal parameters
! units are aJ=10**(-18)J for energies and A=100pm for lengths
      COMMON /PAR0/   V0(1:NI)
      COMMON /PARS/   AR(1:NI),RE(1:NI),FS(1:NI),AS(1:NI),BS(1:NI),         &
     &                E6(1:NI),E8(1:NI),RS(1:NI)
      COMMON /PARB/   AC(1:MI),CE(1:MI),FB1(1:MI),FB2(1:MI),                &
     &                AB1(1:NORD,1:NABORD,1:MI),                            &
     &                AB2(1:NORD,1:NABORD,1:MI),                            &
     &                N1(1:NORD),N2(1:NORD)
      COMMON /PARD/   DL1C(1:NORD,1:MI),DL1O(1:NORD,1:MI),                  &
     &                DL2C(1:NORD,1:MI),DL2O(1:NORD,1:MI),                  &
     &                DX1C(1:NORD,1:MI),DX1O(1:NORD,1:MI),                  &
     &                DX2C(1:NORD,1:MI),DX2O(1:NORD,1:MI),                  &
     &                DLPC(1:MI),DLPO(1:MI),                                &
     &                DXPC(1:MI),DXPO(1:MI),                                &
     &                ND1,ND2,NDCOS
      COMMON /PARP/   DII(1:MI),AII(1:MI),                                  &
     &                RIII,AIII,DIII        
      COMMON /PARSW/  RSWSB,RSWDB,NSWSB,NSWDB
      COMMON /DECOMP/ IFDEC
      IFDEC=0
!
      V0(1)        =   0.00000000              
      AR(1)        =   2.21018960         
      RE(1)        =   1.01279998             
      FS(1)        =   6.42209804             
      AS(1)        =   2.15200162         
      BS(1)        =   0.00000000              
      E6(1)        =   0.844658148E-01          
      E8(1)        =   0.00000000             
      RS(1)        =   2.00000000             
      AC(1)        =   1.47201434         
      CE(1)        =  -0.292883635        
      FB1(1)       =   0.588056500E-01    
      FB2(1)       =   0.128800000              
      AB1(1,1,1)   =   1.73205078             
      AB1(2,1,1)   =  -0.323879000        
      AB1(2,2,1)   =   1.60251900             
      AB1(3,1,1)   =   0.157522000              
      AB1(3,2,1)   =   3.65791400              
      AB1(3,3,1)   =   0.948204500              
      AB1(4,1,1)   =   0.110993300E-01         
      AB1(4,2,1)   =  -3.81903300             
      AB1(4,3,1)   =  -2.21657400             
      AB1(4,4,1)   =  -2.48302800             
      AB2(1,1,1)   =   1.41421354             
      AB2(2,1,1)   =  -1.60000000              
      AB2(2,2,1)   =  -0.800000000        
      AB2(3,1,1)   =  -2.20000000         
      AB2(3,2,1)   =  -2.00000000             
      AB2(3,3,1)   =   0.280000000E-01          
      AB2(4,1,1)   =   3.10000000             
      AB2(4,2,1)   =   4.60000000             
      AB2(4,3,1)   =   2.50000000             
      AB2(4,4,1)   =   0.800000000E-01         
      AB2(4,5,1)   =   0.130000000             
      DII(1)       =   1.08090200               
      AII(1)       =   0.476633900               
      DXPC(1)      =   0.803090200             
      DXPO(1)      =   1.11268400             
      DLPC(1)      =   0.00000000         
      DLPO(1)      =  -0.50000000         
      DX1C(1,1)    =   1.51390000         
      DX1C(2,1)    =   1.51390000              
      DX1C(3,1)    =   1.00000000              
      DX1C(4,1)    =   1.00000000              
      DX1O(1,1)    =   1.40488800         
      DX1O(2,1)    =   1.40488800              
      DX1O(3,1)    =   1.00000000              
      DX1O(4,1)    =   1.00000000              
      DL1C(1,1)    =   0.00000000         
      DL1C(2,1)    =   0.00000000              
      DL1C(3,1)    =   0.00000000              
      DL1C(4,1)    =   0.00000000              
      DL1O(1,1)    =  -1.00000000         
      DL1O(2,1)    =  -1.00000000              
      DL1O(3,1)    =  -1.00000000              
      DL1O(4,1)    =  -1.00000000              
      DX2C(1,1)    =   2.13615600         
      DX2C(2,1)    =   2.13615600              
      DX2C(3,1)    =   1.00000000              
      DX2C(4,1)    =   1.00000000              
      DX2O(1,1)    =  -0.205761900        
      DX2O(2,1)    =  -0.205761900              
      DX2O(3,1)    =   1.00000000              
      DX2O(4,1)    =   1.00000000              
      DL2C(1,1)    =   0.00000000         
      DL2C(2,1)    =   0.00000000             
      DL2C(3,1)    =   0.00000000             
      DL2C(4,1)    =   0.00000000             
      DL2O(1,1)    =  -1.00000000         
      DL2O(2,1)    =  -1.00000000             
      DL2O(3,1)    =  -1.00000000              
      DL2O(4,1)    =  -1.00000000              
      V0(2)        =   0.00000000              
      AR(2)        =   3.06553300              
      RE(2)        =   1.02600002              
      FS(2)        =   6.28749700              
      AS(2)        =   2.24860200              
      BS(2)        =   0.00000000              
      E6(2)        =   0.683101300E-01    
      E8(2)        =   0.00000000              
      RS(2)        =   2.00000000              
      AC(2)        =   1.47201434              
      CE(2)        =  -0.236838147               
      FB1(2)       =   0.137780400               
      FB2(2)       =   0.100000000        
      AB1(1,1,2)   =   1.73205078             
      AB1(2,1,2)   =  -0.378467400              
      AB1(2,2,2)   =  -0.378467400              
      AB1(3,1,2)   =  -0.948867500E-01    
      AB1(3,2,2)   =  -0.948867500E-01    
      AB1(3,3,2)   =  -0.948867500E-01    
      AB1(4,1,2)   =   0.00000000             
      AB1(4,2,2)   =   0.00000000             
      AB1(4,3,2)   =   0.00000000             
      AB1(4,4,2)   =   0.00000000             
      AB2(1,1,2)   =   1.41421354             
      AB2(2,1,2)   =   0.00000000             
      AB2(2,2,2)   =   0.00000000             
      AB2(3,1,2)   =   0.00000000             
      AB2(3,2,2)   =   0.00000000             
      AB2(3,3,2)   =   0.00000000             
      AB2(4,1,2)   =   0.00000000             
      AB2(4,2,2)   =   0.00000000             
      AB2(4,3,2)   =   0.00000000             
      AB2(4,4,2)   =   0.00000000             
      AB2(4,5,2)   =   0.00000000             
      DII(2)       =   1.00000000             
      AII(2)       =   0.501943200              
      DXPC(2)      =   1.00000000             
      DXPO(2)      =   1.00000000             
      DLPC(2)      =   0.00000000             
      DLPO(2)      =  -1.00000000             
      DX1C(1,2)    =   1.00000000         
      DX1C(2,2)    =   1.00000000              
      DX1C(3,2)    =   1.00000000              
      DX1C(4,2)    =   1.00000000              
      DX1O(1,2)    =   1.00000000         
      DX1O(2,2)    =   1.00000000              
      DX1O(3,2)    =   1.00000000              
      DX1O(4,2)    =   1.00000000              
      DL1C(1,2)    =   0.00000000         
      DL1C(2,2)    =   0.00000000              
      DL1C(3,2)    =   0.00000000              
      DL1C(4,2)    =   0.00000000              
      DL1O(1,2)    =  -1.00000000         
      DL1O(2,2)    =  -1.00000000              
      DL1O(3,2)    =  -1.00000000              
      DL1O(4,2)    =  -1.00000000              
      DX2C(1,2)    =   1.00000000         
      DX2C(2,2)    =   1.00000000              
      DX2C(3,2)    =   1.00000000              
      DX2C(4,2)    =   1.00000000              
      DX2O(1,2)    =   1.00000000         
      DX2O(2,2)    =   1.00000000             
      DX2O(3,2)    =   1.00000000             
      DX2O(4,2)    =   1.00000000             
      DL2C(1,2)    =   0.00000000         
      DL2C(2,2)    =   0.00000000             
      DL2C(3,2)    =   0.00000000             
      DL2C(4,2)    =   0.00000000             
      DL2O(1,2)    =  -1.00000000         
      DL2O(2,2)    =  -1.00000000             
      DL2O(3,2)    =  -1.00000000              
      DL2O(4,2)    =  -1.00000000              
      V0(3)        =   0.105836722              
      AR(3)        =   3.06553300             
      RE(3)        =   1.03635001             
      FS(3)        =   6.39284900             
      AS(3)        =   2.15608600             
      BS(3)        =   0.00000000             
      E6(3)        =   0.236717900E-02    
      E8(3)        =  -0.110000000              
      RS(3)        =   2.15000000             
      NSWSB        =   6
      RSWSB        =   2.14571500             
      RIII         =   0.741400003        
      AIII         =   1.94534004             
      NSWDB        =   12             
      RSWDB        =   1.50000000         
      DIII         =   0.760596812        
      ND1          =   1         
      ND2          =   2         
      NDCOS        =   1         
      RETURN            
      END               
!--------------------------------------------------------------------------
      SUBROUTINE DAMPG(N,NI,NNI,NORD,R,RE,C,XA,XL,N1,N2,YD)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      DIMENSION XA(N,NNI,NORD),XL(N,NNI,NORD)
      DIMENSION R(N,NI),RE(N,NI),C(N,NNI)
      DIMENSION YD(N,NNI,NORD)
      PI=ACOS(-1.)
      URF=PI/180.
      DO NO=1,NORD
         DO I=1,NI
            DO J = 1,N 
               R(J,I) = R(J,I) / RE(J,I)
            ENDDO
         ENDDO
         KK=0
         DO I=1,NI
            DO K=I+1,NI
               KK = KK +1
               DO J=1,N
                  YD(J,KK,NO) =                                             &
     &                 EXP(- EXP(XL(J,KK,NO))*                              &
     &                 (R(J,I)-1.)**N1*                                     &
     &                 (R(J,I)-XA(J,KK,NO))**N2)*                           &
     &                 EXP(- EXP(XL(J,KK,NO))*                              &
     &                 (R(J,K)-1.)**N1*                                     &
     &                 (R(J,K)-XA(J,KK,NO))**N2)
               ENDDO
            ENDDO
         ENDDO
         DO I=1,NI
            DO J = 1,N 
               R(J,I) = R(J,I) * RE(J,I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
!-------------------------------------------------------------------------------      
      SUBROUTINE DANG(N,MM,NORD,NDCOS,C,PC,PO,P)
      IMPLICIT REAL(kind=8) (A-H,O-Z)
      DIMENSION PC(N,MM,NORD),PO(N,MM,NORD) 
      DIMENSION P(N,MM,NORD)
      DIMENSION C(N,MM)
      DO NO=1,NORD
         DO KK=1,MM
            DO J=1,N
               P(J,KK,NO)=PC(J,KK,NO)*(0.5*(1.0+C(J,KK)))**NDCOS            &
     &                   +PO(J,KK,NO)*(0.5*(1.0-C(J,KK)))**NDCOS
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
!--------------------------------------------------------------------------
      SUBROUTINE CONDIF(I,PP,F,DF,N)
!
! dummy subroutine
!
      RETURN
      END
#endif
END MODULE POTENTIAL_NH3
!-------------------------------------------------------------------------------


