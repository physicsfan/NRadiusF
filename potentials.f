CCCCCCCCCCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     potentials.f
C
C     A LIBRARY OF FUNCTIONS AND SUBROUTINES THAT CALCULATES CENTRAL
C     ELECTROSTATIC POTENTIALS FOR THREE NUCLEAR CONFIGURATIONS:
C
C     1) POINT NUCLEUS
C     2) HARD-SPHERE NUCLEUS
C     3) FERMI CHARGE DISTRIBUTION (AS GIVEN IN JOHNSON AST)
C
C     THIS LIBRARY IS PRIMARILY A SUPPLIMENT TO THE 'RADIAL' LIBRARY,
C     BUT CAN BE USED INDEPENDENTLY IF REQUIRED.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     INITIALIZATION DATA
      BLOCK DATA PINIT
      INTEGER NMAX
      COMMON NMAX
      DATA NMAX /1E2/
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     POINT NUCLEUS POTENTIAL
      SUBROUTINE POINTP(ZNUC,RSIZE,RVDIM, RV)

      INTEGER I, RSIZE, RVDIM
      DOUBLE PRECISION RV(RSIZE), ZNUC

      DO 12, I = 1, RVDIM
         RV(I) = ZNUC
 12   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     HARD-SPHERE NUCLEAR POTENTIAL
      SUBROUTINE HRDSPH(ZNUC,RADIUS,RSIZE,RVDIM,RGRID, RV)

      INTEGER I, RSIZE, RVDIM
      DOUBLE PRECISION  RADIUS, SPHRAD, ZNUC
      DOUBLE PRECISION  RGRID(RSIZE), RV(RSIZE)
      LOGICAL DGNSTC
      COMMON /TEST/ DGNSTC

C     DIAGNOSTICS
      IF(DGNSTC) OPEN(1,FILE='hsdiagnostics.dat')

C     CONVERT CHARGE RADIUS TO SPHERICAL CHARGE RADIUS
      SPHRAD = SQRT(5.0D0/3.0D0)*RADIUS

C     DIAGNOSTICS
      IF(DGNSTC) WRITE(1,'(F20.15,F25.15/)') RADIUS, SPHRAD

C     CALCULATE HARD-SPHERE POTENTIAL GRID
      DO 22 I=1,RVDIM
         IF(RGRID(I).LE.SPHRAD .OR. RGRID(I).EQ.SPHRAD) THEN
            RV(I)=((ZNUC*RGRID(I))/(2.0D0*SPHRAD))*(3-((RGRID(I)**2)
     1      /(SPHRAD**2)))
         ELSE
            RV(I)=ZNUC
         ENDIF
         
C       DIAGNOSTICS
        IF(DGNSTC)  WRITE(1,'(F20.15,F25.15)') RGRID(I), RV(I)
 22   CONTINUE
      
      IF(DGNSTC) CLOSE(UNIT=1)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     FERMI CHARGE DISTRIBUTION POTENTIAL
      SUBROUTINE FERMIP(ZNUC,RADIUS,RSIZE,RVDIM,RGRID, RV)

      INTEGER I, NMAX, RSIZE, RVDIM
      DOUBLE PRECISION  A,BOHR,C,M,N,PI,RADIUS,R,RTEST,ZNUC
      DOUBLE PRECISION  P, S
      DOUBLE PRECISION RGRID(RSIZE), RV(RSIZE)
      LOGICAL DGNSTC
      COMMON /TEST/ DGNSTC
      COMMON /FERMI/ A, C
      COMMON /MCON/ PI
      COMMON /PCON/ BOHR
      COMMON NMAX
      
      IF(DGNSTC) OPEN(1,FILE='fdiagnostics.dat')

      
C     CALCULATE FERMI POTENTIAL PARAMETERS: A, C
      A = (2.3D0 / (4.0D0*LOG(3.0D0)))*1.0D-15/BOHR
      C = SQRT((5.0D0/3.D0)*RADIUS**2-(7.0D0/3.0D0)*(A*PI)**2)

C     DIAGNOSTIC OUTPUT
      IF(DGNSTC) THEN
         R = RADIUS*BOHR*1.0D15
         WRITE(1,*) BOHR, PI, NMAX
         WRITE(1,*) RADIUS
         WRITE(1,*) R, A, C
      ENDIF
      
C     CALCULATE FERMI POTENTIAL NORMALIZATION TERMS: N, M
      N = 1.0D0 + (PI*A/C)**2 + 6.0D0*(A/C)**3 * S(3)
      M = 1.0D0 + (10.D0/3.D0)*(PI*A/C)**2 + (7.D0/3.D0)*(PI*A/C)**4
     1     +120.0D0*(A/C)**5 * S(5)

C     DIAGNOSTIC OUTPUT
      IF(DGNSTC) THEN
         RTEST = C * SQRT((3.0D0/5.0D0)*(M/N))
         IF (ABS(RTEST-RADIUS).GT. 1E-10) THEN
            WRITE(1,*) "AGREEMENT INSUFFICIENT. STOPPING CALCULATION..."
            RETURN
         ENDIF
      ENDIF

C     CALCULATE FERMI POTENTIAL GRID
      DO 32 I=1,RVDIM
         IF(RGRID(I).LE.C) THEN
            RV(I)=(ZNUC*RGRID(I)/(N*C))*(1.5D0 - 0.5D0*(RGRID(I)/C)**2
     1           + 0.5*(PI*A/C)**2 + 3.0D0*(A/C)**2 * P(2, RGRID(I))
     2           + 6.0D0*(A**3/(C**2*RGRID(I)))*(S(3) - P(3, RGRID(I))))
            IF(RV(I).NE.RV(I)) RV(I) = 0.0D0
         ELSE
            RV(I)=(ZNUC/N)*(1.0D0 + 6.0D0*(A/C)**3*(S(3)-P(3, RGRID(I)))
     1           +(PI*A/C)**2-(3.0D0*RGRID(I)*A**2/C**3)*P(2, RGRID(I)))
         ENDIF
         IF(DGNSTC) WRITE(1,'(F20.15,F25.15)') RGRID(I), RV(I)
 32   CONTINUE
      IF(DGNSTC) CLOSE(UNIT=1)
      RETURN
      END

C*********************************************************************

      DOUBLE PRECISION FUNCTION S(K)

      INTEGER K, N
      DOUBLE PRECISION A, BOHR, C, SSUM
      COMMON /FERMI/ A, C
      COMMON /PARAMS/ BOHR, PI
      COMMON NMAX

      SSUM = 0.0D0
      DO 42, N = 1, NMAX
         SSUM = SSUM + (((-1.0D0)**(N-1))/N**K)*EXP(-N*C/A)
 42   CONTINUE

      S = SSUM

      END

C*********************************************************************

      DOUBLE PRECISION FUNCTION P(K, R)

      INTEGER K, N
      DOUBLE PRECISION A, BOHR, C, PSUM, R
      COMMON /FERMI/ A, C
      COMMON /PARAMS/ BOHR, PI
      COMMON NMAX
      INTRINSIC ABS

      PSUM = 0.0D0
      DO 52, N = 1, NMAX
         PSUM = PSUM + (((-1.0D0)**(N-1))/N**K)*EXP(-N*ABS(R-C)/A)
 52   CONTINUE

      P = PSUM

      END

C*********************************************************************
