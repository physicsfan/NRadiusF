C                                                                     
C         ****************************************                 
C         **  MAIN PROGRAM FOR NRADIUS PACKAGE  **
C         ****************************************                  
C                                                                       
C     THIS IS A COMPLETE FORTRAN77 VERSION OF THE NRADIUS PACKAGE ORIGINALLY WRITTEN
C     IN MATHEMATICA.  ITS PURPOSE IS TO ULTIMATELY EXTRACT THE DEPENDENCE OF THE
C     NUCLEAR RADIUS AS A FUNCTION OF THE TRANSITION ENERGY IN LI-LIKE HEAVY IONS
C     BETWEEN LEAD (Z=82) AND URANIUM (Z=92)
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      PARAMETER (NDIM=8000,NPPG=NDIM+1,NPTG=NDIM+NPPG)                   
      PARAMETER (SL=137.036D0,PI=3.1415926535897932D0)                  
      DIMENSION R0(NDIM),RV0(NDIM), RV0TST(NDIM)
      LOGICAL DGNSTC


C  **** DIAGNOSTICS
      COMMON /TEST/ DGNSTC
C      
C  ****  COULOMB WAVE FUNCTION PARAMETERS.                              
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C      
C  ****  OUTPUT RADIAL FUNCTIONS.                                       
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER

      EXTERNAL POINTP,HRDSPH,FERMIP
C
C  **** INITIALIZE OUTPUT FILES
C      OPEN(8,FILE='RVOUT.DAT')
      
C                                                                       
C  ****  READ FIELD PARAMETERS.                                 
C
      DGNSTC=.TRUE.
 10   CONTINUE

      
      WRITE(6,*) '**** Welcome to the RADIAL Test Program ****'
      WRITE(6,*) '  '
      WRITE(6,*) 'ENTER Ze, Rrms (in fermi''s) ...'                                
      READ(5,*) Z,Rrms                                                                                   
C
C     CONVERT Rrms TO HARTREES
C
      Rrms = Rrms*1.0D-15/5.2917721092D-11

      PRINT *, Rrms
C
C  ****  POTENTIAL GRID.                                                
C
      RATIO=1.15D0
      RNN=50.0D0/83.0D0
C      RNN=0.659999D0                                    
      NV=6000                                                            
      STEP=RNN/(NV-100.0D0)                                             
      CALL GRID(R0,RATIO,STEP,RNN,NV,Z)                                   
C
      print *, size(R0)
      
 20   CONTINUE
      WRITE(6,*) '  '                                                   
      WRITE(6,*) '  SELECT ONE POTENTIAL OPTION ...'                              
      WRITE(6,*) '    1: POINT NUCLEUS POTENTIAL'            
      WRITE(6,*) '    2: HARD-SPHERE POTENTIAL'             
      WRITE(6,*) '    3: FERMI POTENTIAL'
      WRITE(6,*) '    4: EXIT'
      READ(5,*) NUMPOT                                            
      WRITE(6,*) '  '                                                   
      

      IF (NUMPOT .EQ. 1) THEN
         CALL POINTP(Z,NDIM,NV, RV0)
      ELSE IF (NUMPOT .EQ. 2) THEN
         CALL HRDSPH(Z,Rrms,NDIM,NV,R0, RV0)
      ELSE IF (NUMPOT .EQ. 3) THEN
         CALL FERMIP(Z,Rrms,NDIM,NV,R0, RV0)
      ELSE IF (NUMPOT .EQ. 4) THEN
         GOTO 9999
      ELSE
         GO TO 20
      END IF

C     Preliminary check of potential grid
      IF(DGNSTC) THEN
         IF (NUMPOT .EQ. 1) THEN
            DO 21, I = 1, NV
               RV0TST(I) = Z
               IF(ABS(RV0TST(I)-RV0(I)).GE.1.0D-10) THEN
                  PRINT *, "THEY DON'T MATCH AT", RV0(I)
                  GOTO 9999
               ENDIF
 21         CONTINUE
         ENDIF
         
         IF (NUMPOT .EQ. 2) THEN
            
            SPHRAD = SQRT(5.0D0/3.0D0)*Rrms
            
            DO 22 I=1, NV
               IF(R0(I).LE.SPHRAD .OR. R0(I).EQ.SPHRAD) THEN
                  RV0TST(I)=((Z*R0(I))/(2.0D0*SPHRAD))*(3-((R0(I)**2)
     $                 /(SPHRAD**2)))
C     RV0(I)=((Z*R0(I))/(2.0D0*Rm))*(3-((R0(I)**2)/(Rm**2)))
               ELSE
                  RV0TST(I)=Z
               ENDIF
 22         CONTINUE
            
            DO 23 I=1,NV
               IF(ABS(RV0TST(I)-RV0(I)).GE.1.0D-10) THEN
                  PRINT *, "THEY DON'T MATCH STARTING AT", I
                  WRITE(8,3000)
                  DO 24 J=1, NV
                     WRITE(8,3010) R0(J), RV0(J), RV0TST(J)
 24               CONTINUE
                  GOTO 9999
               ENDIF
 23         CONTINUE
            CLOSE(UNIT=8)
 3000       FORMAT(10X,'R0',22X,'RV0',22X,'RV0TST')     
 3010       FORMAT(F20.15,F25.15,F25.15)
         ENDIF
      ENDIF
      
      
C  ----  SPLINE INTERPOLATION TEST.                             
C  (Delete the comments when required. This is time consuming!).        
C     CALL ERRSPL(ERR,R0,RV0,NV)                                        
C     IF(ERR.GT.1.0D-5) THEN                                            
C       WRITE(6,1001) ERR                                               
C1001 FORMAT(1X,'*** ACCUMULATED ERRORS IN THE INTERPOLATING ',         
C    1  'SPLINE. ERROR =',1P,E9.2,'.')                                  
C       STOP                                                            
C     ENDIF                                                             
C                                                                       
      CALL VINT(R0,RV0,NV)                                              
C                                                                       
   30 CONTINUE                                                          
      WRITE(6,*) '  '                                                   
      WRITE(6,*) '  SELECT ONE OPTION ...'                              
      WRITE(6,*) '    1: SCHRODINGER EQUATION. BOUND STATE.'            
      WRITE(6,*) '    2: SCHRODINGER EQUATION. FREE STATE.'             
      WRITE(6,*) '    3: DIRAC EQUATION. BOUND STATE.'                  
      WRITE(6,*) '    4: DIRAC EQUATION. FREE STATE.'
      WRITE(6,*) '    5: GO BACK.'
      READ(5,*) IOPT                                            
      WRITE(6,*) '  '                                                   
C                                                                       
C  ****  SCHRODINGER EQUATION. BOUND STATE.                             
C                                                                       
      IF(IOPT.EQ.1) THEN                                                
        WRITE(6,*) '  ENTER N, L AND EPS ...'                           
        READ(5,*) N,L,EPS                                               
        EPS=DMAX1(EPS,1.0D-15)                                          
        IF(N.LT.1.OR.L.GE.N) GO TO 30                                   
        RN=2.0D3                                                        
        NGP=8000                                                         
        STEP=RN/(NGP-100.0D0)                                           
        CALL GRID(RAD,RATIO,STEP,RN,NGP,Z)                                
        E=-Z**2/(2.0D0*N*N)                                             
        CALL SBOUND(E,EPS,EPS,N,L)                                      
        IF(IER.NE.0) GO TO 30                                           
        WRITE(6,1100) Z,ZS,ALPHA,N,L,EPS,E                              
 1100 FORMAT(/1X,1P,'****  SCHRODINGER EQ. ',                           
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'            
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      
     3  /7X,'BOUND STATE: N=',I4,', L=',I4,'   (EPS=',E8.1,')'          
     4  /7X,'BINDING ENERGY =',E22.15)                                  
C                                                                       
C  ****  SCHRODINGER EQUATION. FREE STATE.                              
C                                                                       
      ELSE IF(IOPT.EQ.2) THEN                                           
        WRITE(6,*) '  ENTER E, L AND EPS ...'                           
        READ(5,*) E,L,EPS                                               
        EPS=DMAX1(EPS,1.0D-15)                                          
        IF(E.LT.0.0D0.OR.L.LT.0) GO TO 30                               
        NGP=6000                                                         
        WAVEL=2.0D0*PI/DSQRT(E+E)                                       
        STEP=0.05D0*WAVEL                                               
        RN=STEP*(NGP-100)                                               
        CALL GRID(RAD,RATIO,STEP,RN,NGP,Z)                                
        CALL SFREE(E,EPS,PHASE,L)                                       
        IF(IER.NE.0) GO TO 30                                           
        WRITE(6,1200) Z,ZS,ALPHA,E,L,EPS,PHASE,DELTA,ETA        
 1200 FORMAT(/1X,1P,'****  SCHRODINGER EQ. ',                           
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'                    
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      
     3  /7X,'FREE STATE: E=',E13.6,', L=',I4,'   (EPS=',E8.1,')'        
     4  /7X,'  INNER PHASE SHIFT=',E22.15,                              
     5  /7X,'COULOMB PHASE SHIFT=',E22.15,'   (ETA=',E13.6,')')          
C                                                                       
C  ****  DIRAC EQUATION. BOUND STATE.                                   
C                                                                       
      ELSE IF(IOPT.EQ.3) THEN                                           
        WRITE(6,*) '  ENTER N, K AND EPS ...'                           
        READ(5,*) N,K,EPS                                               
        EPS=DMAX1(EPS,1.0D-15)                                          
        IF(N.LT.1.OR.K.EQ.0.OR.K.GE.N.OR.K.LT.-N) GO TO 30              
        RN=2.0D3                                                        
        NGP=8000                                                         
        STEP=RN/(NGP-100.0D0)                                           
        CALL GRID(RAD,RATIO,STEP,RN,NGP,Z)                                
        E=-Z**2/(2.0D0*N*N)                                     
        CALL DBOUND(E,EPS,EPS,N,K)                                      
        IF(IER.NE.0) GO TO 30                                           
        WRITE(6,1300) NUMPOT,Z,N,K,EPS,E                              
 1300 FORMAT(/1X,1P,'****  DIRAC EQUATION. ',                           
     1  'POTENTIAL FUNCTION: ',I1,                    
     2  /7X,'Z=',E13.6,                      
     3  /7X,'BOUND STATE: N=',I4,', K=',I4,'   (EPS=',E8.1,')'        10   CONTINUE

      
      WRITE(6,*) '**** Welcome to the RADIAL Test Program ****'
      WRITE(6,*) '  '
      WRITE(6,*) 'ENTER Ze, Rrms (in fermi''s) ...'                                
      READ(5,*) Z,Rrms                                                                                   
C
C     CONVERT Rrms TO HARTREES
C
      Rrms = Rrms*1.0D-15/5.2917721092D-11

      PRINT *, Rrms
C
C  ****  POTENTIAL GRID.                                                
C
      RATIO=1.15D0
      RNN=50.0D0/83.0D0
C      RNN=0.659999D0                                    
      NV=6000                                                            
      STEP=RNN/(NV-100.0D0)                                             
      CALL GRID(R0,RATIO,STEP,RNN,NV,Z)                                   
C
      print *, size(R0)
      
 20   CONTINUE
      WRITE(6,*) '  '                                                   
      WRITE(6,*) '  SELECT ONE POTENTIAL OPTION ...'                              
      WRITE(6,*) '    1: POINT NUCLEUS POTENTIAL'            
      WRITE(6,*) '    2: HARD-SPHERE POTENTIAL'             
      WRITE(6,*) '    3: FERMI POTENTIAL'
      WRITE(6,*) '    4: EXIT'
      READ(5,*) NUMPOT                                            
      WRITE(6,*) '  '                                                   
      

      IF (NUMPOT .EQ. 1) THEN
         CALL POINTP(Z,NDIM,NV, RV0)
      ELSE IF (NUMPOT .EQ. 2) THEN
         CALL HRDSPH(Z,Rrms,NDIM,NV,R0, RV0)
      ELSE IF (NUMPOT .EQ. 3) THEN
         CALL FERMIP(Z,Rrms,NDIM,NV,R0, RV0)
      ELSE IF (NUMPOT .EQ. 4) THEN
         GOTO 9999
      ELSE
         GO TO 20
      END IF

C     Preliminary check of potential grid
      IF(DGNSTC) THEN
         IF (NUMPOT .EQ. 1) THEN
            DO 21, I = 1, NV
               RV0TST(I) = Z
               IF(ABS(RV0TST(I)-RV0(I)).GE.1.0D-10) THEN
                  PRINT *, "THEY DON'T MATCH AT", RV0(I)
                  GOTO 9999
               ENDIF
 21         CONTINUE
         ENDIF
         
         IF (NUMPOT .EQ. 2) THEN
            
            SPHRAD = SQRT(5.0D0/3.0D0)*Rrms
            
            DO 22 I=1, NV
               IF(R0(I).LE.SPHRAD .OR. R0(I).EQ.SPHRAD) THEN
                  RV0TST(I)=((Z*R0(I))/(2.0D0*SPHRAD))*(3-((R0(I)**2)
     $                 /(SPHRAD**2)))
C     RV0(I)=((Z*R0(I))/(2.0D0*Rm))*(3-((R0(I)**2)/(Rm**2)))
               ELSE
                  RV0TST(I)=Z
               ENDIF
 22         CONTINUE
            
            DO 23 I=1,NV
               IF(ABS(RV0TST(I)-RV0(I)).GE.1.0D-10) THEN
                  PRINT *, "THEY DON'T MATCH STARTING AT", I
                  WRITE(8,3000)
                  DO 24 J=1, NV
                     WRITE(8,3010) R0(J), RV0(J), RV0TST(J)
 24               CONTINUE
                  GOTO 9999
               ENDIF
 23         CONTINUE
            CLOSE(UNIT=8)
 3000       FORMAT(10X,'R0',22X,'RV0',22X,'RV0TST')     
 3010       FORMAT(F20.15,F25.15,F25.15)
         ENDIF
      ENDIF
      
      
C  ----  SPLINE INTERPOLATION TEST.                             
C  (Delete the comments when required. This is time consuming!).        
C     CALL ERRSPL(ERR,R0,RV0,NV)                                        
C     IF(ERR.GT.1.0D-5) THEN                                            
C       WRITE(6,1001) ERR                                               
C1001 FORMAT(1X,'*** ACCUMULATED ERRORS IN THE INTERPOLATING ',         
C    1  'SPLINE. ERROR =',1P,E9.2,'.')                                  
C       STOP                                                            
C     ENDIF                                                             
C                                                                       
      CALL VINT(R0,RV0,NV)                                              
C                                                                       
   30 CONTINUE                                                          
      WRITE(6,*) '  '                                                   
      WRITE(6,*) '  SELECT ONE OPTION ...'                              
      WRITE(6,*) '    1: SCHRODINGER EQUATION. BOUND STATE.'            
      WRITE(6,*) '    2: SCHRODINGER EQUATION. FREE STATE.'             
      WRITE(6,*) '    3: DIRAC EQUATION. BOUND STATE.'                  
      WRITE(6,*) '    4: DIRAC EQUATION. FREE STATE.'
      WRITE(6,*) '    5: GO BACK.'
      READ(5,*) IOPT                                            
      WRITE(6,*) '  '                                                   
C                                                                       
C  ****  SCHRODINGER EQUATION. BOUND STATE.                             
C                                                                       
      IF(IOPT.EQ.1) THEN                                                
        WRITE(6,*) '  ENTER N, L AND EPS ...'                           
        READ(5,*) N,L,EPS                                               
        EPS=DMAX1(EPS,1.0D-15)                                          
        IF(N.LT.1.OR.L.GE.N) GO TO 30                                   
        RN=2.0D3                                                        
        NGP=8000                                                         
        STEP=RN/(NGP-100.0D0)                                           
        CALL GRID(RAD,RATIO,STEP,RN,NGP,Z)                                
        E=-Z**2/(2.0D0*N*N)                                             
        CALL SBOUND(E,EPS,EPS,N,L)                                      
        IF(IER.NE.0) GO TO 30                                           
        WRITE(6,1100) Z,ZS,ALPHA,N,L,EPS,E                              
 1100 FORMAT(/1X,1P,'****  SCHRODINGER EQ. ',                           
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'            
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      
     3  /7X,'BOUND STATE: N=',I4,', L=',I4,'   (EPS=',E8.1,')'          
     4  /7X,'BINDING ENERGY =',E22.15)                                  
C                                                                       
C  ****  SCHRODINGER EQUATION. FREE STATE.                              
C                                                                       
      ELSE IF(IOPT.EQ.2) THEN                                           
        WRITE(6,*) '  ENTER E, L AND EPS ...'                           
        READ(5,*) E,L,EPS                                               
        EPS=DMAX1(EPS,1.0D-15)                                          
        IF(E.LT.0.0D0.OR.L.LT.0) GO TO 30                               
        NGP=6000                                                         
        WAVEL=2.0D0*PI/DSQRT(E+E)                                       
        STEP=0.05D0*WAVEL                                               
        RN=STEP*(NGP-100)                                               
        CALL GRID(RAD,RATIO,STEP,RN,NGP,Z)                                
        CALL SFREE(E,EPS,PHASE,L)                                       
        IF(IER.NE.0) GO TO 30                                           
        WRITE(6,1200) Z,ZS,ALPHA,E,L,EPS,PHASE,DELTA,ETA        
 1200 FORMAT(/1X,1P,'****  SCHRODINGER EQ. ',                           
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'                    
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      
     3  /7X,'FREE STATE: E=',E13.6,', L=',I4,'   (EPS=',E8.1,')'        
     4  /7X,'  INNER PHASE SHIFT=',E22.15,                              
     5  /7X,'COULOMB PHASE SHIFT=',E22.15,'   (ETA=',E13.6,')')          
C                                                                       
C  ****  DIRAC EQUATION. BOUND STATE.                                   
C                                                                       
      ELSE IF(IOPT.EQ.3) THEN                                           
        WRITE(6,*) '  ENTER N, K AND EPS ...'                           
        READ(5,*) N,K,EPS                                               
        EPS=DMAX1(EPS,1.0D-15)                                          
        IF(N.LT.1.OR.K.EQ.0.OR.K.GE.N.OR.K.LT.-N) GO TO 30              
        RN=2.0D3                                                        
        NGP=8000                                                         
        STEP=RN/(NGP-100.0D0)                                           
        CALL GRID(RAD,RATIO,STEP,RN,NGP,Z)                                
        E=-Z**2/(2.0D0*N*N)                                     
        CALL DBOUND(E,EPS,EPS,N,K)                                      
        IF(IER.NE.0) GO TO 30                                           
        WRITE(6,1300) NUMPOT,Z,N,K,EPS,E                              
 1300 FORMAT(/1X,1P,'****  DIRAC EQUATION. ',                           
     1  'POTENTIAL FUNCTION: ',I1,                    
     2  /7X,'Z=',E13.6,                      
     3  /7X,'BOUND STATE: N=',I4,', K=',I4,'   (EPS=',E8.1,')'          
     4  /7X,'BINDING ENERGY =',E22.15)                                  
C                                                                       
C  ****  DIRAC EQUATION. FREE STATE.                                   
C                                                                       
      ELSE IF(IOPT.EQ.4) THEN                                           
        WRITE(6,*) '  ENTER E, K AND EPS ...'                           
        READ(5,*) E,K,EPS                                               
        EPS=DMAX1(EPS,1.0D-15)                                          
        IF(E.LT.0.0D0.OR.K.EQ.0) GO TO 30                               
        IF(K.LT.0) THEN                                                 
          L=-K-1                                                        
        ELSE                                                    
          L=K                                                           
        ENDIF                                                           
        NGP=6000                                                         
        WAVEL=2.0D0*PI/DSQRT(E*(2.0D0+E/SL**2))                         
        STEP=0.05D0*WAVEL                                               
        RN=STEP*(NGP-100)                                               
        CALL GRID(RAD,RATIO,STEP,RN,NGP,Z)                                
        CALL DFREE(E,EPS,PHASE,K)                                       
        IF(IER.NE.0) GO TO 30                                           
        WRITE(6,1400) Z,ZS,ALPHA,E,K,EPS,PHASE,DELTA,ETA               
 1400 FORMAT(/1X,1P,'****  DIRAC EQUATION. ',                           
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'                    
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      
     3  /7X,'FREE STATE: E=',E13.6,', K=',I4,'   (EPS=',E8.1,')'        
     4  /7X,'  INNER PHASE SHIFT=',E22.15,                              
     5  /7X,'COULOMB PHASE SHIFT=',E22.15,'   (ETA=',E13.6,')')         
      ELSE                                                              
        GO TO 20                                                        
      ENDIF                                                     
C                                                                       
C  ****  RADIAL WAVE FUNCTIONS WRITTEN ON FILE 'WAVES.DAT'.             
C                                                                       
      OPEN(7,FILE='WAVES.DAT')                                          
      DO 40 I=1,NGP                                                     
      IF(DABS(P(I)).LT.1.0D-35) P(I)=1.0D-35                     
      IF(DABS(Q(I)).LT.1.0D-35) Q(I)=1.0D-35                            
      WRITE(7,'(1X,I4,1P,4E13.5)') I,RAD(I),P(I),Q(I)                   
   40 CONTINUE                                                          
      CLOSE(UNIT=7)                                                     
      GO TO 30   
     4  /7X,'BINDING ENERGY =',E22.15)                                  
C                                                                       
C  ****  DIRAC EQUATION. FREE STATE.                                   
C                                                                       
      ELSE IF(IOPT.EQ.4) THEN                                           
        WRITE(6,*) '  ENTER E, K AND EPS ...'                           
        READ(5,*) E,K,EPS                                               
        EPS=DMAX1(EPS,1.0D-15)                                          
        IF(E.LT.0.0D0.OR.K.EQ.0) GO TO 30                               
        IF(K.LT.0) THEN                                                 
          L=-K-1                                                        
        ELSE                                                    
          L=K                                                           
        ENDIF                                                           
        NGP=6000                                                         
        WAVEL=2.0D0*PI/DSQRT(E*(2.0D0+E/SL**2))                         
        STEP=0.05D0*WAVEL                                               
        RN=STEP*(NGP-100)                                               
        CALL GRID(RAD,RATIO,STEP,RN,NGP,Z)                                
        CALL DFREE(E,EPS,PHASE,K)                                       
        IF(IER.NE.0) GO TO 30                                           
        WRITE(6,1400) Z,ZS,ALPHA,E,K,EPS,PHASE,DELTA,ETA               
 1400 FORMAT(/1X,1P,'****  DIRAC EQUATION. ',                           
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'                    
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      
     3  /7X,'FREE STATE: E=',E13.6,', K=',I4,'   (EPS=',E8.1,')'        
     4  /7X,'  INNER PHASE SHIFT=',E22.15,                              
     5  /7X,'COULOMB PHASE SHIFT=',E22.15,'   (ETA=',E13.6,')')         
      ELSE                                                              
        GO TO 20                                                        
      ENDIF                                                     
C                                                                       
C  ****  RADIAL WAVE FUNCTIONS WRITTEN ON FILE 'WAVES.DAT'.             
C                                                                       
      OPEN(7,FILE='WAVES.DAT')                                          
      DO 40 I=1,NGP                                                     
      IF(DABS(P(I)).LT.1.0D-35) P(I)=1.0D-35                     
      IF(DABS(Q(I)).LT.1.0D-35) Q(I)=1.0D-35                            
      WRITE(7,'(1X,I4,1P,4E13.5)') I,RAD(I),P(I),Q(I)                   
   40 CONTINUE                                                          
      CLOSE(UNIT=7)                                                     
      GO TO 30
 9999 CONTINUE
      END                                                               
C  **************************************************************       
C                        SUBROUTINE GRID                                
C  **************************************************************       
      SUBROUTINE GRID(R,RATIO,STEP,RN,NP,Z)                               
C                                                                       
C     THIS SUBROUTINE SETS UP A RADIAL GRID R(I) (I=1, ..., NP)         
C     SUCH THAT                                                            
C     1) R(1)=0, R(NP)=RN,                                              
C     2) A*R(I)+B*DLOG(R(I))-C=I  (I.GT.0), WITH                        
C        A=1.0/STEP AND B=1.0/DLOG(RATIO).                              
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      PARAMETER (NDIM=8000)                                              
      DIMENSION R(NDIM)                                         
C                                                                       
      IF(NP.LT.50) THEN                                                 
        WRITE(6,1001) NP                                                
 1001 FORMAT(1X,'*** ERROR IN GRID: NP =',I3,                           
     1  ' MUST BE LARGER THAN 50.')                                     
        STOP                                                            
      ENDIF                                                             
      IF(NP.GT.NDIM) THEN                                               
        WRITE(6,1002) NP,NDIM                                           
 1002 FORMAT(1X,'*** ERROR IN GRID: NP =',I5,                           
     1  ' IS LARGER THAN NDIM =',I5,'.')                                
        STOP                                                            
      ENDIF                                                             
      IF(STEP.LT.1.0D-10.OR.(RATIO-1.0D0).LT.1.0D-3) THEN               
        WRITE(6,1003) STEP, RATIO                                       
 1003 FORMAT(1X,'*** ERROR IN GRID: STEP =',1P,E10.3,' OR RATIO',       
     1  ' =',E10.3,' ARE TOO SMALL.')                                   
        STOP                                                            
      ENDIF                                                     
C                                                                       
      A=((1.0D0)*ABS(Z))/STEP                                                    
      B=1.0D0/DLOG(RATIO)                                               
      C=NP-A*(RN/ABS(Z))-B*DLOG(RN/ABS(Z))                                             
C                                                                       
      R(1)=0.0D0                                                        
      RR=1.0D-35                                                        
      FR=A*RR+B*DLOG(RR)+C-1                                            
      IF(FR.GT.0.0D0) THEN                                              
        WRITE(6,1004)                                                   
 1004 FORMAT(1X,'*** ERROR IN GRID: R(2) IS TOO SMALL.')                
        STOP                                                            
      ENDIF                                                             
      DO 3 I=2,NP                                                       
      CI=C-I                                                            
      RL=RR                                                             
      RU=RL                                                             
    1 RU=RU+RU                                                          
      FU=A*RU+B*DLOG(RU)+CI                                     
      IF(FU.LT.0.0D0) GO TO 1                                           
    2 RR=0.5D0*(RU+RL)                                                  
      FR=A*RR+B*DLOG(RR)+CI                                             
      IF(FR.GT.0.0D0) THEN                                              
        RU=RR                                                           
      ELSE                                                              
        RL=RR                                                           
      ENDIF                                                             
      IF(RU-RL.GT.1.0D-15*RR) GO TO 2                                   
      R(I)=RR                                                           
    3 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C  **************************************************************       
C                        SUBROUTINE ERRSPL                              
C  **************************************************************       
      SUBROUTINE ERRSPL(ERR,X,Y,N)                                      
C                                                                       
C     THIS SUBROUTINE ESTIMATES THE ERROR INTRODUCED BY NATURAL 
C  CUBIC SPLINE INTERPOLATION IN A TABLE X(I),Y(I) (I=1,...,N).         
C  THE INTERPOLATION ERROR IN THE VICINITY OF X(K) IS APPROXIMA-        
C  TED BY THE DIFFERENCE BETWEEN Y(K) AND THE VALUE OBTAINED FROM       
C  THE SPLINE THAT INTERPOLATES THE TABLE WITH THE K-TH POINT RE-       
C  MOVED. ERR IS THE LARGEST RELATIVE ERROR ALONG THE TABLE.            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      PARAMETER (NDIM=800,NPPG=NDIM+1,NPTG=NDIM+NPPG)                   
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                
      COMMON/STORE/F(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)              
      DIMENSION X(NDIM),Y(NDIM)                                         
      ERR=0.0D0                                                         
      N1=N-1                                                            
      DO 2 I=2,N1                                                       
      DO 1 J=1,N1                                                       
      IF(J.LT.I) THEN                                                   
        R(J)=X(J)                                                       
        F(J)=Y(J)                                                       
      ELSE                                                              
        R(J)=X(J+1)                                                     
        F(J)=Y(J+1)                                                     
      ENDIF                                                             
    1 CONTINUE                                                          
      CALL SPLINE(R,F,A,B,C,D,0.0D0,0.0D0,N1)                           
      RC=X(I)                                                           
      YI=A(I-1)+RC*(B(I-1)+RC*(C(I-1)+RC*D(I-1)))                       
      IF(DABS(Y(I)).GT.1.0D-3) THEN                                     
        ERRP=1.0D0-YI/Y(I)                                              
      ELSE                                                              
        ERRP=YI-Y(I)                                                    
      ENDIF                                                             
      ERR=DMAX1(ERR,DABS(ERRP))                                         
C     WRITE(6,'(1X,I3,1P,3E18.10,2E9.1)') I,X(I),Y(I),YI,ERRP,ERR       
    2 CONTINUE                                                          
      RETURN                                                            
      END  
