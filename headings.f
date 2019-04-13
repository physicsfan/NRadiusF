CCCCCCCCCCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     heading.f
C
C     A LIBRARY OF FUNCTIONS AND SUBROUTINES THAT THAT WILL PRINT
C     HEADINGS AND MENUS FOR THE NRADIUS PACKAGE.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      SUBROUTINE WELCOM()
      WRITE(6,*) 'WELCOME TO NRADIUS'
      WRITE(6,*) '------------------'
      WRITE(6,'(//)')
      END


      SUBROUTINE MENU()
      WRITE(6,*) '  SELECT ONE OPTION ...'                              
      WRITE(6,*) '    1: DIRAC ENERGY'            
      WRITE(6,*) '    2: ONE-LOOP QED'             
      WRITE(6,*) '    3: TWO-LOOP QED'                              
      WRITE(6,*) '    4: ONE-PHOTON EXCHANGE'                             
      WRITE(6,*) '    5: TWO-PHOTON EXCHANGE'            
      WRITE(6,*) '    6: THREE-PHOTON EXCHANGE'             
      WRITE(6,*) '    7: SCREENED QED'
      WRITE(6,*) '    8: NUCLEAR RECOIL'             
      WRITE(6,*) '    9: NUCLEAR POLARIZATION'
      WRITE(6,*) '    10: FULL CALCULATION'
      WRITE(6,*) '    11: EXIT'
      END


      SUBROUTINE POTMNU()
      WRITE(6,*) '  SELECT ONE POTENTIAL OPTION ...'                              
      WRITE(6,*) '    1: POINT NUCLEUS POTENTIAL'            
      WRITE(6,*) '    2: HARD-SPHERE POTENTIAL'             
      WRITE(6,*) '    3: FERMI POTENTIAL'
      WRITE(6,*) '    4: BACK'
      END
      
      


      
