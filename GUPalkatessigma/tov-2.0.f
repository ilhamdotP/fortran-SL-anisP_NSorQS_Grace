      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     GUP beta EOS + GUP Ctilde
C     modified at 28 September 2021, 17:46
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN, IK, IKK, IKJ
      DIMENSION YA(10), EK(4,10), Y(10)
      DOUBLE PRECISION LAMBDA, LMDCC
      CHARACTER(LEN=10)  FILECT, FILEYDY
      CHARACTER(LEN=50)  FILENAME1, FILENAME2

c---------------------------------------------------------------------
      CT=1.0D7      !CTILDE largest=1.0D7
C      YDY=-1.15D0 !-1.15
C      XC=5.0D-5   
      XC=1.0D-3     !RC
C---------------------------------------------------------------------    

CCC     write parameters to filename WITH YDY
C        WRITE(FILECT, '(D8.3)') CT
C        WRITE(FILEYDY, '(D8.3)') ABS(YDY)
C        FILENAME1 = 'profil_Y=-' // trim(adjustl(FILEYDY)) // 
C      &            '_CT=' // trim(adjustl(FILECT)) //'.dat'
C        FILENAME2 = 'radmass_Y=-' // trim(adjustl(FILEYDY)) // 
C      &            '_CT=' // trim(adjustl(FILECT)) //'.dat'
     
CCC     write parameters to filename WITHOUT YDY
       WRITE(FILECT, '(D8.3)') CT
       FILENAME1 = 'profil_CT=' // trim(adjustl(FILECT)) //'.dat'
       FILENAME2 = 'radmass_CT=' // trim(adjustl(FILECT)) //'.dat'
      


C        OPEN (unit=2,STATUS='unknown',FILE='CatatanB145SC.dat')
       OPEN (unit=3,STATUS='unknown',FILE=trim(FILENAME2))
       OPEN (unit=1,STATUS='unknown',FILE=trim(FILENAME1))
       
       

C     IM = NUMBER OF EQUATIONS Y(I)=Pressure, Y(2)=NS Mass and Y(3)=E density

      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
      IM=5
      IN=IM-1
      
      

       DO 10 IL=1,600,1
       FIXEDIL=300
C        IL=FIXEDIL
       PC=1.D0*IL
C        WRITE(*,*) IL
 
       YA(10)=CT
       YA(9)=YDY
       EDC=FED(PC)
       DEDPC=DEDP(PC)

       PCC=PC-2.D0/3.D0*PI*GS*XC*XC*(PC+EDC)*(3.D0*PC+EDC)
       MCC=4.D0*PI*XC*XC*XC*EDC/(3.D0*MSS)
       ALPCC=-2.D0/3.D0*PI*GS*XC*XC*(3.D0*PC+EDC)
       
       EDENTC=EDC
       PRESSTC=PCC
       ALPHAPC=(GS/(XC*XC))
     &   *(MCC*MSS+4.D0*PI*XC*XC*XC*PRESSTC)
       PRESSPC=-(EDC+PCC)*ALPHAPC
       MASSTPC=4.D0*PI*XC*XC*EDENTC
       CNST2=(EDC+PCC)*(1.D0-DEDPC)
C        SIGMAP=-YDY*(2.D0*GS*MSS*MCC/XC*XC)*PCC
C      &       +YDY*(2.D0*GS*MSS*MCC/XC)*PRESSPC
C      &       +YDY*(2.D0*GS*PCC/XC)*MASSTPC
       SIGMAP=DSDP(PCC)*PRESSPC
       LMDCC=1.D0/3.D0*(-2.D0*CNST2*GS*PI*XC*XC*
     &      (3.D0*PCC*PCC+4.D0*PC*EDC+EDC*EDC)) 
     &       +2.D0*(EDC+PCC)*SIGMAP
     
C        PRINT*, SIG(PCC),DSDP(PCC)
C        STOP
       
       IK=0
       IKK=0
       IKJ=0
       
       
 29   Y(1)=PCC
      Y(2)=MCC
      Y(3)=ALPCC
      Y(4)=LMDCC
      
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      
      P0=Y(1)
 
      Y(5)=FED(P0)
    

     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D0
      NS=8
      !XL=30.0D3

      H=PU/NS
c     XP should be larger than XC in order to avoid the unphysical behavior
c      near center!!    
      XP=1.0D0
      HH=H/(2.0D0)

      IF (XP.LT.XC) THEN
            WRITE(*,*) "XP=",XP," is NOT larger than XC=",XC
            STOP
      END IF

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH

C    COMPUTE K1, L1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB
   
         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K2, L2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
            
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      


         YA(5)=FED(P0)

         XA=XM

         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K3, L3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  
            YA(5)=FED(P0)          

         XA=XM


         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K4, L4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)
         
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

 

            YA(5)=FED(P0)
     

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)
         
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  


          Y(5)=FED(P0)
          
       END DO
       
       XA=XP
       PRESS=Y(1)
       MASST=Y(2)
       ALPHA=Y(3)
       LAMBDA=Y(4)
       EDEN=Y(5)
C        SIGMA=YDY*(2.D0*GS*MSS*MASST/XA)*PRESS
       SIGMA=SIG(PRESS)

C       WRITE(*,*)IL,(XP/1.D3),PRESS,MASST!,
C      &      (ALPHA-0.5D0*LOG(1.D0-2.D0*GS*MSS*MASST/XP)),
C      &      LAMBDA*1.D-6,EDEN,SIGMA
      
      
      IF (IL.EQ.FIXEDIL .AND. IK.EQ.1) THEN 
            WRITE(1,*)(XP/1.D3),PRESS,MASST,
     &      (ALPHA-0.5D0*LOG(1.D0-2.D0*GS*MSS*MASST/XP)),
     &      LAMBDA*1.D-6,EDEN,SIGMA
      ENDIF
     

       PS=Y(1)
       PMIN=1.0D-8

      IF (PS .GT. 1.1*PCC  ) THEN
            IKJ=IKJ+1
            GOTO 12
      ELSE IF (MASST .LT. 0.D0  ) THEN
            IKJ=IKJ+1
            GOTO 13
      ENDIF
      IF (PS .GT. PMIN  ) GOTO 28
      
      IF (IK.EQ.0) THEN
            IKK=IKK+1
C             WRITE(*,*)IL,(ALPHA-0.5D0*LOG(1.D0-2.D0*GS*MSS*Y(2)/XP))
C      &       ,LAMBDA,IKK
            IF (IKK.EQ.100) GOTO 11
      ENDIF
      
      IF (ABS(ALPHA-0.5D0*LOG(1.D0-2.D0*GS*MSS*Y(2)/XP)).GE.1.D-3) THEN
C             PRINT *,"Shooting alphac"
            ALPCC=ALPCC-(ALPHA-0.5D0*LOG(1.D0-2.D0*GS*MSS*Y(2)/XP))
            GOTO 29
      ENDIF
      
      IF (ABS(LAMBDA).GE.1.D4) THEN
C             PRINT *,"Shooting Lambdac"
            LMDCC=LMDCC-LAMBDA
            GOTO 29
      ENDIF
      IF (IK.EQ.0) THEN
C             PRINT *,"Finish both with IKJ=",IKK
            WRITE(*,*)IL,(XP/1.D3),MASST,
     &             2.0D0*GS*Y(2)*MSS/XP,IKK
            WRITE(3,*)IL,(EDC/1.D3),(XP/1.D3),MASST,
     &             2.0D0*GS*Y(2)*MSS/XP,LMDCC*1.D-6,ALPCC
            
            IF (IL.EQ.FIXEDIL) THEN
                  IK=IK+1
C                   PRINT *,"IK=",IK 
                  PRINT *,"Finish print profile data for PCC =",PCC
                  GOTO 29
            ENDIF
      ENDIF

c      WRITE(2,*)IL,(XP/1.D3), (Y(I),I=1,IM)
C         WRITE(*,*)PC,(EDC/1.D3),(XP/1.D3),Y(2),
C      &  2.0D0*GS*Y(2)*MSS/XP  
        
       IF (IKJ.EQ.0) GOTO 14
 
  11   WRITE(*,*) "Lambda not converging to zero. SKIP."
  12   WRITE(*,*) "PCC=",PCC,"P>1.1PC. SKIP"
  13   WRITE(*,*) "PCC=",PCC,"m<0. SKIP"
  
  14   IKJ=0

  10   CONTINUE
      
      
      END
      
      
C=================================================================
C         INCLUDE "EOS.f"


      SUBROUTINE FUNCT(EK,J,YA,XA,H)
C     *********************************************************
C     DEFINES TOV EQUATIONS 
C
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(4,10), YA(10)
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
c--------------------------------------------------------------

      CT=YA(10)
      EDEN=YA(5)
      PRESS=YA(1)
      MASST=YA(2)
      ALPHA=YA(3)
      LAMBDA=YA(4)
      
      YDY=YA(9)
C       SIGMA=YDY*(2.D0*GS*MSS*MASST/XA)*PRESS
      SIGMA=SIG(PRESS)
      
      EDENT=EDEN+2.D0*PI*GS*CT
     &        *(LAMBDA-1.5D0*(EDEN+PRESS)**2
     &        +2.D0*SIGMA*(EDEN+PRESS))

      PRESST=PRESS-2.D0*PI*GS*CT
     &        *(LAMBDA+0.5D0*(EDEN+PRESS)**2
     &        -2.D0*SIGMA*(EDEN+PRESS))
      
      ALPHAP=(GS*MASST*MSS/(XA*XA))
     &   /(1.D0-2.D0*GS*MASST*MSS/XA)
     &   *(1.D0+4.D0*PI*XA*XA*XA/(MASST*MSS)*PRESST)
      
      PRESSP=-(EDEN+PRESS)*ALPHAP
     &       -2.D0*SIGMA/XA
     
      MASSTP=4.D0*PI*XA*XA*EDENT
      
C       SIGMAP=-YDY*(2.D0*GS*MSS*MASST/(XA*XA))*PRESS
C      &       +YDY*(2.D0*GS*MSS*MASST/XA)*PRESSP
C      &       +YDY*(2.D0*GS*PRESS/XA)*MASSTP
      SIGMAP=DSDP(PRESS)*PRESSP
      
      LAMBDAP=PRESSP*((EDEN+PRESS)
     &       *(1.D0-DEDP(PRESS))
     &       +2.D0*SIGMA*(DEDP(PRESS)-1.D0))
     &       +8.D0*SIGMA*(EDEN+PRESS-SIGMA)/XA
     &       +2.D0*SIGMAP*(EDEN+PRESS)
     
      EK(J,1)=H*PRESSP
      
      EK(J,2)=H/MSS*MASSTP
      
      EK(J,3)=H*ALPHAP
      
      EK(J,4)=H*LAMBDAP
     
C       WRITE(*,*)CT,EDEN,PRESS,MASST,ALPHA,LAMBDA,SIGMA
C       WRITE(*,*)EK(J,1),EK(J,2),EK(J,3),EK(J,4)
C       STOP

      RETURN
      END

 
C-------------------------------------------------
C     MIT BAG B=145 MeV^4
c-------------------------------------------------
  
      
C       FUNCTION FED(P0)
C       IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      
C       HC  = 197.327D0
C       B=(145.D0)**4  ! PC until 800
C       !B=(185.D0)**4  ! PC until 1600
C       BMF=B/(HC*HC*HC)
C       FED= 3.D0*P0+4.D0*BMF
   
C       RETURN
C       END
c-----------------------------------------------------------------------  
      
C-----------------------------------------------------------------------
C     DE/DP AS A FUNCTION OF (P0) for NS (ANTO'S VERSION)
C-----------------------------------------------------------------------
      FUNCTION DEDP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL FED
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
   
   
      DEDP = (FED(x4)-8.D0*FED(x2)+8.D0*FED(x1)-FED(x3))/(12.D0*h)
c      WRITE(*,*)xa,DEDP
      RETURN
      END
 
c---------------------------------------------------------------------    

C-----------------------------------------------------------------------
C     DS/DP AS A FUNCTION OF (P0) for NS (ANTO'S VERSION)
C-----------------------------------------------------------------------
      FUNCTION DSDP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL SIG
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      
      DSDP=(SIG(x4)-8.D0*SIG(x2)+8.D0*SIG(x1)-SIG(x3))/(12.D0*h)
      
      RETURN
      END



c---------------------------------------------------------------------       
c     GUP EOS BETA=1.D-8
c---------------------------------------------------------------------

      FUNCTION SIG(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      
      SIG = 0.00893216184253268D0 + 0.0066138230853899135D0*P0 - 
     -  0.00002547856239766327D0*P0**2 + 1.227701424673366D-7*P0**3 - 
     -  2.525768285081147D-10*P0**4 + 2.5092894933230487D-13*P0**5 - 
     -  9.644098826004996D-17*P0**6
   
      RETURN
      END
      
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      
       IF ( P0 .GT. 50.D0 ) THEN
	 
     	FED= 253.2167266972223D0 + 2.042574552844379D0*P0
   
      ELSE IF ( P0 .GT.  0.5658D0 .AND. P0 .LE. 50.D0 ) THEN
 
       FED= 34.996719907638656D0 + 217.38861380763586D0*P0 - 
     -  197.7847221323045D0*P0**2 + 
     -  100.96604597400824D0*P0**3 - 31.000998206090017D0*P0**4 + 
     -  6.201685436984379D0*P0**5 - 0.8545436764597708D0*P0**6 + 
     -  0.08425535863251823D0*P0**7 - 0.006099401533513903D0*P0**8 + 
     -  0.000329567312119461D0*P0**9 - 
     -  0.000013407358939926027D0*P0**10 + 
     -  4.112479560545507D-7*P0**11 - 9.450222969553151D-9*P0**12 + 
     -  1.6006796160515572D-10*P0**13 - 1.9380539806762782D-12*P0**14 + 
     -  1.5862872432492548D-14*P0**15 - 7.861083316562445D-17*P0**16 + 
     -  1.7807392264595506D-19*P0**17
        
       ELSE IF (P0 .GT. 4.99313436D-4 .AND. P0 .LE. 0.5658D0) THEN
   

       FED= 0.25227973223585715D0 + 785.6239013862972D0*P0 - 
     -  7097.193231748922D0*P0**2 + 
     -  44787.78752904974D0*P0**3 - 170674.67647877714D0*P0**4 + 
     -  375798.1030069276D0*P0**5 - 438222.42374165234D0*P0**6 + 
     -  208628.4521716856D0*P0**7

       ELSE
     
       FED= 0.00020663104786823767D0 + 985.7550962048155D0*P0 - 
     -  6.452649548410661D6*P0**2 + 4.045493650683346D10*P0**3 - 
     -  1.2422017897554195D14*P0**4 + 1.765493095354902D17*P0**5 - 
     -  9.1529417012139D19*P0**6
 
      
       END IF
   
      RETURN
      END
      
      
      