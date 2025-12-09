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
      DOUBLE PRECISION  LAMBDA, LMDCC
      CHARACTER(LEN=10)  FILECT, FILEYDY
      CHARACTER(LEN=50)  FILENAME1, FILENAME2

c---------------------------------------------------------------------
      CT=1.D6      !CTILDE largest=1.0D7
      YDY=0.D0 !-1.15
C      XC=5.0D-5   
      XC=1.0D-3     !RC
C---------------------------------------------------------------------    
     
CCC     write parameters to filename WITHOUT YDY
       WRITE(FILECT, '(D8.3)') CT
       FILENAME1 = 'profil_CT=' // trim(adjustl(FILECT)) 
     & //'_aniso.dat'
       FILENAME2 = 'radmass_CT=' // trim(adjustl(FILECT)) 
     & //'_aniso.dat'
      


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
      
      

       DO 10 IL=1,1000,1  !IL=5,500,5
       FIXEDIL=300
C        IL=FIXEDIL
       PC=1.D0*IL
C        WRITE(*,*) IL
 
       YA(10)=CT
       YA(9)=YDY
       PCC=PC-2.D0/3.D0*PI*GS*XC*XC*(PC+EDC)*(3.D0*PC+EDC)
       EDC=FED(PCC)
       ! print*,PC,PCC,EDC
       MCC=4.D0*PI*XC*XC*XC*EDC/(3.D0*MSS)
       ALPCC=-2.D0/3.D0*PI*GS*XC*XC*(3.D0*PC+EDC)
       LMDCC=1.d-12
       
       
C        PRINT*, SIGMAP
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

         CALL FUNCT(EK,J,YA,XA,H,PCC,EDC)

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

         CALL FUNCT(EK,J,YA,XA,H,PCC,EDC)

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


         CALL FUNCT(EK,J,YA,XA,H,PCC,EDC)

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

         CALL FUNCT(EK,J,YA,XA,H,PCC,EDC)

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
       SIGMA=SIG(PRESS,EDEN,PCC,EDC)
      
      
      IF (IL.EQ.FIXEDIL .AND. IK.EQ.1) THEN     ! print data at PC=FIXEDIL
            WRITE(1,*)(XP/1.D3),PRESS,MASST
     &      ,(ALPHA-0.5D0*LOG(1.D0-2.D0*GS*MSS*MASST/XP))
     &      ,LAMBDA*1.D-6
     &      ,(EDEN/1.D3),(SIGMA*1.D0)
      ENDIF
     

       PS=Y(1)
       PMIN=1.0D-8

      IF (PS .GT. 1.1*PCC  ) THEN               !   p > 1.1 pc
            IKJ=IKJ+1
            GOTO 12
      ELSE IF (MASST .LT. 0.D0  ) THEN          !   m<0.
            IKJ=IKJ+1
            GOTO 13
      ENDIF
      IF (PS .GT. PMIN  ) GOTO 28               !   p=0.
      
      IF (IK.EQ.0) THEN
C             PRINT *,"not converging"
            IKK=IKK+1
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
C                   PRINT *,"Finish print profile data for PCC =",PCC
                  GOTO 29
            ENDIF
      ENDIF


        
       IF (IKJ.EQ.0) GOTO 14
 
  11   WRITE(*,*) "Lambda not converging to zero. SKIP."
  12   WRITE(*,*) "PCC=",PCC,"P>1.1PC. SKIP"
  13   WRITE(*,*) "PCC=",PCC,"m<0. SKIP"
  
  14   IKJ=0

  10   CONTINUE
      
      
      END
      
      
C=================================================================
C         INCLUDE "EOS.f"


      SUBROUTINE FUNCT(EK,J,YA,XA,H,PCC,EDC)
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
C       SIGMA=SIG(PRESS)
      SIGMA=SIG(PRESS,EDEN,PCC,EDC) 
      
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
        
C       DPT=1.D0/DEDP(PRESS)          ! ISO
      DPT=DPTDE(EDEN)             ! ANISO
      DSDE=1.D0/DEDP(PRESS)+DPT
      SIGMAP=DSDP(PRESS,EDEN,PCC,EDC)*PRESSP
     &       +DSDE*DEDP(PRESS)*PRESSP
      
      LAMBDAP=PRESSP*((EDEN+PRESS)
     &       *(1.D0-DEDP(PRESS))
     &       +2.D0*SIGMA*(DEDP(PRESS)-1.D0))
     &       +8.D0*SIGMA*(EDEN+PRESS-SIGMA)/XA
     &       +2.D0*SIGMAP*(EDEN+PRESS)
     
     
      EK(J,1)=H*PRESSP
      
      EK(J,2)=H/MSS*MASSTP
      
      EK(J,3)=H*ALPHAP
      
      EK(J,4)=H*LAMBDAP

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
      FUNCTION DSDP(xa,EDEN,PCC,EDC)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL SIG
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      
      DSDP=SIG(x4,EDEN,PCC,EDC)
      DSDP=DSDP-8.D0*SIG(x2,EDEN,PCC,EDC)
      DSDP=DSDP+8.D0*SIG(x1,EDEN,PCC,EDC)
      DSDP=DSDP-SIG(x3,EDEN,PCC,EDC)
      DSDP=DSDP/(12.D0*h)
      
      RETURN
      END
 


      FUNCTION SIG(P0,ED,PCC,EDC)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL FPT
      
      PTAN=FPT(ED,PCC,EDC) 
      SIG=P0-PTAN                       ! ANISO
C       SIG=0.D0                        ! ISO
C       SIG=P0-PCC      !good approx from mathematica
   
      RETURN
      END
      
      FUNCTION FPT(ED,PCC,EDC)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
c     BagCT in MeV fm^-3, a4T dimensionless and ms in MeV. 
c     HC is habarc convertor from MeV into fm
      PI=3.14159265358979D0
      HC  = 197.327D0
      BagCT=57.0D0
      ms=100.0D0
      a4T=0.7D0
      LALPHA=3.D0*PI*DSQRT(a4T)*HC**(1.5D0)
      LBETA=(1.D0-1.D0/a4T)/(12.D0*PI*PI*HC**(3.D0))
      LGAMMA=3.D0*ms**2*DSQRT(a4T)
      ALPHA=ms**2/LALPHA
      BETA=ms**4*LBETA
      GAMMA=8.D0*PI*HC**(1.5D0)/LGAMMA
      KAPPA=3.D0/(1.D0-1.D0/a4T)
     
      FPT=PCC+(ED-4.0D0*BagCT)/3.D0-ALPHA*DSQRT(ED-BagCT)
     &     +BETA*(1.D0+KAPPA*DLOG(GAMMA*DSQRT(ED-BagCT)))
     &     -(EDC-4.0D0*BagCT)/3.D0+ALPHA*DSQRT(EDC-BagCT)
     &     -BETA*(1.D0+KAPPA*DLOG(GAMMA*DSQRT(EDC-BagCT)))
    
      RETURN
      END 
      
      
c--------------------------------------------------------------------------
c  Input to calculate TOV      
c     Energy density MIT BAG
c--------------------------------------------------------------------------      
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
c      BagC=57.0D0
c      FED=3.D0*P0+4.0*BagC
      PI=3.14159265358979D0
      HC =197.327D0
      BagC=57.0D0
      ms=100.0D0
      a4=0.7D0
      LALPHA=3.D0*PI*DSQRT(a4)*HC**(1.5D0)
      LBETA=(1.D0-1.D0/a4)/(12.D0*PI*PI*HC**(3.0D0))
      LGAMMA=3.D0*ms**(2.0D0)*DSQRT(a4)
      ALPHA=ms**2/LALPHA
      BETA=ms**4*LBETA
      GAMMA=8.D0*PI*HC**(1.5D0)/LGAMMA
      KAPPA=3.D0/(1.D0-1.D0/a4)
     
c----------------------------------------------------------
c      ALPHA=0.0D0
c     BETA=0.0
c-----------------------------------------------------------
      DSQ=DSQRT(3.D0*ALPHA**(2.0D0)+4.D0*(BagC+P0))  
      DEPS=1.5D0*DSQRT(3.D0)*ALPHA*DSQ+4.5D0*ALPHA**(2.0D0)
      
      DLO=DLOG(GAMMA*DSQRT(3.D0*(P0+BagC)+DEPS))      
      BPTD=BETA*(1.D0+KAPPA*DLO)
      PTD=P0-BPTD
      
      ALSQ=ALPHA*DSQRT(3.D0*ALPHA**(2.0D0)+4.D0*(BagC+PTD))
     
      FED=1.5D0*(DSQRT(1.5D0)*ALSQ+3.D0*ALPHA**(2.0)+2.D0*(BagC+PTD))
     &    +BagC
      RETURN
      END
      
      

c     function speed of sound squared (analytical derivative tangential  pressure respect to eden)    
      FUNCTION DPTDE(ED)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI=3.14159265358979D0
      HC  = 197.327D0
      BagCT=57.0D0
c       BagCT=90.0D0
      ms=100.0D0
      a4T=0.7D0
C       a4T=0.9D0
      LALPHA=3.D0*PI*DSQRT(a4T)*HC**(1.5D0)
      LBETA=(1.D0-1.D0/a4T)/(12.D0*PI*PI*HC**(3.0D0))
      LGAMMA=3.D0*ms**(2.0D0)*DSQRT(a4T)
      ALPHA=ms**2/LALPHA
      BETA=ms**4*LBETA
      GAMMA=8.D0*PI*HC**(1.5D0)/LGAMMA
      KAPPA=3.D0/(1.D0-1.D0/a4T)
      
     
      DPTDE=1.0D0/3.D0-0.5D0*ALPHA/DSQRT(ED-BagCT)
     &        +0.5D0*BETA*KAPPA/(ED-BagCT)
      RETURN
      END      
      