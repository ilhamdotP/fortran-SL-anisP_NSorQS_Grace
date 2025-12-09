       PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C      TOV, QS MIT Bag Modified
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN    
      DIMENSION YA(10), EK(4,10), Y(10)

c      OPEN (unit=2,STATUS='unknown',FILE='CatatQS.dat')
       OPEN (unit=3,STATUS='unknown',FILE='RadmassQS_aniso-2.dat')
       OPEN (unit=1,STATUS='unknown',FILE='Profiles_anisoQS-2.dat')

C     IM = NUMBER OF EQUATIONS Y(I)=Pressure, Y(2)=NS Mass and Y(3)=E density
      HC  = 197.327D0
      IM=3
      IN=IM-1
c        DO 10 IL=1,20,1
c        PCC=1.D-2*IL         

      DO 10 IL=1,1801,100
      PCC=1.0D0*IL
c
c    Profile/Direct Urca
c 
c      PCC=300.0D0
c       
      Y(1)=PCC

      Y(2)=0.1D-8

      P0=Y(1)
c--------------------------------------------------------------------------
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      


      ED=FED(P0)
c---------------------------------------------------------------------
      EDC=ED
      Y(3)=EDC
    

     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D-1
      NS=32
      XL=30.0D3

      H=PU/NS
      XP=1.0D-3
      HH=H/(2.0D0)

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
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
c---------------------------------------------------------------------

         YA(3)=ED

         XA=XM

         CALL FUNCT(EK,J,YA,XA,H,PCC,EDC)

C    COMPUTE K3, L3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
                 
c---------------------------------------------------------------------- 
            YA(3)=ED          

         XA=XM


         CALL FUNCT(EK,J,YA,XA,H,PCC,EDC)

C    COMPUTE K4, L4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
c---------------------------------------------------------------------
            YA(3)=ED
     

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H,PCC,EDC)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)

c--------------------------------------------------------------------------
c     Call Presure vs energy density parametrization (EOS)
c     SSS speed of sound squared, TSSS tangential speed of sound squared          
c---------------------------------------------------------------------      

        ED=FED(P0)

c---------------------------------------------------------------------  
          Y(3)=ED
c-------------------------------------------------------------------
          PRT=FPT(ED,PCC,EDC)
          PRESS=Y(1)
          SSS=DPDE(PRESS)
          TSSS=DPTDE(ED)
       END DO

c      WRITE(*,*)IL,LI,(XP/1.D3),Y(1),Y(2),Y(3)
C       WRITE(1,*)(XP/1.D3),Y(1),Y(2),Y(3),PRT,SSS,TSSS 
       WRITE(1,*)(XP/1.D3),Y(1),Y(2),Y(3),(Y(1)-PRT)     ! print profil
          

       PS=Y(1)
       PMIN=1.0D-8
c       PMIN=1.0D-4
     
      IF (PS .GT. PMIN  ) GOTO 28
    


c      WRITE(2,*)IL,(XP/1.D3), (Y(I),I=1,IM)
c      WRITE(3,*)IL,PCC,EDC,(XP/1.D3),Y(2)    
      WRITE(3,*)IL,(EDC/1.D3),(XP/1.D3),Y(2),
     &             2.0D0*GS*Y(2)*MSS/XP    
      

c     GS gravitation constant, MSS solar mass      
      GS=1.325D-12
      MSS=1.1155D15      
      PRINT*,PCC,(XP/1.D3),Y(2),2.0D0*GS*Y(2)*MSS/XP

 
   

 10   CONTINUE
      
      STOP
      END

c     function speed of sound squared (numerical derivative  pressure respect to eden)
      
      FUNCTION DPDE(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL FED    
      h=1.D-2
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      DEDP=(FED(x4)-8.D0*FED(x2)+8.D0*FED(x1)-FED(x3))
      DEDP=DEDP/(12.D0*h)
      DPDE=1.0D0/DEDP
      RETURN
      END

c     function speed of sound squared (analytical derivative tangential  pressure respect to eden)    
      FUNCTION DPTDE(ED)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI=3.14159265358979D0
      HC  = 197.327D0
c      BagCT=57.0D0
      BagCT=90.0D0
      ms=100.0D0
c      a4T=0.7D0
      a4T=0.9D0
      ALPHA=ms**2/(3.D0*PI*DSQRT(a4T)*HC**(3/2))
      BETA=ms**4*(1.D0-1.D0/(a4T))/(12.D0*PI**2*HC**(3))
      GAMMA=8.D0*PI/(3.D0*ms**2*DSQRT(a4T)*HC**(3/2))*HC**3
      KAPPA=3.D0/(1.D0-1.D0/(a4T))
      DPTDE=1.0D0/3.D0-0.5D0*ALPHA/DSQRT(ED-BagCT)
     &        +0.5D0*BETA*KAPPA/(ED-BagCT)
      RETURN
      END      
      

      SUBROUTINE FUNCT(EK,J,YA,XA,H,PCC,EDC)
C     *********************************************************
C     DEFINES TOV EQUATIONS 
C
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(4,10), YA(10)
      PI  = 3.14159265358979D0
c----------------------------      
c     GS gravitation constant, MSS solar mass      
      GS=1.325D-12
      MSS=1.1155D15
      
c  EDEN energy density, PRESS pressure, MASST NS mass in solar mass, H dr, XA r in meter
      EDEN=YA(3)
      PRESS=YA(1)
      MASST=YA(2)
      PRESST=FPT(EDEN,PCC,EDC)
c------------------------------------------------------      
c     dP
      
      EK(J,1)=-GS*EDEN*MASST*MSS*H/(XA*XA)
     &        *(1.D0+PRESS/EDEN)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
     &     /(1.D0-2.D0*GS*MASST*MSS/XA)
     &    -2.D0*H*(PRESS-PRESST)/XA 
c     dm expression
      
      EK(J,2)=4.D0*PI*XA*XA*EDEN*H/MSS
      
C       print*,PRESS,(PRESS-PRESST)

      RETURN
      END

c--------------------------------------------------------------------------
c  Input to calculate TOV      
c     Energy density MIT BAG
c--------------------------------------------------------------------------      
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C      BagC=57.0D0
C     FED=3.D0*P0+4.0*BagC
      PI=3.14159265358979D0
      HC =197.327D0
      BagC=57.0D0
      ms=100.0D0
      a4=0.7D0
      ALPHA=ms**2/(3.D0*PI*DSQRT(a4)*HC**(3/2))
      BETA=ms**4*(1.D0-1.D0/(a4))/(12.D0*PI**2*HC**(3))
      GAMMA=8.D0*PI/(3.D0*ms**2*DSQRT(a4)*HC**(3/2))*HC**3
      KAPPA=3.D0/(1.D0-1.D0/(a4))
      DEPS=1.5D0*DSQRT(3.D0)*ALPHA*DSQRT(3.D0*ALPHA**2+4.D0*(BagC+P0))
     &     +4.5D0*ALPHA**2
      PTD=P0
     &     -BETA*(1.D0+KAPPA*DLOG(GAMMA*DSQRT(3.D0*(P0+BagC)+DEPS)))
      FED=1.5D0*(DSQRT(3.0D0)*
     &     ALPHA*DSQRT(3.D0*ALPHA**2+4.D0*(BagC+PTD))+
     &     3.D0*ALPHA**2+2.D0*(BagC+PTD))+BagC
      RETURN
      END

      FUNCTION FPT(ED,PCC,EDC)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
c     BagCT in MeV fm^-3, a4T dimensionless and ms in MeV. HC is habarc
c     convertor from MeV into fm
      PI=3.14159265358979D0
      HC  = 197.327D0
      BagCT=57.0D0
      ms=100.0D0
      a4T=0.7D0
      ALPHA=ms**2/(3.D0*PI*DSQRT(a4T)*HC**(3/2))
      BETA=ms**4*(1.D0-1.D0/a4T)/(12.D0*PI**2*HC**3)
      GAMMA=8.D0*PI/(3.D0*ms**2*DSQRT(a4T)*HC**(3/2))*HC**3
      KAPPA=3.D0/(1.D0-1.D0/a4T)
     
      FPT=PCC+(ED-4.0D0*BagCT)/3.D0-ALPHA*DSQRT(ED-BagCT)
     &     +BETA*(1.D0+KAPPA*DLOG(GAMMA*DSQRT(ED-BagCT)))
     &     -(EDC-4.0D0*BagCT)/3.D0+ALPHA*DSQRT(EDC-BagCT)
     &     -BETA*(1.D0+KAPPA*DLOG(GAMMA*DSQRT(EDC-BagCT)))
    
      RETURN
      END 
 
c---------------------------------------------------------------------
C end of the code
c-------------------------------------------------------------------
