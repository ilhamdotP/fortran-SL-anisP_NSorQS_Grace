      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     TEST1 linear EOS+corr initial boundary 28 Mei 2020
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN    
      DIMENSION YA(10), EK(4,10), Y(10)


c      OPEN (unit=2,STATUS='unknown',FILE='CatatanB145SC.dat')
       OPEN (unit=3,STATUS='unknown',
     & FILE='radmass-TOVaniso.dat')
       OPEN (unit=1,STATUS='unknown',
     & FILE='profil-TOVaniso.dat')

C     IM = NUMBER OF EQUATIONS Y(I)=Pressure, Y(2)=NS Mass and Y(3)=E density

      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
      IM=3
      IN=IM-1
c---------------------------------------------------------------------
c     XL=1.0D-4
C      XL=1.0D-2
C      ALPHA=5.0D-3
C      XC=XL*ALPHA      
      XC=1.0D-3
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      OPEN (unit=2,STATUS='unknown',FILE='fedPN.dat')
      OPEN (unit=4,STATUS='unknown',FILE='fedORI.dat')
      ! DO 100 IJK=1,1000,1
      ! P0=IJK*1.D0
      ! ED=3.07D0*(IJK+5)+4.D0*57.D0
      ! WRITE(2,*) P0,FED(P0),FEDN(P0)
      ! WRITE(4,*) PED(ED),ED
  ! 100 CONTINUE 
      ! STOP
c---------------------------------------------------------------------  
    
       ! DO 10 IL=5,500,5  !IL=5,500,5
       FIXEDIL=1000
       IL=FIXEDIL
       PC=1.D0*IL
 
       YA(10)=XL
       
       EDC=FED(PC)
       
       PCC=PC-2.D0/3.D0*PI*GS*XC*XC*(PC+EDC)*(3.D0*PC+EDC)
       MCC=4.D0*PI*XC*XC*XC*EDC/(3.D0*MSS)
       EDC=FED(PCC)

       
       
      Y(1)=PCC

      Y(2)=MCC
c-----------------------------------------------------------------------
c     Y(1)=PC 
c      Y(2)=1.0D-8
      
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      
      P0=Y(1)
 
      Y(3)=FED(P0)
    

     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D0
      NS=32
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


         YA(3)=FED(P0)

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

  
            YA(3)=FED(P0)          

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

 

            YA(3)=FED(P0)
     

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

  


          Y(3)=FED(P0)
          
       END DO
       
       XA=XP
       PRESS=Y(1)
       MASST=Y(2)
       EDEN=Y(3)
	   SIGMA=SIG(PRESS,EDEN,PCC,EDC)
       
       

      !WRITE(*,*)IL,LI,(XP/1.D3),Y(1),Y(2),Y(3)
      ! IF (IL.EQ.FIXEDIL) WRITE(1,*)IL,LI,(XP/1.D3),Y(1),Y(2),TPA
	  IF (IL.EQ.FIXEDIL) THEN
        WRITE(1,*)(XP/1.D3),PRESS,MASST
C      &      ,(ALPHA-0.5D0*LOG(1.D0-2.D0*GS*MSS*MASST/XP))
C      &      ,LAMBDA*1.D-6
     &      ,(EDEN/1.D3),(SIGMA*1.D0)
        WRITE(2,*) PRESS,FED(PRESS),FEDN(PRESS)
        ED=FED(PRESS)
        BagCT=57.0D0
        W=1.D0/3.D0
        WRITE(4,*) PED(ED),ED,PED(ED)/W+4.D0*BagCT
      ENDIF

     
C       WRITE(*,*)SIGMA,PRESS,EDEN,PCC,EDC
     

       PS=Y(1)
       PMIN=1.0D-8
c       PMIN=2.0D-5
      
      IF (MCC.LT.0.0D0) THEN
            WRITE(*,*) "MCC negative"
            WRITE(*,*) Y(2)
            GOTO 10
      ELSE IF (ABS(PS).GT.10*ABS(PCC)) THEN
            WRITE(*,*) "ABS(P) =",ABS(PS),"> 10 ABS(PCC) =",10*ABS(PCC)
            GOTO 10
      ELSE IF (LI.GE.1000000) THEN
            WRITE(*,*) "Max iteration = 1000000 reached"
            WRITE(*,*)"p=",Y(1),"m=",Y(2),"rho=",Y(3)
            GOTO 10
      END IF

      IF (PS .GT. PMIN  ) GOTO 28
    


c      WRITE(2,*)IL,(XP/1.D3), (Y(I),I=1,IM)
        WRITE(3,*)IL,(EDC/1.D3),(XP/1.D3),MASST,
     &             2.0D0*GS*Y(2)*MSS/XP     
c       write(*,*)X
        WRITE(*,*)IL,(XP/1.D3),MASST,
     &             2.0D0*GS*Y(2)*MSS/XP
 
   

  10   CONTINUE
      
      
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
      GS=1.325D-12
      MSS=1.1155D15
c--------------------------------------------------------------
c     XL=1.0D-4
      XL=YA(10)
      EDEN=YA(3)
      PRESS=YA(1)
      MASST=YA(2)

       XX=(2.D0*GS*MASST*MSS/XA)
     &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
     &   /(1.D0-2.D0*GS*MASST*MSS/XA)
c-------------------------------------------------------------      
c       TPA=0.0D0
c       TMA=1.0D0
c      OPEN (unit=7,STATUS='unknown',FILE='cekTMAXL10mm.dat')
c      write(7,*)TPA,FF,GG,TMA
c-------------------------------------------------------------      
  
      SIGMA=SIG(PRESS,EDEN,PCC,EDC) 
      EK(J,1)=-EDEN*(1.D0+PRESS/EDEN)*XX/(2.D0*XA)*H
     &       -2.D0*SIGMA/XA*H
      
      EK(J,2)=4.D0*PI*XA*XA*EDEN*H/MSS 

      RETURN
      END

 
c---------------------------------------------------------------------       
c     GUP EOS Grace
c----------------------------------------------------------------
  
      
      FUNCTION SIG(P0,ED,PCC,EDC)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      
      PTAN=FPT(ED,PCC,EDC)
C       SIG=P0-PCC      !good approx from mathematica
      SIG=P0-PTAN      !aniso
C       SIG=0.D0          !iso
      
C       PRINT*,"Grace",P0,ED,SIG
   
      RETURN
      END
      
      FUNCTION FEDN(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI=3.14159265358979D0
      HC =197.327D0
      BagC=57.0D0
      ms=100.0D0
      a4=0.7D0
      W=1.D0/3.D0
      
      LALPHA=3.D0*PI*DSQRT(a4)*HC**(1.5D0)
      LBETA=(1.D0-1.D0/a4)/(12.D0*PI*PI*HC**(3.0D0))
      LGAMMA=3.D0*ms**(2.0D0)*DSQRT(a4)
      ALPHA=ms**2/LALPHA
      BETA=ms**4*LBETA
      GAMMA=8.D0*PI*HC**(1.5D0)/LGAMMA
      KAPPA=3.D0/(1.D0-1.D0/a4)
     
      DSQ=-DSQRT(ALPHA**(2.0D0)/W+4.D0*(3.D0*W*BagC+P0))  
      DEPS=DSQRT(W)/(2.D0*W**2.D0)*ALPHA*DSQ
     &     +ALPHA**(2.0D0)/(2.D0*W**2.D0)
      
      DLO=DLOG(GAMMA*DSQRT(P0/W+3.D0*BagC+DEPS))      
      BPTD=BETA*(1.D0+KAPPA*DLO)
      PTD=P0-BPTD
      
      ALSQ=-ALPHA*DSQRT(ALPHA**(2.0D0)/W+4.D0*(3.D0*W*BagC+PTD))
     
      FEDN=1.D0/(2.D0*W)*(ALSQ/DSQRT(W)+(ALPHA**2.D0)/W
     &     +2.D0*(3.D0*W*BagC+PTD))
     &    +BagC
     
C       FEDN=P0/W+4.D0*BagC
      RETURN
      END
      
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI=3.14159265358979D0
      HC =197.327D0
      BagC=57.0D0
      ms=100.0D0
      a4=0.7D0
      W=1.D0/3.D0
      
      LALPHA=3.D0*PI*DSQRT(a4)*HC**(1.5D0)
      LBETA=(1.D0-1.D0/a4)/(12.D0*PI*PI*HC**(3.0D0))
      LGAMMA=3.D0*ms**(2.0D0)*DSQRT(a4)
      ALPHA=ms**2/LALPHA
      BETA=ms**4*LBETA
      GAMMA=8.D0*PI*HC**(1.5D0)/LGAMMA
      KAPPA=3.D0/(1.D0-1.D0/a4)
     
      DSQ=+DSQRT(ALPHA**(2.0D0)/W+4.D0*(3.D0*W*BagC+P0))  
      DEPS=DSQRT(W)/(2.D0*W**2.D0)*ALPHA*DSQ
     &     +ALPHA**(2.0D0)/(2.D0*W**2.D0)
      
      DLO=DLOG(GAMMA*DSQRT(P0/W+3.D0*BagC+DEPS))      
      BPTD=BETA*(1.D0+KAPPA*DLO)
      PTD=P0-BPTD
      
      ALSQ=+ALPHA*DSQRT(ALPHA**(2.0D0)/W+4.D0*(3.D0*W*BagC+PTD))
     
      FED=1.D0/(2.D0*W)*(ALSQ/DSQRT(W)+(ALPHA**2.D0)/W
     &     +2.D0*(3.D0*W*BagC+PTD))
     &    +BagC
     
C       FED=P0/W+4.D0*BagC
      RETURN
      END
      
      FUNCTION PED(ED)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI=3.14159265358979D0
      HC =197.327D0
      BagC=57.0D0
      ms=100.0D0
      a4=0.7D0
      W=1.D0/3.D0
      
      LALPHA=3.D0*PI*DSQRT(a4)*HC**(1.5D0)
      LBETA=(1.D0-1.D0/a4)/(12.D0*PI*PI*HC**(3.D0))
      LGAMMA=3.D0*ms**2*DSQRT(a4)
      ALPHA=ms**2/LALPHA
      BETA=ms**4*LBETA
      GAMMA=8.D0*PI*HC**(1.5D0)/LGAMMA
      KAPPA=3.D0/(1.D0-1.D0/a4)
      
      PED=(ED-4.0D0*BagC)*W-ALPHA*DSQRT(ED-BagC)
     &     +BETA*(1.D0+KAPPA*DLOG(GAMMA*DSQRT(ED-BagC)))
      
      RETURN
      END
      
      
c-----------------------------------------------------------------------     

      FUNCTION FPT(ED,PCC,EDC)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
c     BagCT in MeV fm^-3, a4T dimensionless and ms in MeV. 
c     HC is habarc convertor from MeV into fm
      PI=3.14159265358979D0
      HC  = 197.327D0
      BagCT=57.0D0
      ms=100.0D0
      a4T=0.7D0
      W=1.D0/3.D0
      
      LALPHA=3.D0*PI*DSQRT(a4T)*HC**(1.5D0)
      LBETA=(1.D0-1.D0/a4T)/(12.D0*PI*PI*HC**(3.D0))
      LGAMMA=3.D0*ms**2*DSQRT(a4T)
      ALPHA=ms**2/LALPHA
      BETA=ms**4*LBETA
      GAMMA=8.D0*PI*HC**(1.5D0)/LGAMMA
      KAPPA=3.D0/(1.D0-1.D0/a4T)
     
      FPT=PCC+(ED-4.0D0*BagCT)*W-ALPHA*DSQRT(ED-BagCT)
     &     +BETA*(1.D0+KAPPA*DLOG(GAMMA*DSQRT(ED-BagCT)))
     &     -(EDC-4.0D0*BagCT)*W+ALPHA*DSQRT(EDC-BagCT)
     &     -BETA*(1.D0+KAPPA*DLOG(GAMMA*DSQRT(EDC-BagCT)))
     
C       print*,"Grace",ED,FPT
    
      RETURN
      END  
