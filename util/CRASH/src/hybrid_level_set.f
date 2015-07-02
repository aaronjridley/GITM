C----------------------------------------------------------------------C
C          Extended EULER EQUATIONS IN CONSERVATION FORM               C
C                                                                      C
C                            U=(R,RU,E,R*PSI,P)                        C
C                            -                                         C
C----------------------------------------------------------------------C
            PARAMETER (MX=1600)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(-1:MX+2)
        DIMENSION R(-1:MX+2),U(-1:MX+2)
        DIMENSION P(-1:MX+2),Z(-1:MX+2),G(-1:MX+2)
        DIMENSION IND(MX)
            COMMON /GAMMAS/GAM1,GAM0
            COMMON /GRISIZ/M,MM1,MP1,MP2
        
        OPEN(90,FILE="LEVELSET.in",STATUS="OLD")
            READ(90,*) NOTIST,CFLCOE,M
        CLOSE(90)
            CALL INDATA(M,X,R,U,P,Z,G,DX)
            N=0
            TIME=0.0

            DO 0001 N=1,NOTIST
              
                   CALL CFLTIME(MX,R,U,P,G,DX,DTMIN)
                DT=CFLCOE*DTMIN
                DTDX=DT/DX
               TIME=TIME+DT

               MM1=M-1
               MP1=M+1
               MP2=M+2

C        LOCATE THE INTERFACE

                DO 0002 I=1,M-1
                   IF (Z(I)*Z(I+1).LE.0.0) THEN
                       IND(I)=1
                          ELSE
                 IND(I)=0
                   END IF
 0002      CONTINUE
              
                   CALL ROE1D(N,R,U,P,Z,G,DTDX,IND)
                   PRINT*,N,TIME,DT

 0001        CONTINUE

        DO 0003 I=1,M
                    WRITE(14,102) X(i)
            WRITE(15,102) R(I)
            WRITE(16,102) U(I)
            WRITE(17,102) P(I)
            WRITE(18,102) Z(I)
 0003   CONTINUE

 100    FORMAT(1X,1A,I1,3A)
 101    FORMAT(1X,1A,I2,3A)
 102    FORMAT(1X,f15.7,2X)
        
            END
C----------------------------------------------------------------------C

            SUBROUTINE INDATA(M,X,R,U,P,Z,G,DX)
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (MX=1600) 
        DIMENSION X(-1:MX+2) 
        DIMENSION R(-1:MX+2),U(-1:MX+2)
        DIMENSION P(-1:MX+2),Z(-1:MX+2),G(-1:MX+2)            
                COMMON /GAMMAS/GAM1,GAM0
                
            DATA DOMLEN/1.0D0/
            DX=DOMLEN/REAL(M)

c weak shock hitting an interface 

        gam1=1.4
                gam0=1.67
                
         DO 0001 I=1,M
            R(I)=1.3333
            U(I)=0.3535
            P(I)=1.5
            G(I)=gam1
            IF(I.GT.M/4) THEN
               R(I)=1.0
               U(I)=0.0
               P(I)=1.0
            END IF
            IF(I.GT.M/2) THEN
               G(I)=gam0
               R(I)=0.1379                  
            END IF
 0001   CONTINUE
 
c 2-fluid shock tube
c 
c        gam1=1.4
C        gam0=1.2
c
c        DO 0001 I=1,M
c           R(I)=1.0
c           U(I)=0.0
c           P(I)=1.0
c           G(I)=GAM1
c           IF(I.GT.M/2) THEN
c              R(I)=0.125
c              U(I)=0.0
c              P(I)=0.1
c              G(I)=GAM0
c           END IF
c 0001   CONTINUE

             DO 0002 I=1,M
           X(I)=i*DX
                Z(I)=0.0-(I-M/2-0.5)/(M/2)
                IF (Z(I).GT.0.5) Z(I)= 0.5
                IF (Z(I).LT.-.5) Z(I)=-0.5
 0002   CONTINUE

            RETURN
            END

C----------------------------------------------------------------------C
        SUBROUTINE ROE1D(N,R,U,P,Z,G,DTDX,IND)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
            PARAMETER (MX=1600)
            DIMENSION S(4),EIG(5,4),ALF(4)
            DIMENSION CN(4,-1:MX+2),DWNP(5,4,-1:MX+2)
            DIMENSION NC(4,-1:MX+2),BLIM(4,-1:MX+2)
            DIMENSION W(5,-1:MX+2)
            DIMENSION R(-1:MX+2),U(-1:MX+2),H(-1:MX+2)
            DIMENSION P(-1:MX+2),Z(-1:MX+2),G(-1:MX+2)
            DIMENSION IND(MX)

            COMMON /BEEFUN/B,B1,B2,ACN
            COMMON /GAMMAS/GAM1,GAM0
            COMMON /GRISIZ/M,MM1,MP1,MP2

C        BOUNDARY CONDITIONS (0-order extrapolation)
            R(0)  =R(1)
            R(-1)   =R(0)
            U(0)   =U(1)
            U(-1)  =U(0)
            P(0)   =P(1)
                P(-1)  =P(0)
            Z(0)   =Z(1)
                Z(-1)  =Z(0)
            G(0)   =G(1)
            G(-1)  =G(0)
            R(MP1)  =R(M)
            R(MP2)  =R(MM1)
            U(MP1)  =U(M)
            U(MP2)  =U(MM1)
            P(MP1)  =P(M)
            P(MP2)  =P(MM1)
            Z(MP1)  =Z(M)
            Z(MP2)  =Z(MM1)
            G(MP1)  =G(M)
            G(MP2)  =G(MM1)
        
C        COMPUTE CONSERVED VARIABLES FOR DATA
        DO 0001 I=-1,MP2
               W(1,I)=R(I)
               W(2,I)=R(I)*U(I)
               W(3,I)=0.5*R(I)*U(I)*U(I)+P(I)/(G(I)-1.0)
               W(4,I)=R(I)*Z(I)
                   W(5,I)=P(I)
               H(I)  =(W(3,I)+P(I))/R(I)
 0001   CONTINUE

            DO 3003 I=-1,MP1

C  SOLVE RIGHT RIEMANN PROBLEM RP(I,I+1)

C  LEFT STATE
               RL=R(I)
               UL=U(I)
               PL=P(I)
               HL=H(I)
               ZL=Z(I)
               GL=G(I)

C  RIGHT STATE
               RR=R(I+1)
               UR=U(I+1)
               PR=P(I+1)
               HR=H(I+1)
               ZR=Z(I+1)
               GR=G(I+1)

C  GRADIENTS 
               DR=RR-RL
               DU=UR-UL
               DP=PR-PL
               DZ=ZR-ZL

C  MEAN VALUES ROE AVERAGES
               WT=SQRT(RR/RL)
               WTP1=WT+1.0
               RC=WT*RL
               UC=(UL+WT*UR)/WTP1
                   HC=(HL+WT*HR)/WTP1
               ZC=(ZL+WT*ZR)/WTP1
               QQ=0.5*UC*UC
               IF (DZ.EQ.0.0) THEN 
                   GC=GL
                   XC=0.0
               ELSE
               GC=((ZC-ZL)*GR+(ZR-ZC)*GL)/DZ
               PC=RC*(GC-1)*(HC-QQ)/GC
                   XC=(DP-(GC-1.0)*(PR/(GR-1.0)-PL/(GL-1.0)))/(RC*DZ)
               END IF
               AA=(GC-1.0)*(HC-QQ)
               AC=SQRT(AA)
           
C  WAVE SPEEDS 
               S(1)        =UC-AC
               S(2)        =UC
               S(3) =UC
               S(4)        =UC+AC

C  WAVE STRENGTHS
               QQAA=0.25/AA
               CONST        =RC*AC*DU
               ALF(1)        =QQAA*(DP-CONST)
           ALF(2)        =0.5*(DR-DP/AA)
                ALF(3)   =0.5*RC*DZ
               ALF(4)        =QQAA*(DP+CONST)
           
C  EIGENVECTORS
               UCAC=UC*AC
           
           EIG(1,1)=1.0
               EIG(2,1)=UC-AC
                  EIG(3,1)=HC-UCAC
                  EIG(4,1)=ZC
           EIG(5,1)=AA
           
               EIG(1,2)=1.0
               EIG(2,2)=UC
                  EIG(3,2)=QQ
                  EIG(4,2)=ZC
           EIG(5,2)=0.0

               EIG(1,3)=0.0
               EIG(2,3)=0.0
                  EIG(3,3)=-XC/(GC-1.0)
                  EIG(4,3)=1.0
           EIG(5,3)=0.0
        
               EIG(1,4)=1.0
               EIG(2,4)=UC+AC
                  EIG(3,4)=HC+UCAC
                  EIG(4,4)=ZC
           EIG(5,4)=AA
       
               DO 3002 IW=1,4
C      FIRST-ORDER EFFECTS OF (IW)TH WAVE
                  CN(IW,I)=S(IW)*DTDX
                  NC(IW,I)=SIGN(1.0D0 ,CN(IW,I))
                  ITARG   =I+(1+NC(IW,I))/2
                  DO 3001 IU=1,5
                     WINC        =CN(IW,I)*EIG(IU,IW)*ALF(IW)
                     W(IU,ITARG) =W(IU,ITARG)-2.0*WINC
                     DWNP(IU,IW,I)=(NC(IW,I)-CN(IW,I))*WINC
                     BLIM(IW,I)   =ALF(IW)*CN(IW,I)*(NC(IW,I)-CN(IW,I))
 3001         CONTINUE
 3002      CONTINUE
 3003           CONTINUE

C        SECOND ORDER EFECTS
               DO 4003 I=0,M
               DO 4002 IW=1,4
C  COMPUTE SECOND-ORDER EFFECTS OF (IW)TH WAVE
                   IS=I-NC(IW,I)
                   B1=BLIM(IW,I)
                   B2=BLIM(IW,IS)
               CALL SUPERBEE
c              CALL MINBEE
c                   B=0.0d0
                   DO 4001 IU=1,5
                               W(IU,I+1)=W(IU,I+1)+B*DWNP(IU,IW,I)
                       W(IU,I  )=W(IU,I  )-B*DWNP(IU,IW,I)
 4001          CONTINUE
 4002           CONTINUE
 4003           CONTINUE

C  COMPUTE PHYSICAL VARIABLES

               DO 0002 I=1,M
                  R(I)  =W(1,I)
                  U(I)  =W(2,I)/R(I)
                  Z(I)  =W(4,I)/R(I)
                  IF (Z(I).GE.0.0) THEN
                      G(I)  =GAM1
                  ELSE
                      G(I)  =GAM0
                  END IF
C ACROSS INTERFACE
              IF (IND(I).EQ.1.OR.IND(I-1).EQ.1) THEN
                  P(I)  =W(5,I)
C AWAY FROM INTERFACE
              ELSE
                      P(I)  =(G(I)-1.0)*(W(3,I)-0.5*R(I)*U(I)*U(I))
              END IF
                  H(I)  =(W(3,I)+P(I))/R(I)
 0002           CONTINUE
 
            RETURN
                
 100    format(1X,I3,3X,5(F7.4,2X))
 
        END

C----------------------------------------------------------------------C
       
            SUBROUTINE CFLTIME(M,R,U,P,G,DX,DTMIN)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (MX=1600)  
            DIMENSION R(-1:MX+2),U(-1:MX+2),P(-1:MX+2),G(-1:MX+2)
            COMMON /GAMMAS/GAM1,GAM0
            SMAX=-1.0E+06
            DXMIN=1.0E+06
            DO 0001 I=1,MX
               A=SQRT(G(I)*P(I)/R(I))
               SMUA=ABS(U(I))+A
               IF(SMUA.GT.SMAX)SMAX=SMUA
 0001        CONTINUE
            DTMIN=DX/SMAX
            RETURN
            END
C----------------------------------------------------------------------C

            SUBROUTINE SUPERBEE
                IMPLICIT DOUBLE PRECISION (A-H,O-Z)
            COMMON /BEEFUN/B,B1,B2,ACN
            X=B1*B2/(B1*B1+0.00000001)
            B=2.0
            IF(X.LT.2.0) B=X
            IF(X.LT.1.0) B=1.0
            IF(X.LT.0.5) B=2.0*X
            IF(X.LT.0.0) B=0.0
            RETURN
            END

C----------------------------------------------------------------------C

            SUBROUTINE MINBEE
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
            COMMON /BEEFUN/B,B1,B2,ACN
            X=B1*B2/(B1*B1+0.00000001)
            B=1.0
            IF(X.LT.1.0) B=X
            IF(X.LT.0.0) B=0.0
            RETURN
            END
C
C----------------------------------------------------------------------C

