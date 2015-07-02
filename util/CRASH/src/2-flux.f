C----------------------------------------------------------------------C
C     1D 2-fluid  EULER Equations with level-sets  U=(R,RU,E,RZ);      C
C                                                                      C
C            wave-based second order upwind scheme.                    C
C                                                                      C
C  two different flux functions across  the cell interface             C
C                     separating the two fluids.                       C
C                                                                      C
C----------------------------------------------------------------------C  

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (MDX=3300)
        DIMENSION D(-1:MDX+2),U(-1:MDX+2),P(-1:MDX+2)
        DIMENSION Z(-1:MDX+2),G(-1:MDX+2),PINF(-1:MDX+2),X(1:MDX)
        DIMENSION W(4,-1:MDX+2),DWL(4,-1:MDX+2),DWR(4,-1:MDX+2)
        dimension sum0(3)

        COMMON /GAMMAS/GL,GR,PINFL,PINFR
        COMMON /GRISIZ/M,MM1,MP1,MP2        

        OPEN(90,FILE="LEVELSET.in",STATUS="OLD")
        READ(90,*)NOTIST,CFLCOE,MX
        CLOSE(90)

        CALL INDATA(MX,D,U,P,Z,G,PINF,X,DX,sum0)

        N=0
        TIME=0.0

        DO 0001 N=1,NOTIST

C       Locate the interface between cell i and cell i+1.

        DO 0009 I=1,MX
          IF (Z(I)*Z(I+1).LT.0.0) INT=I
 0009   CONTINUE

            CALL CFLTIME(MX,D,U,P,Z,G,PINF,DX,DTMIN)
            DT=CFLCOE*DTMIN
            DTDX=DT/DX
           DT=DTDX*DX

           DO 0011 I=1,MX
             W(1,I)=D(I)
             W(2,I)=D(I)*U(I)
             W(3,I)=(P(I)+G(I)*PINF(I))/(G(I)-1.0)+0.5*D(I)*U(I)*U(I) 
             W(4,I)=D(I)*Z(I)
 0011      CONTINUE

           call conserve(mx,d,u,p,g,pinf,w,time,n,dtdx,dx,int,sum0)

c          WRITE(6,*)N,TIME,DTDX,DT

              M  =MX
              MM1=MX-1
              MP1=MX+1
              MP2=MX+2

                     CALL ROE(D,U,P,Z,GL,PINFL,DWL,DTDX)
                     CALL ROE(D,U,P,Z,GR,PINFR,DWR,DTDX)

              DO 0004 I=1,M
                 IF (I.LT.INT+1) THEN
                 DO 0002 IU=1,4
                   W(IU,I)=W(IU,I)+DWL(IU,I)
 0002            CONTINUE
                 ELSE
                 DO 0003 IU=1,4
                   W(IU,I)=W(IU,I)+DWR(IU,I)
 0003            CONTINUE
                 ENDIF
 0004         CONTINUE
 
 
C  COMPUTE PHYSICAL VARIABLES
              DO 0005 I=1,M
                 D(I)=W(1,I)
                 U(I)=W(2,I)/W(1,I)
                 P(I)=(G(I)-1.0)*(W(3,I)-0.5*D(I)*U(I)*U(I))
     1                -G(I)*PINF(I)
                 Z(I)=W(4,I)/D(I)
                 IF (Z(I).LT.0.0) THEN
                    G(I)=GL
                    PINF(I)=PINFL
                 ELSE
                    G(I)=GR
                    PINF(I)=PINFR
                 ENDIF
 0005         CONTINUE
 
           TIME=TIME+DT

        

 0001        CONTINUE

         CALL OUTDATA(MX,X,D,U,P,Z)
                 
           write(40,*) time

c          call conserve(mx,d,u,p,g,pinf,w,time+dt,n,dtdx,dx,int,sum0)

        END
C----------------------------------------------------------------------C

        SUBROUTINE INDATA(MX,D,U,P,Z,G,PINF,X,DX,sum0)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (MDX=3300)
        DIMENSION D(-1:MDX+2),U(-1:MDX+2),P(-1:MDX+2)
        DIMENSION Z(-1:MDX+2),G(-1:MDX+2),PINF(-1:MDX+2),X(1:MDX)
        dimension sum0(3)

        COMMON /GAMMAS/GL,GR,PINFL,PINFR

c        DATA dl,ul,pl,gl,pinfl/1.000,0.0,1.000,1.4,0.0/
c        DATA dr,ur,pr,gr,pinfr/0.125,0.0,0.100,1.2,0.0/

c       DATA dl,ul,pl,gl,pinfl/2.000,1.0,1.000,1.4,0.0/
c       DATA dr,ur,pr,gr,pinfr/1.000,1.0,1.000,1.2,0.0/

c       DATA dl,ul,pl,gl,pinfl/1000.0,0.0,1.0d9,4.4,6.0d8/
c       DATA dr,ur,pr,gr,pinfr/50.00 ,0.0,1.0d5,1.4,0.0/

        DX=1.0d0/float(MX)

        DO I=1,MX
           X(I)=(I-0.5)*DX
        end do

C shock-interface data

c
c       helium (light)
c      
c       dratio=0.287/2.08
c       gr=1.67d0

c
c       freon (heavy)
c
c        dratio=0.287/0.091
c        gr=1.249

c
c       air
c
c       dratio=1.0d0
c       gr= 1.4d0

c        read(41,*) dl,ul,pl
c        read(41,*) s

c strong shock wave in air hitting a freon interface

c        xint=0.5

c        gl= 1.4d0
c        pinfl= 0.0d0

c        dl = 5.66981    
c                ul = 9.02990   
c                pl = 100.0000
                
c        dm= 1.0d0
c        um= 0.0d0
c        pm= 1.0d0

c        dratio=0.287/0.091
c        gr=1.249

c        dr= dratio*dm

c        ur= 0.0d0
c        pr= 1.0d0
c        pinfr= 0.0d0

c            do i=1,mx
c           x(i)=(float(i)-0.5)*dx
c           d(i)=dl
c           u(i)=ul
c           p(i)=pl
c           g(i)=gl
c           pinf(i)=pinfl
c           if (i.gt.mx/4) then
c              d(i)=dm
c              u(i)=um
c              p(i)=pm
c           end if
c           if (i.gt.mx/2) then
c              d(i)=dr
c              g(i)=gr
c              pinf(i)=pinfr
c           end if
c        end do

c stiff eos shock tube problem
        gl=4.4
                pinfl=6.0d8
                
                gr=1.4
                pinfr=0.0
                
                xint=0.7
                
        do i=1,mx
          D(I)=1000.0
          U(I)=0.0
          P(I)=1.0d9
          G(I)=gl
          PINF(I)=pinfl
          IF(x(i).GT.xint) THEN
             D(I)=50.0
             U(I)=0.0
             P(I)=1.0d5
             G(I)=gr
             PINF(I)=pinfr
          END IF
                end do
                 
c Initializing level-set function - interface at x=xint
                
        DO I=1,MX
           Z(I)=X(I)-xint
           if (z(i).gt. 0.5) z(i)= 0.5
           if (z(i).lt.-0.5) z(i)=-0.5
        end do

        RETURN
        END

C----------------------------------------------------------------------C

        SUBROUTINE ROE(RHO,U,P,Z,GAMMA,PINFTY,DW,DTDX)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (MD=3300)
        DIMENSION S(4),EIG(4,4),ALF(4)
        DIMENSION RHO(-1:MD+2),U(-1:MD+2),P(-1:MD+2)
        DIMENSION Z(-1:MD+2),H(-1:MD+2)
        DIMENSION DW(4,-1:MD+2),
     1  CN(4,-1:MD+2),DWNP(4,4,-1:MD+2),NC(4,-1:MD+2),BLIM(4,-1:MD+2)
        COMMON /BEEFUN/B,B1,B2
        COMMON /GRISIZ/M,MM1,MP1,MP2

C        COMPUTE ENTHALPY H
        do 0001 I=1,M
              RHOI=RHO(I)
              UI  =U(I)
              PI  =P(I)
              AI2 = GAMMA*(PI+PINFTY)/RHOI
              H(I)= 0.5*UI*UI+AI2/(GAMMA-1.0)
 0001   continue
C        BOUNDARY CONDITIONS
        RHO(-1)=RHO(2)
        RHO(0) =RHO(1)
        U(-1)  =U(2)
        U(0)   =U(1)
        P(-1)  =P(2)
        P(0)   =P(1)
        H(-1)  =H(2)
        H(0)   =H(1)
        Z(-1)  =Z(2)
        Z(0)   =Z(1)
        RHO(MP1)=RHO(M)
        RHO(MP2)=RHO(MM1)
        U(MP1)  =U(M)
        U(MP2)  =U(MM1)
        P(MP1)  =P(M)
        P(MP2)  =P(MM1)
        H(MP1)  =H(M)
        H(MP2)  =H(MM1)
        Z(MP1)  =Z(M)
        Z(MP2)  =Z(MM1)

        DO 0002 I=-1,MP2
        DO 0002 IU=1,4
             DW(IU,I)=0.0
 0002   CONTINUE

        DO 3003 I=-1,MP1
C           SOLVE RIGHT RIEMANN PROBLEM RP(I,I+1)
           RL=RHO(I)
           UL=U(I)
           PL=P(I)
           ZL=Z(I)
           HL=H(I)
           RR=RHO(I+1)
           UR=U(I+1)
           PR=P(I+1)
           ZR=Z(I+1)
           HR=H(I+1)
C  FIND GRADIENTS AND MEAN VALUES
           DU=UR-UL
           DR=RR-RL
           DP=PR-PL
           DZ=ZR-ZL
           WT=SQRT(RR/RL)
           WTP1=WT+1.0
           RC=WT*RL
           UC=(UL+WT*UR)/WTP1
           HC=(HL+WT*HR)/WTP1
           ZC=(ZL+WT*ZR)/WTP1
           QQ=0.5*UC*UC
           AA=(GAMMA-1.0)*(HC-QQ)
           AC=SQRT(AA)
C  FIND WAVESPEEDS AND WAVESTRENGTHS
           S(1)        =UC-AC
           S(2)        =UC
           S(3)        =UC
           S(4)        =UC+AC
           ALF(1)=(DP-RC*AC*DU)/(2.0*AA)
           ALF(2)=(AA*DR-DP)/AA
           ALF(3)= RC*DZ
           ALF(4)=(DP+RC*AC*DU)/(2.0*AA)
C  FIND COMPONENTS OF EIGENVECTORS
           EIG(1,1)=1.0
           EIG(2,1)=UC-AC
           EIG(3,1)=HC-UC*AC
           EIG(4,1)=ZC

           EIG(1,2)=1.0
           EIG(2,2)=UC
           EIG(3,2)=0.5*UC*UC
           EIG(4,2)=ZC

           EIG(1,3)=0.0
           EIG(2,3)=0.0
           EIG(3,3)=0.0
           EIG(4,3)=1.0

           EIG(1,4)=1.0
           EIG(2,4)=UC+AC
           EIG(3,4)=HC+UC*AC
           EIG(4,4)=ZC

C          COMPUTE FIRST-ORDER FLUX Contributions     

           DO 3002 IW=1,4
                 CN(IW,I)=S(IW)*DTDX
                 NC(IW,I)=SIGN(1.0D0,CN(IW,I))
                 ITARG   =I+(1+NC(IW,I))/2
                 DO 3001 IU=1,4
                    WINC        =CN(IW,I)*EIG(IU,IW)*ALF(IW)
                    DW(IU,ITARG) =DW(IU,ITARG)-WINC
                    DWNP(IU,IW,I)=0.5*(NC(IW,I)-CN(IW,I))*WINC
                    BLIM(IW,I)   =ALF(IW)*CN(IW,I)*(NC(IW,I)-CN(IW,I))
 3001           CONTINUE
 3002      CONTINUE
 3003        CONTINUE

C       HIGHER ORDER EFECTS BETWEEN 2 AND M-1

        DO 4003 I=0,M
        DO 4002 IW=1,4
C  COMPUTE HIGHER-ORDER EFFECTS OF (IW)TH WAVE
          IS=I-NC(IW,I)
          B1=BLIM(IW,I)
          B2=BLIM(IW,IS)
          CALL SUPERBEE
c            B=0.0
             DO 4001 IU=1,4
             DW(IU,I+1)=DW(IU,I+1)+B*DWNP(IU,IW,I)
             DW(IU,I  )=DW(IU,I  )-B*DWNP(IU,IW,I)
 4001     CONTINUE
 4002   CONTINUE
 4003   CONTINUE

        RETURN
        END

C----------------------------------------------------------------------C

        SUBROUTINE CFLTIME(MX,D,U,P,Z,G,PINF,DX,DTMIN)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (MDX=3300)
        DIMENSION D(-1:MDX+2),U(-1:MDX+2),P(-1:MDX+2)
        DIMENSION Z(-1:MDX+2),G(-1:MDX+2),PINF(-1:MDX+2)
        COMMON /GAMMAS/GL,GR,PINFL,PINFR
        SMAX=-1.0E+06
        DO 0002 I=1,MX
           A=SQRT(G(I)*(P(I)+PINF(I))/D(I))
           SMUA=ABS(U(I))+A
           IF(SMUA.GT.SMAX)SMAX=SMUA
 0002        CONTINUE
 0001        CONTINUE
        DTMIN=DX/SMAX
        RETURN
        END
C----------------------------------------------------------------------C


        SUBROUTINE SUPERBEE
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMMON /BEEFUN/B,B1,B2
        X=B1*B2/(B1*B1+0.00000001)
        B=2.0
        IF(X.LT.2.0)B=X
        IF(X.LT.1.0)B=1.0
        IF(X.LT.0.5)B=2.0*X
        IF(X.LT.0.0)B=0.0
        RETURN
        END
C
        SUBROUTINE MINBEE
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMMON /BEEFUN/B,B1,B2
        X=B1*B2/(B1*B1+0.00000001)
        B=1.0
        IF(X.LT.1.0)B=X
        IF(X.LT.0.0)B=0.0
        RETURN
        END

C -------------------------------------------------------------------

        SUBROUTINE CONSERVE(MX,D,U,P,G,PINF,W,TIME,N,DTDX,DX,INT,SUM0)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        
        PARAMETER (MDX=3300)
        DIMENSION D(-1:MDX+2),U(-1:MDX+2),P(-1:MDX+2)
        DIMENSION G(-1:MDX+2),PINF(1:MDX)
        DIMENSION W(4,-1:MDX+2)
        DIMENSION SUM(3),SUM0(3)

           DO 0001 J=1,3
           SUM(J)=0.0d0
           DO 0001 I=1,MX
              SUM(J)=SUM(J)+W(J,I)
 0001      CONTINUE

           fr1=d(mx)*u(mx)
           fl1=d( 1)*u( 1)
           fr2=d(mx)*u(mx)*u(mx)+p(mx)
           fl2=d( 1)*u( 1)*u( 1)+p( 1)
           fr3=u(mx)*(g(mx)*p(mx)/(g(mx)-1.0)+0.5*d(mx)*u(mx)*u(mx))
           fl3=u( 1)*(g( 1)*p( 1)/(g( 1)-1.0)+0.5*d( 1)*u( 1)*u( 1))
           sum(1)=dx*sum(1)+time*(fr1-fl1)
           sum(2)=dx*sum(2)+time*(fr2-fl2)
           sum(3)=dx*sum(3)+time*(fr3-fl3)

           err1=(sum(1)-sum0(1))*100.0d0/sum0(1)
           err2=(sum(2)-sum0(2))*100.0d0/sum0(2)
           err3=(sum(3)-sum0(3))*100.0d0/sum0(3)

           WRITE(30,200)err3
           WRITE(31,200)err1

           err3=sum(3)
           err1=sum(1)

           WRITE(32,*)err3
           WRITE(33,*)err1
    
c          WRITE(6,*)N,time,sum(1),sum(2),sum(3),int
c          WRITE(6,*)N,time,sum0(1),sum0(2),sum0(3),int
           
c          WRITE(6,*)N,time,err1,err2,err3,int

           WRITE(6,*)N,time,err3,int
 100       format(1x,i3,4(2x,f10.8),2x,i3)
 200       format(1x,f12.8)

           RETURN
           END

C -------------------------------------------------------------------

        SUBROUTINE OUTDATA(MX,X,D,U,P,Z)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        PARAMETER (MDX=3300)
        DIMENSION D(-1:MDX+2),U(-1:MDX+2),P(-1:MDX+2)
        DIMENSION Z(-1:MDX+2),X(1:MDX)


        DO 0001 I=1,MX
           WRITE(15,100)D(I)
           WRITE(16,100)U(I)
           WRITE(17,*)P(I)
           WRITE(18,100)Z(I)
           WRITE(14,100)X(I)
 0001   CONTINUE

        RETURN
 100    format(1x,f15.10)
        END

C ------------------------------------------------------------------C

