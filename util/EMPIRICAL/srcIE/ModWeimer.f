
      module EIE_ModWeimer

      contains

c
c     The code has been made to implicit real*8 by Mei-Ching Fok on
c     Jan. 30, 2002

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

*

* Subroutines to calculate the electric potentials from the "Weimer 2K" model of

* the polar cap ionospheric electric potentials described in the publication: 

* Weimer, D. R., An improved model of ionospheric electric potentials including

* substorm perturbations and application to the Geospace Environment Modeling

* November 24, 1996 event, Journal of Geophysical Research, Vol. 106, p. 407, 2001.

*

* To use, first call procedure SETMODEL00 with the specified input parameters:

*   angle: IMF Y-Z clock angle in degrees, 0=northward, 180=southward

*   Bt: Magnitude of IMF in Y-Z plane in nT

*   Tilt: dipole tilt angle in degrees.

*   SWVel: solar wind velocity in km/sec

*   SWDen: solar wind density in #/cc

*   ALindex: (optional) AL index in nT

*

* The function EPOTVAL00(gLAT,gMLT) can then be used repeatively to get the

* electric potential in kV at the desired location.

* Input coordinates assume use of 'altitude adjusted' corrected geomagnetic

* coordinates for R=1, also refered to as AACGM0.

*

* The function BOUNDARYLAT00(gMLT) can be used to get the latitude of the boundary

*   where the potential goes to zero.  This boundary is a function of MLT, and

*   varies with the SETMODEL00 parameters.  The potential is zero everywhere below

*   this boundary.

*

* Two data files are provided:

*	'w2klittle.dat' for LITTLE_ENDIAN machines.

*	'w2kbig.dat'    for    BIG_ENDIAN machines.

* You must copy or rename the correct one to the file 'w2k.dat'

*

* This code is protected by copyright and is distributed

* for research or educational use only.

* Commerical use without written permission from Dan Weimer/MRC is prohibited.

*

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

	SUBROUTINE DLEGENDRE00(x,lmax,mmax,Plm,dPlm)

* compute Double Precision Associate Legendre Function P_l^m(x), as well as optional

* derivatives, for all L up to lmax and all M up to mmax.

* Returns results in array Plm, which should be dimensioned as double precision

* with size (0:10,0:10).

* If the first element of dPlm is not zero, then the first derivatives are also

* computed, and put into array dPlm.  To skip the derivatives it is only necessary

* to put a scalar 0.D0 in the dPlm parameter. Otherwise, dPlm should also be a

* double precision array of size (0:10,0:10).

* The recursion formulas keep a count of the exponents of the factor SQRT(1-x^2)

*  in both the numerator and denominator, which may cancel each other in the

*  final evaluation, particularly with the derivatives.  This prevents infinities

*  at x=-1 or +1. 

* If X is out of range ( abs(x)>1 ) then value is returns as if x=1.
       
        implicit real*8 (a-h,o-z)
	real*8           x,xx,Plm(0:10,0:10),P(0:10,0:10,0:1),fact,sfact

	real*8           dPlm(0:10,0:10),dP(0:10,0:10,0:2),anum,term

	LOGICAL DodPlm



	DodPlm=dPlm(0,0) .NE. 0.d0



	DO l=0,lmax

	    DO m=0,mmax

		  Plm(l,m)=0.D0

		  P(l,m,0)=0.D0

		  P(l,m,1)=0.D0

	    ENDDO

	ENDDO

	IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN

	  Print *,'Bad arguments to DLegendre00'

	  RETURN

	ENDIF



* Copy x to xx, and make sure it is in range of -1. to +1.

	xx=MIN(x,1.D0)

	xx=MAX(xx,-1.D0)



	P(0,0,1)=1.D0

	IF(lmax.GT.0) P(1,0,1)=xx

	IF(lmax.GT.1)THEN

	   DO L=2,lmax

	     P(L,0,1)=( (2.D0*L-1)*xx*P(L-1,0,1) - (L-1)*P(L-2,0,1) ) / L

	   ENDDO

	ENDIF



	fact=1.D0-xx**2

	sfact=DSQRT(fact)



	IF(mmax .GT. 0)THEN

		DO M=1,mmax

		  DO L=M,lmax

			L2=MAX( L-2 ,  0 )

			P(L,M,1)= P(L2,M,1) -(2*L-1)*P(L-1,M-1,0)*fact

			P(L,M,0)= P(L2,M,0) -(2*L-1)*P(L-1,M-1,1)

		  ENDDO

	    ENDDO

	ENDIF



	IF(DodPlm)Then !do derivatives

* First zero arrays

		DO l=0,lmax

			DO m=0,mmax

				dPlm(l,m)=0.D0

				dP(l,m,0)=0.D0

				dP(l,m,1)=0.D0

				dP(l,m,2)=0.D0

			ENDDO

		ENDDO



		IF(lmax .GT. 0) dP(1,0,1)=1.D0



		IF(lmax .GT. 1)THEN

			DO L=2,lmax  

				dP(L,0,1)=( (2*L-1)*P(L-1,0,1) + 

     $                  (2*L-1)*xx*dP(L-1,0,1) - 

     $                  (L-1)*dP(L-2,0,1) ) / L

			ENDDO

		ENDIF



		IF(mmax .GT. 0)THEN

		  DO M=1,mmax  

			DO L=M,lmax  

			  L2=MAX( L-2 ,  0 )

			  dP(L,M,1)= dP(L2,M,1) - (2*L-1)*fact*dP(L-1,M-1,0) - 

     $                  (2*L-1)*dP(L-1,M-1,2) + (2*L-1)*xx*P(L-1,M-1,0)

			  dP(L,M,0)= dP(L2,M,0) - (2*L-1)*dP(L-1,M-1,1)

			  dP(L,M,2)=dP(L2,M,2) +(2*L-1)*xx*P(L-1,M-1,1)

			ENDDO

		  ENDDO

		ENDIF



		DO L=0,lmax  

	      mlimit=MIN(mmax,L)

		  DO M=0,mlimit

* Prevent a divide by zero

		    anum=dP(L,M,2) !numerator

			IF(sfact.NE.0.)Then !denominator is OK

			  term=anum/sfact 

			ELSE !denominator is zero

			  IF(DABS(anum).LT.1.D-7)THEN

				term=0.D0 !return 0 in cases where numerator is near zero

			  ELSE !return nearly infinity with same sign as numerator

				term=DSIGN(1.D36,anum) 

			  ENDIF

			ENDIF

			dPlm(L,M)=dP(L,M,1) + dP(L,M,0)*sfact + term

		  ENDDO

		ENDDO



	ENDIF !End doing derivative



	DO L=0,lmax

	    mlimit=MIN(mmax,L)

	    DO M=0,mlimit

		  Plm(L,M)=P(L,M,1) + P(L,M,0)*sfact

	    ENDDO

	ENDDO	



	RETURN

	END SUBROUTINE DLEGENDRE00

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

*	FUNCTION FSVal00(omega,MaxN,FSC)
*
** Fourier Series Value
*
** Return value of Sine/Cosine Fourier series for N terms up to MaxN
*
** at angle omega, given the F.S. coeficients in array FSC
*
*        implicit real*8 (a-h,o-z)
*	REAL*8 omega,FSC(0:1,0:*)
*
*	INTEGER MaxN,n
*
*	REAL*8 Y,theta
*
*	Y=0.
*
*	DO n=0,MaxN
*
*	  theta=omega*n
*
*	  Y=Y + FSC(0,n)*COS(theta) + FSC(1,n)*SIN(theta)
*
*	ENDDO
*
*	FSVal=Y
*
*	RETURN
*
*	END

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

	SUBROUTINE SetModel00(angle,Bt,Tilt,SWVel,SWDen,ALindex,UseAL)

*

* Calculate the complete set of spherical harmonic coeficients,

* given an aribitrary IMF angle (degrees from northward toward +Y),

* magnitude Bt (nT), dipole tilt angle (degrees), 

* solar wind velocity (km/sec), SWDen (#/cc),

* ALindex (nT), and Logical flag to use optional AL index.

*

* Sets the value of Coef and Boundfit in the common block SetW00Coef.

*
        implicit real*8 (a-h,o-z)
	REAL*8 angle,Bt,Tilt,SWVel,SWDen,ALindex

	LOGICAL First,UseAL



	DATA First/.TRUE./

	SAVE First

	INTEGER unit

	CHARACTER*15 cfile

	CHARACTER*30 Copyright

	PARAMETER (MJ=3,ML=4,MM=3,MN=2,MO=2)

	integer :: iMJ,iMO,iML,iMM,iMN

	REAL*4  CS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)

	REAL*4 BCS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:MN)

	REAL*4  SS( 0:1 , 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)

	REAL*4 BSS( 0:1 , 0:1 , 0:MO, 0:1 , 0:MN)

	REAL*8 XA(0:MJ),XB(0:MJ),FSC(0:1,0:4),PSS(0:1)

	REAL*8 Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi

	real*8           dpi

	INTEGER*4 i,j,k,l,m,n,o

	INTEGER*4 maxj,MaxL,MaxM,MaxN,MaxO

	COMMON /EIE_AllW00Coefs/MaxJ,MaxO,CS,BCS,SS,BSS

	COMMON /EIE_SetW00Coef/Coef,BoundFit,pi,dpi,MaxL,MaxM,MaxN



	If(First)Then

	cfile='EIE/w2k.dat' !make sure correct ENDIAN type file is used.

	unit=9

c	OPEN(UNIT=unit,FILE=cfile,STATUS='OLD',form='UNFORMATTED')
	OPEN(UNIT=unit,FILE=cfile,STATUS='old',form='FORMATTED')


	READ(unit,*) Copyright
	
	
	
c       PRINT *,Copyright

	
	read(unit,*) Maxj,MaxL,MaxM,MaxN,MaxO


	If(maxj.NE.MJ .OR. MaxL.NE.ML .OR. MaxM.NE.MM .OR.

     $   MaxN.NE.MN .OR. MaxO.NE.MO)Then

		PRINT *,'Data File Error'

		STOP !Data file did not match sixe expected for arrays

	Endif   

	
	do iMJ=0,MJ
	   do i=0,1
	      do iMO=0,MO
		 do j=0,1
		    do iML=0,ML
		       do iMM=0,MM
			  read(unit,*) CS( iMJ, i , iMO, j , iML, iMM)
		       enddo
		    enddo
		 enddo
	      enddo
	   enddo
	enddo

	
	do iMJ=0,MJ
	   do i=0,1
	      do iMO=0,MO
		 do j=0,1
		    do iMN=0,MN
		       read(unit,*) BCS( iMJ, i , iMO, j , iMN)
		    enddo
		 enddo
	      enddo
	   enddo
	enddo

	
	do i=0,1
	   do j=0,1
	      do iMO=0,MO
		 do k=0,1
		    do iML=0,ML
		       do iMM=0,MM
			  read(unit,*) SS( i , j , iMO, k , iML, iMM)
		       enddo
		    enddo
		 enddo
	      enddo
	   enddo
	enddo
	

	
	do i=0,1
	   do j=0,1
	      do iMO = 0,MO
		 do k=0,1
		    do iMN=0,MN
		       read(unit,*) BSS( i , j , iMO, k , iMN)
		    enddo
		 enddo
	      enddo
	   enddo
	enddo

	CLOSE(unit)
	

	pi=2.*ASIN(1.)

	dpi=2.D0*DASIN(1.D0)

	First=.FALSE.

	Endif



	SinTilt=SIN(Tilt)

	omega=angle*pi/180.

	XA(0)=1.

	XA(1)=Bt**(2./3.) *SWvel

	XA(2)=SinTilt

	XA(3)=SWvel**2 *SWDen

	XB(0)=1.

	XB(1)=Bt

	XB(2)=SinTilt

	XB(3)=SWvel**2 *SWDen



	DO l=0,MaxL

	    mlimit=MIN(l,MaxM)

	    DO m=0,mlimit

		  klimit=MIN(m,1)

		  DO k=0,klimit

			acoef=0. !rezero for summation

			DO j=0,MaxJ

			    DO o=0,MaxO

				  DO i=0,1

					FSC(i,o)=CS(j,i,o,k,l,m)

				  ENDDO

			    ENDDO

     			    acoef=acoef+ XA(j)*FSVAL(omega,MaxO,fsc)

			ENDDO

			IF(UseAL)THEN

			    DO j=0,1

				  DO o=0,MaxO

					DO i=0,1

					    FSC(i,o)=SS(j,i,o,k,l,m)

					ENDDO

				  ENDDO

				  PSS(j)=FSVAL(omega,MaxO,fsc)

			    ENDDO

			    acoef=acoef + PSS(0) + PSS(1)*ALindex

			ENDIF

			Coef(k,l,m)=acoef

		  ENDDO

	    ENDDO

	ENDDO



	DO n=0,MaxN

	    klimit=MIN(n,1)

	    DO k=0,klimit

		  acoef=0. !rezero for summation

		  DO j=0,MaxJ

			DO o=0,MaxO

			    DO i=0,1

				  FSC(i,o)=BCS(j,i,o,k,n)

			    ENDDO

			ENDDO

     			acoef=acoef+ XB(j)*FSVAL(omega,MaxO,fsc)

		  ENDDO

		  IF(UseAL)THEN

			DO j=0,1

			    DO o=0,MaxO

				  DO i=0,1

					FSC(i,o)=BSS(j,i,o,k,n)

				  ENDDO

			    ENDDO

			    PSS(j)=FSVAL(omega,MaxO,fsc)

			ENDDO

			acoef=acoef + PSS(0) + PSS(1)*ALindex

		  ENDIF

		  BoundFit(k,n)=acoef

	    ENDDO

	ENDDO

	RETURN

	END SUBROUTINE SetModel00

****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************

	FUNCTION BoundaryLat00(gmlt)

        implicit real*8 (a-h,o-z)
	REAL*8 gmlt

	REAL*8 Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi

	real*8           dpi

        integer*4 MaxL,MaxM,MaxN
	COMMON /EIE_SetW00Coef/Coef,BoundFit,pi,dpi,MaxL,MaxM,MaxN

	BoundaryLat00 = FSVal(gmlt*pi/12.,MaxN,BoundFit)

	RETURN

	END FUNCTION BoundaryLat00

****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************

	FUNCTION EpotVal00(gLAT,gMLT)

* Return the value of the electric potential in kV at

* corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).

*

* Must first call SetModel00 to set up the model coeficients for

* the desired values of Bt, IMF clock angle, Dipole tilt angle, SW Vel,

* number density, and (optional) AL index.

*
        implicit real*8 (a-h,o-z)
	REAL*8 gLAT,gMLT

	real*8           Phi,Z,O,x,ct,Phim

	real*8            Plm(0:10,0:10),OPlm(0:10,0:10)
	real*8           dPlm(0:10,0:10)


	REAL*8 Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi

	real*8           dpi

        integer*4 MaxL,MaxM,MaxN
	COMMON /EIE_SetW00Coef/Coef,BoundFit,pi,dpi,MaxL,MaxM,MaxN

	blat = BoundaryLat00(gmlt)

	IF(glat .GT. blat)THEN

        Phi=DBLE(gMLT)*dpi/12.D0

	  colat=90.-glat

	  bcolat=90.-blat

	  x=DBLE(colat)*dpi/DBLE(bcolat)

	  DC=DBLE(Coef(0,0,0))

	  Z=DC

	  O=DC

	  ct=DCOS(x)
	  dPlm(:,:)=0.D0
	  CALL DLegendre00(ct,MaxL,MaxM,Plm,dPlm)

!Also find value at outer boundary at angle Pi, or cos(pi)=-1.

	  CALL DLegendre00(-1.D0,MaxL,MaxM,OPlm,dPlm)

	  DO l=1,MaxL

	    Z=Z +  Plm(l,0)*DBLE(Coef(0,l,0))

	    O=O + OPlm(l,0)*DBLE(Coef(0,l,0))

	    mlimit=MIN(l,MaxM)

	    DO m=1,mlimit

	      phim=phi*m

	      Z=Z + Plm(l,m)*(DBLE(Coef(0,l,m))*DCOS(phim) +

     $			    DBLE(Coef(1,l,m))*DSIN(phim) )

	      O=O +OPlm(l,m)*DBLE(Coef(0,l,m))

	    ENDDO

	  ENDDO

	  EpotVal00 = SNGL(Z-O)

	ELSE

	  EpotVal00 = 0.

	ENDIF

	RETURN

	END FUNCTION EpotVal00

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
*
* Subroutines to calculate the electric potentials from the "Weimer 01" model of
* the polar cap ionospheric electric potentials described in the publication: 
* Weimer, D. R., An improved model of ionospheric electric potentials including
* substorm perturbations and application to the Geospace Environment Modeling
* November 24, 1996 event, Journal of Geophysical Research, Vol. 106, p. 407, 2001.
*
* To use, first call procedure SETMODEL01 with the specified input parameters:
*   angle: IMF Y-Z clock angle in degrees, 0=northward, 180=southward
*   Bt: Magnitude of IMF in Y-Z plane in nT
*   Tilt: dipole tilt angle in degrees.
*   SWVel: solar wind velocity in km/sec
*   SWDen: solar wind density in #/cc
*   ALindex: (optional) AL index in nT
*
* The function EPOTVAL01(gLAT,gMLT) can then be used repeatively to get the
* electric potential in kV at the desired location.
* Input coordinates assume use of 'altitude adjusted' corrected geomagnetic
* coordinates for R=1, also refered to as AACGM0.
*
* The function BOUNDARYLAT01(gMLT) can be used to get the latitude of the boundary
*   where the potential goes to zero.  This boundary is a function of MLT, and
*   varies with the SETMODEL01 parameters.  The potential is zero everywhere below
*   this boundary.
*
* This code is protected by copyright and is distributed
* for research or educational use only.
* Commerical use without written permission from Dan Weimer/MRC is prohibited.

CNCAR      Revisions for use at NCAR:
C            (1) Change behavior at minimum magnetic latitude.  When approaching
C                the model equatorial edge (which varies with MLT) the electric
C                potential returned used to go to zero discontinuously; although
C                intended as a flag, it created artificial gradients in the
C                electric field calculation.  Now the potential returned is
C                constant (that of the minimum latitude) for any latitude at
C                or equatorward of the minimum.
C            (2) Accomodate running simultaneously 1996 and 2001 versions.  To
C                avoid name collisions this required: (i) revising names (e.g.,
C                adding '01') for differing subprograms, and (ii) relocating
C                common routines into another file (weicom.f).
C            (3) Pass the coefficients file name and unit number into READCOEF01
C                rather than using hard coded values.
C            (4) Add wrapper subroutines for non-ANSI trig functions which
C                input angles in degrees.
C            (5) Add electric field routine (GECMP01) to deterine the electric
C                potential gradient.
C            (6) Add wrapper routine (WEIEPOT01) for use with AMIE; this is a
C                substitute for calling SETMODEL01 and EPOTVAL01.
C            (7) Remove blanks in some statements to fit in 72 columns to
C                adhere to ANSI Fortran 77.
CNCAR      NCAR changes are delimited by "CNCAR"
*
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
	SUBROUTINE DLEGENDRE01(x,lmax,mmax,Plm)
* compute Double Precision Associate Legendre Function P_l^m(x)
* for all l up to lmax and all m up to mmax.
* Returns results in array Plm.
* The recursion formulas keep a count of the exponents of the factor SQRT(1-x^2)
*  in both the numerator and denominator, which may cancel each other in the
*  final evaluation, particularly with the derivatives.  This prevents infinities
*  at x=-1 or +1. 
* If X is out of range ( abs(x)>1 ) then value is returns as if x=1.
	DOUBLE PRECISION x,xx,Plm(0:10,0:10),P(0:10,0:10,0:1),fact,sfact
	DOUBLE PRECISION dP(0:10,0:10,0:2),anum,term

	DO l=0,lmax
	    DO m=0,mmax
		  Plm(l,m)=0.D0
		  P(l,m,0)=0.D0
		  P(l,m,1)=0.D0
	    ENDDO
	ENDDO
	IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN
	  Print *,'Bad arguments to DLegendre'
	  RETURN
	ENDIF

* Copy x to xx, and make sure it is in range of -1. to +1.
	xx=MIN(x,1.D0)
	xx=MAX(xx,-1.D0)

	P(0,0,1)=1.D0
	IF(lmax.GT.0) P(1,0,1)=xx
	IF(lmax.GT.1)THEN
	   DO L=2,lmax
	    P(L,0,1)=( (2.D0*L-1)*xx*P(L-1,0,1) - (L-1)*P(L-2,0,1) ) / L
	   ENDDO
	ENDIF

	fact=1.D0-xx**2
	sfact=DSQRT(fact)

	IF(mmax .GT. 0)THEN
		DO M=1,mmax
		  DO L=M,lmax
			L2=MAX( L-2 ,  0 )
			P(L,M,1)= P(L2,M,1) -(2*L-1)*P(L-1,M-1,0)*fact
			P(L,M,0)= P(L2,M,0) -(2*L-1)*P(L-1,M-1,1)
		  ENDDO
	    ENDDO
	ENDIF

	DO L=0,lmax
	    mlimit=MIN(mmax,L)
	    DO M=0,mlimit
		  Plm(L,M)=P(L,M,1) + P(L,M,0)*sfact
	    ENDDO
	ENDDO	

	RETURN
	END SUBROUTINE DLEGENDRE01
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

	SUBROUTINE READCOEF01 (udat)
CNCAR      Feb 01: Make a subroutine (a.k.a. ReadCoef) from FIRST=TRUE block
C          of SETMODEL01 (a.k.a. SetModel) so that one can load the coefficients
C          independently of defining geophysical conditions.  This mimics the
C          1996 model.
C          Sep 01: Change argument list to pass both unit number and file name
C          rather than hard coding
CNCAR

	INTEGER udat
	CHARACTER*30 Copyright
	PARAMETER (MJ=3,ML=4,MM=3,MN=2,MO=2)
	REAL  CS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
	REAL BCS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:MN)
	REAL  SS( 0:1 , 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
	REAL BSS( 0:1 , 0:1 , 0:MO, 0:1 , 0:MN)
	REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi
	DOUBLE PRECISION dpi
	INTEGER*4 maxj,MaxL,MaxM,MaxN,MaxO
	COMMON /EIE_AllW01Coefs/MaxJ,MaxO,CS,BCS,SS,BSS
	COMMON /EIE_SetW01Coef/MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi

CNCAR      Feb 01:  Initialize constants used in GECMP01
C          Sep 01:  Omit unneeded min lat variables because of switch to
C                   constant potential at and below min lat (hardcoded
C                   in EPOTVAL01).
      COMMON /EIE_W01CECMP/ ALAMX,STPD,STP2,CSTP,SSTP
C            ALAMX = Absolute max latitude (deg) for normal gradient calc.
C            STPD  = Angular dist (deg) of step @ 300km above earth (r=6371km)
C            STP2  = Denominator in gradient calc
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      STEP = 10.
      STPR = STEP/6671.
      STPD = STPR*R2D
      STP2 = 2.*STEP
      CSTP = COS (STPR)
      SSTP = SQRT (1. - CSTP*CSTP)
      ALAMX = 90. - STPD
CNCAR

CNCAR      Sep 01:  udat is now an argument
C  Assume file has been opened at unit elsewhere
CNCAR
c     OPEN (UNIT=udat,FILE=cfile,STATUS='OLD',form='UNFORMATTED')
      READ (udat) Copyright
      PRINT *,Copyright
      READ (udat) Maxj,MaxL,MaxM,MaxN,MaxO
      If (maxj.NE.MJ .OR. MaxL.NE.ML .OR. MaxM.NE.MM .OR.
     $    MaxN.NE.MN .OR. MaxO.NE.MO) Then
	      PRINT *,'Data File Error'
	      STOP !Data file did not match size expected for arrays
      Endif
      READ (udat) CS
      READ (udat) BCS
      READ (udat) SS
      READ (udat) BSS
      CLOSE (udat)
      pi=2.*ASIN(1.)
      dpi=2.D0*DASIN(1.D0)

      RETURN
      END SUBROUTINE READCOEF01
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
CNCAR      Sep 01: Change name from SetModel to be distinct

	SUBROUTINE SetModel01 (angle,Bt,Tilt,SWVel,SWDen,ALindex,UseAL)
C       SUBROUTINE SetModel   (angle,Bt,Tilt,SWVel,SWDen,ALindex,UseAL)
CNCAR
*
* Calculate the complete set of spherical harmonic coeficients,
* given an aribitrary IMF angle (degrees from northward toward +Y),
* magnitude Bt (nT), dipole tilt angle (degrees), 
* solar wind velocity (km/sec), SWDen (#/cc),
* ALindex (nT), and Logical flag to use optional AL index.
*
* Sets the value of Coef and Boundfit in the common block SetW01Coef.
*
	REAL angle,Bt,Tilt,SWVel,SWDen,ALindex
C       LOGICAL First,UseAL
	LOGICAL UseAL

CNCAR      Feb 01: Move logic block to subroutine RdCoef01
C       DATA First/.TRUE./
C       SAVE First
C       INTEGER unit
C       CHARACTER*15 cfile
C       CHARACTER*30 Copyright
	PARAMETER (MJ=3,ML=4,MM=3,MN=2,MO=2)
	REAL  CS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
	REAL BCS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:MN)
	REAL  SS( 0:1 , 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
	REAL BSS( 0:1 , 0:1 , 0:MO, 0:1 , 0:MN)
	REAL XA(0:MJ),XB(0:MJ),FSC(0:1,0:4),PSS(0:1)
	REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi
	DOUBLE PRECISION dpi
	INTEGER*4 i,j,k,l,m,n,o
	INTEGER*4 maxj,MaxL,MaxM,MaxN,MaxO
	COMMON /EIE_AllW01Coefs/MaxJ,MaxO,CS,BCS,SS,BSS
	COMMON /EIE_SetW01Coef/MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi

CNCAR      Feb 01: Move logic block to subroutine RdCoef01
C  All First=TRUE in new subroutine ReadCoef
C       If(First)Then
C       cfile='w01.dat' !make sure correct ENDIAN type file is used.
C       unit=99
C       OPEN(UNIT=unit,FILE=cfile,STATUS='OLD',form='UNFORMATTED')
C       READ(unit) Copyright
C       PRINT *,Copyright
C       READ(unit) Maxj,MaxL,MaxM,MaxN,MaxO
C       If(maxj.NE.MJ .OR. MaxL.NE.ML .OR. MaxM.NE.MM .OR.
C    $   MaxN.NE.MN .OR. MaxO.NE.MO)Then
C               PRINT *,'Data File Error'
C               STOP !Data file did not match sixe expected for arrays
C       Endif
C       READ(unit) CS
C       READ(unit) BCS
C       READ(unit) SS
C       READ(unit) BSS
C       CLOSE(unit)
C       pi=2.*ASIN(1.)
C       dpi=2.D0*DASIN(1.D0)
C       First=.FALSE.
C       Endif
CNCAR

	SinTilt=SIND(Tilt)
	omega=angle*pi/180.
	XA(0)=1.
	XA(1)=Bt**(2./3.) *SWvel
	XA(2)=SinTilt
	XA(3)=SWvel**2 *SWDen
	XB(0)=1.
	XB(1)=Bt
	XB(2)=SinTilt
	XB(3)=SWvel**2 *SWDen

	DO l=0,MaxL
	    mlimit=MIN(l,MaxM)
	    DO m=0,mlimit
		  klimit=MIN(m,1)
		  DO k=0,klimit
			acoef=0. !rezero for summation
			DO j=0,MaxJ
			    DO o=0,MaxO
				  DO i=0,1
					FSC(i,o)=CS(j,i,o,k,l,m)
				  ENDDO
			    ENDDO
     			    acoef=acoef+ XA(j)*FSVAL(omega,MaxO,fsc)
			ENDDO
			IF(UseAL)THEN
			    DO j=0,1
				  DO o=0,MaxO
					DO i=0,1
					    FSC(i,o)=SS(j,i,o,k,l,m)
					ENDDO
				  ENDDO
				  PSS(j)=FSVAL(omega,MaxO,fsc)
			    ENDDO
			    acoef=acoef + PSS(0) + PSS(1)*ALindex
			ENDIF
			Coef(k,l,m)=acoef
		  ENDDO
	    ENDDO
	ENDDO

	DO n=0,MaxN
	    klimit=MIN(n,1)
	    DO k=0,klimit
		  acoef=0. !rezero for summation
		  DO j=0,MaxJ
			DO o=0,MaxO
			    DO i=0,1
				  FSC(i,o)=BCS(j,i,o,k,n)
			    ENDDO
			ENDDO
     			acoef=acoef+ XB(j)*FSVAL(omega,MaxO,fsc)
		  ENDDO
		  IF(UseAL)THEN
			DO j=0,1
			    DO o=0,MaxO
				  DO i=0,1
					FSC(i,o)=BSS(j,i,o,k,n)
				  ENDDO
			    ENDDO
			    PSS(j)=FSVAL(omega,MaxO,fsc)
			ENDDO
			acoef=acoef + PSS(0) + PSS(1)*ALindex
		  ENDIF
		  BoundFit(k,n)=acoef
	    ENDDO
	ENDDO
	RETURN
	END SUBROUTINE SetModel01
****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************
	FUNCTION BoundaryLat01(gmlt)
	REAL gmlt
	REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi
	DOUBLE PRECISION dpi
	INTEGER MaxL,MaxM,MaxN
	COMMON /EIE_SetW01Coef/MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi
	BoundaryLat01=FSVal(gmlt*pi/12.,MaxN,BoundFit)
	RETURN
	END FUNCTION BoundaryLat01
****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************
CNCAR      Sep 01: Change name from EpotVal to be distinct

	FUNCTION EPOTVAL01 (gLAT,gMLT)
C       FUNCTION EpotVal   (gLAT,gMLT)
CNCAR

* Return the value of the electric potential in kV at
* corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).
*
* Must first call ReadCoef and SetModel to set up the model coeficients for
* the desired values of Bt, IMF clock angle, Dipole tilt angle, and SW Vel.
*
      REAL gLAT,gMLT
      DOUBLE PRECISION Phi,Z,O,x,ct,Phim
      DOUBLE PRECISION Plm(0:10,0:10), OPlm(0:10,0:10)

      REAL Coef(0:1,0:5,0:5), BoundFit(0:1,0:5),pi
      DOUBLE PRECISION dpi
      COMMON /EIE_SetW01Coef/ MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi

      blat = BoundaryLat01(gMLT)
CNCAR      Feb 01:  For latitudes at and equatorward of the model minimum
C          (a function of MLT), limit the latitude used to never be less
C          than the model minimum, thus returning a constant potential for
C          points at and equatorward of BoundaryLat01(gmlt).
C     IF (glat .GT. blat) THEN
      glatlim = amax1 (glat,blat)
CNCAR

      Phi = DBLE (gMLT)*dpi/12.D0

CNCAR
      colat = 90.-glatlim
C     colat = 90.-glat
CNCAR

      bcolat  = 90.-blat
      x  = DBLE (colat)*dpi/DBLE(bcolat)
      DC = DBLE (Coef(0,0,0))
      Z  = DC
      O  = DC
      ct = DCOS(x)
      CALL DLegendre01(ct,   MaxL,MaxM,Plm)
C          Also find value at outer boundary at angle Pi, or cos(pi)=-1.
      CALL DLegendre01(-1.D0,MaxL,MaxM,OPlm)
      DO l=1,MaxL
	Z = Z +  Plm(l,0)*DBLE(Coef(0,l,0))
	O = O + OPlm(l,0)*DBLE(Coef(0,l,0))
	mlimit = MIN(l,MaxM)
	DO m=1,mlimit
	  phim = phi*m
	  Z = Z + Plm(l,m)*(DBLE(Coef(0,l,m))*DCOS(Phim) +
     $                      DBLE(Coef(1,l,m))*DSIN(Phim) )
	  O = O +OPlm(l,m)*DBLE(Coef(0,l,m))
	ENDDO
      ENDDO
CNCAR
      EPOTVAL01 = SNGL(Z-O)
C     EpotVal   = SNGL(Z-O)
C     ELSE
C       EpotVal=0.
C     ENDIF
CNCAR
      RETURN
      END FUNCTION EPOTVAL01

CNCAR
      SUBROUTINE GECMP01 (AMLA,RMLT,ET,EP)
C          Get Electric field components for the 2001 Weimer electrostatic
C          potential model.  Before use, first load coefficients (CALL
C          READCOEF01) and initialize model conditions (CALL SETMODEL01).
C          The electric field, the gradient of the electic potential,
C          is determined using finite differences over a distance of
C          STEP km, defined in READCOEF01 and accessed here via common
C          block CECMP.
C
C          It has been useful to modify STEP.  One application discovered
C          the Weimer electrostatic potential can have minor unrealistic
C          perturbations over short distances.  To avoid these ripples
C          STEP was increased to 5 degrees arc (582 km at 300 km altitude,
C          R=6671 km).
C
C          INPUTS:
C            AMLA = Absolute value of magnetic latitude (deg)
C            RMLT = Magnetic local time (hours).
C          RETURNS:
C            ET = Etheta (magnetic equatorward*) E field component (V/m)
C            EP = Ephi   (magnetic eastward)     E field component (V/m)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when in northern magnetic hemisphere
C              northward when in southern magnetic hemisphere
C          Since the Weimer model uses only NH (positive) latitudes, the
C          sign of ET for SH should be changed outside of this routine.
C
C          HISTORY:
C          Jan 97:  Initial implementation at NCAR for the 1996 version.
C          Feb 01:  Error corrected in determining DPHI.  Old version
C          used wrong spherical right triangle formula (inadvertently
C          computing along a latitude line); more importantly, it
C          neglected to convert mlt from degrees to hour angle input
C          to EPOTVAL01.
C          Sep 01:  Revised equatorward boundary logic to the scheme Barbara
C          uses for AMIE:  As absolute magnetic latitude decreases to the
C          model minimum, the electric potential from EPOTVAL01 now goes to
C          a non-zero constant (rather than zero), thus obviating the need
C          for special handling here.  Therefore, special logic for determining
C          gradients at lower limit has been removed.

      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

C          CECMP contains constants initialized in READCOEF01
      COMMON /EIE_W01CECMP/ ALAMX,STPD,STP2,CSTP,SSTP

      ET = -99999.
      EP = -99999.
      IF (AMLA .LT. 0.) GO TO 100

C          Calculate -(latitude gradient) by stepping along the magnetic
C          latitude line in each direction (flipping coordinates when
C          going over pole to keep lat <= 90).
      KPOL  = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      AMLA1 = AMLA + STPD
      IF (AMLA1 .GT. 90.) THEN
	AMLA1 = 180. - AMLA1
	XMLT1 = XMLT1 + 12.
      ENDIF
      P1 = EPOTVAL01 (AMLA1    ,XMLT1)
      P2 = EPOTVAL01 (AMLA-STPD,XMLT )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

C          Calculate -(lon gradient).  Step along the magnetic meridion
C          in both directions to obtain the electric potential
      IF (AMLA .LT. ALAMX) THEN
	AMLA1 = ASIN(SIN(AMLA*D2R)*CSTP)
	DPHI  = ASIN (SSTP/COS(AMLA1))*R2D/15.      ! 15 converts from degrees to hours
	AMLA1 = AMLA1*R2D
	P1 = EPOTVAL01 (AMLA1,XMLT+DPHI)
	P2 = EPOTVAL01 (AMLA1,XMLT-DPHI)
      ELSE
	AMLA = 90.
	XMLT = XMLT + 6.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

  100 RETURN
      END    SUBROUTINE GECMP01

      SUBROUTINE WEIEPOT01 (IYR,IMO,IDA,IHR,IMN,SWS,SWD,BY,BZ,AL,USEAL,
     +                     IHEM,ISETM,RMLA,RMLT, ET,EP,EPOT)
C          Interface to Weimer-01 (a.k.a w01) for AMIE's combined electrostatic
C          potential models.  This replaces two calls (SETMODEL01 and
C          EPOTVAL01) and automates sign changes for southern hemisphere.
C          INPUTS:
C            IYR   = UT year ((4-digit)
C            IMO   = month of IYR
C            IDA   = day of IMO
C            IHR   = hour of day
C            IMN   = min of hour
C            SWS   = Solar wind speed (km/s)
C            SWD   = Solar wind density (#/cm3)
C            BY    = IMF By component in GSM coordinates (nt)
C            BZ    = IMF Bz component in GSM coordinates (nt)
C            AL    = AL index (nT)
C            USEAL = logical (T or F) to use AL index or not
C            IHEM  = Hemisphere flag: (-1) southern, (1) northern
C            ISETM = Model conditions change flag: (0) no-change
C                                                  (1) time, IMF or SW changed
C            RMLA  = Magnetic latitude  (deg) of point to determine Potential
C                    (RMLA should be NH positive latitudes only, since SH
C                    values found by changing sign of By (or ANGL) and tilt.)
C            RMLT  = Magnetic longitude (hrs) of point to determine Potential
C          RETURNS:
C            ET   = Etheta (magnetic equatorward*) E field component (V/m)
C            EP   = Ephi   (magnetic eastward)     E field component (V/m)
C            EPOT = Electrostatic potential (kV)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when in northern magnetic hemisphere
C              northward when in southern magnetic hemisphere
C
C          Since the AMIE model assumes a NH solution although the
C          latitudes are negative for SH data, the sign of ET
C          for the SH should be changed outside of this routine.

C          HISTORY:
C          Jan 97:  Initial implementation. B. Emery.
C          Feb 01:  Remove bug (HR was not defined)
C          Sep 01:  Add common block s.t. computed values are available to
C          the calling routine without changing the argument list
      COMMON /EIE_W01AT/ ANGL, TILT
      PARAMETER (R2D=57.2957795130823208767981548147)
      LOGICAL USEAL

      H = REAL (IHEM)

      IF (ISETM .EQ. 1) THEN
	ANGL = ATAN2 (BY,BZ)*R2D
	BT   = SQRT (BY*BY + BZ*BZ)
	HR = FLOAT(IHR) + FLOAT(IMN)/60.
	TILT = GET_TILT (IYR,IMO,IDA,HR)

	HANGL = H * ANGL
	HTILT = H * TILT
	CALL SETMODEL01 (HANGL,BT,HTILT,SWS,SWD,AL,USEAL)
C     WRITE (6,'('' WEIMER-2001: '',I4,5I3,8F7.2,L)') IYR,IMO,IDA,IHR,
C    +                    IMN,IHEM,SWS,SWD,AL,BY,BZ,BT,HANGL,HTILT,USEAL
      ENDIF

C          NH assumed latitudes only
      AMLA  = ABS (RMLA)
      EPOT = EPOTVAL01 (AMLA,RMLT)
      CALL GECMP01 (AMLA,RMLT,ET,EP)

      RETURN
      END SUBROUTINE WEIEPOT01
CNCAR
************************ Copyright 1996, Dan Weimer/MRC ***********************
*
* Subroutines to calculate the electric potentials from the Weimer '96 model of
* the polar cap ionospheric electric potentials.
*
* To use, first call subroutine ReadCoef once.
* Next, call SetModel with the specified input parameters.
* The function EpotVal96(gLAT,gMLT) can then be used repeatively to get the
* electric potential at the desired location in geomagnetic coordinates.
*
* This code is protected by copyright and is
* distributed for research or educational use only.
* Commerical use without written permission from Dan Weimer/MRC is prohibited.
*

CNCAR      Revisions for use at NCAR:
C            (1) Change behavior at minimum magnetic latitude.  When approaching
C                the model minimum absolute latitude (45 degrees) the electric
C                potential returned used to go to zero discontinuously; although
C                intended as a flag, it created artificial gradients in the
C                electric field calculation.  Now the potential returned is
C                constant (that of the minimum latitude) for any latitude at
C                or equatorward of the minimum.
C            (2) Accommodate running sumultaneously 1996 and 2001 versions.  To
C                avoid name collisions this required: (i) revising names (e.g.,
C                adding '96') for differing subprograms, and (ii) relocating
C                common routines into another file (weicom.f).
C            (3) Pass the coefficients file name and unit number into RDCOEF96
C                rather than using hard coded values.
C            (4) Add wrapper subroutines for non-ANSI trig functions which
C                input angles in degrees.
C            (5) Add electric field routine (GECMP96) to deterine the electric
C                potential gradient.
C            (6) Add wrapper routine (WEIEPOT96) for use with AMIE; this is a
C                substitute for calling SETMODEL96 and EPOTVAL96.
CNCAR      NCAR changes are delimited by "CNCAR"

************************ Copyright 1996, Dan Weimer/MRC ***********************

CNCAR      Oct 01:  Distinguish from the 2001 version

	FUNCTION EPOTVAL96(gLAT,gMLT)
CNCAR
* Return the value of the electric potential in kV at
* corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).
*
* Must first call ReadCoef and SetModel to set up the model coeficients for
* the desired values of Bt, IMF clock angle, Dipole tilt angle, and SW Vel.
*
	REAL gLAT,gMLT
	Real Theta,Phi,Z,ct,Phim
	REAL Plm(0:20,0:20)

	REAL Coef(0:1,0:8,0:3),pi
	INTEGER ML,MM
	COMMON /EIE_SetW96Coef/ML,MM,Coef,pi

	r = 90.-gLAT

CNCAR      Sep 01: Limit the colatitude (r) such that an input absolute
C          magnetic latitude of 45 degrees or less always calculates the
C          potential at 45.
	R = AMIN1 (R,45.)
C       IF(r .LT. 45.)THEN
CNCAR
	  Theta=r*pi/45.
          Phi=gMLT*pi/12.
	  Z=Coef(0,0,0)
	  ct=COS(Theta)
	  CALL Legendre96(ct,ML,MM,Plm)
	  DO l=1,ML
	    Z=Z + Coef(0,l,0)*Plm(l,0)
	    IF(l.LT.MM)THEN
	      limit=l
	    ELSE
	      limit=MM
	    ENDIF
	    DO m=1,limit
	      phim=phi*m
	      Z=Z + Coef(0,l,m)*Plm(l,m)*COS(phim) +
     $		   Coef(1,l,m)*Plm(l,m)*SIN(phim)
	    ENDDO
	  ENDDO
CNCAR
C       ELSE
C         Z=0.
C       ENDIF
	EPOTVAL96 = Z
CNCAR
	RETURN
	END FUNCTION EPOTVAL96

************************ Copyright 1996, Dan Weimer/MRC ***********************

CNCAR      Sep 01:  Pass in the unit no. and file name rather than hard coding
C          inside.  Also change name to distinguish 2001 version.

	SUBROUTINE READCOEF96 (udat)
CNCAR

*
* Read in the data file with the model coefficients
*
	INTEGER udat
CNCAR
	CHARACTER    skip*15
C       CHARACTER*15 cfile,skip
CNCAR
	REAL C(0:3)
	REAL Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
	INTEGER MaxL,MaxM,MaxN
	COMMON /EIE_AllW96Coefs/MaxL,MaxM,MaxN,Cn


CNCAR      Jan 97:  Initialize constants used in GECMP96
C          Sep 01:  Omit unneeded min lat variables because of switch to
C                   constant potential at and below min lat (hardcoded
C                   in EPOTVAL96).
      COMMON /EIE_W96CECMP/ ALAMX,STPD,STP2,CSTP,SSTP
C            ALAMX = Absolute max latitude (deg) for normal gradient calc.
C            STPD  = Angular dist (deg) of step @ 300km above earth (r=6371km)
C            STP2  = Denominator in gradient calc
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      STEP = 10.
      STPR = STEP/6671.
      STPD = STPR*R2D
      STP2 = 2.*STEP
      CSTP = COS (STPR)
      SSTP = SQRT (1. - CSTP*CSTP)
      ALAMX = 90. - STPD
CNCAR

CNCAR      Sep 01:  udat,cfile are now a formal arguments
C       cfile='wei96.cofcnts'
C       udat=99
CNCAR
C	OPEN(udat,FILE=cfile,STATUS='OLD')
  900   FORMAT(A15)
 1000	FORMAT(3I8)
 2000	FORMAT(3I2)
 3000	FORMAT(2I2,4E15.6)

	READ(udat,900) skip
	READ(udat,1000) MaxL,MaxM,MaxN
	DO l=0,MaxL
	  IF(l.LT.MaxM)THEN
	    mlimit=l
	  ELSE
	    mlimit=MaxM
	  ENDIF
	  DO m=0,mlimit
	    IF(m.LT.1)THEN
	      klimit=0
	    ELSE
	      klimit=1
	    ENDIF
	    DO k=0,klimit
	      READ(udat,2000) ll,mm,kk
	      IF(ll.NE.l .OR. mm.NE.m .OR. kk.NE.k)THEN
		PRINT *,'Data File Format Error'
		STOP
	      ENDIF
	      DO n=0,MaxN
	        IF(n.LT.1)THEN
	          ilimit=0
	        ELSE
	          ilimit=1
	        ENDIF
		DO i=0,ilimit
		  READ(udat,3000) nn,ii,C
	          IF(nn.NE.n .OR. ii.NE.i)THEN
		    PRINT *,'Data File Format Error'
		    STOP
		  ENDIF
		  Cn(0,i,n,k,l,m)=C(0)
		  Cn(1,i,n,k,l,m)=C(1)
		  Cn(2,i,n,k,l,m)=C(2)
		  Cn(3,i,n,k,l,m)=C(3)
		ENDDO
	      ENDDO
	    ENDDO
	  ENDDO
	ENDDO

	CLOSE(udat)
	RETURN
	END SUBROUTINE READCOEF96
************************ Copyright 1996, Dan Weimer/MRC ***********************
CNCAR      Oct 01:  Change name from SetModel to be distinct

	SUBROUTINE SETMODEL96 (angle,Bt,Tilt,SWVel)
CNCAR
*
* Calculate the complete set of spherical harmonic coefficients,
* given an arbitrary IMF angle (degrees from northward toward +Y),
* magnitude Bt (nT), dipole tilt angle (degrees),
* and solar wind velocity (km/sec).
* Returns the Coef in the common block SetCoef.
*
	REAL angle,Bt,Tilt,SWVel
	REAL FSC(0:1,0:4)
	REAL Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
	INTEGER MaxL,MaxM,MaxN
	COMMON /EIE_AllW96Coefs/MaxL,MaxM,MaxN,Cn

	REAL Coef(0:1,0:8,0:3),pi
	INTEGER ML,MM
	COMMON /EIE_SetW96Coef/ML,MM,Coef,pi

	pi=2.*ASIN(1.)
	ML=MaxL
	MM=MaxM
	SinTilt=SIND(Tilt)

	omega=angle*pi/180.
	DO l=0,MaxL
	  IF(l.LT.MaxM)THEN
	    mlimit=l
	  ELSE
	    mlimit=MaxM
	  ENDIF
	  DO m=0,mlimit
	    IF(m.LT.1)THEN
	      klimit=0
	    ELSE
	      klimit=1
	    ENDIF
	    DO k=0,klimit
* Retrieve the regression coefficients and evaluate the function
* as a function of Bt,Tilt,and SWVel to get each Fourier coefficient.
	      DO n=0,MaxN
	        IF(n.LT.1)THEN
	          ilimit=0
* in n=0 case the FSC(1,0) element was never set but used in FSVAL
                  FSC(1,0)=0.0
	        ELSE
	          ilimit=1
	        ENDIF
		DO i=0,ilimit
		  FSC(i,n)=Cn(0,i,n,k,l,m) + Bt*Cn(1,i,n,k,l,m) +
     $		   SinTilt*Cn(2,i,n,k,l,m) + SWVel*Cn(3,i,n,k,l,m)
		ENDDO
	      ENDDO
* Next evaluate the Fourier series as a function of angle.
      	      Coef(k,l,m)=FSVal(omega,MaxN,FSC)
	    ENDDO
	  ENDDO
	ENDDO
	RETURN
	END SUBROUTINE SETMODEL96
************************ Copyright 1996, Dan Weimer/MRC ***********************
	SUBROUTINE LEGENDRE96(x,lmax,mmax,Plm)
* compute Associate Legendre Function P_l^m(x)
* for all l up to lmax and all m up to mmax.
* returns results in array Plm
* if X is out of range ( abs(x)>1 ) then value is returned as if x=1.
	DIMENSION Plm(0:20,0:20)
	  DO l=0,20
	    DO m=0,20
		Plm(l,m)=0.
	    ENDDO
	  ENDDO
	xx=MIN(x,1.)
	xx=MAX(xx,-1.)
	IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN
	  Print *,'Bad arguments to Legendre96'
	  RETURN
	ENDIF
* First calculate all Pl0 for l=0 to l
	Plm(0,0)=1.
	IF(lmax.GT.0)Plm(1,0)=xx
	IF (lmax .GT. 1 )THEN
	  DO L=2,lmax
	    Plm(L,0)=( (2.*L-1)*xx*Plm(L-1,0) - (L-1)*Plm(L-2,0) )/L
	  ENDDO
	ENDIF
	IF (mmax .EQ. 0 )RETURN
	fact=SQRT( (1.-xx)*(1.+xx) )
	DO M=1,mmax
	  DO L=m,lmax
	    lm2=MAX(L-2,0)
	    Plm(L,M)=Plm(lm2,M) - ( 2*L-1)*fact*Plm(L-1,M-1)
	  ENDDO
	ENDDO
	RETURN
	END SUBROUTINE LEGENDRE96

CNCAR
      SUBROUTINE GECMP96 (AMLA,RMLT,ET,EP)
C          Get Electric field components for the 1996 Weimer electrostatic
C          potential model.  Before use, first load coefficients (CALL
C          READCOEF96) and initialize model conditions (CALL SETMODEL96).
C          The electric field, the gradient of the electic potential,
C          is determined using finite differences.
C
C          It has been useful to modify STEP.  One application discovered
C          the Weimer electrostatic potential can have minor unrealistic
C          perturbations over short distances.  To avoid these ripples
C          STEP was increased to 5 degrees arc (582 km at 300 km altitude,
C          R=6671 km).
C
C          INPUTS:
C            AMLA = Absolute value of magnetic latitude (deg)
C            RMLT = Magnetic local time (hours).
C          RETURNS:
C            ET = Etheta (magnetic equatorward*) E field component (V/m)
C            EP = Ephi   (magnetic eastward)     E field component (V/m)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when in northern magnetic hemisphere
C              northward when in southern magnetic hemisphere
C
C          HISTORY:
C          Jan 97:  Initial implementation at NCAR.  R.Barnes
C          Feb 01:  Error corrected in determining DPHI.  Old version
C          used wrong spherical right triangle formula (inadvertently
C          computing along a latitude line); more importantly, it
C          neglected to convert mlt from degrees to hour angle input
C          to EPOTVAL96.
C          Sep 01:  Revised equatorward boundary logic to the scheme Barbara
C          uses for AMIE:  As absolute magnetic latitude decreases to the
C          model minimum, the electric potential from EPOTVAL96 now goes to
C          a non-zero constant (rather than zero), thus obviating the need
C          for special handling here.  Therefore, special logic for determining
C          gradients at lower limit has been removed.

      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

C          CECMP contains constants initialized in READCOEF
      COMMON /EIE_W96CECMP/ ALAMX,STPD,STP2,CSTP,SSTP

      ET = -99999.
      EP = -99999.
      IF (AMLA .LT. 0.) GO TO 100

C          Calculate -(latitude gradient) by stepping along the magnetic
C          latitude line in each direction (flipping coordinates when
C          going over pole to keep lat <= 90).
      KPOL  = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      AMLA1 = AMLA + STPD
      IF (AMLA1 .GT. 90.) THEN
	AMLA1 = 180. - AMLA1
	XMLT1 = XMLT1 + 12.
      ENDIF
      P1 = EPOTVAL96(AMLA1    ,XMLT1)
      P2 = EPOTVAL96(AMLA-STPD,XMLT )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

C          Calculate -(lon gradient).  Step along the magnetic meridion
C          in both directions to obtain the electric potential
      IF (AMLA .LT. ALAMX) THEN
	AMLA1 = ASIN(SIN(AMLA*D2R)*CSTP)

C       DPHI  = ASIN (SSTP/SIN(AMLA1))*R2D          ! old and wrong! (changed Feb 01)
	DPHI  = ASIN (SSTP/COS(AMLA1))*R2D/15.      ! 15 converts from degrees to hours

	AMLA1 = AMLA1*R2D
	P1 = EPOTVAL96(AMLA1,XMLT+DPHI)
	P2 = EPOTVAL96(AMLA1,XMLT-DPHI)
      ELSE
C          At the pole.  Avoid a divide by zero by using Art's trick
C          where Ephi(90,lon) = Etheta(90,lon+90).
	AMLA = 90.
	XMLT = XMLT + 6.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

  100 RETURN
      END    SUBROUTINE GECMP96

      SUBROUTINE WEIEPOT96 (IYR,IMO,IDA,IHR,IMN,SWS,BY,BZ,IHEM,ISETM,
     +                     RMLA,RMLT, ET,EP,EPOT)
C          Interface to Weimer-96 for AMIE's combined electrostatic
C          potential models.  This replaces two calls (SETMODEL96 and
C          EPOTVAL96) and automates sign changes for southern hemisphere.
C          INPUTS:
C            IYR   = UT year ((4-digit)
C            IMO   = month of IYR
C            IDA   = day of IMO
C            IHR   = hour of day
C            IMN   = min of hour
C            SWS   = Solar wind speed (km/s)
C            BY    = IMF By component in GSM coordinates (nt)
C            BZ    = IMF Bz component in GSM coordinates (nt)
C            IHEM  = Hemisphere flag: (-1) southern, (1) northern
C            ISETM = Model conditions change flag: (0) no-change
C                                                  (1) time, IMF or SW changed
C            RMLA  = Magnetic latitude (deg) of point to determine Potential.
C                    RMLA should be positive (NH) only, since SH values are
C                    obtained by changing sign of By (or ANGL) and TILT.
C            RMLT  = Magnetic longitude (hrs) of point to determine Potential
C          RETURNS:
C            ET    = Etheta (magnetic equatorward*) E field component (V/m)
C            EP    = Ephi   (magnetic eastward)     E field component (V/m)
C            EPOT  = Electric potential (kV)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when in northern magnetic hemisphere
C              northward when in southern magnetic hemisphere
C
C          Since the AMIE model assumes a NH solution although the
C          latitudes are negative for SH data, the sign of ET
C          for the SH should be changed outside of this routine.

C          HISTORY:
C          Jan 97:  Initial implementation. B. Emery.
C          Feb 01:  Remove bug (HR was not defined)
C          Sep 01:  Add common block s.t. computed values are available to
C          the calling routine without changing the argument list
C          Oct 01:  Revise routine names s.t. this can run with the 2001 version
      COMMON /EIE_W96AT/ ANGL, TILT
      PARAMETER (R2D=57.2957795130823208767981548147)

      H = REAL (IHEM)

      IF (ISETM .EQ. 1) THEN
	ANGL = ATAN2 (BY,BZ)*R2D
	BT   = SQRT (BY*BY + BZ*BZ)
	HR   = FLOAT(IHR) + FLOAT(IMN)/60.       ! define HR (bug fix Feb 01)
	TILT = GET_TILT (IYR,IMO,IDA,HR)

	HANGL = H * ANGL
	HTILT = H * TILT
	CALL SETMODEL96 (HANGL,BT,HTILT,SWS)
C       WRITE (6,'('' WEIMER-1996: sw by bz bt angle tilt='',6f7.2)')
C    +                            SWS,BY,BZ,BT,ANGL,TILT
      ENDIF

C          NH assumed latitudes only
      AMLA = ABS(RMLA)
      EPOT = EPOTVAL96(AMLA, RMLT)
      CALL GECMP96(AMLA,RMLT,ET,EP)

      RETURN
      END SUBROUTINE WEIEPOT96
CNCAR
C          Routines common to both 1996 and 2001 versions of Weimer's
C          Ionospheric Electrostatic Potential Model

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
	FUNCTION FSVal(omega,MaxN,FSC)
* Fourier Series Value
* Return value of Sine/Cosine Fourier series for N terms up to MaxN
* at angle omega, given the F.S. coeficients in array FSC
	REAL omega,FSC(0:1,0:*)
	INTEGER MaxN,n
	REAL Y,theta
	Y=0.
	DO n=0,MaxN
	  theta=omega*n
	  Y=Y + FSC(0,n)*COS(theta) + FSC(1,n)*SIN(theta)
	ENDDO
	FSVal=Y
	RETURN
	END FUNCTION FSVal
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
* COORDINATE TRANSFORMATION UTILITIES

CNCAR      Feb 01:  Changed TRANS to GET_TILT s.t. the dipole tilt angle is
C          returned.

	FUNCTION GET_TILT (YEAR,MONTH,DAY,HOUR)
C       SUBROUTINE TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)
CNCAR

	INTEGER YEAR,MONTH,DAY,IDBUG
	REAL HOUR
C         
C      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
C      TRANSFORMATIONS, IDENTIFIED BY K.
C          K=1 TRANSFORMS GSE to GEO
C          K=2     "      GEO to MAG
C          K=3     "      GSE to MAG
C          K=4     "      GSE to GSM
C          K=5     "      GEO to GSM
C          K=6     "      GSM to MAG
C          K=7     "      GSE to GEI
C          K=8     "      GEI to GEO
C          K=9     "      GSM to SM 
C	   K=10    "      GEO to SM 
C	   K=11    "      MAG to SM 
C
C      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
C      FILE UNIT=IDBUG
C	
	INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
	INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
	INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
	INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
	INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

	PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
	PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
	PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
	PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
	PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
	PARAMETER (MAGSM =11,SMMAG =-11)
C
C      The formal names of the coordinate systems are:
C	GSE - Geocentric Solar Ecliptic
C	GEO - Geographic
C	MAG - Geomagnetic
C	GSM - Geocentric Solar Magnetospheric
C	SM  - Solar Magnetic
C	
C      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
C      ST(I) AND CT(I) ARE SINES & COSINES.       
C
C      Program author:  D. R. Weimer
C
C      Some of this code has been copied from subroutines which had been
C      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
C      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
C
C      The formulas for the calculation of Greenwich mean sidereal time (GMST)
C      and the sun's location are from "Almanac for Computers 1990",
C      U.S. Naval Observatory.
C
	REAL UT,T0,GMSTD,GMSTH,ECLIP,MA,LAMD,SUNLON

	COMMON /EIE_WEIMFIELD/EPOCH,TH0,PH0,DIPOLE
	COMMON /EIE_WEITRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
C         
C UM - Changed Data statement to the following:

	EPOCH = 1980. 
	TH0 = 11.19
	PH0 = -70.76
	DIPOLE = .30574

C UM

CNCAR      Feb 01:  Prevent debug printing to a file
	IDBUG = 0
CNCAR

	IF(YEAR.LT.1900)THEN
	  IYR=1900+YEAR
	ELSE
	  IYR=YEAR
	ENDIF
	UT=HOUR
	JD=JULDAY(MONTH,DAY,IYR)
	MJD=JD-2400001
	T0=(FLOAT(MJD)-51544.5)/36525.0
	GMSTD=100.4606184 + 36000.770*T0 + 3.87933E-4*T0*T0 +
     $        15.0410686*UT
	CALL ADJUST(GMSTD)
	GMSTH=GMSTD*24./360.
	ECLIP=23.439 - 0.013*T0
	MA=357.528 + 35999.050*T0 + 0.041066678*UT
	CALL ADJUST(MA)
	LAMD=280.460 + 36000.772*T0 + 0.041068642*UT
	CALL ADJUST(LAMD)
	SUNLON=LAMD + (1.915-0.0048*T0)*SIND(MA) + 0.020*SIND(2.*MA)
	CALL ADJUST(SUNLON)
	IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) YEAR,MONTH,DAY,HOUR
	  WRITE(IDBUG,*) 'MJD=',MJD
	  WRITE(IDBUG,*) 'T0=',T0
	  WRITE(IDBUG,*) 'GMSTH=',GMSTH
	  WRITE(IDBUG,*) 'ECLIPTIC OBLIQUITY=',ECLIP
	  WRITE(IDBUG,*) 'MEAN ANOMALY=',MA
	  WRITE(IDBUG,*) 'MEAN LONGITUDE=',LAMD
	  WRITE(IDBUG,*) 'TRUE LONGITUDE=',SUNLON
	ENDIF

	CX(1)= GMSTD
	CX(2) = ECLIP
	CX(3) = SUNLON
	CX(4) = TH0
	CX(5) = PH0
c Derived later:
c       CX(6) = Dipole tilt angle  
c       CX(7) = Angle between sun and magnetic pole
c       CX(8) = Subsolar point latitude
c       CX(9) = Subsolar point longitude

	DO I=1,5
	  ST(I) = SIND(CX(I))
	  CT(I) = COSD(CX(I))
	ENDDO
C         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0.         
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
C         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0.         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0.         
      AM(3,1,GEIGEO) = 0.         
      AM(3,2,GEIGEO) = 0.         
      AM(3,3,GEIGEO) = 1.         
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) +
     $    AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +  AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
C         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0.
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
C         
      DO I=1,3   
      DO J=1,3   
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) +
     $   AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +  AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
C         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0.) B3 = -B3    
C         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1.
      AM(1,2,GSEGSM) = 0.
      AM(1,3,GSEGSM) = 0.
      AM(2,1,GSEGSM) = 0.
      AM(3,1,GSEGSM) = 0.
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) +
     $    AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) +
     $   AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO
C
	ST(6) = AM(3,1,GSEMAG)       
	CT(6) = SQRT(1.-ST(6)*ST(6))      
	CX(6) = ASIND(ST(6))     

        AM(1,1,GSMSM) = CT(6)
        AM(1,2,GSMSM) = 0.
        AM(1,3,GSMSM) = -ST(6)
        AM(2,1,GSMSM) = 0.
        AM(2,2,GSMSM) = 1.
        AM(2,3,GSMSM) = 0.
        AM(3,1,GSMSM) = ST(6)
        AM(3,2,GSMSM) = 0.
        AM(3,3,GSMSM) = CT(6)
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) +
     $    AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) +
     $   AM(I,2,GSMSM)*AM(J,2,GSMMAG) + AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
C
      CX(7)=ATAN2D( AM(2,1,11) , AM(1,1,11) )
      CX(8)=ASIND( AM(3,1,1) )
      CX(9)=ATAN2D( AM(2,1,1) , AM(1,1,1) )

      IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) 'Dipole tilt angle=',CX(6)
	  WRITE(IDBUG,*) 'Angle between sun and magnetic pole=',CX(7)
	  WRITE(IDBUG,*) 'Subsolar point latitude=',CX(8)
	  WRITE(IDBUG,*) 'Subsolar point longitude=',CX(9)

        DO K=1,11
         WRITE(IDBUG,1001) K
         DO I=1,3
           WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

CNCAR      Mar 96: return the dipole tilt from this function call.
      GET_TILT = CX(6)
CNCAR

      RETURN
      END  FUNCTION GET_TILT
******************************************************************************
CNCAR      Feb 01:  Eliminate unused routines from translib.for: ROTATE,
C          ROTATEV, FROMCART, TOCART, MLT, MAGLONG, SUNLOC.  Remaining
C          are ADJUST and JULDAY
CNCAR
	SUBROUTINE ADJUST(ANGLE)
C	ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
 10	CONTINUE
	IF(ANGLE.LT.0.)THEN
	  ANGLE=ANGLE+360.
	  GOTO 10
	ENDIF
 20	CONTINUE
 	IF(ANGLE.GE.360.)THEN
	  ANGLE=ANGLE-360.
	  GOTO 20
	ENDIF
	RETURN
	END SUBROUTINE ADJUST
******************************************************************************
      FUNCTION JULDAY(MM,ID,IYYY)
      PARAMETER (IGREG=15+31*(10+12*1582))
      IF (IYYY.EQ.0) PAUSE 'There is no Year Zero.'
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF
      RETURN
      END FUNCTION JULDAY

CNCAR      Routines added to work around non-ANSI trig functions which
C          input degrees instead of radians:  SIND, COSD, ASIND, ATAN2D

      FUNCTION SIND (DEG)
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      SIND = SIN (DEG * D2R)
      RETURN
      END FUNCTION SIND

      FUNCTION COSD (DEG)
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

      COSD = COS (DEG * D2R)
      RETURN
      END FUNCTION COSD

      FUNCTION ASIND (RNUM)
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      ASIND = R2D * ASIN (RNUM)
      RETURN
      END FUNCTION ASIND

      FUNCTION ATAN2D (RNUM1,RNUM2)
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      ATAN2D = R2D * ATAN2 (RNUM1,RNUM2)
      RETURN
      END FUNCTION ATAN2D
CNCAR


      end module EIE_ModWeimer
