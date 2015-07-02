C
C     file fish.f
C
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C

!     this module is used by all fishpack solvers to allocate
!     real and complex work space
      MODULE fish
	TYPE fishworkspace
	  REAL, allocatable, DIMENSION(:) :: rew
	  COMPLEX, allocatable ,DIMENSION(:) :: cxw
          INTEGER:: IFAC(15)
          INTEGER:: IFAC2(15)
	END TYPE fishworkspace

	CONTAINS

	SUBROUTINE allocatfish(irwk,icwk,wsave,ierror)
	IMPLICIT NONE
	TYPE (fishworkspace) :: wsave
!       irwk is the required real work space length
!       icwk is the required integer work space length
	INTEGER, INTENT(IN) :: irwk,icwk
!       ierror is set to 20 if the dynamic allocation is unsuccessful
!       (e.g., this would happen if m,n are too large for the computers memory
	INTEGER, INTENT(INOUT) :: ierror
	INTEGER :: istatus
!       first deallocate to avoid memory leakage

 	if(allocated(wsave%rew))DEALLOCATE(wsave%rew)
 	if(allocated(wsave%cxw))DEALLOCATE(wsave%cxw)

!       allocate irwk words of real work space
	if (irwk > 0) then
	     ALLOCATE(wsave%rew(irwk),STAT = istatus)
	end if
!       allocate icwk words of complex work space
	if (icwk > 0) then
	     ALLOCATE(wsave%cxw(icwk),STAT = istatus)
	end if
	ierror = 0
!       flag fatal error if allocation fails
c       IF (istatus /= 0) THEN
	if (istatus .ne. 0 ) then
	  ierror = 20
	END IF
	RETURN
	END SUBROUTINE allocatfish

	SUBROUTINE BLK_space(N,M,irwk,icwk)
!       this subroutine computes the real and complex work space
!       requirements (generous estimate) of blktri for N,M values
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: N,M
	INTEGER,INTENT(OUT) :: irwk,icwk
	INTEGER :: L,log2n
!       compute nearest integer greater than or equal to
!       log base 2 of n+1, i.e., log2n is smallest integer
!       such that 2**log2n >= n+1
	log2n = 1
	do
	   log2n = log2n+1
	   if (n+1 <= 2**log2n) EXIT
	end do
	L = 2**(log2n+1)
	irwk = (log2n-2)*L+5+MAX0(2*N,6*M)+log2n+2*n
	icwk = ((log2n-2)*L+5+log2n)/2+3*M+N
	RETURN
	END SUBROUTINE BLK_space

	SUBROUTINE GEN_space(N,M,irwk)
!       this subroutine computes the real work space
!       requirement (generously) of genbun for the current N,M
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: N,M
	INTEGER,INTENT(OUT) :: irwk
	INTEGER :: log2n
!       compute nearest integer greater than or equal to
!       log base 2 of n+1, i.e., log2n is smallest integer
!       such that 2**log2n >= n+1
	log2n = 1
	do
	   log2n = log2n+1
	   if (n+1 <= 2**log2n) EXIT
	end do
	irwk = 4*N + (10 + log2n)*M
	RETURN
	END SUBROUTINE GEN_space

	SUBROUTINE fishfin(wsave)
!       this subroutine releases allocated work space
!       fishfin should be called after a fishpack solver has finished
!       TYPE (fishworkspace) variable wsave.
	IMPLICIT NONE
	TYPE (fishworkspace) :: wsave
	INTEGER :: istatus

 	if(allocated(wsave%rew))DEALLOCATE(wsave%rew)
 	if(allocated(wsave%cxw))DEALLOCATE(wsave%cxw)

	RETURN
	END SUBROUTINE fishfin

C
C     file blktri.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C     SUBROUTINE BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,
C    +                   IERROR,W)
C
C
C                                                                       
C DIMENSION OF           AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N),
C ARGUMENTS
C                                                                       
C LATEST REVISION        JUNE 2004
C                                                                       
C USAGE                  CALL BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,    
C                                     CM,IDIMY,Y,IERROR,W)              
C                                                                       
C PURPOSE                BLKTRI SOLVES A SYSTEM OF LINEAR EQUATIONS     
C                        OF THE FORM                                    
C                                                                       
C                        AN(J)*X(I,J-1) + AM(I)*X(I-1,J) +              
C                        (BN(J)+BM(I))*X(I,J) + CN(J)*X(I,J+1) +        
C                        CM(I)*X(I+1,J) = Y(I,J)                        
C                                                                       
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.         
C                                                                       
C                        I+1 AND I-1 ARE EVALUATED MODULO M AND         
C                        J+1 AND J-1 MODULO N, I.E.,                    
C                                                                       
C                        X(I,0) = X(I,N),  X(I,N+1) = X(I,1),           
C                        X(0,J) = X(M,J),  X(M+1,J) = X(1,J).           
C                                                                       
C                        THESE EQUATIONS USUALLY RESULT FROM THE        
C                        DISCRETIZATION OF SEPARABLE ELLIPTIC           
C                        EQUATIONS.  BOUNDARY CONDITIONS MAY BE         
C                        DIRICHLET, NEUMANN, OR PERIODIC.               
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT               IFLG                                           
C                                                                       
C                          = 0  INITIALIZATION ONLY.                    
C                               CERTAIN QUANTITIES THAT DEPEND ON NP,   
C                               N, AN, BN, AND CN ARE COMPUTED AND      
C                               STORED IN DERIVED data type w (see
c                               description of w below)
C                                                                       
C                          = 1  THE QUANTITIES THAT WERE COMPUTED       
C                               IN THE INITIALIZATION ARE USED          
C                               TO OBTAIN THE SOLUTION X(I,J).          
C                                                                       
C                               NOTE:                                   
C                               A CALL WITH IFLG=0 TAKES                
C                               APPROXIMATELY ONE HALF THE TIME         
C                               AS A CALL WITH IFLG = 1.                
C                               HOWEVER, THE INITIALIZATION DOES        
C                               NOT HAVE TO BE REPEATED UNLESS NP,      
C                               N, AN, BN, OR CN CHANGE.                
C                                                                       
C                        NP                                             
C                          = 0  IF AN(1) AND CN(N) ARE NOT ZERO,        
C                               WHICH CORRESPONDS TO PERIODIC           
C                               BOUNARY CONDITIONS.                     
C                                                                       
C                          = 1  IF AN(1) AND CN(N) ARE ZERO.            
C                                                                       
C                        N                                              
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.   
C                          N MUST BE GREATER THAN 4.                    
C                          THE OPERATION COUNT IS PROPORTIONAL TO       
C                          MNLOG2(N), HENCE N SHOULD BE SELECTED        
C                          LESS THAN OR EQUAL TO M.                     
C                                                                       
C                        AN,BN,CN                                       
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH N           
C                          THAT SPECIFY THE COEFFICIENTS IN THE         
C                          LINEAR EQUATIONS GIVEN ABOVE.                
C                                                                       
C                        MP                                             
C                          = 0  IF AM(1) AND CM(M) ARE NOT ZERO,        
C                               WHICH CORRESPONDS TO PERIODIC           
C                               BOUNDARY CONDITIONS.                    
C                                                                       
C                          = 1  IF AM(1) = CM(M) = 0  .                 
C                                                                       
C                        M                                              
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.   
C                           M MUST BE GREATER THAN 4.                   
C                                                                       
C                        AM,BM,CM                                       
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT      
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR       
C                          EQUATIONS GIVEN ABOVE.                       
C                                                                       
C                        IDIMY                                          
C                          THE ROW (OR FIRST) DIMENSION OF THE          
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS        
C                          IN THE PROGRAM CALLING BLKTRI.               
C                          THIS PARAMETER IS USED TO SPECIFY THE        
C                          VARIABLE DIMENSION OF Y.                     
C                          IDIMY MUST BE AT LEAST M.                    
C                                                                       
C                        Y                                              
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES       
C                          THE VALUES OF THE RIGHT SIDE OF THE LINEAR   
C                          SYSTEM OF EQUATIONS GIVEN ABOVE.             
C                          Y MUST BE DIMENSIONED AT LEAST M*N.          
C                                                                       
C                        W
c                          A fortran 90 derived TYPE (fishworkspace) variable
c                          that must be declared by the user.  The first
c                          two declarative statements in the user program
c                          calling BLKTRI must be:
c
c                               USE fish
c                               TYPE (fishworkspace) :: W
c
c                          The first statement makes the fishpack module
c                          defined in the file "fish.f" available to the
c                          user program calling BLKTRI.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in BLKTRI to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (IFLG=0) call to BLKTRI.  These
c                          must be preserved on non-initial calls (IFLG=1)
c                          to BLKTRI.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling BLKTRI should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          BLKTRI.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
C
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON OUTPUT              Y                                              
C                          CONTAINS THE SOLUTION X.                     
C                                                                       
C                        IERROR                                         
C                          AN ERROR FLAG THAT INDICATES INVALID         
C                          INPUT PARAMETERS.  EXCEPT FOR NUMBER ZER0,   
C                          A SOLUTION IS NOT ATTEMPTED.                 
C                                                                       
C                        = 0  NO ERROR.                                 
C                        = 1  M IS LESS THAN 5                          
C                        = 2  N IS LESS THAN 5                          
C                        = 3  IDIMY IS LESS THAN M.                     
C                        = 4  BLKTRI FAILED WHILE COMPUTING RESULTS     
C                             THAT DEPEND ON THE COEFFICIENT ARRAYS     
C                             AN, BN, CN.  CHECK THESE ARRAYS.          
C                        = 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.  
C                                                                       
C                             POSSIBLE REASONS FOR THIS CONDITION ARE   
C                             1. THE ARRAYS AN AND CN ARE NOT CORRECT   
C                             2. TOO LARGE A GRID SPACING WAS USED      
C                                IN THE DISCRETIZATION OF THE ELLIPTIC  
C                                EQUATION.                              
C                             3. THE LINEAR EQUATIONS RESULTED FROM A   
C                                PARTIAL DIFFERENTIAL EQUATION WHICH    
C                                WAS NOT ELLIPTIC.                      
C                                                                       
C                        = 20 If the dynamic allocation of real and
C                             complex work space in the derived type
C                             (fishworkspace) variable W fails (e.g.,
c                             if N,M are too large for the platform used)
C
C                                                                       
C                        W
c                             The derived type (fishworkspace) variable W
c                             contains real and complex values that must not
C                             be destroyed if BLKTRI is called again with
C                             IFLG=1.
C                                                                       
C                                                                       
C SPECIAL CONDITIONS     THE ALGORITHM MAY FAIL IF ABS(BM(I)+BN(J))     
C                        IS LESS THAN ABS(AM(I))+ABS(AN(J))+            
C                        ABS(CM(I))+ABS(CN(J))                          
C                        FOR SOME I AND J. THE ALGORITHM WILL ALSO      
C                        FAIL IF AN(J)*CN(J-1) IS LESS THAN ZERO FOR    
C                        SOME J.                                        
C                        SEE THE DESCRIPTION OF THE OUTPUT PARAMETER    
C                        IERROR.                                        
C                                                                       
C I/O                    NONE                                           
C                                                                       
C PRECISION              SINGLE                                         
C                                                                       
C REQUIRED FILES         fish.f,comf.f
C                                                                       
C LANGUAGE               FORTRAN 90
C                                                                       
C HISTORY                WRITTEN BY PAUL SWARZTRAUBER AT NCAR IN THE    
C                        EARLY 1970'S.  REWRITTEN AND RELEASED IN       
C                        LIBRARIES IN JANUARY 1980. Revised in June
C                        2004 using Fortan 90 dynamically allocated work
c                        space and derived data types to eliminate mixed
c                        mode conflicts in the earlier versions.
C                                                                       
C ALGORITHM              GENERALIZED CYCLIC REDUCTION                   
C                                                                       
C PORTABILITY            FORTRAN 90.  APPROXIMATE MACHINE ACCURACY
C                        IS COMPUTED IN FUNCTION EPMACH.                
C                                                                       
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, 'EFFICIENT       
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF        
C                        ELLIPTIC EQUATIONS'                            
C                        NCAR TN/IA-109, JULY, 1975, 138 PP.            
C                                                                       
C                        SWARZTRAUBER P. N.,A DIRECT METHOD FOR         
C                        THE DISCRETE SOLUTION OF SEPARABLE             
C                        ELLIPTIC EQUATIONS, S.I.A.M.                   
C                        J. NUMER. ANAL.,11(1974) PP. 1136-1150.        
C***********************************************************************
      SUBROUTINE BLKTRI(IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, 
     1   IDIMY, Y, IERROR, W)

      Implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IFLG
      INTEGER  :: NP
      INTEGER  :: N
      INTEGER  :: MP
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL  :: Y(IDIMY,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
      IF (M < 5) THEN
         IERROR = 1
         RETURN 
      ENDIF
      IF (N < 3) THEN
         IERROR = 2
         RETURN 
      ENDIF
      IF (IDIMY < M) THEN
         IERROR = 3
         RETURN 
      ENDIF
      IF (IFLG == 0) THEN
!     compute and allocate real and complex required work space
         CALL BLK_SPACE (N, M, IRWK, ICWK)
         CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
         IF (IERROR == 20) RETURN 
      ENDIF
      call blktrii (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,
     +             IERROR,w%rew,w%cxw)
      RETURN 
      END SUBROUTINE BLKTRI


 
      SUBROUTINE BLKTRII(IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, 
     1   IDIMY, Y, IERROR, W, WC)
      Implicit none
!!!      external prod,prodp,cprod,cprodp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IFLG
      INTEGER , INTENT(IN) :: NP
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: MP
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL  :: Y(IDIMY,*)
      REAL  :: W(*)
      COMPLEX  :: WC(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IW1, IW2, IW3, IWW, IWU, IWD, NL, NH, IWAH, IWBH
      integer :: iwc1=1, iwc2=1, iwc3=1

      SAVE IW1, IW2, IW3, IWW, IWU, IWD, NL
C-----------------------------------------------
 
C
C TEST M AND N FOR THE PROPER FORM
C
      NM = N
!     check again for solvers which call blktrii directly
      IF (M < 5) THEN
         IERROR = 1
         RETURN 
      ENDIF
      IF (NM < 3) THEN
         IERROR = 2
         RETURN 
      ENDIF
      IF (IDIMY < M) THEN
         IERROR = 3
         RETURN 
      ENDIF
 
      IF (IFLG == 0) THEN
         NH = N
         NPP = NP
         IF (NPP /= 0) THEN
            NH = NH + 1
         ENDIF
         IK = 2
         K = 1
         IK = IK + IK
         K = K + 1
         DO WHILE(NH - IK > 0)
            IK = IK + IK
            K = K + 1
         END DO
         NL = IK
         IK = IK + IK
         NL = NL - 1
         IWAH = (K - 2)*IK + K + 5
         IF (NPP == 0) THEN
            IWBH = IWAH + NM + NM
            IW1 = IWBH
            NM = N - 1
         ELSE
            IW1 = IWAH
            IWBH = IW1 + NM
         ENDIF
!     set pointers in real,complex work space
         IW2 = IW1 + M
         IW3 = IW2 + M
         IWD = IW3 + M
         IWW = IWD + M
         IWU = IWW + M
         CALL COMPB (NL, IERROR, AN, BN, CN, W, WC, W(IWAH), W(IWBH))
         RETURN 
      ENDIF
! *** important to reset nm for np = 0
      IF (NPP == 0) NM = N - 1

      if(NCMPLX==0)then
         IWC1 = 1; IWC2 = 1; IWC3 = 1
      else
         IWC1 = IW1; IWC2 = IW2; IWC3 = IW3
      endif
 
      IF (MP /= 0) THEN
         CALL BLKTR1 (NL, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, W, WC,
     1        W(IW1), W(IW2), W(IW3), W(IWD), W(IWW), W(IWU), 
     2        WC(IWC1), WC(IWC2), WC(IWC3), PROD, CPROD)

      ELSE
         CALL BLKTR1 (NL, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, W, WC, 
     1      W(IW1), W(IW2), W(IW3), W(IWD), W(IWW), W(IWU), 
     2      WC(IWC1), WC(IWC2), WC(IWC3), PRODP, CPRODP)
      ENDIF
      RETURN 
      END SUBROUTINE BLKTRII


      SUBROUTINE BLKTR1(N, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, B, BC, 
     1   W1, W2, W3, WD, WW, WU, CW1, CW2, CW3, PRDCT, CPRDCT)
      Implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      INTEGER  :: M
      INTEGER , INTENT(IN) :: IDIMY
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL  :: Y(IDIMY,*)
      REAL  :: B(*)
      REAL  :: W1(*)
      REAL  :: W2(*)
      REAL  :: W3(*)
      REAL  :: WD(*)
      REAL  :: WW(*)
      REAL  :: WU(*)
      COMPLEX  :: BC(*)
      COMPLEX  :: CW1(*)
      COMPLEX  :: CW2(*)
      COMPLEX  :: CW3(*)

      EXTERNAL :: PRDCT, CPRDCT

C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: KDO, L, IR, I2, I1, I3, I4, IRM1, IM2, NM2, IM3, NM3, 
     1   IM1, NM1, I0, IF, I, IPI1, IPI2, IPI3, IDXC, NC, IDXA, NA, IP2
     2   , NP2, IP1, NP1, IP3, NP3, J, IZ, NZ, IZR, LL, IFD, IP, NP, 
     3   IMI1, IMI2
      REAL :: DUM(1)
C-----------------------------------------------
C
C BLKTR1 SOLVES THE LINEAR SYSTEM
C
C B  CONTAINS THE ROOTS OF ALL THE B POLYNOMIALS
C W1,W2,W3,WD,WW,WU  ARE ALL WORKING ARRAYS
C PRDCT  IS EITHER PRODP OR PROD DEPENDING ON WHETHER THE BOUNDARY
C CONDITIONS IN THE M DIRECTION ARE PERIODIC OR NOT
C CPRDCT IS EITHER CPRODP OR CPROD WHICH ARE THE COMPLEX VERSIONS
C OF PRODP AND PROD. THESE ARE CALLED IN THE EVENT THAT SOME
C OF THE ROOTS OF THE B SUB P POLYNOMIAL ARE COMPLEX
C
C
C
C BEGIN REDUCTION PHASE
C
      KDO = K - 1
      DO L = 1, KDO
         IR = L - 1
         I2 = 2**IR
         I1 = I2/2
         I3 = I2 + I1
         I4 = I2 + I2
         IRM1 = IR - 1
         CALL INDXB (I2, IR, IM2, NM2)
         CALL INDXB (I1, IRM1, IM3, NM3)
         CALL INDXB (I3, IRM1, IM1, NM1)
         I0 = 0

         if(NM1==0) IM1=1
         if(NM2==0) IM2=1
         if(NM3==0) IM3=1

         CALL PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), I0, DUM, 
     1      Y(1,I2), W3, M, AM, BM, CM, WD, WW, WU)
         IF = 2**K
         DO I = I4, IF, I4
            IF (I - NM > 0) CYCLE 
            IPI1 = I + I1
            IPI2 = I + I2
            IPI3 = I + I3
            CALL INDXC (I, IR, IDXC, NC)
            IF (I - IF >= 0) CYCLE 
            CALL INDXA (I, IR, IDXA, NA)
            CALL INDXB (I - I1, IRM1, IM1, NM1)
            CALL INDXB (IPI2, IR, IP2, NP2)
            CALL INDXB (IPI1, IRM1, IP1, NP1)
            CALL INDXB (IPI3, IRM1, IP3, NP3)

            if(NM1==0) IM1 = 1
            CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W3, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            IF (IPI2 - NM > 0) THEN
               W3(:M) = 0.
               W2(:M) = 0.
            ELSE
               if(NP1==0) IP1=1
               if(NP2==0) IP2=1
               if(NP3==0) IP3=1

               CALL PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM
     1            , Y(1,IPI2), W3, M, AM, BM, CM, WD, WW, WU)
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W3
     1            , W2, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            Y(:M,I) = W1(:M) + W2(:M) + Y(:M,I)
         END DO
      END DO
      IF (NPP == 0) THEN
         IF = 2**K
         I = IF/2
         I1 = I/2
         CALL INDXB (I - I1, K - 2, IM1, NM1)
         CALL INDXB (I + I1, K - 2, IP1, NP1)
         CALL INDXB (I, K - 1, IZ, NZ)
         CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, Y(1,I)
     1      , W1, M, AM, BM, CM, WD, WW, WU)
         IZR = I
         W2(:M) = W1(:M)
         DO LL = 2, K
            L = K - LL + 1
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I = I2
            CALL INDXC (I, IR, IDXC, NC)
            CALL INDXB (I, IR, IZ, NZ)
            CALL INDXB (I - I1, IR - 1, IM1, NM1)
            CALL INDXB (I + I1, IR - 1, IP1, NP1)
            CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W1, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            W1(:M) = Y(:M,I) + W1(:M)
            CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1
     1         , W1, M, AM, BM, CM, WD, WW, WU)
         END DO
         L118: DO LL = 2, K
            L = K - LL + 1
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I4 = I2 + I2
            IFD = IF - I2
            DO I = I2, IFD, I4
               IF (I - I2 - IZR /= 0) CYCLE 
               IF (I - NM > 0) CYCLE  L118
               CALL INDXA (I, IR, IDXA, NA)
               CALL INDXB (I, IR, IZ, NZ)
               CALL INDXB (I - I1, IR - 1, IM1, NM1)
               CALL INDXB (I + I1, IR - 1, IP1, NP1)
               CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W2
     1            , W2, M, AM, BM, CM, WD, WW, WU)
               W2(:M) = Y(:M,I) + W2(:M)
               CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, 
     1            W2, W2, M, AM, BM, CM, WD, WW, WU)
               IZR = I
               IF (I - NM == 0) EXIT  L118
            END DO
         END DO L118
  119    CONTINUE
         Y(:M,NM+1) = Y(:M,NM+1) - CN(NM+1)*W1(:M) - AN(NM+1)*W2(:M)
         CALL INDXB (IF/2, K - 1, IM1, NM1)
         CALL INDXB (IF, K - 1, IP, NP)
         IF (NCMPLX /= 0) THEN
            CALL CPRDCT (NM + 1, BC(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(
     1         1,NM+1), Y(1,NM+1), M, AM, BM, CM, CW1, CW2, CW3)
         ELSE
            CALL PRDCT (NM + 1, B(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(1,
     1         NM+1), Y(1,NM+1), M, AM, BM, CM, WD, WW, WU)
         ENDIF
         W1(:M) = AN(1)*Y(:M,NM+1)
         W2(:M) = CN(NM)*Y(:M,NM+1)
         Y(:M,1) = Y(:M,1) - W1(:M)
         Y(:M,NM) = Y(:M,NM) - W2(:M)
         DO L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I4 = I2 + I2
            I1 = I2/2
            I = I4
            CALL INDXA (I, IR, IDXA, NA)
            CALL INDXB (I - I2, IR, IM2, NM2)
            CALL INDXB (I - I2 - I1, IR - 1, IM3, NM3)
            CALL INDXB (I - I1, IR - 1, IM1, NM1)
            CALL PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), 0, DUM, 
     1         W1, W1, M, AM, BM, CM, WD, WW, WU)
            CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W1, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            Y(:M,I) = Y(:M,I) - W1(:M)
         END DO
C
         IZR = NM
         L131: DO L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I3 = I2 + I1
            I4 = I2 + I2
            IRM1 = IR - 1
            DO I = I4, IF, I4
               IPI1 = I + I1
               IPI2 = I + I2
               IPI3 = I + I3
               IF (IPI2 - IZR /= 0) THEN
                  IF (I - IZR /= 0) CYCLE 
                  CYCLE  L131
               ENDIF
               CALL INDXC (I, IR, IDXC, NC)
               CALL INDXB (IPI2, IR, IP2, NP2)
               CALL INDXB (IPI1, IRM1, IP1, NP1)
               CALL INDXB (IPI3, IRM1, IP3, NP3)
               CALL PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM
     1            , W2, W2, M, AM, BM, CM, WD, WW, WU)
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W2
     1            , W2, M, AM, BM, CM, WD, WW, WU)
               Y(:M,I) = Y(:M,I) - W2(:M)
               IZR = I
               CYCLE  L131
            END DO
         END DO L131
      ENDIF
C
C BEGIN BACK SUBSTITUTION PHASE
C
      DO LL = 1, K
         L = K - LL + 1
         IR = L - 1
         IRM1 = IR - 1
         I2 = 2**IR
         I1 = I2/2
         I4 = I2 + I2
         IFD = IF - I2
         DO I = I2, IFD, I4
            IF (I - NM > 0) CYCLE 
            IMI1 = I - I1
            IMI2 = I - I2
            IPI1 = I + I1
            IPI2 = I + I2
            CALL INDXA (I, IR, IDXA, NA)
            CALL INDXC (I, IR, IDXC, NC)
            CALL INDXB (I, IR, IZ, NZ)
            CALL INDXB (IMI1, IRM1, IM1, NM1)
            CALL INDXB (IPI1, IRM1, IP1, NP1)
            IF (I - I2 <= 0) THEN
               W1(:M) = 0.
            ELSE
               if(NM1 == 0) IM1 = 1
               CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), Y(
     1            1,IMI2), W1, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            IF (IPI2 - NM > 0) THEN
               W2(:M) = 0.
            ELSE
               if(NP1 == 0) IP1 = 1
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), Y(
     1            1,IPI2), W2, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            W1(:M) = Y(:M,I) + W1(:M) + W2(:M)
            if( NZ == 0) IZ = 1
            if(NM1 == 0) IM1 = 1
            if(NP1 == 0) IP1 = 1
            CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1
     1         , Y(1,I), M, AM, BM, CM, WD, WW, WU)
         END DO
      END DO
      RETURN 
      END SUBROUTINE BLKTR1


      REAL FUNCTION BSRH (XLL, XRR, IZ, C, A, BH, F, SGN)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IZ
      REAL , INTENT(IN) :: XLL
      REAL , INTENT(IN) :: XRR
      REAL  :: F
      REAL , INTENT(IN) :: SGN
      REAL  :: C(*)
      REAL  :: A(*)
      REAL  :: BH(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: R1, XL, XR, DX, X
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
      XL = XLL
      XR = XRR
      DX = 0.5*ABS(XR - XL)
      X = 0.5*(XL + XR)
      R1 = SGN*F(X,IZ,C,A,BH)
      IF (R1 >= 0.) THEN
         IF (R1 == 0.) GO TO 105
         XR = X
      ELSE
         XL = X
      ENDIF
      DX = 0.5*DX
      DO WHILE(DX - CNV > 0.)
         X = 0.5*(XL + XR)
         R1 = SGN*F(X,IZ,C,A,BH)
         IF (R1 >= 0.) THEN
            IF (R1 == 0.) GO TO 105
            XR = X
         ELSE
            XL = X
         ENDIF
         DX = 0.5*DX
      END DO
  105 CONTINUE
      BSRH = 0.5*(XL + XR)
      RETURN 
      END FUNCTION BSRH


      SUBROUTINE COMPB(N, IERROR, AN, BN, CN, B, BC, AH, BH)
      implicit none
!!!      real :: epmach
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      INTEGER  :: IERROR
      REAL  :: AN(*)
      REAL , INTENT(IN) :: BN(*)
      REAL  :: CN(*)
      REAL  :: B(*)
      REAL  :: AH(*)
      REAL  :: BH(*)
      COMPLEX  :: BC(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, IF, KDO, L, IR, I2, I4, IPL, IFD, I, IB, NB, JS, JF
     1   , LS, LH, NMP, L1, L2, J2, J1, N2M2
      REAL :: DUM, BNORM, ARG, D1, D2, D3
C-----------------------------------------------
C
C     COMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS USING SUBROUTINE
C     TEVLS WHICH IS A MODIFICATION THE EISPACK PROGRAM TQLRAT.
C     IERROR IS SET TO 4 IF EITHER TEVLS FAILS OR IF A(J+1)*C(J) IS
C     LESS THAN ZERO FOR SOME J.  AH,BH ARE TEMPORARY WORK ARRAYS.
C
      EPS = EPMACH(DUM)
      BNORM = ABS(BN(1))
      DO J = 2, NM
         BNORM = AMAX1(BNORM,ABS(BN(J)))
         ARG = AN(J)*CN(J-1)
         IF (ARG < 0.) GO TO 119
         B(J) = SIGN(SQRT(ARG),AN(J))
      END DO
      CNV = EPS*BNORM
      IF = 2**K
      KDO = K - 1
      L108: DO L = 1, KDO
         IR = L - 1
         I2 = 2**IR
         I4 = I2 + I2
         IPL = I4 - 1
         IFD = IF - I4
         DO I = I4, IFD, I4
            CALL INDXB (I, L, IB, NB)
            IF (NB <= 0) CYCLE  L108
            JS = I - IPL
            JF = JS + NB - 1
            LS = 0
            BH(:JF-JS+1) = BN(JS:JF)
            AH(:JF-JS+1) = B(JS:JF)
            CALL TEVLS (NB, BH, AH, IERROR)
            IF (IERROR /= 0) GO TO 118
            LH = IB - 1
            IF (NB > 0) THEN
               B(LH+1:NB+LH) = -BH(:NB)
               LH = NB + LH
            ENDIF
         END DO
      END DO L108
      B(:NM) = -BN(:NM)
      IF (NPP == 0) THEN
         NMP = NM + 1
         NB = NM + NMP
         DO J = 1, NB
            L1 = MOD(J - 1,NMP) + 1
            L2 = MOD(J + NM - 1,NMP) + 1
            ARG = AN(L1)*CN(L2)
            IF (ARG < 0.) GO TO 119
            BH(J) = SIGN(SQRT(ARG),(-AN(L1)))
            AH(J) = -BN(L1)
         END DO
         CALL TEVLS (NB, AH, BH, IERROR)
         IF (IERROR /= 0) GO TO 118
         CALL INDXB (IF, K - 1, J2, LH)
         CALL INDXB (IF/2, K - 1, J1, LH)
         J2 = J2 + 1
         LH = J2
         N2M2 = J2 + NM + NM - 2
  114    CONTINUE
         D1 = ABS(B(J1)-B(J2-1))
         D2 = ABS(B(J1)-B(J2))
         D3 = ABS(B(J1)-B(J2+1))
         IF (D2>=D1 .OR. D2>=D3) THEN
            B(LH) = B(J2)
            J2 = J2 + 1
            LH = LH + 1
            IF (J2 - N2M2 <= 0) GO TO 114
         ELSE
            J2 = J2 + 1
            J1 = J1 + 1
            IF (J2 - N2M2 <= 0) GO TO 114
         ENDIF
         B(LH) = B(N2M2+1)
         CALL INDXB (IF, K - 1, J1, J2)
         J2 = J1 + NMP + NMP
         CALL PPADD (NM + 1, IERROR, AN, CN, BC(J1), B(J1), B(J2))
      ENDIF
      RETURN 
  118 CONTINUE
      IERROR = 4
      RETURN 
  119 CONTINUE
      IERROR = 5
      RETURN 
      END SUBROUTINE COMPB


      SUBROUTINE CPROD(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,YY,M,A,B,C,D,W,Y)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      REAL , INTENT(IN) :: X(*)
      REAL , INTENT(OUT) :: YY(*)
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      COMPLEX , INTENT(IN) :: BD(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, ID, M1, M2, IA, IFLG, K
      REAL :: RT
      COMPLEX :: CRT, DEN, Y1, Y2
C-----------------------------------------------
C
C PROD APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN YY           (COMPLEX CASE)
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C NA IS THE LENGTH OF THE ARRAY AA
C X,YY THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS YY
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,W,Y ARE WORKING ARRAYS
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      DO J = 1, M
         Y(J) = CMPLX(X(J),0.)
      END DO
      MM = M - 1
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IFLG = 0
      IF (ID > 0) THEN
         CRT = BD(ID)
         ID = ID - 1
C
C BEGIN SOLUTION TO SYSTEM
C
         D(M) = A(M)/(B(M)-CRT)
         W(M) = Y(M)/(B(M)-CRT)
         DO J = 2, MM
            K = M - J
            DEN = B(K+1) - CRT - C(K+1)*D(K+2)
            D(K+1) = A(K+1)/DEN
            W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
         END DO
         DEN = B(1) - CRT - C(1)*D(2)
         IF (CABS(DEN) /= 0.) THEN
            Y(1) = (Y(1)-C(1)*W(2))/DEN
         ELSE
            Y(1) = (1.,0.)
         ENDIF
         DO J = 2, M
            Y(J) = W(J) - D(J)*Y(J-1)
         END DO
      ENDIF
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 121
         RT = BM2(M2)
         M2 = M2 - 1
      ELSE
         IF (M2 <= 0) THEN
            RT = BM1(M1)
            M1 = M1 - 1
         ELSE
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) > 0.) THEN
               RT = BM1(M1)
               M1 = M1 - 1
            ELSE
               RT = BM2(M2)
               M2 = M2 - 1
            ENDIF
         ENDIF
      ENDIF
      Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2)
      IF (MM - 2 >= 0) THEN
         DO J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
         END DO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M)
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  121 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IA = IA - 1
         IFLG = 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*Y(:M)
      ENDIF
      IF (IFLG > 0) GO TO 102
      DO J = 1, M
         YY(J) = REAL(Y(J))
      END DO
      RETURN 
      END SUBROUTINE CPROD


      SUBROUTINE CPRODP(ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, YY, M, A
     1   , B, C, D, U, Y)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      REAL , INTENT(IN) :: X(*)
      REAL , INTENT(OUT) :: YY(*)
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      COMPLEX , INTENT(IN) :: BD(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: U(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, MM2, ID, M1, M2, IA, IFLG, K
      REAL :: RT
      COMPLEX :: V, DEN, BH, YM, AM, Y1, Y2, YH, CRT
C-----------------------------------------------
C
C PRODP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN YY       PERIODIC BOUNDARY CONDITIONS
C AND  COMPLEX  CASE
C
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,YY THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS YY
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,U,Y ARE WORKING ARRAYS
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      DO J = 1, M
         Y(J) = CMPLX(X(J),0.)
      END DO
      MM = M - 1
      MM2 = M - 2
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IFLG = 0
      IF (ID > 0) THEN
         CRT = BD(ID)
         ID = ID - 1
         IFLG = 1
C
C BEGIN SOLUTION TO SYSTEM
C
         BH = B(M) - CRT
         YM = Y(M)
         DEN = B(1) - CRT
         D(1) = C(1)/DEN
         U(1) = A(1)/DEN
         Y(1) = Y(1)/DEN
         V = CMPLX(C(M),0.)
         IF (MM2 - 2 >= 0) THEN
            DO J = 2, MM2
               DEN = B(J) - CRT - A(J)*D(J-1)
               D(J) = C(J)/DEN
               U(J) = -A(J)*U(J-1)/DEN
               Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
               BH = BH - V*U(J-1)
               YM = YM - V*Y(J-1)
               V = -V*D(J-1)
            END DO
         ENDIF
         DEN = B(M-1) - CRT - A(M-1)*D(M-2)
         D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
         Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN
         AM = A(M) - V*D(M-2)
         BH = BH - V*U(M-2)
         YM = YM - V*Y(M-2)
         DEN = BH - AM*D(M-1)
         IF (CABS(DEN) /= 0.) THEN
            Y(M) = (YM - AM*Y(M-1))/DEN
         ELSE
            Y(M) = (1.,0.)
         ENDIF
         Y(M-1) = Y(M-1) - D(M-1)*Y(M)
         DO J = 2, MM
            K = M - J
            Y(K) = Y(K) - D(K)*Y(K+1) - U(K)*Y(M)
         END DO
      ENDIF
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 123
         RT = BM2(M2)
         M2 = M2 - 1
      ELSE
         IF (M2 <= 0) THEN
            RT = BM1(M1)
            M1 = M1 - 1
         ELSE
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) > 0.) THEN
               RT = BM1(M1)
               M1 = M1 - 1
            ELSE
               RT = BM2(M2)
               M2 = M2 - 1
C
C MATRIX MULTIPLICATION
C
            ENDIF
         ENDIF
      ENDIF
      YH = Y(1)
      Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2) + A(1)*Y(M)
      IF (MM - 2 >= 0) THEN
         DO J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
         END DO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M) + C(M)*YH
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  123 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IA = IA - 1
         IFLG = 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*Y(:M)
      ENDIF
      IF (IFLG > 0) GO TO 102
      DO J = 1, M
         YY(J) = REAL(Y(J))
      END DO
      RETURN 
      END SUBROUTINE CPRODP


      SUBROUTINE INDXA(I, IR, IDXA, NA)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDXA
      INTEGER , INTENT(OUT) :: NA
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
      NA = 2**IR
      IDXA = I - NA + 1
      IF (I - NM > 0) THEN
         NA = 0
      ENDIF
      RETURN 
      END SUBROUTINE INDXA


      SUBROUTINE INDXB(I, IR, IDX, IDP)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDX
      INTEGER , INTENT(OUT) :: IDP
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IZH, ID, IPL
C-----------------------------------------------
C
C B(IDX) IS THE LOCATION OF THE FIRST ROOT OF THE B(I,IR) POLYNOMIAL
C
      IDP = 0
      IF (IR >= 0) THEN
         IF (IR <= 0) THEN
            IF (I - NM > 0) GO TO 107
            IDX = I
            IDP = 1
            RETURN 
         ENDIF
         IZH = 2**IR
         ID = I - IZH - IZH
         IDX = ID + ID + (IR - 1)*IK + IR + (IK - I)/IZH + 4
         IPL = IZH - 1
         IDP = IZH + IZH - 1
         IF (I - IPL - NM > 0) THEN
            IDP = 0
            RETURN 
         ENDIF
         IF (I + IPL - NM > 0) THEN
            IDP = NM + IPL - I + 1
         ENDIF
      ENDIF
  107 CONTINUE
      RETURN 
      END SUBROUTINE INDXB


      SUBROUTINE INDXC(I, IR, IDXC, NC)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDXC
      INTEGER , INTENT(OUT) :: NC
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
      NC = 2**IR
      IDXC = I
      IF (IDXC + NC - 1 - NM > 0) THEN
         NC = 0
      ENDIF
      RETURN 
      END SUBROUTINE INDXC


      SUBROUTINE PPADD(N, IERROR, A, C, CBP, BP, BH)
      implicit none
!!!      external psgf,ppspf,ppsgf,bsrh
!!!      real psgf,ppspf,ppsgf,bsrh
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: IERROR
      REAL  :: A(*)
      REAL  :: C(*)
      REAL , INTENT(INOUT) :: BP(*)
      REAL  :: BH(*)
      COMPLEX , INTENT(INOUT) :: CBP(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::IZ,IZM,IZM2,J,NT,MODIZ,IS,IF,IG,IT,ICV,I3,I2,NHALF
      REAL :: R4, R5, R6, SCNV, XL, DB, SGN, XR, XM, PSG
      COMPLEX :: CF, CX, FSG, HSG, DD, F, FP, FPP, CDIS, R1, R2, R3
C-----------------------------------------------
C
C     PPADD COMPUTES THE EIGENVALUES OF THE PERIODIC TRIDIAGONAL MATRIX
C     WITH COEFFICIENTS AN,BN,CN
C
C N IS THE ORDER OF THE BH AND BP POLYNOMIALS
C     ON OUTPUT BP CONTIANS THE EIGENVALUES
C CBP IS THE SAME AS BP EXCEPT TYPE COMPLEX
C BH IS USED TO TEMPORARILY STORE THE ROOTS OF THE B HAT POLYNOMIAL
C WHICH ENTERS THROUGH BP
C
      SCNV = SQRT(CNV)
      IZ = N
      IZM = IZ - 1
      IZM2 = IZ - 2
      IF (BP(N) - BP(1) <= 0.) THEN
         IF (BP(N) - BP(1) == 0.) GO TO 142
         BH(:N) = BP(N:1:(-1))
      ELSE
         BH(:N) = BP(:N)
      ENDIF
      NCMPLX = 0
      MODIZ = MOD(IZ,2)
      IS = 1
      IF (MODIZ /= 0) THEN
         IF (A(1) < 0.) GO TO 110
         IF (A(1) == 0.) GO TO 142
      ENDIF
      XL = BH(1)
      DB = BH(3) - BH(1)
      XL = XL - DB
      R4 = PSGF(XL,IZ,C,A,BH)
      DO WHILE(R4 <= 0.)
         XL = XL - DB
         R4 = PSGF(XL,IZ,C,A,BH)
      END DO
      SGN = -1.
      CBP(1) = CMPLX(BSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.)
      BP(1) = REAL(CBP(1))
      IS = 2
  110 CONTINUE
      IF = IZ - 1
      IF (MODIZ /= 0) THEN
         IF (A(1) > 0.) GO TO 115
         IF (A(1) == 0.) GO TO 142
      ENDIF
      XR = BH(IZ)
      DB = BH(IZ) - BH(IZ-2)
      XR = XR + DB
      R5 = PSGF(XR,IZ,C,A,BH)
      DO WHILE(R5 < 0.)
         XR = XR + DB
         R5 = PSGF(XR,IZ,C,A,BH)
      END DO
      SGN = 1.
      CBP(IZ) = CMPLX(BSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.)
      IF = IZ - 2
  115 CONTINUE
      DO IG = IS, IF, 2
         XL = BH(IG)
         XR = BH(IG+1)
         SGN = -1.
         XM = BSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)
         PSG = PSGF(XM,IZ,C,A,BH)
         IF (ABS(PSG) - EPS <= 0.) GO TO 118
         R6 = PSG*PPSGF(XM,IZ,C,A,BH)
         IF (R6 > 0.) GO TO 119
         IF (R6 == 0.) GO TO 118
         SGN = 1.
         CBP(IG) = CMPLX(BSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.)
c        bp(ig) = real(cbp(ig))
         SGN = -1.
         CBP(IG+1) = CMPLX(BSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.)
c        bp(ig) = real(cbp(ig))
c        bp(ig+1) = real(cbp(ig+1))
         CYCLE 
C
C     CASE OF A MULTIPLE ZERO
C
  118    CONTINUE
         CBP(IG) = CMPLX(XM,0.)
         CBP(IG+1) = CMPLX(XM,0.)
c        bp(ig) = real(cbp(ig))
c        bp(ig+1) = real(cbp(ig+1))
         CYCLE 
C
C     CASE OF A COMPLEX ZERO
C
  119    CONTINUE
         IT = 0
         ICV = 0
         CX = CMPLX(XM,0.)
  120    CONTINUE
         FSG = (1.,0.)
         HSG = (1.,0.)
         FP = (0.,0.)
         FPP = (0.,0.)
         DO J = 1, IZ
            DD = 1./(CX - BH(J))
            FSG = FSG*A(J)*DD
            HSG = HSG*C(J)*DD
            FP = FP + DD
            FPP = FPP - DD*DD
         END DO
         IF (MODIZ == 0) THEN
            F = (1.,0.) - FSG - HSG
         ELSE
            F = (1.,0.) + FSG + HSG
         ENDIF
         I3 = 0
         IF (CABS(FP) > 0.) THEN
            I3 = 1
            R3 = -F/FP
         ENDIF
         I2 = 0
         IF (CABS(FPP) > 0.) THEN
            I2 = 1
            CDIS = CSQRT(FP**2 - 2.*F*FPP)
            R1 = CDIS - FP
            R2 = (-FP) - CDIS
            IF (CABS(R1) - CABS(R2) > 0.) THEN
               R1 = R1/FPP
            ELSE
               R1 = R2/FPP
            ENDIF
            R2 = 2.*F/FPP/R1
            IF (CABS(R2) < CABS(R1)) R1 = R2
            IF (I3 <= 0) GO TO 133
            IF (CABS(R3) < CABS(R1)) R1 = R3
            GO TO 133
         ENDIF
         R1 = R3
  133    CONTINUE
         CX = CX + R1
         IT = IT + 1
         IF (IT > 50) GO TO 142
         IF (CABS(R1) > SCNV) GO TO 120
         IF (ICV > 0) GO TO 135
         ICV = 1
         GO TO 120
  135    CONTINUE
         CBP(IG) = CX
         CBP(IG+1) = CONJG(CX)
      END DO
      IF (CABS(CBP(N)) - CABS(CBP(1)) <= 0.) THEN
         IF (CABS(CBP(N)) - CABS(CBP(1)) == 0.) GO TO 142
         NHALF = N/2
         DO J = 1, NHALF
            NT = N - J
            CX = CBP(J)
            CBP(J) = CBP(NT+1)
            CBP(NT+1) = CX
         END DO
      ENDIF
      NCMPLX = 1
      DO J = 2, IZ
         IF (AIMAG(CBP(J)) /= 0.) GO TO 143
      END DO
      NCMPLX = 0
      BP(1) = REAL(CBP(1))
      DO J = 2, IZ
         BP(J) = REAL(CBP(J))
      END DO
      GO TO 143
  142 CONTINUE
      IERROR = 4
  143 CONTINUE
      RETURN 
      END SUBROUTINE PPADD


      SUBROUTINE PROD(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,U)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BD(*)
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      REAL , INTENT(IN) :: X(*)
      REAL , INTENT(INOUT) :: Y(*)
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL , INTENT(INOUT) :: D(*)
      REAL , INTENT(INOUT) :: W(*)
      REAL  :: U(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, ID, IBR, M1, M2, IA, K
      REAL :: RT, DEN
C-----------------------------------------------
C
C PROD APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,W,U ARE WORKING ARRAYS
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      W(:M) = X(:M)
      Y(:M) = W(:M)
      MM = M - 1
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IF (ND == 0) RT = -RT
         IA = IA - 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*W(:M)
      ENDIF
      IF (ID <= 0) GO TO 125
      RT = BD(ID)
      ID = ID - 1
      IF (ID == 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      D(M) = A(M)/(B(M)-RT)
      W(M) = Y(M)/(B(M)-RT)
      DO J = 2, MM
         K = M - J
         DEN = B(K+1) - RT - C(K+1)*D(K+2)
         D(K+1) = A(K+1)/DEN
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
      END DO
      DEN = B(1) - RT - C(1)*D(2)
      W(1) = 1.
      IF (DEN /= 0.) THEN
         W(1) = (Y(1)-C(1)*W(2))/DEN
      ENDIF
      DO J = 2, M
         W(J) = W(J) - D(J)*W(J-1)
      END DO
      IF (NA > 0) GO TO 102
      GO TO 113
  111 CONTINUE
      Y(:M) = W(:M)
      IBR = 1
      GO TO 102
  113 CONTINUE
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 111
      ELSE
         IF (M2 > 0) THEN
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) <= 0.) GO TO 120
         ENDIF
         IF (IBR <= 0) THEN
            IF (ABS(BM1(M1)-BD(ID)) - ABS(BM1(M1)-RT) < 0.) GO TO 111
         ENDIF
         RT = RT - BM1(M1)
         M1 = M1 - 1
         GO TO 123
      ENDIF
  120 CONTINUE
      IF (IBR <= 0) THEN
         IF (ABS(BM2(M2)-BD(ID)) - ABS(BM2(M2)-RT) < 0.) GO TO 111
      ENDIF
      RT = RT - BM2(M2)
      M2 = M2 - 1
  123 CONTINUE
      Y(:M) = Y(:M) + RT*W(:M)
      GO TO 102
  125 CONTINUE
      RETURN 
      END SUBROUTINE PROD


      SUBROUTINE PRODP(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BD(*)
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      REAL , INTENT(IN) :: X(*)
      REAL , INTENT(INOUT) :: Y(*)
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL , INTENT(INOUT) :: D(*)
      REAL , INTENT(INOUT) :: U(*)
      REAL , INTENT(INOUT) :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, MM2, ID, IBR, M1, M2, IA, K
      REAL :: RT, BH, YM, DEN, V, AM
C-----------------------------------------------
C
C PRODP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y        PERIODIC BOUNDARY CONDITIONS
C
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,U,W ARE WORKING ARRAYS
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
 
      Y(:M) = X(:M)
      W(:M) = Y(:M)
      MM = M - 1
      MM2 = M - 2
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IF (ND == 0) RT = -RT
         IA = IA - 1
         Y(:M) = RT*W(:M)
      ENDIF
      IF (ID <= 0) GO TO 128
      RT = BD(ID)
      ID = ID - 1
      IF (ID == 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      BH = B(M) - RT
      YM = Y(M)
      DEN = B(1) - RT
      D(1) = C(1)/DEN
      U(1) = A(1)/DEN
      W(1) = Y(1)/DEN
      V = C(M)
      IF (MM2 - 2 >= 0) THEN
         DO J = 2, MM2
            DEN = B(J) - RT - A(J)*D(J-1)
            D(J) = C(J)/DEN
            U(J) = -A(J)*U(J-1)/DEN
            W(J) = (Y(J)-A(J)*W(J-1))/DEN
            BH = BH - V*U(J-1)
            YM = YM - V*W(J-1)
            V = -V*D(J-1)
         END DO
      ENDIF
      DEN = B(M-1) - RT - A(M-1)*D(M-2)
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
      W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN
      AM = A(M) - V*D(M-2)
      BH = BH - V*U(M-2)
      YM = YM - V*W(M-2)
      DEN = BH - AM*D(M-1)
      IF (DEN /= 0.) THEN
         W(M) = (YM - AM*W(M-1))/DEN
      ELSE
         W(M) = 1.
      ENDIF
      W(M-1) = W(M-1) - D(M-1)*W(M)
      DO J = 2, MM
         K = M - J
         W(K) = W(K) - D(K)*W(K+1) - U(K)*W(M)
      END DO
      IF (NA > 0) GO TO 102
      GO TO 116
  114 CONTINUE
      Y(:M) = W(:M)
      IBR = 1
      GO TO 102
  116 CONTINUE
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 114
      ELSE
         IF (M2 > 0) THEN
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) <= 0.) GO TO 123
         ENDIF
         IF (IBR <= 0) THEN
            IF (ABS(BM1(M1)-BD(ID)) - ABS(BM1(M1)-RT) < 0.) GO TO 114
         ENDIF
         RT = RT - BM1(M1)
         M1 = M1 - 1
         GO TO 126
      ENDIF
  123 CONTINUE
      IF (IBR <= 0) THEN
         IF (ABS(BM2(M2)-BD(ID)) - ABS(BM2(M2)-RT) < 0.) GO TO 114
      ENDIF
      RT = RT - BM2(M2)
      M2 = M2 - 1
  126 CONTINUE
      Y(:M) = Y(:M) + RT*W(:M)
      GO TO 102
  128 CONTINUE
      RETURN 
      END SUBROUTINE PRODP


      SUBROUTINE TEVLS(N, D, E2, IERR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: IERR
      REAL , INTENT(INOUT) :: D(N)
      REAL , INTENT(INOUT) :: E2(N)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CBLKT/
      COMMON /FISH_CBLKT/ NPP, K, MACHEP, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   MACHEP, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J, L, M, II, L1, MML, NHALF, NTOP
      REAL :: B, C, F, G, H, P, R, S, DHOLD
C-----------------------------------------------
C
C
C     REAL SQRT,ABS,SIGN
C
C
C     THIS SUBROUTINE IS A MODIFICATION OF THE EISPACK SUBROUTINE TQLRAT
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
C
      IERR = 0
      IF (N /= 1) THEN
C
         E2(:N-1) = E2(2:N)*E2(2:N)
C
         F = 0.0
         B = 0.0
         E2(N) = 0.0
C
         DO L = 1, N
            J = 0
            H = MACHEP*(ABS(D(L))+SQRT(E2(L)))
            IF (B <= H) THEN
               B = H
               C = B*B
            ENDIF
C
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
C
            DO M = L, N
               IF (E2(M) > C) CYCLE 
               EXIT 
C
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
C
            END DO
C
            IF (M /= L) THEN
  105          CONTINUE
               IF (J == 30) GO TO 114
               J = J + 1
C
C     ********** FORM SHIFT **********
C
               L1 = L + 1
               S = SQRT(E2(L))
               G = D(L)
               P = (D(L1)-G)/(2.0*S)
               R = SQRT(P*P + 1.0)
               D(L) = S/(P + SIGN(R,P))
               H = G - D(L)
C
               D(L1:N) = D(L1:N) - H
C
               F = F + H
C
C     ********** RATIONAL QL TRANSFORMATION **********
C
               G = D(M)
               IF (G == 0.0) G = B
               H = G
               S = 0.0
               MML = M - L
C
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
C
               DO II = 1, MML
                  I = M - II
                  P = G*H
                  R = P + E2(I)
                  E2(I+1) = S*R
                  S = E2(I)/R
                  D(I+1) = H + S*(H + D(I))
                  G = D(I) - E2(I)/G
                  IF (G == 0.0) G = B
                  H = G*P/R
               END DO
C
               E2(L) = S*G
               D(L) = H
C
C     ********** GUARD AGAINST UNDERFLOWED H **********
C
               IF (H == 0.0) GO TO 108
               IF (ABS(E2(L)) <= ABS(C/H)) GO TO 108
               E2(L) = H*E2(L)
               IF (E2(L) /= 0.0) GO TO 105
            ENDIF
  108       CONTINUE
            P = D(L) + F
C
C     ********** ORDER EIGENVALUES **********
C
            IF (L /= 1) THEN
C
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
C
               DO II = 2, L
                  I = L + 2 - II
                  IF (P >= D(I-1)) GO TO 111
                  D(I) = D(I-1)
               END DO
            ENDIF
C
            I = 1
  111       CONTINUE
            D(I) = P
         END DO
C
         IF (ABS(D(N)) >= ABS(D(1))) GO TO 115
         NHALF = N/2
         DO I = 1, NHALF
            NTOP = N - I
            DHOLD = D(I)
            D(I) = D(NTOP+1)
            D(NTOP+1) = DHOLD
         END DO
         GO TO 115
C
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
C
  114    CONTINUE
         IERR = L
      ENDIF
  115 CONTINUE
      RETURN 
C
C     ********** LAST CARD OF TQLRAT **********
C
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE TEVLS
C
C     file cblktri.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE CBLKTR (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,
C    +                   IERROR)
C
C                                                                       
C DIMENSION OF           AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N)
C ARGUMENTS
C                                                                       
C LATEST REVISION        JUNE 2004
C                                                                       
C PURPOSE                CBLKTR SOLVES A SYSTEM OF LINEAR EQUATIONS     
C                        OF THE FORM                                    
C                                                                       
C                        AN(J)*X(I,J-1) + AM(I)*X(I-1,J) +              
C                        (BN(J)+BM(I))*X(I,J) + CN(J)*X(I,J+1) +        
C                        CM(I)*X(I+1,J) = Y(I,J)                        
C                                                                       
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.         
C                                                                       
C                        I+1 AND I-1 ARE EVALUATED MODULO M AND         
C                        J+1 AND J-1 MODULO N, I.E.,                    
C                                                                       
C                        X(I,0) = X(I,N),  X(I,N+1) = X(I,1),           
C                        X(0,J) = X(M,J),  X(M+1,J) = X(1,J).           
C                                                                       
C                        THESE EQUATIONS USUALLY RESULT FROM THE        
C                        DISCRETIZATION OF SEPARABLE ELLIPTIC           
C                        EQUATIONS.  BOUNDARY CONDITIONS MAY BE         
C                        DIRICHLET, NEUMANN, OR PERIODIC.               
C                                                                       
C                        CBLKTRI IS A COMPLEX VERSION OF PACKAGE        
C                        BLKTRI ON ULIB.                                
C                                                                       
C USAGE                  CALL CBLKTR (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,    
C                                     CM,IDIMY,Y,IERROR,W)              
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT               IFLG                                           
C                                                                       
C                          = 0  INITIALIZATION ONLY.                    
C                               CERTAIN QUANTITIES THAT DEPEND ON NP,   
C                               N, AN, BN, AND CN ARE COMPUTED AND      
C                               STORED IN THE DERIVED DATA TYPE W
C                                                                       
C                          = 1  THE QUANTITIES THAT WERE COMPUTED       
C                               IN THE INITIALIZATION ARE USED          
C                               TO OBTAIN THE SOLUTION X(I,J).          
C                                                                       
C                               NOTE:                                   
C                               A CALL WITH IFLG=0 TAKES                
C                               APPROXIMATELY ONE HALF THE TIME         
C                               AS A CALL WITH IFLG = 1.                
C                               HOWEVER, THE INITIALIZATION DOES        
C                               NOT HAVE TO BE REPEATED UNLESS NP,      
C                               N, AN, BN, OR CN CHANGE.                
C                                                                       
C                        NP                                             
C                          = 0  IF AN(1) AND CN(N) ARE NOT ZERO,        
C                               WHICH CORRESPONDS TO PERIODIC           
C                               BOUNARY CONDITIONS.                     
C                                                                       
C                          = 1  IF AN(1) AND CN(N) ARE ZERO.            
C                                                                       
C                        N                                              
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.   
C                          N MUST BE GREATER THAN 4.                    
C                          THE OPERATION COUNT IS PROPORTIONAL TO       
C                          MNLOG2(N), HENCE N SHOULD BE SELECTED        
C                          LESS THAN OR EQUAL TO M.                     
C                                                                       
C                        AN,BN,CN                                       
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH N           
C                          THAT SPECIFY THE COEFFICIENTS IN THE         
C                          LINEAR EQUATIONS GIVEN ABOVE.                
C                                                                       
C                        MP                                             
C                          = 0  IF AM(1) AND CM(M) ARE NOT ZERO,        
C                               WHICH CORRESPONDS TO PERIODIC           
C                               BOUNDARY CONDITIONS.                    
C                                                                       
C                          = 1  IF AM(1) = CM(M) = 0  .                 
C                                                                       
C                        M                                              
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.   
C                           M MUST BE GREATER THAN 4.                   
C                                                                       
C                        AM,BM,CM                                       
C                          COMPLEX ONE-DIMENSIONAL ARRAYS OF LENGTH M   
C                          THAT SPECIFY THE COEFFICIENTS IN THE LINEAR  
C                          EQUATIONS GIVEN ABOVE.                       
C                                                                       
C                        IDIMY                                          
C                          THE ROW (OR FIRST) DIMENSION OF THE          
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS        
C                          IN THE PROGRAM CALLING CBLKTR.               
C                          THIS PARAMETER IS USED TO SPECIFY THE        
C                          VARIABLE DIMENSION OF Y.                     
C                          IDIMY MUST BE AT LEAST M.                    
C                                                                       
C                        Y                                              
C                          A COMPLEX TWO-DIMENSIONAL ARRAY THAT         
C                          SPECIFIES THE VALUES OF THE RIGHT SIDE OF    
C                          THE LINEAR SYSTEM OF EQUATIONS GIVEN ABOVE.  
C                          Y MUST BE DIMENSIONED Y(IDIMY,N) WITH        
C                          IDIMY .GE. M.                                
C                                                                       
C                        W
c                          A fortran 90 derived TYPE (fishworkspace) variable
c                          that must be declared by the user.  The first
c                          two declarative statements in the user program
c                          calling CBLKTRI must be:
c
c                               USE fish
c                               TYPE (fishworkspace) :: W
c
c                          The first statement makes the fishpack module
c                          defined in the file "fish.f" available to the
c                          user program calling CBLKTRI.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in CBLKTRI to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (IFLG=0) call to CBLKTRI.  These
c                          must be preserved on non-initial calls (IFLG=1)
c                          to CBLKTRI.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling CBLKTRI should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          CBLKTRI.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
c
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON OUTPUT              Y                                              
C                          CONTAINS THE SOLUTION X.                     
C                                                                       
C                        IERROR                                         
C                          AN ERROR FLAG THAT INDICATES INVALID         
C                          INPUT PARAMETERS.  EXCEPT FOR NUMBER ZER0,   
C                          A SOLUTION IS NOT ATTEMPTED.                 
C                                                                       
C                        = 0  NO ERROR.                                 
C                        = 1  M IS LESS THAN 5                          
C                        = 2  N IS LESS THAN 5                          
C                        = 3  IDIMY IS LESS THAN M.                     
C                        = 4  CBLKTR FAILED WHILE COMPUTING RESULTS     
C                             THAT DEPEND ON THE COEFFICIENT ARRAYS     
C                             AN, BN, CN.  CHECK THESE ARRAYS.          
C                        = 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.  
C                                                                       
C                             POSSIBLE REASONS FOR THIS CONDITION ARE   
C                             1. THE ARRAYS AN AND CN ARE NOT CORRECT   
C                             2. TOO LARGE A GRID SPACING WAS USED      
C                                IN THE DISCRETIZATION OF THE ELLIPTIC  
C                                EQUATION.                              
C                             3. THE LINEAR EQUATIONS RESULTED FROM A   
C                                PARTIAL DIFFERENTIAL EQUATION WHICH    
C                                WAS NOT ELLIPTIC.                      
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space in the derived type
C                               (fishworkspace) variable W fails (e.g.,
c                               if N,M are too large for the platform used)
c
C                                                                       
C                                                                       
C SPECIAL CONDITIONS     THE ALGORITHM MAY FAIL IF ABS(BM(I)+BN(J))     
C                        IS LESS THAN ABS(AM(I))+ABS(AN(J))+            
C                        ABS(CM(I))+ABS(CN(J))                          
C                        FOR SOME I AND J. THE ALGORITHM WILL ALSO      
C                        FAIL IF AN(J)*CN(J-1) IS LESS THAN ZERO FOR    
C                        SOME J.                                        
C                        SEE THE DESCRIPTION OF THE OUTPUT PARAMETER    
C                        IERROR.                                        
C                                                                       
C                                                                       
C I/O                    NONE                                           
C                                                                       
C PRECISION              SINGLE                                         
C                                                                       
C REQUIRED LIBRARY       comf.f,fish.f
C FILES
C                                                                       
C LANGUAGE               FORTRAN 90
C                                                                       
C HISTORY                WRITTEN BY PAUL SWARZTRAUBER AT NCAR IN        
C                        THE EARLY 1970'S.  REWRITTEN AN RELEASED       
C                        ON NCAR'S PUBLIC SOFTWARE LIBRARIES IN         
C                        JANUARY, 1980. Revised in June 2004 by John
C                        Adams using Fortan 90 dynamically allocated
c                        space and derived data types to eliminate mixed
c                        mode conflicts in the earlier versions.
C                                                                       
C ALGORITHM              GENERALIZED CYCLIC REDUCTION                   
C                        (SEE REFERENCE BELOW)                          
C                                                                       
C PORTABILITY
C                        THE APPROXIMATE MACHINE ACCURACY IS COMPUTED   
C                        IN FUNCTION EPMACH                             
C                                                                       
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, 'EFFICIENT       
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF        
C                        ELLIPTIC EQUATIONS'                            
C                        NCAR TN/IA-109, JULY, 1975, 138 PP.            
C                                                                       
C                        SWARZTRAUBER P. N.,A DIRECT METHOD FOR         
C                        THE DISCRETE SOLUTION OF SEPARABLE             
C                        ELLIPTIC EQUATIONS, S.I.A.M.                   
C                        J. NUMER. ANAL.,11(1974) PP. 1136-1150.        
C                                                                       
C***********************************************************************
      SUBROUTINE CBLKTRI(IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, 
     1   IDIMY, Y, IERROR, W)

      Implicit none
      TYPE (fishworkspace) :: w
!!!      EXTERNAL        PROC       ,PROCP      ,CPROC      ,CPROCP        
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IFLG
      INTEGER , INTENT(IN) :: NP
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: MP
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      COMPLEX  :: AM(*)
      COMPLEX  :: BM(*)
      COMPLEX  :: CM(*)
      COMPLEX  :: Y(IDIMY,*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::M2,NH,NL,IWAH,IW1,IWBH,IW2,IW3,IWD,IWW,IWU,IRWK,ICWK
C-----------------------------------------------
C
C TEST M AND N FOR THE PROPER FORM
C
      NM = N
      M2 = M + M
      IERROR = 0
      IF (M - 5 < 0) THEN
         IERROR = 1
      ELSE
         IF (NM - 3 < 0) THEN
            IERROR = 2
         ELSE
            IF (IDIMY - M < 0) THEN
               IERROR = 3
            ELSE
               NH = N
               NPP = NP
               IF (NPP /= 0) THEN
                  NH = NH + 1
               ENDIF
               IK = 2
               K = 1
               IK = IK + IK
               K = K + 1
               DO WHILE(NH - IK > 0)
                  IK = IK + IK
                  K = K + 1
               END DO
               NL = IK
               IK = IK + IK
               NL = NL - 1
               IWAH = (K - 2)*IK + K + 6
               IF (NPP /= 0) THEN
                  IW1 = IWAH
                  IWBH = IW1 + NM
               ELSE
                  IWBH = IWAH + NM + NM
                  IW1 = IWBH
                  NM = NM - 1
               ENDIF
C
C SUBROUTINE COMP B COMPUTES THE ROOTS OF THE B POLYNOMIALS
C
               IF (IERROR == 0) THEN
                  IW2 = IW1 + M
                  IW3 = IW2 + M
                  IWD = IW3 + M
                  IWW = IWD + M
                  IWU = IWW + M
                  IF (IFLG == 0) THEN
                     IRWK = IW1 + 2*N
                     ICWK = IW1 + 6*M
                     CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
                     IF (IERROR /= 0) RETURN 
!     COMPUTE b poly roots (real and complex)
      call ccompb(NL,ierror,an,bn,cn,w%rew,w%cxw,w%rew(iwah),
     +            w%rew(iwbh))
                  ELSE
                     IF (MP /= 0) THEN
      CALL CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,w%rew,w%cxw,
     +w%cxw(iw1),w%cxw(iw2),w%cxw(iw3),w%cxw(iwd),w%cxw(iww),
     +w%cxw(iwu),PROC,CPROC)
                     ELSE
      CALL CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,w%rew,w%cxw,
     +w%cxw(iw1),w%cxw(iw2),w%cxw(iw3),w%cxw(iwd),w%cxw(iww),
     +w%cxw(iwu),PROCP,CPROCP)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      RETURN 
      END SUBROUTINE CBLKTRI


      SUBROUTINE CBLKT1(N, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, B, BC, 
     1   W1, W2, W3, WD, WW, WU, PRDCT, CPRDCT)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      INTEGER  :: M
      INTEGER , INTENT(IN) :: IDIMY
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL  :: B(*)
      COMPLEX  :: AM(*)
      COMPLEX  :: BM(*)
      COMPLEX  :: CM(*)
      COMPLEX  :: Y(IDIMY,*)
      COMPLEX  :: BC(*)
      COMPLEX  :: W1(*)
      COMPLEX  :: W2(*)
      COMPLEX  :: W3(*)
      COMPLEX  :: WD(*)
      COMPLEX  :: WW(*)
      COMPLEX  :: WU(*)

      EXTERNAL :: PRDCT, CPRDCT
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: KDO, L, IR, I2, I1, I3, I4, IRM1, IM2, NM2, IM3, NM3, 
     1   IM1, NM1, IF, I, IPI1, IPI2, IPI3, IDXC, NC, IDXA, NA, IP2, NP2
     2   , IP1, NP1, IP3, NP3, J, IZ, NZ, IZR, LL, IFD, IP, NP, IMI1, 
     3   IMI2
      REAL :: DUM
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C CBLKT1 SOLVES THE LINEAR SYSTEM
C
C B  CONTAINS THE ROOTS OF ALL THE B POLYNOMIALS
C W1,W2,W3,WD,WW,WU  ARE ALL WORKING ARRAYS
C PRDCT IS EITHER PROCP OR PROC DEPENDING ON WHETHER THE BOUNDARY
C CONDITIONS IN THE M DIRECTION ARE PERIODIC OR NOT
C CPRDCT IS EITHER CPROCP OR CPROC WHICH ARE CALLED IF SOME OF THE ZEROS
C OF THE B POLYNOMIALS ARE COMPLEX.
C
C
C
C BEGIN REDUCTION PHASE
C
      KDO = K - 1
      DO L = 1, KDO
         IR = L - 1
         I2 = 2**IR
         I1 = I2/2
         I3 = I2 + I1
         I4 = I2 + I2
         IRM1 = IR - 1
         CALL CINDXB (I2, IR, IM2, NM2)
         CALL CINDXB (I1, IRM1, IM3, NM3)
         CALL CINDXB (I3, IRM1, IM1, NM1)
         CALL PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), 0, DUM, Y(1,
     1      I2), W3, M, AM, BM, CM, WD, WW, WU)
         IF = 2**K
         DO I = I4, IF, I4
            IF (I - NM > 0) CYCLE 
            IPI1 = I + I1
            IPI2 = I + I2
            IPI3 = I + I3
            CALL CINDXC (I, IR, IDXC, NC)
            IF (I - IF >= 0) CYCLE 
            CALL CINDXA (I, IR, IDXA, NA)
            CALL CINDXB (I - I1, IRM1, IM1, NM1)
            CALL CINDXB (IPI2, IR, IP2, NP2)
            CALL CINDXB (IPI1, IRM1, IP1, NP1)
            CALL CINDXB (IPI3, IRM1, IP3, NP3)
            CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W3, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            IF (IPI2 - NM > 0) THEN
               W3(:M) = (0.,0.)
               W2(:M) = (0.,0.)
            ELSE
               CALL PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM
     1            , Y(1,IPI2), W3, M, AM, BM, CM, WD, WW, WU)
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W3
     1            , W2, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            Y(:M,I) = W1(:M) + W2(:M) + Y(:M,I)
         END DO
      END DO
      IF (NPP == 0) THEN
         IF = 2**K
         I = IF/2
         I1 = I/2
         CALL CINDXB (I - I1, K - 2, IM1, NM1)
         CALL CINDXB (I + I1, K - 2, IP1, NP1)
         CALL CINDXB (I, K - 1, IZ, NZ)
         CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, Y(1,I)
     1      , W1, M, AM, BM, CM, WD, WW, WU)
         IZR = I
         W2(:M) = W1(:M)
         DO LL = 2, K
            L = K - LL + 1
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I = I2
            CALL CINDXC (I, IR, IDXC, NC)
            CALL CINDXB (I, IR, IZ, NZ)
            CALL CINDXB (I - I1, IR - 1, IM1, NM1)
            CALL CINDXB (I + I1, IR - 1, IP1, NP1)
            CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W1, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            W1(:M) = Y(:M,I) + W1(:M)
            CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1
     1         , W1, M, AM, BM, CM, WD, WW, WU)
         END DO
         L118: DO LL = 2, K
            L = K - LL + 1
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I4 = I2 + I2
            IFD = IF - I2
            DO I = I2, IFD, I4
               IF (I - I2 - IZR /= 0) CYCLE 
               IF (I - NM > 0) CYCLE  L118
               CALL CINDXA (I, IR, IDXA, NA)
               CALL CINDXB (I, IR, IZ, NZ)
               CALL CINDXB (I - I1, IR - 1, IM1, NM1)
               CALL CINDXB (I + I1, IR - 1, IP1, NP1)
               CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W2
     1            , W2, M, AM, BM, CM, WD, WW, WU)
               W2(:M) = Y(:M,I) + W2(:M)
               CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, 
     1            W2, W2, M, AM, BM, CM, WD, WW, WU)
               IZR = I
               IF (I - NM == 0) EXIT  L118
            END DO
         END DO L118
  119    CONTINUE
         Y(:M,NM+1) = Y(:M,NM+1) - CN(NM+1)*W1(:M) - AN(NM+1)*W2(:M)
         CALL CINDXB (IF/2, K - 1, IM1, NM1)
         CALL CINDXB (IF, K - 1, IP, NP)
         IF (NCMPLX /= 0) THEN
            CALL CPRDCT (NM + 1, BC(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(
     1         1,NM+1), Y(1,NM+1), M, AM, BM, CM, W1, W3, WW)
         ELSE
            CALL PRDCT (NM + 1, B(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(1,
     1         NM+1), Y(1,NM+1), M, AM, BM, CM, WD, WW, WU)
         ENDIF
         W1(:M) = AN(1)*Y(:M,NM+1)
         W2(:M) = CN(NM)*Y(:M,NM+1)
         Y(:M,1) = Y(:M,1) - W1(:M)
         Y(:M,NM) = Y(:M,NM) - W2(:M)
         DO L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I4 = I2 + I2
            I1 = I2/2
            I = I4
            CALL CINDXA (I, IR, IDXA, NA)
            CALL CINDXB (I - I2, IR, IM2, NM2)
            CALL CINDXB (I - I2 - I1, IR - 1, IM3, NM3)
            CALL CINDXB (I - I1, IR - 1, IM1, NM1)
            CALL PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), 0, DUM, 
     1         W1, W1, M, AM, BM, CM, WD, WW, WU)
            CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W1, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            Y(:M,I) = Y(:M,I) - W1(:M)
         END DO
C
         IZR = NM
         L131: DO L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I3 = I2 + I1
            I4 = I2 + I2
            IRM1 = IR - 1
            DO I = I4, IF, I4
               IPI1 = I + I1
               IPI2 = I + I2
               IPI3 = I + I3
               IF (IPI2 - IZR /= 0) THEN
                  IF (I - IZR /= 0) CYCLE 
                  CYCLE  L131
               ENDIF
               CALL CINDXC (I, IR, IDXC, NC)
               CALL CINDXB (IPI2, IR, IP2, NP2)
               CALL CINDXB (IPI1, IRM1, IP1, NP1)
               CALL CINDXB (IPI3, IRM1, IP3, NP3)
               CALL PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM
     1            , W2, W2, M, AM, BM, CM, WD, WW, WU)
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W2
     1            , W2, M, AM, BM, CM, WD, WW, WU)
               Y(:M,I) = Y(:M,I) - W2(:M)
               IZR = I
               CYCLE  L131
            END DO
         END DO L131
      ENDIF
C
C BEGIN BACK SUBSTITUTION PHASE
C
      DO LL = 1, K
         L = K - LL + 1
         IR = L - 1
         IRM1 = IR - 1
         I2 = 2**IR
         I1 = I2/2
         I4 = I2 + I2
         IFD = IF - I2
         DO I = I2, IFD, I4
            IF (I - NM > 0) CYCLE 
            IMI1 = I - I1
            IMI2 = I - I2
            IPI1 = I + I1
            IPI2 = I + I2
            CALL CINDXA (I, IR, IDXA, NA)
            CALL CINDXC (I, IR, IDXC, NC)
            CALL CINDXB (I, IR, IZ, NZ)
            CALL CINDXB (IMI1, IRM1, IM1, NM1)
            CALL CINDXB (IPI1, IRM1, IP1, NP1)
            IF (I - I2 <= 0) THEN
               W1(:M) = (0.,0.)
            ELSE
               CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), Y(
     1            1,IMI2), W1, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            IF (IPI2 - NM > 0) THEN
               W2(:M) = (0.,0.)
            ELSE
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), Y(
     1            1,IPI2), W2, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            W1(:M) = Y(:M,I) + W1(:M) + W2(:M)
            CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1
     1         , Y(1,I), M, AM, BM, CM, WD, WW, WU)
         END DO
      END DO
      RETURN 
      END SUBROUTINE CBLKT1


      REAL FUNCTION CBSRH (XLL, XRR, IZ, C, A, BH, F, SGN)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IZ
      REAL , INTENT(IN) :: XLL
      REAL , INTENT(IN) :: XRR
      REAL  :: F
      REAL , INTENT(IN) :: SGN
      REAL  :: C(*)
      REAL  :: A(*)
      REAL  :: BH(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: R1, XL, XR, DX, X
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
      XL = XLL
      XR = XRR
      DX = 0.5*ABS(XR - XL)
      X = 0.5*(XL + XR)
      R1 = SGN*F(X,IZ,C,A,BH)
      IF (R1 >= 0.) THEN
         IF (R1 == 0.) GO TO 105
         XR = X
      ELSE
         XL = X
      ENDIF
      DX = 0.5*DX
      DO WHILE(DX - CNV > 0.)
         X = 0.5*(XL + XR)
         R1 = SGN*F(X,IZ,C,A,BH)
         IF (R1 >= 0.) THEN
            IF (R1 == 0.) GO TO 105
            XR = X
         ELSE
            XL = X
         ENDIF
         DX = 0.5*DX
      END DO
  105 CONTINUE
      CBSRH = 0.5*(XL + XR)
      RETURN 
      END FUNCTION CBSRH


      SUBROUTINE CCOMPB(N, IERROR, AN, BN, CN, B, BC, AH, BH)
      IMPLICIT NONE
!!!      Real epmach
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      INTEGER  :: IERROR
      REAL  :: AN(*)
      REAL , INTENT(IN) :: BN(*)
      REAL  :: CN(*)
      REAL  :: B(*)
      REAL  :: AH(*)
      REAL  :: BH(*)
      COMPLEX  :: BC(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, IF, KDO, L, IR, I2, I4, IPL, IFD, I, IB, NB, JS, JF
     1   , LS, LH, NMP, L1, L2, J2, J1, N2M2
      REAL :: DUM, BNORM, ARG, D1, D2, D3
C-----------------------------------------------
C
C     CCOMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS USING SUBROUTINE
C     CTEVLS WHICH IS A MODIFICATION THE EISPACK PROGRAM TQLRAT.
C     IERROR IS SET TO 4 IF EITHER CTEVLS FAILS OR IF A(J+1)*C(J) IS
C     LESS THAN ZERO FOR SOME J.  AH,BH ARE TEMPORARY WORK ARRAYS.
C
 
      EPS = EPMACH(DUM)
      BNORM = ABS(BN(1))
      DO J = 2, NM
         BNORM = AMAX1(BNORM,ABS(BN(J)))
         ARG = AN(J)*CN(J-1)
         IF (ARG < 0.) GO TO 119
         B(J) = SIGN(SQRT(ARG),AN(J))
      END DO
      CNV = EPS*BNORM
      IF = 2**K
      KDO = K - 1
      L108: DO L = 1, KDO
         IR = L - 1
         I2 = 2**IR
         I4 = I2 + I2
         IPL = I4 - 1
         IFD = IF - I4
         DO I = I4, IFD, I4
            CALL CINDXB (I, L, IB, NB)
            IF (NB <= 0) CYCLE  L108
            JS = I - IPL
            JF = JS + NB - 1
            LS = 0
            BH(:JF-JS+1) = BN(JS:JF)
            AH(:JF-JS+1) = B(JS:JF)
            CALL CTEVLS (NB, BH, AH, IERROR)
            IF (IERROR /= 0) GO TO 118
            LH = IB - 1
            IF (NB > 0) THEN
               B(LH+1:NB+LH) = -BH(:NB)
               LH = NB + LH
            ENDIF
         END DO
      END DO L108
      B(:NM) = -BN(:NM)
      IF (NPP == 0) THEN
         NMP = NM + 1
         NB = NM + NMP
         DO J = 1, NB
            L1 = MOD(J - 1,NMP) + 1
            L2 = MOD(J + NM - 1,NMP) + 1
            ARG = AN(L1)*CN(L2)
            IF (ARG < 0.) GO TO 119
            BH(J) = SIGN(SQRT(ARG),(-AN(L1)))
            AH(J) = -BN(L1)
         END DO
         CALL CTEVLS (NB, AH, BH, IERROR)
         IF (IERROR /= 0) GO TO 118
         CALL CINDXB (IF, K - 1, J2, LH)
         CALL CINDXB (IF/2, K - 1, J1, LH)
         J2 = J2 + 1
         LH = J2
         N2M2 = J2 + NM + NM - 2
  114    CONTINUE
         D1 = ABS(B(J1)-B(J2-1))
         D2 = ABS(B(J1)-B(J2))
         D3 = ABS(B(J1)-B(J2+1))
         IF (D2>=D1 .OR. D2>=D3) THEN
            B(LH) = B(J2)
            J2 = J2 + 1
            LH = LH + 1
            IF (J2 - N2M2 <= 0) GO TO 114
         ELSE
            J2 = J2 + 1
            J1 = J1 + 1
            IF (J2 - N2M2 <= 0) GO TO 114
         ENDIF
         B(LH) = B(N2M2+1)
         CALL CINDXB (IF, K - 1, J1, J2)
         J2 = J1 + NMP + NMP
         CALL CPPADD (NM + 1, IERROR, AN, CN, BC(J1), B(J1), B(J2))
      ENDIF
      RETURN 
  118 CONTINUE
      IERROR = 4
      RETURN 
  119 CONTINUE
      IERROR = 5
      RETURN 
      END SUBROUTINE CCOMPB


      SUBROUTINE CPROC(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,YY)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      COMPLEX , INTENT(IN) :: BD(*)
      COMPLEX , INTENT(IN) :: X(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W(*)
      COMPLEX  :: YY(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, ID, M1, M2, IA, IFLG, K
      REAL :: RT
      COMPLEX :: CRT, DEN, Y1, Y2
C-----------------------------------------------
C
C PROC APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,W ARE WORK ARRAYS
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      Y(:M) = X(:M)
      MM = M - 1
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IFLG = 0
      IF (ID > 0) THEN
         CRT = BD(ID)
         ID = ID - 1
C
C BEGIN SOLUTION TO SYSTEM
C
         D(M) = A(M)/(B(M)-CRT)
         W(M) = Y(M)/(B(M)-CRT)
         DO J = 2, MM
            K = M - J
            DEN = B(K+1) - CRT - C(K+1)*D(K+2)
            D(K+1) = A(K+1)/DEN
            W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
         END DO
         DEN = B(1) - CRT - C(1)*D(2)
         IF (CABS(DEN) /= 0.) THEN
            Y(1) = (Y(1)-C(1)*W(2))/DEN
         ELSE
            Y(1) = (1.,0.)
         ENDIF
         DO J = 2, M
            Y(J) = W(J) - D(J)*Y(J-1)
         END DO
      ENDIF
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 121
         RT = BM2(M2)
         M2 = M2 - 1
      ELSE
         IF (M2 <= 0) THEN
            RT = BM1(M1)
            M1 = M1 - 1
         ELSE
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) > 0.) THEN
               RT = BM1(M1)
               M1 = M1 - 1
            ELSE
               RT = BM2(M2)
               M2 = M2 - 1
            ENDIF
         ENDIF
      ENDIF
      Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2)
      IF (MM - 2 >= 0) THEN
         DO J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
         END DO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M)
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  121 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IA = IA - 1
         IFLG = 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*Y(:M)
      ENDIF
      IF (IFLG > 0) GO TO 102
      RETURN 
      END SUBROUTINE CPROC


      SUBROUTINE CPROCP(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,YY)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      REAL  :: YY(*)
      COMPLEX , INTENT(IN) :: BD(*)
      COMPLEX , INTENT(IN) :: X(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: U(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, MM2, ID, M1, M2, IA, IFLG, K
      REAL :: RT
      COMPLEX :: V, DEN, BH, YM, AM, Y1, Y2, YH, CRT
C-----------------------------------------------
C
C CPROCP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y
C
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,U ARE WORK ARRAYS
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      Y(:M) = X(:M)
      MM = M - 1
      MM2 = M - 2
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IFLG = 0
      IF (ID > 0) THEN
         CRT = BD(ID)
         ID = ID - 1
         IFLG = 1
C
C BEGIN SOLUTION TO SYSTEM
C
         BH = B(M) - CRT
         YM = Y(M)
         DEN = B(1) - CRT
         D(1) = C(1)/DEN
         U(1) = A(1)/DEN
         Y(1) = Y(1)/DEN
         V = C(M)
         IF (MM2 - 2 >= 0) THEN
            DO J = 2, MM2
               DEN = B(J) - CRT - A(J)*D(J-1)
               D(J) = C(J)/DEN
               U(J) = -A(J)*U(J-1)/DEN
               Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
               BH = BH - V*U(J-1)
               YM = YM - V*Y(J-1)
               V = -V*D(J-1)
            END DO
         ENDIF
         DEN = B(M-1) - CRT - A(M-1)*D(M-2)
         D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
         Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN
         AM = A(M) - V*D(M-2)
         BH = BH - V*U(M-2)
         YM = YM - V*Y(M-2)
         DEN = BH - AM*D(M-1)
         IF (CABS(DEN) /= 0.) THEN
            Y(M) = (YM - AM*Y(M-1))/DEN
         ELSE
            Y(M) = (1.,0.)
         ENDIF
         Y(M-1) = Y(M-1) - D(M-1)*Y(M)
         DO J = 2, MM
            K = M - J
            Y(K) = Y(K) - D(K)*Y(K+1) - U(K)*Y(M)
         END DO
      ENDIF
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 123
         RT = BM2(M2)
         M2 = M2 - 1
      ELSE
         IF (M2 <= 0) THEN
            RT = BM1(M1)
            M1 = M1 - 1
         ELSE
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) > 0.) THEN
               RT = BM1(M1)
               M1 = M1 - 1
            ELSE
               RT = BM2(M2)
               M2 = M2 - 1
C
C MATRIX MULTIPLICATION
C
            ENDIF
         ENDIF
      ENDIF
      YH = Y(1)
      Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2) + A(1)*Y(M)
      IF (MM - 2 >= 0) THEN
         DO J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
         END DO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M) + C(M)*YH
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  123 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IA = IA - 1
         IFLG = 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*Y(:M)
      ENDIF
      IF (IFLG > 0) GO TO 102
      RETURN 
      END SUBROUTINE CPROCP


      SUBROUTINE CINDXA(I, IR, IDXA, NA)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDXA
      INTEGER , INTENT(OUT) :: NA
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
      NA = 2**IR
      IDXA = I - NA + 1
      IF (I - NM > 0) THEN
         NA = 0
      ENDIF
      RETURN 
      END SUBROUTINE CINDXA


      SUBROUTINE CINDXB(I, IR, IDX, IDP)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDX
      INTEGER , INTENT(OUT) :: IDP
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IZH, ID, IPL
C-----------------------------------------------
C
C B(IDX) IS THE LOCATION OF THE FIRST ROOT OF THE B(I,IR) POLYNOMIAL
C
      IDP = 0
      IF (IR >= 0) THEN
         IF (IR <= 0) THEN
            IF (I - NM > 0) GO TO 107
            IDX = I
            IDP = 1
            RETURN 
         ENDIF
         IZH = 2**IR
         ID = I - IZH - IZH
         IDX = ID + ID + (IR - 1)*IK + IR + (IK - I)/IZH + 4
         IPL = IZH - 1
         IDP = IZH + IZH - 1
         IF (I - IPL - NM > 0) THEN
            IDP = 0
            RETURN 
         ENDIF
         IF (I + IPL - NM > 0) THEN
            IDP = NM + IPL - I + 1
         ENDIF
      ENDIF
  107 CONTINUE
      RETURN 
      END SUBROUTINE CINDXB


      SUBROUTINE CINDXC(I, IR, IDXC, NC)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDXC
      INTEGER , INTENT(OUT) :: NC
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
      NC = 2**IR
      IDXC = I
      IF (IDXC + NC - 1 - NM > 0) THEN
         NC = 0
      ENDIF
      RETURN 
      END SUBROUTINE CINDXC


      SUBROUTINE CPPADD(N, IERROR, A, C, CBP, BP, BH)
      IMPLICIT NONE
!!!      real psgf,ppspf,ppsgf,cbsrh
!!!      EXTERNAL        PSGF       ,PPSPF      ,PPSGF                     
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: IERROR
      REAL  :: A(*)
      REAL  :: C(*)
      REAL , INTENT(INOUT) :: BP(*)
      REAL  :: BH(*)
      COMPLEX , INTENT(INOUT) :: CBP(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::IZ,IZM,IZM2,J,NT,MODIZ,IS,IF,IG,IT,ICV,I3,I2,NHALF
      REAL :: R4, R5, R6, SCNV, XL, DB, SGN, XR, XM, PSG
      COMPLEX :: CF, CX, FSG, HSG, DD, F, FP, FPP, CDIS, R1, R2, R3
C-----------------------------------------------
C
C     CPPADD COMPUTES THE EIGENVALUES OF THE PERIODIC TRIDIAGONAL
C     MATRIX WITH COEFFICIENTS AN,BN,CN
C
C     N IS THE ORDER OF THE BH AND BP POLYNOMIALS
C     ON OUTPUT BP CONTAINS THE EIGENVALUES
C     CBP IS THE SAME AS BP EXCEPT TYPE COMPLEX
C     BH IS USED TO TEMPORARILY STORE THE ROOTS OF THE B HAT POLYNOMIAL
C       WHICH ENTERS THROUGH BP
C
      SCNV = SQRT(CNV)
      IZ = N
      IZM = IZ - 1
      IZM2 = IZ - 2
      IF (BP(N) - BP(1) <= 0.) THEN
         IF (BP(N) - BP(1) == 0.) GO TO 142
         BH(:N) = BP(N:1:(-1))
      ELSE
         BH(:N) = BP(:N)
      ENDIF
      NCMPLX = 0
      MODIZ = MOD(IZ,2)
      IS = 1
      IF (MODIZ /= 0) THEN
         IF (A(1) < 0.) GO TO 110
         IF (A(1) == 0.) GO TO 142
      ENDIF
      XL = BH(1)
      DB = BH(3) - BH(1)
      XL = XL - DB
      R4 = PSGF(XL,IZ,C,A,BH)
      DO WHILE(R4 <= 0.)
         XL = XL - DB
         R4 = PSGF(XL,IZ,C,A,BH)
      END DO
      SGN = -1.
      CBP(1) = CMPLX(CBSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.)
      IS = 2
  110 CONTINUE
      IF = IZ - 1
      IF (MODIZ /= 0) THEN
         IF (A(1) > 0.) GO TO 115
         IF (A(1) == 0.) GO TO 142
      ENDIF
      XR = BH(IZ)
      DB = BH(IZ) - BH(IZ-2)
      XR = XR + DB
      R5 = PSGF(XR,IZ,C,A,BH)
      DO WHILE(R5 < 0.)
         XR = XR + DB
         R5 = PSGF(XR,IZ,C,A,BH)
      END DO
      SGN = 1.
      CBP(IZ) = CMPLX(CBSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.)
      IF = IZ - 2
  115 CONTINUE
      DO IG = IS, IF, 2
         XL = BH(IG)
         XR = BH(IG+1)
         SGN = -1.
         XM = CBSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)
         PSG = PSGF(XM,IZ,C,A,BH)
         IF (ABS(PSG) - EPS <= 0.) GO TO 118
         R6 = PSG*PPSGF(XM,IZ,C,A,BH)
         IF (R6 > 0.) GO TO 119
         IF (R6 == 0.) GO TO 118
         SGN = 1.
         CBP(IG) = CMPLX(CBSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.)
         SGN = -1.
         CBP(IG+1) = CMPLX(CBSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.)
         CYCLE 
C
C     CASE OF A MULTIPLE ZERO
C
  118    CONTINUE
         CBP(IG) = CMPLX(XM,0.)
         CBP(IG+1) = CMPLX(XM,0.)
         CYCLE 
C
C     CASE OF A COMPLEX ZERO
C
  119    CONTINUE
         IT = 0
         ICV = 0
         CX = CMPLX(XM,0.)
  120    CONTINUE
         FSG = (1.,0.)
         HSG = (1.,0.)
         FP = (0.,0.)
         FPP = (0.,0.)
         DO J = 1, IZ
            DD = 1./(CX - BH(J))
            FSG = FSG*A(J)*DD
            HSG = HSG*C(J)*DD
            FP = FP + DD
            FPP = FPP - DD*DD
         END DO
         IF (MODIZ == 0) THEN
            F = (1.,0.) - FSG - HSG
         ELSE
            F = (1.,0.) + FSG + HSG
         ENDIF
         I3 = 0
         IF (CABS(FP) > 0.) THEN
            I3 = 1
            R3 = -F/FP
         ENDIF
         I2 = 0
         IF (CABS(FPP) > 0.) THEN
            I2 = 1
            CDIS = CSQRT(FP**2 - 2.*F*FPP)
            R1 = CDIS - FP
            R2 = (-FP) - CDIS
            IF (CABS(R1) - CABS(R2) > 0.) THEN
               R1 = R1/FPP
            ELSE
               R1 = R2/FPP
            ENDIF
            R2 = 2.*F/FPP/R1
            IF (CABS(R2) < CABS(R1)) R1 = R2
            IF (I3 <= 0) GO TO 133
            IF (CABS(R3) < CABS(R1)) R1 = R3
            GO TO 133
         ENDIF
         R1 = R3
  133    CONTINUE
         CX = CX + R1
         IT = IT + 1
         IF (IT > 50) GO TO 142
         IF (CABS(R1) > SCNV) GO TO 120
         IF (ICV > 0) GO TO 135
         ICV = 1
         GO TO 120
  135    CONTINUE
         CBP(IG) = CX
         CBP(IG+1) = CONJG(CX)
      END DO
      IF (CABS(CBP(N)) - CABS(CBP(1)) <= 0.) THEN
         IF (CABS(CBP(N)) - CABS(CBP(1)) == 0.) GO TO 142
         NHALF = N/2
         DO J = 1, NHALF
            NT = N - J
            CX = CBP(J)
            CBP(J) = CBP(NT+1)
            CBP(NT+1) = CX
         END DO
      ENDIF
      NCMPLX = 1
      DO J = 2, IZ
         IF (AIMAG(CBP(J)) /= 0.) GO TO 143
      END DO
      NCMPLX = 0
      DO J = 2, IZ
         BP(J) = REAL(CBP(J))
      END DO
      GO TO 143
  142 CONTINUE
      IERROR = 4
  143 CONTINUE
      RETURN 
      END SUBROUTINE CPPADD


      SUBROUTINE PROC(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,U)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BD(*)
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      COMPLEX , INTENT(IN) :: X(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W(*)
      COMPLEX  :: U(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, ID, IBR, M1, M2, IA, K
      REAL :: RT
      COMPLEX :: DEN
C-----------------------------------------------
C
C PROC APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,W,U ARE WORKING ARRAYS
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      W(:M) = X(:M)
      Y(:M) = W(:M)
      MM = M - 1
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IF (ND == 0) RT = -RT
         IA = IA - 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*W(:M)
      ENDIF
      IF (ID <= 0) GO TO 125
      RT = BD(ID)
      ID = ID - 1
      IF (ID == 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      D(M) = A(M)/(B(M)-RT)
      W(M) = Y(M)/(B(M)-RT)
      DO J = 2, MM
         K = M - J
         DEN = B(K+1) - RT - C(K+1)*D(K+2)
         D(K+1) = A(K+1)/DEN
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
      END DO
      DEN = B(1) - RT - C(1)*D(2)
      W(1) = (1.,0.)
      IF (CABS(DEN) /= 0.) THEN
         W(1) = (Y(1)-C(1)*W(2))/DEN
      ENDIF
      DO J = 2, M
         W(J) = W(J) - D(J)*W(J-1)
      END DO
      IF (NA > 0) GO TO 102
      GO TO 113
  111 CONTINUE
      Y(:M) = W(:M)
      IBR = 1
      GO TO 102
  113 CONTINUE
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 111
      ELSE
         IF (M2 > 0) THEN
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) <= 0.) GO TO 120
         ENDIF
         IF (IBR <= 0) THEN
            IF (ABS(BM1(M1)-BD(ID)) - ABS(BM1(M1)-RT) < 0.) GO TO 111
         ENDIF
         RT = RT - BM1(M1)
         M1 = M1 - 1
         GO TO 123
      ENDIF
  120 CONTINUE
      IF (IBR <= 0) THEN
         IF (ABS(BM2(M2)-BD(ID)) - ABS(BM2(M2)-RT) < 0.) GO TO 111
      ENDIF
      RT = RT - BM2(M2)
      M2 = M2 - 1
  123 CONTINUE
      Y(:M) = Y(:M) + RT*W(:M)
      GO TO 102
  125 CONTINUE
      RETURN 
      END SUBROUTINE PROC


      SUBROUTINE PROCP(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,W)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BD(*)
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      COMPLEX , INTENT(IN) :: X(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: U(*)
      COMPLEX , INTENT(INOUT) :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, MM2, ID, IBR, M1, M2, IA, K
      REAL :: RT
      COMPLEX :: DEN, YM, V, BH, AM
C-----------------------------------------------
C
C PROCP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y        PERIODIC BOUNDARY CONDITIONS
C
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,U,W ARE WORKING ARRAYS
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      Y(:M) = X(:M)
      W(:M) = Y(:M)
      MM = M - 1
      MM2 = M - 2
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IF (ND == 0) RT = -RT
         IA = IA - 1
         Y(:M) = RT*W(:M)
      ENDIF
      IF (ID <= 0) GO TO 128
      RT = BD(ID)
      ID = ID - 1
      IF (ID == 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      BH = B(M) - RT
      YM = Y(M)
      DEN = B(1) - RT
      D(1) = C(1)/DEN
      U(1) = A(1)/DEN
      W(1) = Y(1)/DEN
      V = C(M)
      IF (MM2 - 2 >= 0) THEN
         DO J = 2, MM2
            DEN = B(J) - RT - A(J)*D(J-1)
            D(J) = C(J)/DEN
            U(J) = -A(J)*U(J-1)/DEN
            W(J) = (Y(J)-A(J)*W(J-1))/DEN
            BH = BH - V*U(J-1)
            YM = YM - V*W(J-1)
            V = -V*D(J-1)
         END DO
      ENDIF
      DEN = B(M-1) - RT - A(M-1)*D(M-2)
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
      W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN
      AM = A(M) - V*D(M-2)
      BH = BH - V*U(M-2)
      YM = YM - V*W(M-2)
      DEN = BH - AM*D(M-1)
      IF (CABS(DEN) /= 0.) THEN
         W(M) = (YM - AM*W(M-1))/DEN
      ELSE
         W(M) = (1.,0.)
      ENDIF
      W(M-1) = W(M-1) - D(M-1)*W(M)
      DO J = 2, MM
         K = M - J
         W(K) = W(K) - D(K)*W(K+1) - U(K)*W(M)
      END DO
      IF (NA > 0) GO TO 102
      GO TO 116
  114 CONTINUE
      Y(:M) = W(:M)
      IBR = 1
      GO TO 102
  116 CONTINUE
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 114
      ELSE
         IF (M2 > 0) THEN
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) <= 0.) GO TO 123
         ENDIF
         IF (IBR <= 0) THEN
            IF (ABS(BM1(M1)-BD(ID)) - ABS(BM1(M1)-RT) < 0.) GO TO 114
         ENDIF
         RT = RT - BM1(M1)
         M1 = M1 - 1
         GO TO 126
      ENDIF
  123 CONTINUE
      IF (IBR <= 0) THEN
         IF (ABS(BM2(M2)-BD(ID)) - ABS(BM2(M2)-RT) < 0.) GO TO 114
      ENDIF
      RT = RT - BM2(M2)
      M2 = M2 - 1
  126 CONTINUE
      Y(:M) = Y(:M) + RT*W(:M)
      GO TO 102
  128 CONTINUE
      RETURN 
      END SUBROUTINE PROCP


      SUBROUTINE CTEVLS(N, D, E2, IERR)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: IERR
      REAL , INTENT(INOUT) :: D(N)
      REAL , INTENT(INOUT) :: E2(N)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /FISH_CCBLK/ NPP, K, MACHEP, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   MACHEP, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J, L, M, II, L1, MML, NHALF, NTOP
      REAL :: B, C, F, G, H, P, R, S, DHOLD
C-----------------------------------------------
C
C
C     REAL SQRT,ABS,SIGN
C
C
C     THIS SUBROUTINE IS A MODIFICATION OF THE EISPACK SUBROUTINE TQLRAT
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
C
      IERR = 0
      IF (N /= 1) THEN
C
         E2(:N-1) = E2(2:N)*E2(2:N)
C
         F = 0.0
         B = 0.0
         E2(N) = 0.0
C
         DO L = 1, N
            J = 0
            H = MACHEP*(ABS(D(L))+SQRT(E2(L)))
            IF (B <= H) THEN
               B = H
               C = B*B
            ENDIF
C
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
C
            DO M = L, N
               IF (E2(M) > C) CYCLE 
               EXIT 
C
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
C
            END DO
C
            IF (M /= L) THEN
  105          CONTINUE
               IF (J == 30) GO TO 114
               J = J + 1
C
C     ********** FORM SHIFT **********
C
               L1 = L + 1
               S = SQRT(E2(L))
               G = D(L)
               P = (D(L1)-G)/(2.0*S)
               R = SQRT(P*P + 1.0)
               D(L) = S/(P + SIGN(R,P))
               H = G - D(L)
C
               D(L1:N) = D(L1:N) - H
C
               F = F + H
C
C     ********** RATIONAL QL TRANSFORMATION **********
C
               G = D(M)
               IF (G == 0.0) G = B
               H = G
               S = 0.0
               MML = M - L
C
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
C
               DO II = 1, MML
                  I = M - II
                  P = G*H
                  R = P + E2(I)
                  E2(I+1) = S*R
                  S = E2(I)/R
                  D(I+1) = H + S*(H + D(I))
                  G = D(I) - E2(I)/G
                  IF (G == 0.0) G = B
                  H = G*P/R
               END DO
C
               E2(L) = S*G
               D(L) = H
C
C     ********** GUARD AGAINST UNDERFLOWED H **********
C
               IF (H == 0.0) GO TO 108
               IF (ABS(E2(L)) <= ABS(C/H)) GO TO 108
               E2(L) = H*E2(L)
               IF (E2(L) /= 0.0) GO TO 105
            ENDIF
  108       CONTINUE
            P = D(L) + F
C
C     ********** ORDER EIGENVALUES **********
C
            IF (L /= 1) THEN
C
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
C
               DO II = 2, L
                  I = L + 2 - II
                  IF (P >= D(I-1)) GO TO 111
                  D(I) = D(I-1)
               END DO
            ENDIF
C
            I = 1
  111       CONTINUE
            D(I) = P
         END DO
C
         IF (ABS(D(N)) >= ABS(D(1))) GO TO 115
         NHALF = N/2
         DO I = 1, NHALF
            NTOP = N - I
            DHOLD = D(I)
            D(I) = D(NTOP+1)
            D(NTOP+1) = DHOLD
         END DO
         GO TO 115
C
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
C
  114    CONTINUE
         IERR = L
      ENDIF
  115 CONTINUE
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE CTEVLS
C
C     file cmgnbn.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE CMGNBN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR)
C
C
C DIMENSION OF           A(M),B(M),C(M),Y(IDIMY,N)
C ARGUMENTS
C
C LATEST REVISION        NOVEMBER 2004
C
C PURPOSE                THE NAME OF THIS PACKAGE IS A MNEMONIC FOR THE
C                        COMPLEX GENERALIZED BUNEMAN ALGORITHM.
C                        IT SOLVES THE COMPLEX LINEAR SYSTEM OF EQUATION
C
C                        A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C                        + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
C
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.
C
C                        INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
C                        I.E., X(0,J) = X(M,J) AND X(M+1,J) = X(1,J),
C                        AND X(I,0) MAY EQUAL 0, X(I,2), OR X(I,N),
C                        AND X(I,N+1) MAY EQUAL 0, X(I,N-1), OR X(I,1)
C                        DEPENDING ON AN INPUT PARAMETER.
C
C USAGE                  CALL CMGNBN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
C                                     IERROR)
C
C ARGUMENTS
C
C ON INPUT               NPEROD
C
C                          INDICATES THE VALUES THAT X(I,0) AND
C                          X(I,N+1) ARE ASSUMED TO HAVE.
C
C                          = 0  IF X(I,0) = X(I,N) AND X(I,N+1) =
C                               X(I,1).
C                          = 1  IF X(I,0) = X(I,N+1) = 0  .
C                          = 2  IF X(I,0) = 0 AND X(I,N+1) = X(I,N-1).
C                          = 3  IF X(I,0) = X(I,2) AND X(I,N+1) =
C                               X(I,N-1).
C                          = 4  IF X(I,0) = X(I,2) AND X(I,N+1) = 0.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        MPEROD
C                          = 0 IF A(1) AND C(M) ARE NOT ZERO
C                          = 1 IF A(1) = C(M) = 0
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        A,B,C
C                          ONE-DIMENSIONAL COMPLEX ARRAYS OF LENGTH M
C                          THAT SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0
C                          THE ARRAY ELEMENTS MUST NOT DEPEND UPON
C                          THE INDEX I, BUT MUST BE CONSTANT.
C                          SPECIFICALLY, THE SUBROUTINE CHECKS THE
C                          FOLLOWING CONDITION .
C
C                            A(I) = C(1)
C                            C(I) = C(1)
C                            B(I) = B(1)
C
C                          FOR I=1,2,...,M.
C
C                        IDIMY
C                          THE ROW (OR FIRST) DIMENSION OF THE
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS
C                          IN THE PROGRAM CALLING CMGNBN.
C                          THIS PARAMETER IS USED TO SPECIFY THE
C                          VARIABLE DIMENSION OF Y.
C                          IDIMY MUST BE AT LEAST M.
C
C                        Y
C                          A TWO-DIMENSIONAL COMPLEX ARRAY THAT
C                          SPECIFIES THE VALUES OF THE RIGHT SIDE
C                          OF THE LINEAR SYSTEM OF EQUATIONS GIVEN
C                          ABOVE.
C                          Y MUST BE DIMENSIONED AT LEAST M*N.
C
C
C  ON OUTPUT             Y
C
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG WHICH INDICATES INVALID
C                          INPUT PARAMETERS  EXCEPT FOR NUMBER
C                          ZERO, A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR.
C                          = 1  M .LE. 2  .
C                          = 2  N .LE. 2
C                          = 3  IDIMY .LT. M
C                          = 4  NPEROD .LT. 0 OR NPEROD .GT. 4
C                          = 5  MPEROD .LT. 0 OR MPEROD .GT. 1
C                          = 6  A(I) .NE. C(1) OR C(I) .NE. C(1) OR
C                               B(I) .NE. B(1) FOR
C                               SOME I=1,2,...,M.
C                          = 7  A(1) .NE. 0 OR C(M) .NE. 0 AND
C                                 MPEROD = 1
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C SPECIAL CONDITONS      NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       comf.f,fish.f
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
C                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C ALGORITHM              THE LINEAR SYSTEM IS SOLVED BY A CYCLIC
C                        REDUCTION ALGORITHM DESCRIBED IN THE
C                        REFERENCE BELOW.
C
C PORTABILITY            FORTRAN 90.  ALL MACHINE DEPENDENT CONSTANTS
C                        ARE DEFINED IN FUNCTION P1MACH.
C
C REFERENCES             SWEET, R., 'A CYCLIC REDUCTION ALGORITHM FOR
C                        SOLVING BLOCK TRIDIAGONAL SYSTEMS OF ARBITRARY
C                        DIMENSIONS,' SIAM J. ON NUMER. ANAL.,
C                          14(SEPT., 1977), PP. 706-720.
C
C ACCURACY               THIS TEST WAS PERFORMED ON A Platform with
c                        64 bit floating point arithmetic.
C                        A UNIFORM RANDOM NUMBER GENERATOR WAS USED
C                        TO CREATE A SOLUTION ARRAY X FOR THE SYSTEM
C                        GIVEN IN THE 'PURPOSE' DESCRIPTION ABOVE
C                        WITH
C                          A(I) = C(I) = -0.5*B(I) = 1, I=1,2,...,M
C
C                        AND, WHEN MPEROD = 1
C
C                          A(1) = C(M) = 0
C                          A(M) = C(1) = 2.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE
C                        GIVEN SYSTEM  AND A RIGHT SIDE Y WAS
C                        COMPUTED.  USING THIS ARRAY Y, SUBROUTINE
C                        CMGNBN WAS CALLED TO PRODUCE APPROXIMATE
C                        SOLUTION Z.  THEN RELATIVE ERROR
C                          E = MAX(CABS(Z(I,J)-X(I,J)))/
C                              MAX(CABS(X(I,J)))
C                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
C                        OVER I=1,2,...,M AND J=1,...,N.
C
C                        THE VALUE OF E IS GIVEN IN THE TABLE
C                        BELOW FOR SOME TYPICAL VALUES OF M AND N.
C
C                   M (=N)    MPEROD    NPEROD       E
C                   ------    ------    ------     ------
C
C                     31        0         0        1.E-12
C                     31        1         1        4.E-13
C                     31        1         3        2.E-12
C                     32        0         0        7.E-14
C                     32        1         1        5.E-13
C                     32        1         3        2.E-13
C                     33        0         0        6.E-13
C                     33        1         1        5.E-13
C                     33        1         3        3.E-12
C                     63        0         0        5.E-12
C                     63        1         1        6.E-13
C                     63        1         3        1.E-11
C                     64        0         0        1.E-13
C                     64        1         1        3.E-12
C                     64        1         3        3.E-13
C                     65        0         0        2.E-12
C                     65        1         1        5.E-13
C                     65        1         3        1.E-11
C
C***********************************************************************
      SUBROUTINE CMGNBN(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, IERROR)

      implicit none
      TYPE(fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      COMPLEX  :: A(*)
      COMPLEX  :: B(*)
      COMPLEX  :: C(*)
      COMPLEX  :: Y(IDIMY,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, ICWK, IRWK
      COMPLEX :: A1
C-----------------------------------------------
      IERROR = 0
      IF (M <= 2) IERROR = 1
      IF (N <= 2) IERROR = 2
      IF (IDIMY < M) IERROR = 3
      IF (NPEROD<0 .OR. NPEROD>4) IERROR = 4
      IF (MPEROD<0 .OR. MPEROD>1) IERROR = 5
      IF (MPEROD /= 1) THEN
         DO I = 2, M
            IF (CABS(A(I)-C(1)) /= 0.) GO TO 103
            IF (CABS(C(I)-C(1)) /= 0.) GO TO 103
            IF (CABS(B(I)-B(1)) /= 0.) GO TO 103
         END DO
         GO TO 104
      ENDIF
      IF (CABS(A(1))/=0. .AND. CABS(C(M))/=0.) IERROR = 7
      GO TO 104
  103 CONTINUE
      IERROR = 6
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
!     allocate required complex work space
      ICWK = (10 + INT(ALOG(FLOAT(N))/ALOG(2.0)))*M + 4*N
      IRWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed
      IF (IERROR == 20) RETURN 
      call cmgnbnn(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,w%cxw)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE CMGNBN


 
      SUBROUTINE CMGNBNN(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX  :: Y(IDIMY,*)
      COMPLEX  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IWBA, IWBB, IWBC, IWB2, IWB3, IWW1, IWW2, IWW3, IWD, 
     1   IWTCOS, IWP, I, K, J, MP, NP, IPSTOR, IREV, MH, MHM1, MODD, 
     2   MHPI, MHMI, NBY2, MSKIP
      COMPLEX :: A1
C-----------------------------------------------
      IWBA = M + 1
      IWBB = IWBA + M
      IWBC = IWBB + M
      IWB2 = IWBC + M
      IWB3 = IWB2 + M
      IWW1 = IWB3 + M
      IWW2 = IWW1 + M
      IWW3 = IWW2 + M
      IWD = IWW3 + M
      IWTCOS = IWD + M
      IWP = IWTCOS + 4*N
      DO I = 1, M
         K = IWBA + I - 1
         W(K) = -A(I)
         K = IWBC + I - 1
         W(K) = -C(I)
         K = IWBB + I - 1
         W(K) = 2. - B(I)
         Y(I,:N) = -Y(I,:N)
      END DO
      MP = MPEROD + 1
      NP = NPEROD + 1
      GO TO (114,107) MP
  107 CONTINUE
      GO TO (108,109,110,111,123) NP
  108 CONTINUE
      CALL CMPOSP (M, N, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W(IWB2)
     1   , W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), W(IWP)
     2   )
      GO TO 112
  109 CONTINUE
      CALL CMPOSD (M, N, 1, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W(
     1   IWW1), W(IWD), W(IWTCOS), W(IWP))
      GO TO 112
  110 CONTINUE
      CALL CMPOSN (M, N, 1, 2, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W
     1   (IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), 
     2   W(IWP))
      GO TO 112
  111 CONTINUE
      CALL CMPOSN (M, N, 1, 1, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W
     1   (IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), 
     2   W(IWP))
  112 CONTINUE
      IPSTOR = REAL(W(IWW1))
      IREV = 2
      IF (NPEROD == 4) GO TO 124
  113 CONTINUE
      GO TO (127,133) MP
  114 CONTINUE
      MH = (M + 1)/2
      MHM1 = MH - 1
      MODD = 1
      IF (MH*2 == M) MODD = 2
      DO J = 1, N
         DO I = 1, MHM1
            W(I) = Y(MH-I,J) - Y(I+MH,J)
            W(I+MH) = Y(MH-I,J) + Y(I+MH,J)
         END DO
         W(MH) = 2.*Y(MH,J)
         GO TO (117,116) MODD
  116    CONTINUE
         W(M) = 2.*Y(M,J)
  117    CONTINUE
         Y(:M,J) = W(:M)
      END DO
      K = IWBC + MHM1 - 1
      I = IWBA + MHM1
      W(K) = (0.,0.)
      W(I) = (0.,0.)
      W(K+1) = 2.*W(K+1)
      SELECT CASE (MODD) 
      CASE DEFAULT
         K = IWBB + MHM1 - 1
         W(K) = W(K) - W(I-1)
         W(IWBC-1) = W(IWBC-1) + W(IWBB-1)
      CASE (2) 
         W(IWBB-1) = W(K+1)
      END SELECT
  122 CONTINUE
      GO TO 107
C
C     REVERSE COLUMNS WHEN NPEROD = 4
C
  123 CONTINUE
      IREV = 1
      NBY2 = N/2
  124 CONTINUE
      DO J = 1, NBY2
         MSKIP = N + 1 - J
         DO I = 1, M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
         END DO
      END DO
      GO TO (110,113) IREV
  127 CONTINUE
      DO J = 1, N
         W(MH-1:MH-MHM1:(-1)) = 0.5*(Y(MH+1:MHM1+MH,J)+Y(:MHM1,J))
         W(MH+1:MHM1+MH) = 0.5*(Y(MH+1:MHM1+MH,J)-Y(:MHM1,J))
         W(MH) = 0.5*Y(MH,J)
         GO TO (130,129) MODD
  129    CONTINUE
         W(M) = 0.5*Y(M,J)
  130    CONTINUE
         Y(:M,J) = W(:M)
      END DO
  133 CONTINUE
      W(1) = CMPLX(FLOAT(IPSTOR + IWP - 1),0.)
      RETURN 
      END SUBROUTINE CMGNBNN


      SUBROUTINE CMPOSD(MR,NR,ISTAG,BA,BB,BC,Q,IDIMQ,B,W,D,TCOS,P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MR
      INTEGER , INTENT(IN) :: NR
      INTEGER , INTENT(IN) :: ISTAG
      INTEGER , INTENT(IN) :: IDIMQ
      COMPLEX  :: BA(*)
      COMPLEX  :: BB(*)
      COMPLEX  :: BC(*)
      COMPLEX , INTENT(INOUT) :: Q(IDIMQ,1)
      COMPLEX  :: B(*)
      COMPLEX  :: W(*)
      COMPLEX  :: D(*)
      COMPLEX  :: TCOS(*)
      COMPLEX , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: M, N, IP, IPSTOR, JSH, KR, IRREG, JSTSAV, I, LR, NUN, 
     1   JST, JSP, L, NODD, J, JM1, JP1, JM2, JP2, JM3, JP3, NODDPR, IP1
     2   , KRPI, IDEG, JDEG
      REAL :: FI
      COMPLEX :: T
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION FOR DIRICHLET BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A+I.
C
      M = MR
      N = NR
      FI = 1./FLOAT(ISTAG)
      IP = -M
      IPSTOR = 0
      JSH = 0
      SELECT CASE (ISTAG) 
      CASE DEFAULT
         KR = 0
         IRREG = 1
         IF (N > 1) GO TO 106
         TCOS(1) = (0.,0.)
      CASE (2) 
         KR = 1
         JSTSAV = 1
         IRREG = 2
         IF (N > 1) GO TO 106
         TCOS(1) = CMPLX(-1.,0.)
      END SELECT
  103 CONTINUE
      B(:M) = Q(:M,1)
      CALL CMPTRX (1, 0, M, BA, BB, BC, B, TCOS, D, W)
      Q(:M,1) = B(:M)
      GO TO 183
  106 CONTINUE
      LR = 0
      DO I = 1, M
         P(I) = CMPLX(0.,0.)
      END DO
      NUN = N
      JST = 1
      JSP = N
C
C     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
C
  108 CONTINUE
      L = 2*JST
      NODD = 2 - 2*((NUN + 1)/2) + NUN
C
C     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
C
      SELECT CASE (NODD) 
      CASE DEFAULT
         JSP = JSP - L
      CASE (1) 
         JSP = JSP - JST
         IF (IRREG /= 1) JSP = JSP - L
      END SELECT
  111 CONTINUE
      CALL CMPCSG (JST, 1, 0.5, 0.0, TCOS)
      IF (L <= JSP) THEN
         DO J = L, JSP, L
            JM1 = J - JSH
            JP1 = J + JSH
            JM2 = J - JST
            JP2 = J + JST
            JM3 = JM2 - JSH
            JP3 = JP2 + JSH
            IF (JST == 1) THEN
               B(:M) = 2.*Q(:M,J)
               Q(:M,J) = Q(:M,JM2) + Q(:M,JP2)
            ELSE
               DO I = 1, M
                  T = Q(I,J) - Q(I,JM1) - Q(I,JP1) + Q(I,JM2) + Q(I,JP2)
                  B(I) = T + Q(I,J) - Q(I,JM3) - Q(I,JP3)
                  Q(I,J) = T
               END DO
            ENDIF
            CALL CMPTRX (JST, 0, M, BA, BB, BC, B, TCOS, D, W)
            Q(:M,J) = Q(:M,J) + B(:M)
         END DO
      ENDIF
C
C     REDUCTION FOR LAST UNKNOWN
C
      SELECT CASE (NODD) 
      CASE DEFAULT
         GO TO (152,120) IRREG
C
C     ODD NUMBER OF UNKNOWNS
C
  120    CONTINUE
         JSP = JSP + L
         J = JSP
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         JM3 = JM2 - JSH
         GO TO (123,121) ISTAG
  121    CONTINUE
         IF (JST /= 1) GO TO 123
         DO I = 1, M
            B(I) = Q(I,J)
            Q(I,J) = CMPLX(0.,0.)
         END DO
         GO TO 130
  123    CONTINUE
         SELECT CASE (NODDPR) 
         CASE DEFAULT
            B(:M) = 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3)) + P(IP+1:M+IP)
     1          + Q(:M,J)
         CASE (2) 
            B(:M) = 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3)) + Q(:M,JP2) - Q(
     1         :M,JP1) + Q(:M,J)
         END SELECT
  128    CONTINUE
         Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1))
  130    CONTINUE
         CALL CMPTRX (JST, 0, M, BA, BB, BC, B, TCOS, D, W)
         IP = IP + M
         IPSTOR = MAX0(IPSTOR,IP + M)
         P(IP+1:M+IP) = Q(:M,J) + B(:M)
         B(:M) = Q(:M,JP2) + P(IP+1:M+IP)
         IF (LR == 0) THEN
            DO I = 1, JST
               KRPI = KR + I
               TCOS(KRPI) = TCOS(I)
            END DO
         ELSE
            CALL CMPCSG (LR, JSTSAV, 0., FI, TCOS(JST+1))
            CALL CMPMRG (TCOS, 0, JST, JST, LR, KR)
         ENDIF
         CALL CMPCSG (KR, JSTSAV, 0.0, FI, TCOS)
         CALL CMPTRX (KR, KR, M, BA, BB, BC, B, TCOS, D, W)
         Q(:M,J) = Q(:M,JM2) + B(:M) + P(IP+1:M+IP)
         LR = KR
         KR = KR + L
C
C     EVEN NUMBER OF UNKNOWNS
C
      CASE (2) 
         JSP = JSP + L
         J = JSP
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         JM3 = JM2 - JSH
         SELECT CASE (IRREG) 
         CASE DEFAULT
            JSTSAV = JST
            IDEG = JST
            KR = L
         CASE (2) 
            CALL CMPCSG (KR, JSTSAV, 0.0, FI, TCOS)
            CALL CMPCSG (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
            IDEG = KR
            KR = KR + JST
         END SELECT
  139    CONTINUE
         IF (JST == 1) THEN
            IRREG = 2
            B(:M) = Q(:M,J)
            Q(:M,J) = Q(:M,JM2)
         ELSE
            B(:M) = Q(:M,J) + 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3))
            SELECT CASE (IRREG) 
            CASE DEFAULT
               Q(:M,J) = Q(:M,JM2) + 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1))
               IRREG = 2
            CASE (2) 
               SELECT CASE (NODDPR) 
               CASE DEFAULT
                  Q(:M,J) = Q(:M,JM2) + P(IP+1:M+IP)
                  IP = IP - M
               CASE (2) 
                  Q(:M,J) = Q(:M,JM2) + Q(:M,J) - Q(:M,JM1)
               END SELECT
            END SELECT
         ENDIF
  150    CONTINUE
         CALL CMPTRX (IDEG, LR, M, BA, BB, BC, B, TCOS, D, W)
         Q(:M,J) = Q(:M,J) + B(:M)
      END SELECT
  152 CONTINUE
      NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN >= 2) GO TO 108
C
C     START SOLUTION.
C
      J = JSP
      B(:M) = Q(:M,J)
      SELECT CASE (IRREG) 
      CASE DEFAULT
         CALL CMPCSG (JST, 1, 0.5, 0.0, TCOS)
         IDEG = JST
      CASE (2) 
         KR = LR + JST
         CALL CMPCSG (KR, JSTSAV, 0.0, FI, TCOS)
         CALL CMPCSG (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
         IDEG = KR
      END SELECT
  156 CONTINUE
      CALL CMPTRX (IDEG, LR, M, BA, BB, BC, B, TCOS, D, W)
      JM1 = J - JSH
      JP1 = J + JSH
      SELECT CASE (IRREG) 
      CASE DEFAULT
         Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1)) + B(:M)
      CASE (2) 
         SELECT CASE (NODDPR) 
         CASE DEFAULT
            Q(:M,J) = P(IP+1:M+IP) + B(:M)
            IP = IP - M
         CASE (2) 
            Q(:M,J) = Q(:M,J) - Q(:M,JM1) + B(:M)
         END SELECT
      END SELECT
  164 CONTINUE
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN > N) GO TO 183
      DO J = JST, N, L
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         IF (J <= JST) THEN
            B(:M) = Q(:M,J) + Q(:M,JP2)
         ELSE
            IF (JP2 <= N) GO TO 168
            B(:M) = Q(:M,J) + Q(:M,JM2)
            IF (JST < JSTSAV) IRREG = 1
            GO TO (170,171) IRREG
  168       CONTINUE
            B(:M) = Q(:M,J) + Q(:M,JM2) + Q(:M,JP2)
         ENDIF
  170    CONTINUE
         CALL CMPCSG (JST, 1, 0.5, 0.0, TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    CONTINUE
         IF (J + L > N) LR = LR - JST
         KR = JST + LR
         CALL CMPCSG (KR, JSTSAV, 0.0, FI, TCOS)
         CALL CMPCSG (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL CMPTRX (IDEG, JDEG, M, BA, BB, BC, B, TCOS, D, W)
         IF (JST <= 1) THEN
            Q(:M,J) = B(:M)
         ELSE
            IF (JP2 > N) GO TO 177
  175       CONTINUE
            Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1)) + B(:M)
            CYCLE 
  177       CONTINUE
            GO TO (175,178) IRREG
  178       CONTINUE
            IF (J + JSH <= N) THEN
               Q(:M,J) = B(:M) + P(IP+1:M+IP)
               IP = IP - M
            ELSE
               Q(:M,J) = B(:M) + Q(:M,J) - Q(:M,JM1)
            ENDIF
         ENDIF
      END DO
      L = L/2
      GO TO 164
  183 CONTINUE
      W(1) = CMPLX(FLOAT(IPSTOR),0.)
      RETURN 
      END SUBROUTINE CMPOSD


      SUBROUTINE CMPOSN(M, N, ISTAG, MIXBND, A, BB, C, Q, IDIMQ, B, B2, 
     1   B3, W, W2, W3, D, TCOS, P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: ISTAG
      INTEGER , INTENT(IN) :: MIXBND
      INTEGER , INTENT(IN) :: IDIMQ
      COMPLEX  :: A(*)
      COMPLEX  :: BB(*)
      COMPLEX  :: C(*)
      COMPLEX , INTENT(INOUT) :: Q(IDIMQ,*)
      COMPLEX  :: B(*)
      COMPLEX  :: B2(*)
      COMPLEX  :: B3(*)
      COMPLEX  :: W(*)
      COMPLEX  :: W2(*)
      COMPLEX  :: W3(*)
      COMPLEX  :: D(*)
      COMPLEX  :: TCOS(*)
      COMPLEX , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER , DIMENSION(4) :: K
      INTEGER :: K1, K2, K3, K4, MR, IP, IPSTOR, I2R, JR, NR, NLAST, KR
     1   , LR, I, NROD, JSTART, JSTOP, I2RBY2, J, JP1, JP2, JP3, JM1, 
     2   JM2, JM3, NRODPR, II, I1, I2, JR2, NLASTP, JSTEP
      REAL :: FISTAG, FNUM, FDEN
      COMPLEX :: FI, T
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION WITH NEUMANN BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS A-I.
C     MIXBND = 1 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTH BOUNDARIES.
C     MIXBND = 2 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTTOM AND
C     DIRICHLET CONDITION AT TOP.  (FOR THIS CASE, MUST HAVE ISTAG = 1.)
C
      EQUIVALENCE (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
      FISTAG = 3 - ISTAG
      FNUM = 1./FLOAT(ISTAG)
      FDEN = 0.5*FLOAT(ISTAG - 1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103) ISTAG
  101 CONTINUE
      Q(:MR,N) = 0.5*Q(:MR,N)
      GO TO (103,104) MIXBND
  103 CONTINUE
      IF (N <= 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 == NR) NROD = 0
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         JSTART = 1
      CASE (2) 
         JSTART = JR
         NROD = 1 - NROD
      END SELECT
  107 CONTINUE
      JSTOP = NLAST - JR
      IF (NROD == 0) JSTOP = JSTOP - I2R
      CALL CMPCSG (I2R, 1, 0.5, 0.0, TCOS)
      I2RBY2 = I2R/2
      IF (JSTOP < JSTART) THEN
         J = JR
      ELSE
         DO J = JSTART, JSTOP, JR
            JP1 = J + I2RBY2
            JP2 = J + I2R
            JP3 = JP2 + I2RBY2
            JM1 = J - I2RBY2
            JM2 = J - I2R
            JM3 = JM2 - I2RBY2
            IF (J == 1) THEN
               JM1 = JP1
               JM2 = JP2
               JM3 = JP3
            ENDIF
            IF (I2R == 1) THEN
               IF (J == 1) JM2 = JP2
               B(:MR) = 2.*Q(:MR,J)
               Q(:MR,J) = Q(:MR,JM2) + Q(:MR,JP2)
            ELSE
               DO I = 1, MR
                  FI = Q(I,J)
                  Q(I,J)=Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
                  B(I) = FI + Q(I,J) - Q(I,JM3) - Q(I,JP3)
               END DO
            ENDIF
            CALL CMPTRX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,J) = Q(:MR,J) + B(:MR)
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
         END DO
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
         J = JSTOP + JR
      ENDIF
      NLAST = J
      JM1 = J - I2RBY2
      JM2 = J - I2R
      JM3 = JM2 - I2RBY2
      IF (NROD /= 0) THEN
C
C     ODD NUMBER OF UNKNOWNS
C
         IF (I2R == 1) THEN
            B(:MR) = FISTAG*Q(:MR,J)
            Q(:MR,J) = Q(:MR,JM2)
         ELSE
            B(:MR) = Q(:MR,J) + 0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
            IF (NRODPR == 0) THEN
               Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP)
               IP = IP - MR
            ELSE
               Q(:MR,J) = Q(:MR,J) - Q(:MR,JM1) + Q(:MR,JM2)
            ENDIF
            IF (LR /= 0) THEN
               CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(KR+1))
            ELSE
               B(:MR) = FISTAG*B(:MR)
            ENDIF
         ENDIF
         CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS)
         CALL CMPTRX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         KR = KR + I2R
      ELSE
         JP1 = J + I2RBY2
         JP2 = J + I2R
         IF (I2R == 1) THEN
            B(:MR) = Q(:MR,J)
            CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            IP = 0
            IPSTOR = MR
            SELECT CASE (ISTAG) 
            CASE DEFAULT
               P(:MR) = B(:MR)
               B(:MR) = B(:MR) + Q(:MR,N)
               TCOS(1) = CMPLX(1.,0.)
               TCOS(2) = CMPLX(0.,0.)
               CALL CMPTRX (1, 1, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,JM2) + P(:MR) + B(:MR)
               GO TO 150
            CASE (1) 
               P(:MR) = B(:MR)
               Q(:MR,J) = Q(:MR,JM2) + 2.*Q(:MR,JP2) + 3.*B(:MR)
               GO TO 150
            END SELECT
         ENDIF
         B(:MR) = Q(:MR,J) + 0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
         IF (NRODPR == 0) THEN
            B(:MR) = B(:MR) + P(IP+1:MR+IP)
         ELSE
            B(:MR) = B(:MR) + Q(:MR,JP2) - Q(:MR,JP1)
         ENDIF
         CALL CMPTRX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
         IP = IP + MR
         IPSTOR = MAX0(IPSTOR,IP + MR)
         P(IP+1:MR+IP) = B(:MR) + 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         B(:MR) = P(IP+1:MR+IP) + Q(:MR,JP2)
         IF (LR /= 0) THEN
            CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(I2R+1))
            CALL CMPMRG (TCOS, 0, I2R, I2R, LR, KR)
         ELSE
            DO I = 1, I2R
               II = KR + I
               TCOS(II) = TCOS(I)
            END DO
         ENDIF
         CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS)
         IF (LR == 0) THEN
            GO TO (146,145) ISTAG
         ENDIF
  145    CONTINUE
         CALL CMPTRX (KR, KR, MR, A, BB, C, B, TCOS, D, W)
         GO TO 148
  146    CONTINUE
         B(:MR) = FISTAG*B(:MR)
  148    CONTINUE
         Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP) + B(:MR)
  150    CONTINUE
         LR = KR
         KR = KR + JR
      ENDIF
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         NR = (NLAST - 1)/JR + 1
         IF (NR <= 3) GO TO 155
      CASE (2) 
         NR = NLAST/JR
         IF (NR <= 1) GO TO 192
      END SELECT
  154 CONTINUE
      I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
      J = 1 + JR
      JM1 = J - I2R
      JP1 = J + I2R
      JM2 = NLAST - I2R
      IF (NR /= 2) THEN
         IF (LR /= 0) GO TO 170
         IF (N == 3) THEN
C
C     CASE N = 3.
C
            GO TO (156,168) ISTAG
  156       CONTINUE
            B(:MR) = Q(:MR,2)
            TCOS(1) = CMPLX(0.,0.)
            CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,2) = B(:MR)
            B(:MR) = 4.*B(:MR) + Q(:MR,1) + 2.*Q(:MR,3)
            TCOS(1) = CMPLX(-2.,0.)
            TCOS(2) = CMPLX(2.,0.)
            I1 = 2
            I2 = 0
            CALL CMPTRX (I1, I2, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,2) = Q(:MR,2) + B(:MR)
            B(:MR) = Q(:MR,1) + 2.*Q(:MR,2)
            TCOS(1) = (0.,0.)
            CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,1) = B(:MR)
            JR = 1
            I2R = 0
            GO TO 194
         ENDIF
C
C     CASE N = 2**P+1
C
         GO TO (162,170) ISTAG
  162    CONTINUE
         B(:MR) = Q(:MR,J) + 0.5*Q(:MR,1) - Q(:MR,JM1) + Q(:MR,NLAST) - 
     1      Q(:MR,JM2)
         CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS)
         CALL CMPTRX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1)) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,NLAST) + 4.*Q(:MR,J)
         JR2 = 2*JR
         CALL CMPCSG (JR, 1, 0.0, 0.0, TCOS)
         TCOS(JR+1:JR*2) = -TCOS(JR:1:(-1))
         CALL CMPTRX (JR2, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,J)
         CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS)
         CALL CMPTRX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168    CONTINUE
         B(:MR) = Q(:MR,2)
         Q(:MR,2) = (0.,0.)
         B2(:MR) = Q(:MR,3)
         B3(:MR) = Q(:MR,1)
         JR = 1
         I2R = 0
         J = 2
         GO TO 177
  170    CONTINUE
         B(:MR) = 0.5*Q(:MR,1) - Q(:MR,JM1) + Q(:MR,J)
         IF (NROD == 0) THEN
            B(:MR) = B(:MR) + P(IP+1:MR+IP)
         ELSE
            B(:MR) = B(:MR) + Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
         DO I = 1, MR
            T = 0.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
            Q(I,J) = T
            B2(I) = Q(I,NLAST) + T
            B3(I) = Q(I,1) + 2.*T
         END DO
  177    CONTINUE
         K1 = KR + 2*JR - 1
         K2 = KR + JR
         TCOS(K1+1) = (-2.,0.)
         K4 = K1 + 3 - ISTAG
         CALL CMPCSG (K2 + ISTAG - 2, 1, 0.0, FNUM, TCOS(K4))
         K4 = K1 + K2 + 1
         CALL CMPCSG (JR - 1, 1, 0.0, 1.0, TCOS(K4))
         CALL CMPMRG (TCOS, K1, K2, K1 + K2, JR - 1, 0)
         K3 = K1 + K2 + LR
         CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS(K3+1))
         K4 = K3 + JR + 1
         CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS(K4))
         CALL CMPMRG (TCOS, K3, JR, K3 + JR, KR, K1)
         IF (LR /= 0) THEN
            CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(K4))
            CALL CMPMRG (TCOS, K3, JR, K3 + JR, LR, K3 - LR)
            CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS(K4))
         ENDIF
         K3 = KR
         K4 = KR
         CALL CMPTR3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
         B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
         TCOS(1) = (2.,0.)
         CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,J)
         CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS)
         CALL CMPTRX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         IF (JR == 1) THEN
            Q(:MR,1) = B(:MR)
            GO TO 194
         ENDIF
         Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 194
      ENDIF
      IF (N == 2) THEN
C
C     CASE  N = 2
C
         B(:MR) = Q(:MR,1)
         TCOS(1) = (0.,0.)
         CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = B(:MR)
         B(:MR) = 2.*(Q(:MR,2)+B(:MR))*FISTAG
         TCOS(1) = CMPLX((-FISTAG),0.)
         TCOS(2) = CMPLX(2.,0.)
         CALL CMPTRX (2, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = Q(:MR,1) + B(:MR)
         JR = 1
         I2R = 0
         GO TO 194
      ENDIF
      B3(:MR) = (0.,0.)
      B(:MR) = Q(:MR,1) + 2.*P(IP+1:MR+IP)
      Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1)
      B2(:MR) = 2.*(Q(:MR,1)+Q(:MR,NLAST))
      K1 = KR + JR - 1
      TCOS(K1+1) = (-2.,0.)
      K4 = K1 + 3 - ISTAG
      CALL CMPCSG (KR + ISTAG - 2, 1, 0.0, FNUM, TCOS(K4))
      K4 = K1 + KR + 1
      CALL CMPCSG (JR - 1, 1, 0.0, 1.0, TCOS(K4))
      CALL CMPMRG (TCOS, K1, KR, K1 + KR, JR - 1, 0)
      CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS(K1+1))
      K2 = KR
      K4 = K1 + K2 + 1
      CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(K4))
      K3 = LR
      K4 = 0
      CALL CMPTR3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
      B(:MR) = B(:MR) + B2(:MR)
      TCOS(1) = (2.,0.)
      CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
      Q(:MR,1) = Q(:MR,1) + B(:MR)
      GO TO 194
  192 CONTINUE
      B(:MR) = Q(:MR,NLAST)
      GO TO 196
  194 CONTINUE
      J = NLAST - JR
      B(:MR) = Q(:MR,NLAST) + Q(:MR,J)
  196 CONTINUE
      JM2 = NLAST - I2R
      IF (JR == 1) THEN
         Q(:MR,NLAST) = (0.,0.)
      ELSE
         IF (NROD == 0) THEN
            Q(:MR,NLAST) = P(IP+1:MR+IP)
            IP = IP - MR
         ELSE
            Q(:MR,NLAST) = Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
      ENDIF
      CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS)
      CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(KR+1))
      IF (LR == 0) THEN
         B(:MR) = FISTAG*B(:MR)
      ENDIF
      CALL CMPTRX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
      Q(:MR,NLAST) = Q(:MR,NLAST) + B(:MR)
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR == 0) GO TO 222
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         JSTART = 1 + JR
      CASE (2) 
         JSTART = JR
      END SELECT
  209 CONTINUE
      KR = KR - JR
      IF (NLAST + JR <= N) THEN
         KR = KR - JR
         NLAST = NLAST + JR
         JSTOP = NLAST - JSTEP
      ELSE
         JSTOP = NLAST - JR
      ENDIF
      LR = KR - JR
      CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS)
      DO J = JSTART, JSTOP, JSTEP
         JM2 = J - JR
         JP2 = J + JR
         IF (J == JR) THEN
            B(:MR) = Q(:MR,J) + Q(:MR,JP2)
         ELSE
            B(:MR) = Q(:MR,J) + Q(:MR,JM2) + Q(:MR,JP2)
         ENDIF
         IF (JR == 1) THEN
            Q(:MR,J) = (0.,0.)
         ELSE
            JM1 = J - I2R
            JP1 = J + I2R
            Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         ENDIF
         CALL CMPTRX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
      END DO
      NROD = 1
      IF (NLAST + I2R <= N) NROD = 0
      IF (NLASTP /= NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
      W(1) = CMPLX(FLOAT(IPSTOR),0.)
      RETURN 
      END SUBROUTINE CMPOSN


      SUBROUTINE CMPOSP(M,N,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,W3,D,TCOS,P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER  :: IDIMQ
      COMPLEX  :: A(*)
      COMPLEX  :: BB(*)
      COMPLEX  :: C(*)
      COMPLEX  :: Q(IDIMQ,1)
      COMPLEX  :: B(*)
      COMPLEX  :: B2(*)
      COMPLEX  :: B3(*)
      COMPLEX  :: W(*)
      COMPLEX  :: W2(*)
      COMPLEX  :: W3(*)
      COMPLEX  :: D(*)
      COMPLEX  :: TCOS(*)
      COMPLEX  :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MR, NR, NRM1, J, NRMJ, NRPJ, I, IPSTOR, LH
      COMPLEX :: S, T
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON EQUATION WITH PERIODIC BOUNDARY
C     CONDITIONS.
C
      MR = M
      NR = (N + 1)/2
      NRM1 = NR - 1
      IF (2*NR == N) THEN
C
C     EVEN NUMBER OF UNKNOWNS
C
         DO J = 1, NRM1
            NRMJ = NR - J
            NRPJ = NR + J
            DO I = 1, MR
               S = Q(I,NRMJ) - Q(I,NRPJ)
               T = Q(I,NRMJ) + Q(I,NRPJ)
               Q(I,NRMJ) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 2.*Q(:MR,NR)
         Q(:MR,N) = 2.*Q(:MR,N)
         CALL CMPOSD (MR, NRM1, 1, A, BB, C, Q, IDIMQ, B, W, D, TCOS, P)
         IPSTOR = REAL(W(1))
         CALL CMPOSN (MR, NR + 1, 1, 1, A, BB, C, Q(1,NR), IDIMQ, B, B2
     1      , B3, W, W2, W3, D, TCOS, P)
         IPSTOR = MAX0(IPSTOR,INT(REAL(W(1))))
         DO J = 1, NRM1
            NRMJ = NR - J
            NRPJ = NR + J
            DO I = 1, MR
               S = 0.5*(Q(I,NRPJ)+Q(I,NRMJ))
               T = 0.5*(Q(I,NRPJ)-Q(I,NRMJ))
               Q(I,NRMJ) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 0.5*Q(:MR,NR)
         Q(:MR,N) = 0.5*Q(:MR,N)
      ELSE
         DO J = 1, NRM1
            NRPJ = N + 1 - J
            DO I = 1, MR
               S = Q(I,J) - Q(I,NRPJ)
               T = Q(I,J) + Q(I,NRPJ)
               Q(I,J) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 2.*Q(:MR,NR)
         LH = NRM1/2
         DO J = 1, LH
            NRMJ = NR - J
            DO I = 1, MR
               S = Q(I,J)
               Q(I,J) = Q(I,NRMJ)
               Q(I,NRMJ) = S
            END DO
         END DO
         CALL CMPOSD (MR, NRM1, 2, A, BB, C, Q, IDIMQ, B, W, D, TCOS, P)
         IPSTOR = REAL(W(1))
         CALL CMPOSN (MR, NR, 2, 1, A, BB, C, Q(1,NR), IDIMQ, B, B2, B3
     1      , W, W2, W3, D, TCOS, P)
         IPSTOR = MAX0(IPSTOR,INT(REAL(W(1))))
         DO J = 1, NRM1
            NRPJ = NR + J
            DO I = 1, MR
               S = 0.5*(Q(I,NRPJ)+Q(I,J))
               T = 0.5*(Q(I,NRPJ)-Q(I,J))
               Q(I,NRPJ) = T
               Q(I,J) = S
            END DO
         END DO
         Q(:MR,NR) = 0.5*Q(:MR,NR)
         DO J = 1, LH
            NRMJ = NR - J
            DO I = 1, MR
               S = Q(I,J)
               Q(I,J) = Q(I,NRMJ)
               Q(I,NRMJ) = S
            END DO
         END DO
      ENDIF
      W(1) = CMPLX(FLOAT(IPSTOR),0.)
      RETURN 
      END SUBROUTINE CMPOSP


      SUBROUTINE CMPCSG(N, IJUMP, FNUM, FDEN, A)
      implicit none
      REAL PIMACH
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IJUMP
      REAL , INTENT(IN) :: FNUM
      REAL , INTENT(IN) :: FDEN
      COMPLEX , INTENT(OUT) :: A(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K3, K4, K, K1, K5, I, K2, NP1
      REAL :: PI, DUM, PIBYN, X, Y
C-----------------------------------------------
C
C
C     THIS SUBROUTINE COMPUTES REQUIRED COSINE VALUES IN ASCENDING
C     ORDER.  WHEN IJUMP .GT. 1 THE ROUTINE COMPUTES VALUES
C
C        2*COS(J*PI/L) , J=1,2,...,L AND J .NE. 0(MOD N/IJUMP+1)
C
C     WHERE L = IJUMP*(N/IJUMP+1).
C
C
C     WHEN IJUMP = 1 IT COMPUTES
C
C            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
C
C     WHERE
C        FNUM = 0.5, FDEN = 0.0,  FOR REGULAR REDUCTION VALUES
C        FNUM = 0.0, FDEN = 1.0, FOR B-R AND C-R WHEN ISTAG = 1
C        FNUM = 0.0, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
C        FNUM = 0.5, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
C                                IN CMPOSN ONLY.
C
C
      PI = 4.0*ATAN(1.0)
      IF (N /= 0) THEN
         IF (IJUMP /= 1) THEN
            K3 = N/IJUMP + 1
            K4 = K3 - 1
            PIBYN = PI/FLOAT(N + IJUMP)
            DO K = 1, IJUMP
               K1 = (K - 1)*K3
               K5 = (K - 1)*K4
               DO I = 1, K4
                  X = K1 + I
                  K2 = K5 + I
                  A(K2) = CMPLX((-2.*COS(X*PIBYN)),0.)
               END DO
            END DO
         ELSE
            NP1 = N + 1
            Y = PI/(FLOAT(N) + FDEN)
            DO I = 1, N
               X = FLOAT(NP1 - I) - FNUM
               A(I) = CMPLX(2.*COS(X*Y),0.)
            END DO
         ENDIF
      ENDIF
      RETURN 
      END SUBROUTINE CMPCSG


      SUBROUTINE CMPMRG(TCOS, I1, M1, I2, M2, I3)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: M1
      INTEGER , INTENT(IN) :: I2
      INTEGER , INTENT(IN) :: M2
      INTEGER , INTENT(IN) :: I3
      COMPLEX , INTENT(INOUT) :: TCOS(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J11, J3, J1, J2, J, L, K, M
      COMPLEX :: X, Y
C-----------------------------------------------
C
C
C     THIS SUBROUTINE MERGES TWO ASCENDING STRINGS OF NUMBERS IN THE
C     ARRAY TCOS.  THE FIRST STRING IS OF LENGTH M1 AND STARTS AT
C     TCOS(I1+1).  THE SECOND STRING IS OF LENGTH M2 AND STARTS AT
C     TCOS(I2+1).  THE MERGED STRING GOES INTO TCOS(I3+1).
C
C
      J1 = 1
      J2 = 1
      J = I3
      IF (M1 == 0) GO TO 107
      IF (M2 == 0) GO TO 104
  101 CONTINUE
      J11 = J1
      J3 = MAX(M1,J11)
      DO J1 = J11, J3
         J = J + 1
         L = J1 + I1
         X = TCOS(L)
         L = J2 + I2
         Y = TCOS(L)
         IF (REAL(X - Y) > 0.) GO TO 103
         TCOS(J) = X
      END DO
      GO TO 106
  103 CONTINUE
      TCOS(J) = Y
      J2 = J2 + 1
      IF (J2 <= M2) GO TO 101
      IF (J1 > M1) GO TO 109
  104 CONTINUE
      K = J - J1 + 1
      DO J = J1, M1
         M = K + J
         L = J + I1
         TCOS(M) = TCOS(L)
      END DO
      GO TO 109
  106 CONTINUE
      IF (J2 > M2) GO TO 109
  107 CONTINUE
      K = J - J2 + 1
      DO J = J2, M2
         M = K + J
         L = J + I2
         TCOS(M) = TCOS(L)
      END DO
  109 CONTINUE
      RETURN 
      END SUBROUTINE CMPMRG


      SUBROUTINE CMPTRX(IDEGBR, IDEGCR, M, A, B, C, Y, TCOS, D, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDEGBR
      INTEGER , INTENT(IN) :: IDEGCR
      INTEGER , INTENT(IN) :: M
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: TCOS(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MM1, IFB, IFC, L, LINT, K, I, IP
      COMPLEX :: X, XX, Z
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE A SYSTEM OF LINEAR EQUATIONS WHERE THE
C     COEFFICIENT MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
C     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ).
C
      MM1 = M - 1
      IFB = IDEGBR + 1
      IFC = IDEGCR + 1
      L = IFB/IFC
      LINT = 1
      DO K = 1, IDEGBR
         X = TCOS(K)
         IF (K == L) THEN
            I = IDEGBR + LINT
            XX = X - TCOS(I)
            W(:M) = Y(:M)
            Y(:M) = XX*Y(:M)
         ENDIF
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y(1) = Y(1)*Z
         DO I = 2, MM1
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y(I) = (Y(I)-A(I)*Y(I-1))*Z
         END DO
         Z = B(M) - X - A(M)*D(MM1)
         IF (CABS(Z) == 0.) THEN
            Y(M) = (0.,0.)
         ELSE
            Y(M) = (Y(M)-A(M)*Y(MM1))/Z
         ENDIF
         DO IP = 1, MM1
            Y(M-IP) = Y(M-IP) - D(M-IP)*Y(M+1-IP)
         END DO
         IF (K /= L) CYCLE 
         Y(:M) = Y(:M) + W(:M)
         LINT = LINT + 1
         L = (LINT*IFB)/IFC
      END DO
      RETURN 
      END SUBROUTINE CMPTRX


      SUBROUTINE CMPTR3(M, A, B, C, K, Y1, Y2, Y3, TCOS, D, W1, W2, W3)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: K(4)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: Y1(*)
      COMPLEX , INTENT(INOUT) :: Y2(*)
      COMPLEX , INTENT(INOUT) :: Y3(*)
      COMPLEX , INTENT(IN) :: TCOS(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W1(*)
      COMPLEX , INTENT(INOUT) :: W2(*)
      COMPLEX , INTENT(INOUT) :: W3(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MM1, K1, K2, K3, K4, IF1, IF2, IF3, IF4, K2K3K4, L1, L2
     1   , L3, LINT1, LINT2, LINT3, KINT1, KINT2, KINT3, N, I, IP
      COMPLEX :: X, XX, Z
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE TRIDIAGONAL SYSTEMS
C
      MM1 = M - 1
      K1 = K(1)
      K2 = K(2)
      K3 = K(3)
      K4 = K(4)
      IF1 = K1 + 1
      IF2 = K2 + 1
      IF3 = K3 + 1
      IF4 = K4 + 1
      K2K3K4 = K2 + K3 + K4
      IF (K2K3K4 /= 0) THEN
         L1 = IF1/IF2
         L2 = IF1/IF3
         L3 = IF1/IF4
         LINT1 = 1
         LINT2 = 1
         LINT3 = 1
         KINT1 = K1
         KINT2 = KINT1 + K2
         KINT3 = KINT2 + K3
      ENDIF
      DO N = 1, K1
         X = TCOS(N)
         IF (K2K3K4 /= 0) THEN
            IF (N == L1) THEN
               W1(:M) = Y1(:M)
            ENDIF
            IF (N == L2) THEN
               W2(:M) = Y2(:M)
            ENDIF
            IF (N == L3) THEN
               W3(:M) = Y3(:M)
            ENDIF
         ENDIF
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y1(1) = Y1(1)*Z
         Y2(1) = Y2(1)*Z
         Y3(1) = Y3(1)*Z
         DO I = 2, M
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y1(I) = (Y1(I)-A(I)*Y1(I-1))*Z
            Y2(I) = (Y2(I)-A(I)*Y2(I-1))*Z
            Y3(I) = (Y3(I)-A(I)*Y3(I-1))*Z
         END DO
         DO IP = 1, MM1
            Y1(M-IP) = Y1(M-IP) - D(M-IP)*Y1(M+1-IP)
            Y2(M-IP) = Y2(M-IP) - D(M-IP)*Y2(M+1-IP)
            Y3(M-IP) = Y3(M-IP) - D(M-IP)*Y3(M+1-IP)
         END DO
         IF (K2K3K4 == 0) CYCLE 
         IF (N == L1) THEN
            I = LINT1 + KINT1
            XX = X - TCOS(I)
            Y1(:M) = XX*Y1(:M) + W1(:M)
            LINT1 = LINT1 + 1
            L1 = (LINT1*IF1)/IF2
         ENDIF
         IF (N == L2) THEN
            I = LINT2 + KINT2
            XX = X - TCOS(I)
            Y2(:M) = XX*Y2(:M) + W2(:M)
            LINT2 = LINT2 + 1
            L2 = (LINT2*IF1)/IF3
         ENDIF
         IF (N /= L3) CYCLE 
         I = LINT3 + KINT3
         XX = X - TCOS(I)
         Y3(:M) = XX*Y3(:M) + W3(:M)
         LINT3 = LINT3 + 1
         L3 = (LINT3*IF1)/IF4
      END DO
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE CMPTR3
C
C     file comf.f
C
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C PACKAGE COMF           THE ENTRIES IN THIS PACKAGE ARE LOWLEVEL
C                        ENTRIES, SUPPORTING FISHPACK ENTRIES BLKTRI
C                        AND CBLKTRI. THAT IS, THESE ROUTINES ARE
C                        NOT CALLED DIRECTLY BY USERS, BUT RATHER
C                        BY ENTRIES WITHIN BLKTRI AND CBLKTRI.
C                        DESCRIPTION OF ENTRIES EPMACH AND PIMACH
C                        FOLLOW BELOW.
C
C LATEST REVISION        JUNE 2004
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       NONE
C FILES
C
C LANGUAGE               FORTRAN 90
C ********************************************************************
C
C FUNCTION EPMACH (DUM)
C
C PURPOSE                TO COMPUTE AN APPROXIMATE MACHINE ACCURACY
C                        EPSILON ACCORDING TO THE FOLLOWING DEFINITION:
C                        EPSILON IS THE SMALLEST NUMBER SUCH THAT
C                        (1.+EPSILON).GT.1.)
C
C USAGE                  EPS = EPMACH (DUM)
C
C ARGUMENTS
C ON INPUT               DUM
C                          DUMMY VALUE
C
C ARGUMENTS
C ON OUTPUT              NONE
C
C HISTORY                THE ORIGINAL VERSION, WRITTEN WHEN THE
C                        BLKTRI PACKAGE WAS CONVERTED FROM THE
C                        CDC 7600 TO RUN ON THE CRAY-1, CALCULATED
C                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS
C                        BY 10.  USE OF THIS CONSTANT CAUSED BLKTRI
C                        TO COMPUTE SOLUTIONS ON THE CRAY-1 WITH FOUR
C                        FEWER PLACES OF ACCURACY THAN THE VERSION
C                        ON THE 7600.  IT WAS FOUND THAT COMPUTING
C                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS
C                        OF 2 PRODUCED A MACHINE ACCURACY 29% LESS
C                        THAN THE VALUE GENERATED BY SUCCESSIVE
C                        DIVISIONS BY 10, AND THAT USE OF THIS
C                        MACHINE CONSTANT IN THE BLKTRI PACKAGE
C                        RECOVERED THE ACCURACY THAT APPEARED TO
C                        BE LOST ON CONVERSION.
C
C ALGORITHM              COMPUTES MACHINE ACCURACY BY SUCCESSIVE
C                        DIVISIONS OF TWO.
C
C PORTABILITY            THIS CODE WILL EXECUTE ON MACHINES OTHER
C                        THAN THE CRAY1, BUT THE RETURNED VALUE MAY
C                        BE UNSATISFACTORY.  SEE HISTORY ABOVE.
C ********************************************************************
C
C FUNCTION PIMACH (DUM)
C
C PURPOSE                TO SUPPLY THE VALUE OF THE CONSTANT PI
C                        CORRECT TO MACHINE PRECISION WHERE
C                        PI=3.141592653589793238462643383279502884197
C                             1693993751058209749446
C
C USAGE                  PI = PIMACH (DUM)
C
C ARGUMENTS
C ON INPUT               DUM
C                          DUMMY VALUE
C
C ARGUMENTS
C ON OUTPUT              NONE
C
C ALGORITHM              THE VALUE OF PI IS SET TO 4.*ATAN(1.0)
C
C PORTABILITY            THIS ENTRY IS PORTABLE, BUT USERS SHOULD
C                        CHECK TO SEE WHETHER GREATER ACCURACY IS
C                        REQUIRED.
C
C***********************************************************************
      REAL FUNCTION EPMACH (DUM) 
      IMPLICIT NONE
      SAVE ALL
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL  :: DUM 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: ALL, EPS,V
      COMMON /FISH_VALUE/  V                                                 
      EPS = 1. 
      EPS = EPS/2. 
      CALL STRWRD (EPS + 1.) 
      DO WHILE(V - 1. > 0.) 
         EPS = EPS/2. 
         CALL STRWRD (EPS + 1.) 
      END DO 
      EPMACH = 100.*EPS 
      RETURN  
      END FUNCTION EPMACH 


      SUBROUTINE STRWRD(X) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL , INTENT(IN) :: X 
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /VALUE/ 
      COMMON /FISH_VALUE/ V 
      REAL   V 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: ALL 

      SAVE ALL 
C-----------------------------------------------
      V = X 
      RETURN  
      END SUBROUTINE STRWRD 


      REAL FUNCTION PIMACH (DUM) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL  :: DUM 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
      PIMACH = 4.*ATAN(1.0) 
      RETURN  
      END FUNCTION PIMACH 


      REAL FUNCTION PPSGF (X, IZ, C, A, BH) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL , INTENT(IN) :: X 
      REAL  :: C(*) 
      REAL  :: A(*) 
      REAL , INTENT(IN) :: BH(*) 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J 
      REAL :: ALL, SUM 

      SAVE ALL 
C-----------------------------------------------
      SUM = 0. 
      DO J = 1, IZ 
         SUM = SUM - 1./(X - BH(J))**2 
      END DO 
      PPSGF = SUM 
      RETURN  
      END FUNCTION PPSGF 


      REAL FUNCTION PPSPF (X, IZ, C, A, BH) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL , INTENT(IN) :: X 
      REAL  :: C(*) 
      REAL  :: A(*) 
      REAL , INTENT(IN) :: BH(*) 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J 
      REAL :: ALL, SUM 

      SAVE ALL 
C-----------------------------------------------
      SUM = 0. 
      DO J = 1, IZ 
         SUM = SUM + 1./(X - BH(J)) 
      END DO 
      PPSPF = SUM 
      RETURN  
      END FUNCTION PPSPF 


      REAL FUNCTION PSGF (X, IZ, C, A, BH) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL , INTENT(IN) :: X 
      REAL , INTENT(IN) :: C(*) 
      REAL , INTENT(IN) :: A(*) 
      REAL , INTENT(IN) :: BH(*) 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J 
      REAL :: ALL, FSG, HSG, DD 

      SAVE ALL 
C-----------------------------------------------
      FSG = 1. 
      HSG = 1. 
      DO J = 1, IZ 
         DD = 1./(X - BH(J)) 
         FSG = FSG*A(J)*DD 
         HSG = HSG*C(J)*DD 
      END DO 
      IF (MOD(IZ,2) == 0) THEN 
         PSGF = 1. - FSG - HSG 
         RETURN  
      ENDIF 
      PSGF = 1. + FSG + HSG 
      RETURN  
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
C June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END FUNCTION PSGF 
C
C     file fftpack.f
C
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C LATEST REVISION
C ---------------
C     June 2004      (VERSION 5.0) FORTRAN 90 CHANGES
C
C PURPOSE
C -------
C     THIS PACKAGE CONSISTS OF PROGRAMS WHICH PERFORM FAST FOURIER
C     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND
C     CERTAIN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW.
C
C USAGE
C -----
C     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB
C     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
C     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
C
C     4.   EZFFTI    INITIALIZE EZFFTF AND EZFFTB
C     5.   EZFFTF    A SIMPLIFIED REAL PERIODIC FORWARD TRANSFORM
C     6.   EZFFTB    A SIMPLIFIED REAL PERIODIC BACKWARD TRANSFORM
C
C     7.   SINTI     INITIALIZE SINT
C     8.   SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE
C
C     9.   COSTI     INITIALIZE COST
C     10.  COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE
C
C     11.  SINQI     INITIALIZE SINQF AND SINQB
C     12.  SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS
C     13.  SINQB     UNNORMALIZED INVERSE OF SINQF
C
C     14.  COSQI     INITIALIZE COSQF AND COSQB
C     15.  COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS
C     16.  COSQB     UNNORMALIZED INVERSE OF COSQF
C
C     17.  CFFTI     INITIALIZE CFFTF AND CFFTB
C     18.  CFFTF     FORWARD TRANSFORM OF A COMPLEX PERIODIC SEQUENCE
C     19.  CFFTB     UNNORMALIZED INVERSE OF CFFTF
C
C SPECIAL CONDITIONS
C ------------------
C     BEFORE CALLING ROUTINES EZFFTB AND EZFFTF FOR THE FIRST TIME,
C     OR BEFORE CALLING EZFFTB AND EZFFTF WITH A DIFFERENT LENGTH,
C     USERS MUST INITIALIZE BY CALLING ROUTINE EZFFTI.
C
C I/O
C ---
C     NONE
C
C PRECISION
C ---------
C     NONE
C
C REQUIRED LIBRARY FILES
C ----------------------
C     NONE
C
C LANGUAGE
C --------
C     FORTRAN
C
C HISTORY
C -------
C     DEVELOPED AT NCAR IN BOULDER, COLORADO BY PAUL N. SWARZTRAUBER
C     OF THE SCIENTIFIC COMPUTING DIVISION.  RELEASED ON NCAR'S PUBLIC
C     SOFTWARE LIBRARIES IN JANUARY 1980.  MODIFIED MAY 29, 1985 TO
C     INCREASE EFFICIENCY.
C
C PORTABILITY
C -----------
C     FORTRAN 77
C
C **********************************************************************
C
C     SUBROUTINE RFFTI(N,WSAVE)
C
C     SUBROUTINE RFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     BOTH RFFTF AND RFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
C             THE SAME WORK ARRAY CAN BE USED FOR BOTH RFFTF AND RFFTB
C             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
C             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
C             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF RFFTF OR RFFTB.
C
C **********************************************************************
C
C     SUBROUTINE RFFTF(N,R,WSAVE)
C
C     SUBROUTINE RFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL
C     PERODIC SEQUENCE (FOURIER ANALYSIS). THE TRANSFORM IS DEFINED
C     BELOW AT OUTPUT PARAMETER R.
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
C
C     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C             TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
C             IN THE PROGRAM THAT CALLS RFFTF. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.
C
C
C     OUTPUT PARAMETERS
C
C     R       R(1) = THE SUM FROM I=1 TO I=N OF R(I)
C
C             IF N IS EVEN SET L =N/2   , IF N IS ODD SET L = (N+1)/2
C
C               THEN FOR K = 2,...,L
C
C                  R(2*K-2) = THE SUM FROM I = 1 TO I = N OF
C
C                       R(I)*COS((K-1)*(I-1)*2*PI/N)
C
C                  R(2*K-1) = THE SUM FROM I = 1 TO I = N OF
C
C                      -R(I)*SIN((K-1)*(I-1)*2*PI/N)
C
C             IF N IS EVEN
C
C                  R(N) = THE SUM FROM I = 1 TO I = N OF
C
C                       (-1)**(I-1)*R(I)
C
C      *****  NOTE
C                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
C                  FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
C                  SEQUENCE BY N.
C
C     WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
C             CALLS OF RFFTF OR RFFTB.
C
C
C **********************************************************************
C
C     SUBROUTINE RFFTB(N,R,WSAVE)
C
C     SUBROUTINE RFFTB COMPUTES THE REAL PERODIC SEQUENCE FROM ITS
C     FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS DEFINED
C     BELOW AT OUTPUT PARAMETER R.
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
C
C     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C             TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
C             IN THE PROGRAM THAT CALLS RFFTB. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.
C
C
C     OUTPUT PARAMETERS
C
C     R       FOR N EVEN AND FOR I = 1,...,N
C
C                  R(I) = R(1)+(-1)**(I-1)*R(N)
C
C                       PLUS THE SUM FROM K=2 TO K=N/2 OF
C
C                        2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                       -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C             FOR N ODD AND FOR I = 1,...,N
C
C                  R(I) = R(1) PLUS THE SUM FROM K=2 TO K=(N+1)/2 OF
C
C                       2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                      -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C      *****  NOTE
C                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
C                  FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
C                  SEQUENCE BY N.
C
C     WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
C             CALLS OF RFFTB OR RFFTF.
C
C
C **********************************************************************
C
C     SUBROUTINE EZFFTI(N,WSAVE)
C
C     SUBROUTINE EZFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     BOTH EZFFTF AND EZFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C             THE SAME WORK ARRAY CAN BE USED FOR BOTH EZFFTF AND EZFFTB
C             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
C             ARE REQUIRED FOR DIFFERENT VALUES OF N.
C
C
C **********************************************************************
C
C     SUBROUTINE EZFFTF(N,R,AZERO,A,B,WSAVE)
C
C     SUBROUTINE EZFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL
C     PERODIC SEQUENCE (FOURIER ANALYSIS). THE TRANSFORM IS DEFINED
C     BELOW AT OUTPUT PARAMETERS AZERO,A AND B. EZFFTF IS A SIMPLIFIED
C     BUT SLOWER VERSION OF RFFTF.
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
C             IS MUST EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
C
C     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C             TO BE TRANSFORMED. R IS NOT DESTROYED.
C
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C             IN THE PROGRAM THAT CALLS EZFFTF. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE EZFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY EZFFTF AND EZFFTB.
C
C     OUTPUT PARAMETERS
C
C     AZERO   THE SUM FROM I=1 TO I=N OF R(I)/N
C
C     A,B     FOR N EVEN B(N/2)=0. AND A(N/2) IS THE SUM FROM I=1 TO
C             I=N OF (-1)**(I-1)*R(I)/N
C
C             FOR N EVEN DEFINE KMAX=N/2-1
C             FOR N ODD  DEFINE KMAX=(N-1)/2
C
C             THEN FOR  K=1,...,KMAX
C
C                  A(K) EQUALS THE SUM FROM I=1 TO I=N OF
C
C                       2./N*R(I)*COS(K*(I-1)*2*PI/N)
C
C                  B(K) EQUALS THE SUM FROM I=1 TO I=N OF
C
C                       2./N*R(I)*SIN(K*(I-1)*2*PI/N)
C
C
C **********************************************************************
C
C     SUBROUTINE EZFFTB(N,R,AZERO,A,B,WSAVE)
C
C     SUBROUTINE EZFFTB COMPUTES A REAL PERODIC SEQUENCE FROM ITS
C     FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS
C     DEFINED BELOW AT OUTPUT PARAMETER R. EZFFTB IS A SIMPLIFIED
C     BUT SLOWER VERSION OF RFFTB.
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE OUTPUT ARRAY R.  THE METHOD IS MOST
C             EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
C
C     AZERO   THE CONSTANT FOURIER COEFFICIENT
C
C     A,B     ARRAYS WHICH CONTAIN THE REMAINING FOURIER COEFFICIENTS
C             THESE ARRAYS ARE NOT DESTROYED.
C
C             THE LENGTH OF THESE ARRAYS DEPENDS ON WHETHER N IS EVEN OR
C             ODD.
C
C             IF N IS EVEN N/2    LOCATIONS ARE REQUIRED
C             IF N IS ODD (N-1)/2 LOCATIONS ARE REQUIRED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C             IN THE PROGRAM THAT CALLS EZFFTB. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE EZFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY EZFFTF AND EZFFTB.
C
C
C     OUTPUT PARAMETERS
C
C     R       IF N IS EVEN DEFINE KMAX=N/2
C             IF N IS ODD  DEFINE KMAX=(N-1)/2
C
C             THEN FOR I=1,...,N
C
C                  R(I)=AZERO PLUS THE SUM FROM K=1 TO K=KMAX OF
C
C                  A(K)*COS(K*(I-1)*2*PI/N)+B(K)*SIN(K*(I-1)*2*PI/N)
C
C     ********************* COMPLEX NOTATION **************************
C
C             FOR J=1,...,N
C
C             R(J) EQUALS THE SUM FROM K=-KMAX TO K=KMAX OF
C
C                  C(K)*EXP(I*K*(J-1)*2*PI/N)
C
C             WHERE
C
C                  C(K) = .5*CMPLX(A(K),-B(K))   FOR K=1,...,KMAX
C
C                  C(-K) = CONJG(C(K))
C
C                  C(0) = AZERO
C
C                       AND I=SQRT(-1)
C
C     *************** AMPLITUDE - PHASE NOTATION ***********************
C
C             FOR I=1,...,N
C
C             R(I) EQUALS AZERO PLUS THE SUM FROM K=1 TO K=KMAX OF
C
C                  ALPHA(K)*COS(K*(I-1)*2*PI/N+BETA(K))
C
C             WHERE
C
C                  ALPHA(K) = SQRT(A(K)*A(K)+B(K)*B(K))
C
C                  COS(BETA(K))=A(K)/ALPHA(K)
C
C                  SIN(BETA(K))=-B(K)/ALPHA(K)
C
C **********************************************************************
C
C     SUBROUTINE SINTI(N,WSAVE)
C
C     SUBROUTINE SINTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     SUBROUTINE SINT. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N+1 IS A PRODUCT OF SMALL PRIMES.
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WITH AT LEAST INT(2.5*N+15) LOCATIONS.
C             DIFFERENT WSAVE ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
C             OF N. THE CONTENTS OF WSAVE MUST NOT BE CHANGED BETWEEN
C             CALLS OF SINT.
C
C **********************************************************************
C
C     SUBROUTINE SINT(N,X,WSAVE)
C
C     SUBROUTINE SINT COMPUTES THE DISCRETE FOURIER SINE TRANSFORM
C     OF AN ODD SEQUENCE X(I). THE TRANSFORM IS DEFINED BELOW AT
C     OUTPUT PARAMETER X.
C
C     SINT IS THE UNNORMALIZED INVERSE OF ITSELF SINCE A CALL OF SINT
C     FOLLOWED BY ANOTHER CALL OF SINT WILL MULTIPLY THE INPUT SEQUENCE
C     X BY 2*(N+1).
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINT MUST BE
C     INITIALIZED BY CALLING SUBROUTINE SINTI(N,WSAVE).
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N+1 IS THE PRODUCT OF SMALL PRIMES.
C
C     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C
C     WSAVE   A WORK ARRAY WITH DIMENSION AT LEAST INT(2.5*N+15)
C             IN THE PROGRAM THAT CALLS SINT. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE SINTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C     OUTPUT PARAMETERS
C
C     X       FOR I=1,...,N
C
C                  X(I)= THE SUM FROM K=1 TO K=N
C
C                       2*X(K)*SIN(K*I*PI/(N+1))
C
C                  A CALL OF SINT FOLLOWED BY ANOTHER CALL OF
C                  SINT WILL MULTIPLY THE SEQUENCE X BY 2*(N+1).
C                  HENCE SINT IS THE UNNORMALIZED INVERSE
C                  OF ITSELF.
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
C             DESTROYED BETWEEN CALLS OF SINT.
C
C **********************************************************************
C
C     SUBROUTINE COSTI(N,WSAVE)
C
C     SUBROUTINE COSTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     SUBROUTINE COST. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF SMALL PRIMES.
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C             DIFFERENT WSAVE ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
C             OF N. THE CONTENTS OF WSAVE MUST NOT BE CHANGED BETWEEN
C             CALLS OF COST.
C
C **********************************************************************
C
C     SUBROUTINE COST(N,X,WSAVE)
C
C     SUBROUTINE COST COMPUTES THE DISCRETE FOURIER COSINE TRANSFORM
C     OF AN EVEN SEQUENCE X(I). THE TRANSFORM IS DEFINED BELOW AT OUTPUT
C     PARAMETER X.
C
C     COST IS THE UNNORMALIZED INVERSE OF ITSELF SINCE A CALL OF COST
C     FOLLOWED BY ANOTHER CALL OF COST WILL MULTIPLY THE INPUT SEQUENCE
C     X BY 2*(N-1). THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COST MUST BE
C     INITIALIZED BY CALLING SUBROUTINE COSTI(N,WSAVE).
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE SEQUENCE X. N MUST BE GREATER THAN 1.
C             THE METHOD IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF
C             SMALL PRIMES.
C
C     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
C             IN THE PROGRAM THAT CALLS COST. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE COSTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C     OUTPUT PARAMETERS
C
C     X       FOR I=1,...,N
C
C                 X(I) = X(1)+(-1)**(I-1)*X(N)
C
C                  + THE SUM FROM K=2 TO K=N-1
C
C                      2*X(K)*COS((K-1)*(I-1)*PI/(N-1))
C
C                  A CALL OF COST FOLLOWED BY ANOTHER CALL OF
C                  COST WILL MULTIPLY THE SEQUENCE X BY 2*(N-1)
C                  HENCE COST IS THE UNNORMALIZED INVERSE
C                  OF ITSELF.
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
C             DESTROYED BETWEEN CALLS OF COST.
C
C **********************************************************************
C
C     SUBROUTINE SINQI(N,WSAVE)
C
C     SUBROUTINE SINQI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     BOTH SINQF AND SINQB. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED. THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C             THE SAME WORK ARRAY CAN BE USED FOR BOTH SINQF AND SINQB
C             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
C             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
C             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF SINQF OR SINQB.
C
C **********************************************************************
C
C     SUBROUTINE SINQF(N,X,WSAVE)
C
C     SUBROUTINE SINQF COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
C     WAVE DATA. THAT IS , SINQF COMPUTES THE COEFFICIENTS IN A SINE
C     SERIES REPRESENTATION WITH ONLY ODD WAVE NUMBERS. THE TRANSFORM
C     IS DEFINED BELOW AT OUTPUT PARAMETER X.
C
C     SINQB IS THE UNNORMALIZED INVERSE OF SINQF SINCE A CALL OF SINQF
C     FOLLOWED BY A CALL OF SINQB WILL MULTIPLY THE INPUT SEQUENCE X
C     BY 4*N.
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINQF MUST BE
C     INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE).
C
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C             IN THE PROGRAM THAT CALLS SINQF. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C     OUTPUT PARAMETERS
C
C     X       FOR I=1,...,N
C
C                  X(I) = (-1)**(I-1)*X(N)
C
C                     + THE SUM FROM K=1 TO K=N-1 OF
C
C                     2*X(K)*SIN((2*I-1)*K*PI/(2*N))
C
C                  A CALL OF SINQF FOLLOWED BY A CALL OF
C                  SINQB WILL MULTIPLY THE SEQUENCE X BY 4*N.
C                  THEREFORE SINQB IS THE UNNORMALIZED INVERSE
C                  OF SINQF.
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
C             BE DESTROYED BETWEEN CALLS OF SINQF OR SINQB.
C
C **********************************************************************
C
C     SUBROUTINE SINQB(N,X,WSAVE)
C
C     SUBROUTINE SINQB COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
C     WAVE DATA. THAT IS , SINQB COMPUTES A SEQUENCE FROM ITS
C     REPRESENTATION IN TERMS OF A SINE SERIES WITH ODD WAVE NUMBERS.
C     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
C
C     SINQF IS THE UNNORMALIZED INVERSE OF SINQB SINCE A CALL OF SINQB
C     FOLLOWED BY A CALL OF SINQF WILL MULTIPLY THE INPUT SEQUENCE X
C     BY 4*N.
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINQB MUST BE
C     INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE).
C
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C             IN THE PROGRAM THAT CALLS SINQB. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C     OUTPUT PARAMETERS
C
C     X       FOR I=1,...,N
C
C                  X(I)= THE SUM FROM K=1 TO K=N OF
C
C                    4*X(K)*SIN((2K-1)*I*PI/(2*N))
C
C                  A CALL OF SINQB FOLLOWED BY A CALL OF
C                  SINQF WILL MULTIPLY THE SEQUENCE X BY 4*N.
C                  THEREFORE SINQF IS THE UNNORMALIZED INVERSE
C                  OF SINQB.
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
C             BE DESTROYED BETWEEN CALLS OF SINQB OR SINQF.
C
C **********************************************************************
C
C     SUBROUTINE COSQI(N,WSAVE)
C
C     SUBROUTINE COSQI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     BOTH COSQF AND COSQB. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE ARRAY TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
C             THE SAME WORK ARRAY CAN BE USED FOR BOTH COSQF AND COSQB
C             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
C             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
C             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF COSQF OR COSQB.
C
C **********************************************************************
C
C     SUBROUTINE COSQF(N,X,WSAVE)
C
C     SUBROUTINE COSQF COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
C     WAVE DATA. THAT IS , COSQF COMPUTES THE COEFFICIENTS IN A COSINE
C     SERIES REPRESENTATION WITH ONLY ODD WAVE NUMBERS. THE TRANSFORM
C     IS DEFINED BELOW AT OUTPUT PARAMETER X
C
C     COSQF IS THE UNNORMALIZED INVERSE OF COSQB SINCE A CALL OF COSQF
C     FOLLOWED BY A CALL OF COSQB WILL MULTIPLY THE INPUT SEQUENCE X
C     BY 4*N.
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COSQF MUST BE
C     INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE).
C
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
C             IN THE PROGRAM THAT CALLS COSQF. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C     OUTPUT PARAMETERS
C
C     X       FOR I=1,...,N
C
C                  X(I) = X(1) PLUS THE SUM FROM K=2 TO K=N OF
C
C                     2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
C
C                  A CALL OF COSQF FOLLOWED BY A CALL OF
C                  COSQB WILL MULTIPLY THE SEQUENCE X BY 4*N.
C                  THEREFORE COSQB IS THE UNNORMALIZED INVERSE
C                  OF COSQF.
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
C             BE DESTROYED BETWEEN CALLS OF COSQF OR COSQB.
C
C **********************************************************************
C
C     SUBROUTINE COSQB(N,X,WSAVE)
C
C     SUBROUTINE COSQB COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
C     WAVE DATA. THAT IS , COSQB COMPUTES A SEQUENCE FROM ITS
C     REPRESENTATION IN TERMS OF A COSINE SERIES WITH ODD WAVE NUMBERS.
C     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
C
C     COSQB IS THE UNNORMALIZED INVERSE OF COSQF SINCE A CALL OF COSQB
C     FOLLOWED BY A CALL OF COSQF WILL MULTIPLY THE INPUT SEQUENCE X
C     BY 4*N.
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COSQB MUST BE
C     INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE).
C
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C
C     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY THAT MUST BE DIMENSIONED AT LEAST 3*N+15
C             IN THE PROGRAM THAT CALLS COSQB. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C
C     OUTPUT PARAMETERS
C
C     X       FOR I=1,...,N
C
C                  X(I)= THE SUM FROM K=1 TO K=N OF
C
C                    4*X(K)*COS((2*K-1)*(I-1)*PI/(2*N))
C
C                  A CALL OF COSQB FOLLOWED BY A CALL OF
C                  COSQF WILL MULTIPLY THE SEQUENCE X BY 4*N.
C                  THEREFORE COSQF IS THE UNNORMALIZED INVERSE
C                  OF COSQB.
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
C             BE DESTROYED BETWEEN CALLS OF COSQB OR COSQF.
C
C **********************************************************************
C
C     SUBROUTINE CFFTI(N,WSAVE)
C
C     SUBROUTINE CFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     BOTH CFFTF AND CFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4*N+15
C             THE SAME WORK ARRAY CAN BE USED FOR BOTH CFFTF AND CFFTB
C             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
C             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
C             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF CFFTF OR CFFTB.
C
C **********************************************************************
C
C     SUBROUTINE CFFTF(N,C,WSAVE)
C
C     SUBROUTINE CFFTF COMPUTES THE FORWARD COMPLEX DISCRETE FOURIER
C     TRANSFORM (THE FOURIER ANALYSIS). EQUIVALENTLY , CFFTF COMPUTES
C     THE FOURIER COEFFICIENTS OF A COMPLEX PERIODIC SEQUENCE.
C     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
C
C     THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED TRANSFORM
C     THE OUTPUT MUST BE DIVIDED BY N. OTHERWISE A CALL OF CFFTF
C     FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE SEQUENCE BY N.
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE CFFTF MUST BE
C     INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).
C
C     INPUT PARAMETERS
C
C
C     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
C            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES. N
C
C     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C
C     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
C             IN THE PROGRAM THAT CALLS CFFTF. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY CFFTF AND CFFTB.
C
C     OUTPUT PARAMETERS
C
C     C      FOR J=1,...,N
C
C                C(J)=THE SUM FROM K=1,...,N OF
C
C                      C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)
C
C                            WHERE I=SQRT(-1)
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
C             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB
C
C **********************************************************************
C
C     SUBROUTINE CFFTB(N,C,WSAVE)
C
C     SUBROUTINE CFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER
C     TRANSFORM (THE FOURIER SYNTHESIS). EQUIVALENTLY , CFFTB COMPUTES
C     A COMPLEX PERIODIC SEQUENCE FROM ITS FOURIER COEFFICIENTS.
C     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
C
C     A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE
C     SEQUENCE BY N.
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE CFFTB MUST BE
C     INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).
C
C     INPUT PARAMETERS
C
C
C     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
C            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
C
C     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C
C     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
C             IN THE PROGRAM THAT CALLS CFFTB. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY CFFTF AND CFFTB.
C
C     OUTPUT PARAMETERS
C
C     C      FOR J=1,...,N
C
C                C(J)=THE SUM FROM K=1,...,N OF
C
C                      C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)
C
C                            WHERE I=SQRT(-1)
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
C             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB
C **********************************************************************
      SUBROUTINE EZFFTF(N, R, AZERO, A, B, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL , INTENT(OUT) :: AZERO
      REAL , INTENT(IN) :: R(*)
      REAL , INTENT(OUT) :: A(*)
      REAL , INTENT(OUT) :: B(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, NS2, NS2M
      REAL :: CF, CFM
C-----------------------------------------------
C
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            AZERO = R(1)
            RETURN 
         ENDIF
         AZERO = 0.5*(R(1)+R(2))
         A(1) = 0.5*(R(1)-R(2))
         RETURN 
      ENDIF
      WSAVE(:N) = R(:N)
      CALL RFFTF (N, WSAVE, WSAVE(N+1), IFAC)
      CF = 2./FLOAT(N)
      CFM = -CF
      AZERO = 0.5*CF*WSAVE(1)
      NS2 = (N + 1)/2
      NS2M = NS2 - 1
      A(:NS2M) = CF*WSAVE(2:NS2M*2:2)
      B(:NS2M) = CFM*WSAVE(3:NS2M*2+1:2)
      IF (MOD(N,2) == 1) RETURN 
      A(NS2) = 0.5*CF*WSAVE(N)
      B(NS2) = 0.
      RETURN 
      END SUBROUTINE EZFFTF


      SUBROUTINE EZFFTB(N, R, AZERO, A, B, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL , INTENT(IN) :: AZERO
      REAL  :: R(*)
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NS2, I
C-----------------------------------------------
C
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            R(1) = AZERO
            RETURN 
         ENDIF
         R(1) = AZERO + A(1)
         R(2) = AZERO - A(1)
         RETURN 
      ENDIF
      NS2 = (N - 1)/2
      R(2:NS2*2:2) = 0.5*A(:NS2)
      R(3:NS2*2+1:2) = -0.5*B(:NS2)
      R(1) = AZERO
      IF (MOD(N,2) == 0) R(N) = A(NS2+1)
      CALL RFFTB (N, R, WSAVE(N+1), IFAC)
      RETURN 
      END SUBROUTINE EZFFTB


      SUBROUTINE EZFFTI(N, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C
      IF (N == 1) RETURN 
      CALL EZFFT1 (N, WSAVE(2*N+1), IFAC)
      RETURN 
      END SUBROUTINE EZFFTI


      SUBROUTINE EZFFT1(N, WA, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: IFAC(*)
      REAL , INTENT(INOUT) :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER::NL,NF,J,NTRY,NQ,NR,I,IB,IS,NFM1,L1,K1,IP,L2,IDO,IPM,II
      REAL :: TPI, DUM, ARGH, ARG1, CH1, SH1, DCH1, DSH1, CH1H
C-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/ 
      TPI = 8.0*ATAN(1.0)
      NL = N
      NF = 0
      J = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF - 1
      L1 = 1
      IF (NFM1 == 0) RETURN 
      DO K1 = 1, NFM1
         IP = IFAC(K1+2)
         L2 = L1*IP
         IDO = N/L2
         IPM = IP - 1
         ARG1 = FLOAT(L1)*ARGH
         CH1 = 1.
         SH1 = 0.
         DCH1 = COS(ARG1)
         DSH1 = SIN(ARG1)
         DO J = 1, IPM
            CH1H = DCH1*CH1 - DSH1*SH1
            SH1 = DCH1*SH1 + DSH1*CH1
            CH1 = CH1H
            I = IS + 2
            WA(I-1) = CH1
            WA(I) = SH1
            IF (IDO >= 5) THEN
               DO II = 5, IDO, 2
                  I = I + 2
                  WA(I-1) = CH1*WA(I-3) - SH1*WA(I-2)
                  WA(I) = CH1*WA(I-2) + SH1*WA(I-3)
               END DO
            ENDIF
            IS = IS + IDO
         END DO
         L1 = L2
      END DO
      RETURN 
      END SUBROUTINE EZFFT1


      SUBROUTINE COSTI(N, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N
      REAL                :: WSAVE(*)
      INTEGER             :: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NM1, NP1, NS2, K, KC
      REAL :: PI, DUM, DT, FK
C-----------------------------------------------
C
      PI = 4.0*ATAN(1.0)
      IF (N <= 3) RETURN 
      NM1 = N - 1
      NP1 = N + 1
      NS2 = N/2
      DT = PI/FLOAT(NM1)
      FK = 0.
      DO K = 2, NS2
         KC = NP1 - K
         FK = FK + 1.
         WSAVE(K) = 2.*SIN(FK*DT)
         WSAVE(KC) = 2.*COS(FK*DT)
      END DO
      CALL RFFTI (NM1, WSAVE(N+1), IFAC)
      RETURN 
      END SUBROUTINE COSTI


      SUBROUTINE COST(N, X, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL  :: X(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NM1, NP1, NS2, K, KC, MODN, I
      REAL :: X1H, X1P3, TX2, C1, T1, T2, XIM2, XI
C-----------------------------------------------
C
      NM1 = N - 1
      NP1 = N + 1
      NS2 = N/2
      IF (N - 2 >= 0) THEN
         IF (N - 2 <= 0) THEN
            X1H = X(1) + X(2)
            X(2) = X(1) - X(2)
            X(1) = X1H
            RETURN 
         ENDIF
         IF (N <= 3) THEN
            X1P3 = X(1) + X(3)
            TX2 = X(2) + X(2)
            X(2) = X(1) - X(3)
            X(1) = X1P3 + TX2
            X(3) = X1P3 - TX2
            RETURN 
         ENDIF
         C1 = X(1) - X(N)
         X(1) = X(1) + X(N)
         DO K = 2, NS2
            KC = NP1 - K
            T1 = X(K) + X(KC)
            T2 = X(K) - X(KC)
            C1 = C1 + WSAVE(KC)*T2
            T2 = WSAVE(K)*T2
            X(K) = T1 - T2
            X(KC) = T1 + T2
         END DO
         MODN = MOD(N,2)
         IF (MODN /= 0) X(NS2+1) = X(NS2+1) + X(NS2+1)
         CALL RFFTF (NM1, X, WSAVE(N+1), IFAC)
         XIM2 = X(2)
         X(2) = C1
         DO I = 4, N, 2
            XI = X(I)
            X(I) = X(I-2) - X(I-1)
            X(I-1) = XIM2
            XIM2 = XI
         END DO
         IF (MODN /= 0) X(N) = XIM2
      ENDIF
      RETURN 
      END SUBROUTINE COST


      SUBROUTINE SINTI(N, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN):: N
      REAL                :: WSAVE(*)
      INTEGER             :: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NS2, NP1, K
      REAL :: PI, DUM, DT
C-----------------------------------------------
C
      PI = 4.0*ATAN(1.0)
      IF (N <= 1) RETURN 
      NS2 = N/2
      NP1 = N + 1
      DT = PI/FLOAT(NP1)
      DO K = 1, NS2
         WSAVE(K) = 2.*SIN(K*DT)
      END DO
      CALL RFFTI (NP1, WSAVE(NS2+1), IFAC)
      RETURN 
      END SUBROUTINE SINTI


      SUBROUTINE SINT(N, X, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: X(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NP1, IW1, IW2
C-----------------------------------------------
C
      NP1 = N + 1
      IW1 = N/2 + 1
      IW2 = IW1 + NP1
      CALL SINT1 (N, X, WSAVE, WSAVE(IW1), WSAVE(IW2), IFAC)
      RETURN 
      END SUBROUTINE SINT


      SUBROUTINE SINT1(N, WAR, WAS, XH, X, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER  :: IFAC(*)
      REAL  :: WAR(*)
      REAL , INTENT(IN) :: WAS(*)
      REAL  :: XH(*)
      REAL  :: X(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, NP1, NS2, K, KC, MODN
      REAL :: SQRT3, XHOLD, T1, T2
C-----------------------------------------------
      DATA SQRT3/ 1.73205080756888/ 
      XH(:N) = WAR(:N)
      WAR(:N) = X(:N)
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            XH(1) = XH(1) + XH(1)
            GO TO 106
         ENDIF
         XHOLD = SQRT3*(XH(1)+XH(2))
         XH(2) = SQRT3*(XH(1)-XH(2))
         XH(1) = XHOLD
         GO TO 106
      ENDIF
      NP1 = N + 1
      NS2 = N/2
      X(1) = 0.
      DO K = 1, NS2
         KC = NP1 - K
         T1 = XH(K) - XH(KC)
         T2 = WAS(K)*(XH(K)+XH(KC))
         X(K+1) = T1 + T2
         X(KC+1) = T2 - T1
      END DO
      MODN = MOD(N,2)
      IF (MODN /= 0) X(NS2+2) = 4.*XH(NS2+1)
      CALL RFFTF1 (NP1, X, XH, WAR, IFAC)
      XH(1) = 0.5*X(1)
      DO I = 3, N, 2
         XH(I-1) = -X(I)
         XH(I) = XH(I-2) + X(I-1)
      END DO
      IF (MODN == 0) XH(N) = -X(N+1)
  106 CONTINUE
      X(:N) = WAR(:N)
      WAR(:N) = XH(:N)
      RETURN 
      END SUBROUTINE SINT1


      SUBROUTINE COSQI(N, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K
      REAL :: PIH, DUM, DT, FK
C-----------------------------------------------
C
      PIH = 2.0*ATAN(1.0)
      DT = PIH/FLOAT(N)
      FK = 0.
      DO K = 1, N
         FK = FK + 1.
         WSAVE(K) = COS(FK*DT)
      END DO
      CALL RFFTI (N, WSAVE(N+1), IFAC)
      RETURN 
      END SUBROUTINE COSQI


      SUBROUTINE COSQF(N, X, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: X(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: SQRT2, TSQX
C-----------------------------------------------
      DATA SQRT2/ 1.4142135623731/ 
C
      IF (N - 2 >= 0) THEN
         IF (N - 2 > 0) GO TO 103
         TSQX = SQRT2*X(2)
         X(2) = X(1) - TSQX
         X(1) = X(1) + TSQX
      ENDIF
      RETURN 
  103 CONTINUE
      CALL COSQF1 (N, X, WSAVE, WSAVE(N+1), IFAC)
      RETURN 
      END SUBROUTINE COSQF


      SUBROUTINE COSQF1(N, X, W, XH, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: X(*)
      REAL , INTENT(IN) :: W(*)
      REAL  :: XH(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NS2, NP2, K, KC, MODN, I
      REAL :: XIM1
C-----------------------------------------------
      NS2 = (N + 1)/2
      NP2 = N + 2
      DO K = 2, NS2
         KC = NP2 - K
         XH(K) = X(K) + X(KC)
         XH(KC) = X(K) - X(KC)
      END DO
      MODN = MOD(N,2)
      IF (MODN == 0) XH(NS2+1) = X(NS2+1) + X(NS2+1)
      DO K = 2, NS2
         KC = NP2 - K
         X(K) = W(K-1)*XH(KC) + W(KC-1)*XH(K)
         X(KC) = W(K-1)*XH(K) - W(KC-1)*XH(KC)
      END DO
      IF (MODN == 0) X(NS2+1) = W(NS2)*XH(NS2+1)
      CALL RFFTF (N, X, XH, IFAC)
      DO I = 3, N, 2
         XIM1 = X(I-1) - X(I)
         X(I) = X(I-1) + X(I)
         X(I-1) = XIM1
      END DO
      RETURN 
      END SUBROUTINE COSQF1


      SUBROUTINE COSQB(N, X, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: X(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: TSQRT2, X1
C-----------------------------------------------
      DATA TSQRT2/ 2.82842712474619/ 
C
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            X(1) = 4.*X(1)
            RETURN 
         ENDIF
         X1 = 4.*(X(1)+X(2))
         X(2) = TSQRT2*(X(1)-X(2))
         X(1) = X1
         RETURN 
      ENDIF
      CALL COSQB1 (N, X, WSAVE, WSAVE(N+1), IFAC)
      RETURN 
      END SUBROUTINE COSQB


      SUBROUTINE COSQB1(N, X, W, XH, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: X(*)
      REAL , INTENT(IN) :: W(*)
      REAL  :: XH(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NS2, NP2, I, MODN, K, KC
      REAL :: XIM1
C-----------------------------------------------
      NS2 = (N + 1)/2
      NP2 = N + 2
      DO I = 3, N, 2
         XIM1 = X(I-1) + X(I)
         X(I) = X(I) - X(I-1)
         X(I-1) = XIM1
      END DO
      X(1) = X(1) + X(1)
      MODN = MOD(N,2)
      IF (MODN == 0) X(N) = X(N) + X(N)
      CALL RFFTB (N, X, XH, IFAC)
      DO K = 2, NS2
         KC = NP2 - K
         XH(K) = W(K-1)*X(KC) + W(KC-1)*X(K)
         XH(KC) = W(K-1)*X(K) - W(KC-1)*X(KC)
      END DO
      IF (MODN == 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
      DO K = 2, NS2
         KC = NP2 - K
         X(K) = XH(K) + XH(KC)
         X(KC) = XH(K) - XH(KC)
      END DO
      X(1) = X(1) + X(1)
      RETURN 
      END SUBROUTINE COSQB1


      SUBROUTINE SINQI(N, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C
      CALL COSQI (N, WSAVE, IFAC)
      RETURN 
      END SUBROUTINE SINQI


      SUBROUTINE SINQF(N, X, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: X(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NS2, K, KC
      REAL :: XHOLD
C-----------------------------------------------
C
      IF (N == 1) RETURN 
      NS2 = N/2
      DO K = 1, NS2
         KC = N - K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
      END DO
      CALL COSQF (N, X, WSAVE, IFAC)
      X(2:N:2) = -X(2:N:2)
      RETURN 
      END SUBROUTINE SINQF


      SUBROUTINE SINQB(N, X, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: X(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NS2, K, KC
      REAL :: XHOLD
C-----------------------------------------------
C
      IF (N <= 1) THEN
         X(1) = 4.*X(1)
         RETURN 
      ENDIF
      NS2 = N/2
      X(2:N:2) = -X(2:N:2)
      CALL COSQB (N, X, WSAVE, IFAC)
      DO K = 1, NS2
         KC = N - K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
      END DO
      RETURN 
      END SUBROUTINE SINQB


      SUBROUTINE CFFTI(N, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IW1
C-----------------------------------------------
C
      IF (N == 1) RETURN 
      IW1 = 2*N + 1
      CALL CFFTI1 (N, WSAVE(IW1), IFAC)
      RETURN 
      END SUBROUTINE CFFTI

      SUBROUTINE CFFTI1(N, WA, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: IFAC(*)
      REAL , INTENT(INOUT) :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER :: NL, NF, J, NTRY, NQ, NR, I, IB, L1, K1, IP, LD, L2, IDO
     1   , IDOT, IPM, I1, II
      REAL :: TPI, DUM, ARGH, FI, ARGLD, ARG
C-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 3, 4, 2, 5/ 
      NL = N
      NF = 0
      J = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.0*ATAN(1.0)
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO + IDO + 2
         IPM = IP - 1
         DO J = 1, IPM
            I1 = I
            WA(I-1) = 1.
            WA(I) = 0.
            LD = LD + L1
            FI = 0.
            ARGLD = FLOAT(LD)*ARGH
            DO II = 4, IDOT, 2
               I = I + 2
               FI = FI + 1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            END DO
            IF (IP <= 5) CYCLE 
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
         END DO
         L1 = L2
      END DO
      RETURN 
      END SUBROUTINE CFFTI1


      SUBROUTINE CFFTB(N, C, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: C(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C
      IF (N == 1) RETURN 
      CALL CFFTB1 (N, C, WSAVE, WSAVE(2*N+1), IFAC)
      RETURN 
      END SUBROUTINE CFFTB

      SUBROUTINE CFFTB1(N, C, CH, WA, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      REAL  :: C(*)
      REAL  :: CH(*)
      REAL  :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::NF,NA,L1,IW,K1,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,NAC,N2,I
C-----------------------------------------------
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO + IDO
         IDL1 = IDOT*L1
         IF (IP == 4) THEN
            IX2 = IW + IDOT
            IX3 = IX2 + IDOT
            IF (NA == 0) THEN
               CALL PASSB4 (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL PASSB4 (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL PASSB2 (IDOT, L1, C, CH, WA(IW))
               ELSE
                  CALL PASSB2 (IDOT, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDOT
                  IF (NA == 0) THEN
                     CALL PASSB3 (IDOT, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL PASSB3 (IDOT, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDOT
                     IX3 = IX2 + IDOT
                     IX4 = IX3 + IDOT
                     IF (NA == 0) THEN
                        CALL PASSB5 (IDOT, L1, C, CH, WA(IW), WA(IX2), 
     1                     WA(IX3), WA(IX4))
                     ELSE
                        CALL PASSB5 (IDOT, L1, CH, C, WA(IW), WA(IX2), 
     1                     WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL PASSB (NAC, IDOT, IP, L1, IDL1, C, C, C, CH
     1                     , CH, WA(IW))
                     ELSE
                        CALL PASSB (NAC, IDOT, IP, L1, IDL1, CH, CH, CH
     1                     , C, C, WA(IW))
                     ENDIF
                     IF (NAC /= 0) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDOT
      END DO
      IF (NA == 0) RETURN 
      N2 = N + N
      C(:N2) = CH(:N2)
      RETURN 
      END SUBROUTINE CFFTB1


      SUBROUTINE PASSB2(IDO, L1, CC, CH, WA1)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,2,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,2)
      REAL , INTENT(IN) :: WA1(1)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, I
      REAL :: TR2, TI2
C-----------------------------------------------
      IF (IDO <= 2) THEN
         CH(1,:,1) = CC(1,1,:) + CC(1,2,:)
         CH(1,:,2) = CC(1,1,:) - CC(1,2,:)
         CH(2,:,1) = CC(2,1,:) + CC(2,2,:)
         CH(2,:,2) = CC(2,1,:) - CC(2,2,:)
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            CH(I-1,K,1) = CC(I-1,1,K) + CC(I-1,2,K)
            TR2 = CC(I-1,1,K) - CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K) + CC(I,2,K)
            TI2 = CC(I,1,K) - CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2 + WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2 - WA1(I)*TI2
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB2


      SUBROUTINE PASSB3(IDO, L1, CC, CH, WA1, WA2)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,3,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,3)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, I
      REAL::TAUR,TAUI,TR2,CR2,TI2,CI2,CR3,CI3,DR2,DR3,DI2,DI3
C-----------------------------------------------
      DATA TAUR, TAUI/ -.5, 0.866025403784439/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TR2 = CC(1,2,K) + CC(1,3,K)
            CR2 = CC(1,1,K) + TAUR*TR2
            CH(1,K,1) = CC(1,1,K) + TR2
            TI2 = CC(2,2,K) + CC(2,3,K)
            CI2 = CC(2,1,K) + TAUR*TI2
            CH(2,K,1) = CC(2,1,K) + TI2
            CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
            CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
            CH(1,K,2) = CR2 - CI3
            CH(1,K,3) = CR2 + CI3
            CH(2,K,2) = CI2 + CR3
            CH(2,K,3) = CI2 - CR3
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TR2 = CC(I-1,2,K) + CC(I-1,3,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,2,K) + CC(I,3,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I,K,2) = WA1(I-1)*DI2 + WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2 - WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3 + WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3 - WA2(I)*DI3
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB3


      SUBROUTINE PASSB4(IDO, L1, CC, CH, WA1, WA2, WA3)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,4,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,4)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
      REAL , INTENT(IN) :: WA3(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, I
      REAL::TI1,TI2,TR4,TI3,TR1,TR2,TI4,TR3,CR3,CI3,CR2,CR4,CI2,CI4
C-----------------------------------------------
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI1 = CC(2,1,K) - CC(2,3,K)
            TI2 = CC(2,1,K) + CC(2,3,K)
            TR4 = CC(2,4,K) - CC(2,2,K)
            TI3 = CC(2,2,K) + CC(2,4,K)
            TR1 = CC(1,1,K) - CC(1,3,K)
            TR2 = CC(1,1,K) + CC(1,3,K)
            TI4 = CC(1,2,K) - CC(1,4,K)
            TR3 = CC(1,2,K) + CC(1,4,K)
            CH(1,K,1) = TR2 + TR3
            CH(1,K,3) = TR2 - TR3
            CH(2,K,1) = TI2 + TI3
            CH(2,K,3) = TI2 - TI3
            CH(1,K,2) = TR1 + TR4
            CH(1,K,4) = TR1 - TR4
            CH(2,K,2) = TI1 + TI4
            CH(2,K,4) = TI1 - TI4
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI1 = CC(I,1,K) - CC(I,3,K)
            TI2 = CC(I,1,K) + CC(I,3,K)
            TI3 = CC(I,2,K) + CC(I,4,K)
            TR4 = CC(I,4,K) - CC(I,2,K)
            TR1 = CC(I-1,1,K) - CC(I-1,3,K)
            TR2 = CC(I-1,1,K) + CC(I-1,3,K)
            TI4 = CC(I-1,2,K) - CC(I-1,4,K)
            TR3 = CC(I-1,2,K) + CC(I-1,4,K)
            CH(I-1,K,1) = TR2 + TR3
            CR3 = TR2 - TR3
            CH(I,K,1) = TI2 + TI3
            CI3 = TI2 - TI3
            CR2 = TR1 + TR4
            CR4 = TR1 - TR4
            CI2 = TI1 + TI4
            CI4 = TI1 - TI4
            CH(I-1,K,2) = WA1(I-1)*CR2 - WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2 + WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3 - WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3 + WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4 - WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4 + WA3(I)*CR4
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB4


      SUBROUTINE PASSB5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,5,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,5)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
      REAL , INTENT(IN) :: WA3(*)
      REAL , INTENT(IN) :: WA4(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, I
      REAL :: TR11, TI11, TR12, TI12, TI5, TI2, TI4, TI3, TR5, TR2, TR4
     1   , TR3, CR2, CI2, CR3, CI3, CR5, CI5, CR4, CI4, DR3, DR4, DI3, 
     2   DI4, DR5, DR2, DI5, DI2
C-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, 0.951056516295154
     1   , -.809016994374947, 0.587785252292473/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI5 = CC(2,2,K) - CC(2,5,K)
            TI2 = CC(2,2,K) + CC(2,5,K)
            TI4 = CC(2,3,K) - CC(2,4,K)
            TI3 = CC(2,3,K) + CC(2,4,K)
            TR5 = CC(1,2,K) - CC(1,5,K)
            TR2 = CC(1,2,K) + CC(1,5,K)
            TR4 = CC(1,3,K) - CC(1,4,K)
            TR3 = CC(1,3,K) + CC(1,4,K)
            CH(1,K,1) = CC(1,1,K) + TR2 + TR3
            CH(2,K,1) = CC(2,1,K) + TI2 + TI3
            CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(2,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(2,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            CH(1,K,2) = CR2 - CI5
            CH(1,K,5) = CR2 + CI5
            CH(2,K,2) = CI2 + CR5
            CH(2,K,3) = CI3 + CR4
            CH(1,K,3) = CR3 - CI4
            CH(1,K,4) = CR3 + CI4
            CH(2,K,4) = CI3 - CR4
            CH(2,K,5) = CI2 - CR5
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI5 = CC(I,2,K) - CC(I,5,K)
            TI2 = CC(I,2,K) + CC(I,5,K)
            TI4 = CC(I,3,K) - CC(I,4,K)
            TI3 = CC(I,3,K) + CC(I,4,K)
            TR5 = CC(I-1,2,K) - CC(I-1,5,K)
            TR2 = CC(I-1,2,K) + CC(I-1,5,K)
            TR4 = CC(I-1,3,K) - CC(I-1,4,K)
            TR3 = CC(I-1,3,K) + CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-1)*DR2 - WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2 + WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3 - WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3 + WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4 - WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4 + WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5 - WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5 + WA4(I)*DR5
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB5


      SUBROUTINE PASSB(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(OUT) :: NAC
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL , INTENT(IN) :: CC(IDO,IP,L1)
      REAL , INTENT(OUT) :: C1(IDO,L1,IP)
      REAL , INTENT(INOUT) :: C2(IDL1,IP)
      REAL , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL , INTENT(IN) :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDOT, NT, IPP2, IPPH, IDP, J, JC, K, I, IDL, INC, L, LC
     1   , IK, IDLJ, IDIJ, IDJ
      REAL :: WAR, WAI
C-----------------------------------------------
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IDP = IP*IDO
C
      IF (IDO >= L1) THEN
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ELSE
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      IDL = 2 - IDO
      INC = 0
      DO L = 2, IPPH
         LC = IPP2 - L
         IDL = IDL + IDO
         C2(:,L) = CH2(:,1) + WA(IDL-1)*CH2(:,2)
         C2(:,LC) = WA(IDL)*CH2(:,IP)
         IDLJ = IDL
         INC = INC + IDO
         DO J = 3, IPPH
            JC = IPP2 - J
            IDLJ = IDLJ + INC
            IF (IDLJ > IDP) IDLJ = IDLJ - IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            C2(:,L) = C2(:,L) + WAR*CH2(:,J)
            C2(:,LC) = C2(:,LC) + WAI*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH2(:IDL1-1:2,J) = C2(:IDL1-1:2,J) - C2(2:IDL1:2,JC)
         CH2(:IDL1-1:2,JC) = C2(:IDL1-1:2,J) + C2(2:IDL1:2,JC)
         CH2(2:IDL1:2,J) = C2(2:IDL1:2,J) + C2(:IDL1-1:2,JC)
         CH2(2:IDL1:2,JC) = C2(2:IDL1:2,J) - C2(:IDL1-1:2,JC)
      END DO
      NAC = 1
      IF (IDO == 2) RETURN 
      NAC = 0
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      C1(2,:,2:IP) = CH(2,:,2:IP)
      IF (IDOT <= L1) THEN
         IDIJ = 0
         DO J = 2, IP
            IDIJ = IDIJ + 2
            DO I = 4, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) - WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) + WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
         RETURN 
      ENDIF
      IDJ = 2 - IDO
      DO J = 2, IP
         IDJ = IDJ + IDO
         DO K = 1, L1
            IDIJ = IDJ
            C1(3:IDO-1:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(3:IDO-1:2,K,J
     1         ) - WA(IDIJ+2:IDO-2+IDIJ:2)*CH(4:IDO:2,K,J)
            C1(4:IDO:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(4:IDO:2,K,J) + 
     1         WA(IDIJ+2:IDO-2+IDIJ:2)*CH(3:IDO-1:2,K,J)
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB


      SUBROUTINE CFFTF(N, C, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: C(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C
      IF (N == 1) RETURN 
      CALL CFFTF1 (N, C, WSAVE, WSAVE(2*N+1), IFAC)
      RETURN 
      END SUBROUTINE CFFTF

      SUBROUTINE CFFTF1(N, C, CH, WA, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      REAL  :: C(*)
      REAL  :: CH(*)
      REAL  :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::NF,NA,L1,IW,K1,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,NAC,N2,I
C-----------------------------------------------
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO + IDO
         IDL1 = IDOT*L1
         IF (IP == 4) THEN
            IX2 = IW + IDOT
            IX3 = IX2 + IDOT
            IF (NA == 0) THEN
               CALL PASSF4 (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL PASSF4 (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL PASSF2 (IDOT, L1, C, CH, WA(IW))
               ELSE
                  CALL PASSF2 (IDOT, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDOT
                  IF (NA == 0) THEN
                     CALL PASSF3 (IDOT, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL PASSF3 (IDOT, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDOT
                     IX3 = IX2 + IDOT
                     IX4 = IX3 + IDOT
                     IF (NA == 0) THEN
                        CALL PASSF5 (IDOT, L1, C, CH, WA(IW), WA(IX2), 
     1                     WA(IX3), WA(IX4))
                     ELSE
                        CALL PASSF5 (IDOT, L1, CH, C, WA(IW), WA(IX2), 
     1                     WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL PASSF (NAC, IDOT, IP, L1, IDL1, C, C, C, CH
     1                     , CH, WA(IW))
                     ELSE
                        CALL PASSF (NAC, IDOT, IP, L1, IDL1, CH, CH, CH
     1                     , C, C, WA(IW))
                     ENDIF
                     IF (NAC /= 0) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDOT
      END DO
      IF (NA == 0) RETURN 
      N2 = N + N
      C(:N2) = CH(:N2)
      RETURN 
      END SUBROUTINE CFFTF1


      SUBROUTINE PASSF2(IDO, L1, CC, CH, WA1)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,2,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,2)
      REAL , INTENT(IN) :: WA1(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, I
      REAL :: TR2, TI2
C-----------------------------------------------
      IF (IDO <= 2) THEN
         CH(1,:,1) = CC(1,1,:) + CC(1,2,:)
         CH(1,:,2) = CC(1,1,:) - CC(1,2,:)
         CH(2,:,1) = CC(2,1,:) + CC(2,2,:)
         CH(2,:,2) = CC(2,1,:) - CC(2,2,:)
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            CH(I-1,K,1) = CC(I-1,1,K) + CC(I-1,2,K)
            TR2 = CC(I-1,1,K) - CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K) + CC(I,2,K)
            TI2 = CC(I,1,K) - CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2 - WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2 + WA1(I)*TI2
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF2


      SUBROUTINE PASSF3(IDO, L1, CC, CH, WA1, WA2)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,3,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,3)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, I
      REAL::TAUR,TAUI,TR2,CR2,TI2,CI2,CR3,CI3,DR2,DR3,DI2,DI3
C-----------------------------------------------
      DATA TAUR, TAUI/ -.5,  - 0.866025403784439/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TR2 = CC(1,2,K) + CC(1,3,K)
            CR2 = CC(1,1,K) + TAUR*TR2
            CH(1,K,1) = CC(1,1,K) + TR2
            TI2 = CC(2,2,K) + CC(2,3,K)
            CI2 = CC(2,1,K) + TAUR*TI2
            CH(2,K,1) = CC(2,1,K) + TI2
            CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
            CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
            CH(1,K,2) = CR2 - CI3
            CH(1,K,3) = CR2 + CI3
            CH(2,K,2) = CI2 + CR3
            CH(2,K,3) = CI2 - CR3
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TR2 = CC(I-1,2,K) + CC(I-1,3,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,2,K) + CC(I,3,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I,K,2) = WA1(I-1)*DI2 - WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2 + WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3 - WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3 + WA2(I)*DI3
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF3


      SUBROUTINE PASSF4(IDO, L1, CC, CH, WA1, WA2, WA3)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,4,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,4)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
      REAL , INTENT(IN) :: WA3(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, I
      REAL::TI1,TI2,TR4,TI3,TR1,TR2,TI4,TR3,CR3,CI3,CR2,CR4,CI2,CI4
C-----------------------------------------------
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI1 = CC(2,1,K) - CC(2,3,K)
            TI2 = CC(2,1,K) + CC(2,3,K)
            TR4 = CC(2,2,K) - CC(2,4,K)
            TI3 = CC(2,2,K) + CC(2,4,K)
            TR1 = CC(1,1,K) - CC(1,3,K)
            TR2 = CC(1,1,K) + CC(1,3,K)
            TI4 = CC(1,4,K) - CC(1,2,K)
            TR3 = CC(1,2,K) + CC(1,4,K)
            CH(1,K,1) = TR2 + TR3
            CH(1,K,3) = TR2 - TR3
            CH(2,K,1) = TI2 + TI3
            CH(2,K,3) = TI2 - TI3
            CH(1,K,2) = TR1 + TR4
            CH(1,K,4) = TR1 - TR4
            CH(2,K,2) = TI1 + TI4
            CH(2,K,4) = TI1 - TI4
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI1 = CC(I,1,K) - CC(I,3,K)
            TI2 = CC(I,1,K) + CC(I,3,K)
            TI3 = CC(I,2,K) + CC(I,4,K)
            TR4 = CC(I,2,K) - CC(I,4,K)
            TR1 = CC(I-1,1,K) - CC(I-1,3,K)
            TR2 = CC(I-1,1,K) + CC(I-1,3,K)
            TI4 = CC(I-1,4,K) - CC(I-1,2,K)
            TR3 = CC(I-1,2,K) + CC(I-1,4,K)
            CH(I-1,K,1) = TR2 + TR3
            CR3 = TR2 - TR3
            CH(I,K,1) = TI2 + TI3
            CI3 = TI2 - TI3
            CR2 = TR1 + TR4
            CR4 = TR1 - TR4
            CI2 = TI1 + TI4
            CI4 = TI1 - TI4
            CH(I-1,K,2) = WA1(I-1)*CR2 + WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2 - WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3 + WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3 - WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4 + WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4 - WA3(I)*CR4
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF4


      SUBROUTINE PASSF5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,5,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,5)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
      REAL , INTENT(IN) :: WA3(*)
      REAL , INTENT(IN) :: WA4(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, I
      REAL :: TR11, TI11, TR12, TI12, TI5, TI2, TI4, TI3, TR5, TR2, TR4
     1   , TR3, CR2, CI2, CR3, CI3, CR5, CI5, CR4, CI4, DR3, DR4, DI3, 
     2   DI4, DR5, DR2, DI5, DI2
C-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, -.951056516295154
     1   , -.809016994374947,  - 0.587785252292473/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI5 = CC(2,2,K) - CC(2,5,K)
            TI2 = CC(2,2,K) + CC(2,5,K)
            TI4 = CC(2,3,K) - CC(2,4,K)
            TI3 = CC(2,3,K) + CC(2,4,K)
            TR5 = CC(1,2,K) - CC(1,5,K)
            TR2 = CC(1,2,K) + CC(1,5,K)
            TR4 = CC(1,3,K) - CC(1,4,K)
            TR3 = CC(1,3,K) + CC(1,4,K)
            CH(1,K,1) = CC(1,1,K) + TR2 + TR3
            CH(2,K,1) = CC(2,1,K) + TI2 + TI3
            CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(2,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(2,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            CH(1,K,2) = CR2 - CI5
            CH(1,K,5) = CR2 + CI5
            CH(2,K,2) = CI2 + CR5
            CH(2,K,3) = CI3 + CR4
            CH(1,K,3) = CR3 - CI4
            CH(1,K,4) = CR3 + CI4
            CH(2,K,4) = CI3 - CR4
            CH(2,K,5) = CI2 - CR5
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI5 = CC(I,2,K) - CC(I,5,K)
            TI2 = CC(I,2,K) + CC(I,5,K)
            TI4 = CC(I,3,K) - CC(I,4,K)
            TI3 = CC(I,3,K) + CC(I,4,K)
            TR5 = CC(I-1,2,K) - CC(I-1,5,K)
            TR2 = CC(I-1,2,K) + CC(I-1,5,K)
            TR4 = CC(I-1,3,K) - CC(I-1,4,K)
            TR3 = CC(I-1,3,K) + CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-1)*DR2 + WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2 - WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3 + WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3 - WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4 + WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4 - WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5 + WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5 - WA4(I)*DR5
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF5


      SUBROUTINE PASSF(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(OUT) :: NAC
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL , INTENT(IN) :: CC(IDO,IP,L1)
      REAL , INTENT(OUT) :: C1(IDO,L1,IP)
      REAL , INTENT(INOUT) :: C2(IDL1,IP)
      REAL , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL , INTENT(IN) :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDOT, NT, IPP2, IPPH, IDP, J, JC, K, I, IDL, INC, L, LC
     1   , IK, IDLJ, IDIJ, IDJ
      REAL :: WAR, WAI
C-----------------------------------------------
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IDP = IP*IDO
C
      IF (IDO >= L1) THEN
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ELSE
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      IDL = 2 - IDO
      INC = 0
      DO L = 2, IPPH
         LC = IPP2 - L
         IDL = IDL + IDO
         C2(:,L) = CH2(:,1) + WA(IDL-1)*CH2(:,2)
         C2(:,LC) = -WA(IDL)*CH2(:,IP)
         IDLJ = IDL
         INC = INC + IDO
         DO J = 3, IPPH
            JC = IPP2 - J
            IDLJ = IDLJ + INC
            IF (IDLJ > IDP) IDLJ = IDLJ - IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            C2(:,L) = C2(:,L) + WAR*CH2(:,J)
            C2(:,LC) = C2(:,LC) - WAI*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH2(:IDL1-1:2,J) = C2(:IDL1-1:2,J) - C2(2:IDL1:2,JC)
         CH2(:IDL1-1:2,JC) = C2(:IDL1-1:2,J) + C2(2:IDL1:2,JC)
         CH2(2:IDL1:2,J) = C2(2:IDL1:2,J) + C2(:IDL1-1:2,JC)
         CH2(2:IDL1:2,JC) = C2(2:IDL1:2,J) - C2(:IDL1-1:2,JC)
      END DO
      NAC = 1
      IF (IDO == 2) RETURN 
      NAC = 0
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      C1(2,:,2:IP) = CH(2,:,2:IP)
      IF (IDOT <= L1) THEN
         IDIJ = 0
         DO J = 2, IP
            IDIJ = IDIJ + 2
            DO I = 4, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) + WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) - WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
         RETURN 
      ENDIF
      IDJ = 2 - IDO
      DO J = 2, IP
         IDJ = IDJ + IDO
         DO K = 1, L1
            IDIJ = IDJ
            C1(3:IDO-1:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(3:IDO-1:2,K,J
     1         ) + WA(IDIJ+2:IDO-2+IDIJ:2)*CH(4:IDO:2,K,J)
            C1(4:IDO:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(4:IDO:2,K,J) - 
     1         WA(IDIJ+2:IDO-2+IDIJ:2)*CH(3:IDO-1:2,K,J)
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF


      SUBROUTINE RFFTI(N, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER:: N
      REAL   :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C
      IF (N == 1) RETURN 
      CALL RFFTI1 (N, WSAVE(N+1), IFAC)
      RETURN 
      END SUBROUTINE RFFTI


      SUBROUTINE RFFTI1(N, WA, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: IFAC(*)
      REAL , INTENT(OUT) :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER :: NL, NF, J, NTRY, NQ, NR, I, IB, IS, NFM1, L1, K1, IP, 
     1   LD, L2, IDO, IPM, II
      REAL :: TPI, DUM, ARGH, ARGLD, FI, ARG
C-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/ 
      NL = N
      NF = 0
      J = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.0*ATAN(1.0)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF - 1
      L1 = 1
      IF (NFM1 == 0) RETURN 
      DO K1 = 1, NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP - 1
         DO J = 1, IPM
            LD = LD + L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            DO II = 3, IDO, 2
               I = I + 2
               FI = FI + 1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            END DO
            IS = IS + IDO
         END DO
         L1 = L2
      END DO
      RETURN 
      END SUBROUTINE RFFTI1


      SUBROUTINE RFFTB(N, R, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: R(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C
      IF (N == 1) RETURN 
      CALL RFFTB1 (N, R, WSAVE, WSAVE(N+1), IFAC)
      RETURN 
      END SUBROUTINE RFFTB

      SUBROUTINE RFFTB1(N, C, CH, WA, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      REAL  :: C(*)
      REAL  :: CH(*)
      REAL  :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NF, NA, L1, IW, K1, IP, L2, IDO, IDL1, IX2, IX3, IX4, I
C-----------------------------------------------
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP == 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA == 0) THEN
               CALL RADB4 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL RADB4 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL RADB2 (IDO, L1, C, CH, WA(IW))
               ELSE
                  CALL RADB2 (IDO, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDO
                  IF (NA == 0) THEN
                     CALL RADB3 (IDO, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL RADB3 (IDO, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDO
                     IX3 = IX2 + IDO
                     IX4 = IX3 + IDO
                     IF (NA == 0) THEN
                        CALL RADB5 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(
     1                     IX3), WA(IX4))
                     ELSE
                        CALL RADB5 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(
     1                     IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL RADBG(IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
                     ELSE
                        CALL RADBG(IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
                     ENDIF
                     IF (IDO == 1) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDO
      END DO
      IF (NA == 0) RETURN 
      C(:N) = CH(:N)
      RETURN 
      END SUBROUTINE RFFTB1


      SUBROUTINE RADB2(IDO, L1, CC, CH, WA1)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,2,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,2)
      REAL , INTENT(IN) :: WA1(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: TR2, TI2
C-----------------------------------------------
      CH(1,:,1) = CC(1,1,:) + CC(IDO,2,:)
      CH(1,:,2) = CC(1,1,:) - CC(IDO,2,:)
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  CH(I-1,K,1) = CC(I-1,1,K) + CC(IC-1,2,K)
                  TR2 = CC(I-1,1,K) - CC(IC-1,2,K)
                  CH(I,K,1) = CC(I,1,K) - CC(IC,2,K)
                  TI2 = CC(I,1,K) + CC(IC,2,K)
                  CH(I-1,K,2) = WA1(I-2)*TR2 - WA1(I-1)*TI2
                  CH(I,K,2) = WA1(I-2)*TI2 + WA1(I-1)*TR2
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         CH(IDO,:,1) = CC(IDO,1,:) + CC(IDO,1,:)
         CH(IDO,:,2) = -(CC(1,2,:)+CC(1,2,:))
      ENDIF
      RETURN 
      END SUBROUTINE RADB2


      SUBROUTINE RADB3(IDO, L1, CC, CH, WA1, WA2)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,3,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,3)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL::TAUR,TAUI,TR2,CR2,CI3,TI2,CI2,CR3,DR2,DR3,DI2,DI3
C-----------------------------------------------
      DATA TAUR, TAUI/ -.5, 0.866025403784439/ 
      DO K = 1, L1
         TR2 = CC(IDO,2,K) + CC(IDO,2,K)
         CR2 = CC(1,1,K) + TAUR*TR2
         CH(1,K,1) = CC(1,1,K) + TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2 - CI3
         CH(1,K,3) = CR2 + CI3
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            TR2 = CC(I-1,3,K) + CC(IC-1,2,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,3,K) - CC(IC,2,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I-1,K,2) = WA1(I-2)*DR2 - WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2 + WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3 - WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3 + WA2(I-1)*DR3
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADB3


      SUBROUTINE RADB4(IDO, L1, CC, CH, WA1, WA2, WA3)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,4,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,4)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
      REAL , INTENT(IN) :: WA3(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: SQRT2, TR1, TR2, TR3, TR4, TI1, TI2, TI3, TI4, CR3, CI3, 
     1   CR2, CR4, CI2, CI4
C-----------------------------------------------
      DATA SQRT2/ 1.414213562373095/ 
      DO K = 1, L1
         TR1 = CC(1,1,K) - CC(IDO,4,K)
         TR2 = CC(1,1,K) + CC(IDO,4,K)
         TR3 = CC(IDO,2,K) + CC(IDO,2,K)
         TR4 = CC(1,3,K) + CC(1,3,K)
         CH(1,K,1) = TR2 + TR3
         CH(1,K,2) = TR1 - TR4
         CH(1,K,3) = TR2 - TR3
         CH(1,K,4) = TR1 + TR4
      END DO
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  TI1 = CC(I,1,K) + CC(IC,4,K)
                  TI2 = CC(I,1,K) - CC(IC,4,K)
                  TI3 = CC(I,3,K) - CC(IC,2,K)
                  TR4 = CC(I,3,K) + CC(IC,2,K)
                  TR1 = CC(I-1,1,K) - CC(IC-1,4,K)
                  TR2 = CC(I-1,1,K) + CC(IC-1,4,K)
                  TI4 = CC(I-1,3,K) - CC(IC-1,2,K)
                  TR3 = CC(I-1,3,K) + CC(IC-1,2,K)
                  CH(I-1,K,1) = TR2 + TR3
                  CR3 = TR2 - TR3
                  CH(I,K,1) = TI2 + TI3
                  CI3 = TI2 - TI3
                  CR2 = TR1 - TR4
                  CR4 = TR1 + TR4
                  CI2 = TI1 + TI4
                  CI4 = TI1 - TI4
                  CH(I-1,K,2) = WA1(I-2)*CR2 - WA1(I-1)*CI2
                  CH(I,K,2) = WA1(I-2)*CI2 + WA1(I-1)*CR2
                  CH(I-1,K,3) = WA2(I-2)*CR3 - WA2(I-1)*CI3
                  CH(I,K,3) = WA2(I-2)*CI3 + WA2(I-1)*CR3
                  CH(I-1,K,4) = WA3(I-2)*CR4 - WA3(I-1)*CI4
                  CH(I,K,4) = WA3(I-2)*CI4 + WA3(I-1)*CR4
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         DO K = 1, L1
            TI1 = CC(1,2,K) + CC(1,4,K)
            TI2 = CC(1,4,K) - CC(1,2,K)
            TR1 = CC(IDO,1,K) - CC(IDO,3,K)
            TR2 = CC(IDO,1,K) + CC(IDO,3,K)
            CH(IDO,K,1) = TR2 + TR2
            CH(IDO,K,2) = SQRT2*(TR1 - TI1)
            CH(IDO,K,3) = TI2 + TI2
            CH(IDO,K,4) = -SQRT2*(TR1 + TI1)
         END DO
      ENDIF
      RETURN 
      END SUBROUTINE RADB4


      SUBROUTINE RADB5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,5,L1)
      REAL , INTENT(OUT) :: CH(IDO,L1,5)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
      REAL , INTENT(IN) :: WA3(*)
      REAL , INTENT(IN) :: WA4(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: TR11, TI11, TR12, TI12, TI5, TI4, TR2, TR3, CR2, CR3, CI5
     1   , CI4, TI2, TI3, TR5, TR4, CI2, CI3, CR5, CR4, DR3, DR4, DI3, 
     2   DI4, DR5, DR2, DI5, DI2
C-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, 0.951056516295154
     1   , -.809016994374947, 0.587785252292473/ 
      DO K = 1, L1
         TI5 = CC(1,3,K) + CC(1,3,K)
         TI4 = CC(1,5,K) + CC(1,5,K)
         TR2 = CC(IDO,2,K) + CC(IDO,2,K)
         TR3 = CC(IDO,4,K) + CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K) + TR2 + TR3
         CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
         CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
         CI5 = TI11*TI5 + TI12*TI4
         CI4 = TI12*TI5 - TI11*TI4
         CH(1,K,2) = CR2 - CI5
         CH(1,K,3) = CR3 - CI4
         CH(1,K,4) = CR3 + CI4
         CH(1,K,5) = CR2 + CI5
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            TI5 = CC(I,3,K) + CC(IC,2,K)
            TI2 = CC(I,3,K) - CC(IC,2,K)
            TI4 = CC(I,5,K) + CC(IC,4,K)
            TI3 = CC(I,5,K) - CC(IC,4,K)
            TR5 = CC(I-1,3,K) - CC(IC-1,2,K)
            TR2 = CC(I-1,3,K) + CC(IC-1,2,K)
            TR4 = CC(I-1,5,K) - CC(IC-1,4,K)
            TR3 = CC(I-1,5,K) + CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-2)*DR2 - WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2 + WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3 - WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3 + WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4 - WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4 + WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5 - WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5 + WA4(I-1)*DR5
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADB5


      SUBROUTINE RADBG(IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL , INTENT(IN) :: CC(IDO,IP,L1)
      REAL , INTENT(INOUT) :: C1(IDO,L1,IP)
      REAL , INTENT(INOUT) :: C2(IDL1,IP)
      REAL , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL , INTENT(IN) :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::IDP2,NBD,IPP2,IPPH,K,I,J,JC,J2,IC,L,LC,IK,IS,IDIJ
      REAL::TPI,DUM,ARG,DCP,DSP,AR1,AI1,AR1H,DC2,DS2,AR2,AI2,AR2H
C-----------------------------------------------
      TPI = 8.0*ATAN(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO + 2
      NBD = (IDO - 1)/2
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IF (IDO >= L1) THEN
         CH(:,:,1) = CC(:,1,:)
      ELSE
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      DO J = 2, IPPH
         JC = IPP2 - J
         J2 = J + J
         CH(1,:,J) = CC(IDO,J2-2,:) + CC(IDO,J2-2,:)
         CH(1,:,JC) = CC(1,J2-1,:) + CC(1,J2-1,:)
      END DO
      IF (IDO /= 1) THEN
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = CC(2:IDO-1:2,2*J-1,:) + CC(IDP2-4:
     1            IDP2-1-IDO:(-2),2*J-2,:)
               CH(2:IDO-1:2,:,JC) = CC(2:IDO-1:2,2*J-1,:) - CC(IDP2-4:
     1            IDP2-1-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,J) = CC(3:IDO:2,2*J-1,:) - CC(IDP2-3:IDP2-
     1            IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,JC) = CC(3:IDO:2,2*J-1,:) + CC(IDP2-3:IDP2-
     1            IDO:(-2),2*J-2,:)
            END DO
         ELSE
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = CC(2:IDO-1:2,2*J-1,:) + CC(IDP2-4:
     1            IDP2-1-IDO:(-2),2*J-2,:)
               CH(2:IDO-1:2,:,JC) = CC(2:IDO-1:2,2*J-1,:) - CC(IDP2-4:
     1            IDP2-1-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,J) = CC(3:IDO:2,2*J-1,:) - CC(IDP2-3:IDP2-
     1            IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,JC) = CC(3:IDO:2,2*J-1,:) + CC(IDP2-3:IDP2-
     1            IDO:(-2),2*J-2,:)
            END DO
         ENDIF
      ENDIF
      AR1 = 1.
      AI1 = 0.
      DO L = 2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         C2(:,L) = CH2(:,1) + AR1*CH2(:,2)
         C2(:,LC) = AI1*CH2(:,IP)
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO J = 3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            C2(:,L) = C2(:,L) + AR2*CH2(:,J)
            C2(:,LC) = C2(:,LC) + AI2*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH(1,:,J) = C1(1,:,J) - C1(1,:,JC)
         CH(1,:,JC) = C1(1,:,J) + C1(1,:,JC)
      END DO
      IF (IDO /= 1) THEN
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = C1(2:IDO-1:2,:,J) - C1(3:IDO:2,:,JC)
               CH(2:IDO-1:2,:,JC) = C1(2:IDO-1:2,:,J) + C1(3:IDO:2,:,JC)
               CH(3:IDO:2,:,J) = C1(3:IDO:2,:,J) + C1(2:IDO-1:2,:,JC)
               CH(3:IDO:2,:,JC) = C1(3:IDO:2,:,J) - C1(2:IDO-1:2,:,JC)
            END DO
         ELSE
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = C1(2:IDO-1:2,:,J) - C1(3:IDO:2,:,JC)
               CH(2:IDO-1:2,:,JC) = C1(2:IDO-1:2,:,J) + C1(3:IDO:2,:,JC)
               CH(3:IDO:2,:,J) = C1(3:IDO:2,:,J) + C1(2:IDO-1:2,:,JC)
               CH(3:IDO:2,:,JC) = C1(3:IDO:2,:,J) - C1(2:IDO-1:2,:,JC)
            END DO
         ENDIF
      ENDIF
      IF (IDO == 1) RETURN 
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      IF (NBD <= L1) THEN
         IS = -IDO
         DO J = 2, IP
            IS = IS + IDO
            IDIJ = IS
            DO I = 3, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) - WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) + WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
      ELSE
         IS = -IDO
         DO J = 2, IP
            IS = IS + IDO
            DO K = 1, L1
               IDIJ = IS
               C1(2:IDO-1:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*CH(2:IDO-1:2,
     1            K,J) - WA(IDIJ+2:IDO-1+IDIJ:2)*CH(3:IDO:2,K,J)
               C1(3:IDO:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*CH(3:IDO:2,K,J)
     1             + WA(IDIJ+2:IDO-1+IDIJ:2)*CH(2:IDO-1:2,K,J)
            END DO
         END DO
      ENDIF
      RETURN 
      END SUBROUTINE RADBG


      SUBROUTINE RFFTF(N, R, WSAVE, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      REAL  :: R(*)
      REAL  :: WSAVE(*)
      INTEGER:: IFAC(*)
C-----------------------------------------------
C
      IF (N == 1) RETURN 
      CALL RFFTF1 (N, R, WSAVE, WSAVE(N+1), IFAC)
      RETURN 
      END SUBROUTINE RFFTF


      SUBROUTINE RFFTF1(N, C, CH, WA, IFAC)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      REAL  :: C(*)
      REAL  :: CH(*)
      REAL  :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::NF,NA,L2,IW,K1,KH,IP,L1,IDO,IDL1,IX2,IX3,IX4,I
C-----------------------------------------------
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO K1 = 1, NF
         KH = NF - K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW - (IP - 1)*IDO
         NA = 1 - NA
         IF (IP == 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA == 0) THEN
               CALL RADF4 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
               GO TO 110
            ENDIF
            CALL RADF4 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            GO TO 110
         ENDIF
         IF (IP == 2) THEN
            IF (NA == 0) THEN
               CALL RADF2 (IDO, L1, C, CH, WA(IW))
               GO TO 110
            ENDIF
            CALL RADF2 (IDO, L1, CH, C, WA(IW))
            GO TO 110
         ENDIF
  104    CONTINUE
         IF (IP == 3) THEN
            IX2 = IW + IDO
            IF (NA == 0) THEN
               CALL RADF3 (IDO, L1, C, CH, WA(IW), WA(IX2))
               GO TO 110
            ENDIF
            CALL RADF3 (IDO, L1, CH, C, WA(IW), WA(IX2))
            GO TO 110
         ENDIF
  106    CONTINUE
         IF (IP == 5) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IX4 = IX3 + IDO
            IF (NA == 0) THEN
               CALL RADF5(IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
               GO TO 110
            ENDIF
            CALL RADF5(IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            GO TO 110
         ENDIF
  108    CONTINUE
         IF (IDO == 1) NA = 1 - NA
         IF (NA == 0) THEN
            CALL RADFG (IDO, IP, L1, IDL1, C, C, C, CH, CH, WA(IW))
            NA = 1
         ELSE
            CALL RADFG (IDO, IP, L1, IDL1, CH, CH, CH, C, C, WA(IW))
            NA = 0
         ENDIF
  110    CONTINUE
         L2 = L1
      END DO
      IF (NA == 1) RETURN 
      C(:N) = CH(:N)
      RETURN 
      END SUBROUTINE RFFTF1


      SUBROUTINE RADF2(IDO, L1, CC, CH, WA1)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,L1,2)
      REAL , INTENT(OUT) :: CH(IDO,2,L1)
      REAL , INTENT(IN) :: WA1(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: TR2, TI2
C-----------------------------------------------
      CH(1,1,:) = CC(1,:,1) + CC(1,:,2)
      CH(IDO,2,:) = CC(1,:,1) - CC(1,:,2)
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  TR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
                  TI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
                  CH(I,1,K) = CC(I,K,1) + TI2
                  CH(IC,2,K) = TI2 - CC(I,K,1)
                  CH(I-1,1,K) = CC(I-1,K,1) + TR2
                  CH(IC-1,2,K) = CC(I-1,K,1) - TR2
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         CH(1,2,:) = -CC(IDO,:,2)
         CH(IDO,1,:) = CC(IDO,:,1)
      ENDIF
      RETURN 
      END SUBROUTINE RADF2


      SUBROUTINE RADF3(IDO, L1, CC, CH, WA1, WA2)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,L1,3)
      REAL , INTENT(OUT) :: CH(IDO,3,L1)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL::TAUR,TAUI,CR2,DR2,DI2,DR3,DI3,CI2,TR2,TI2,TR3,TI3
C-----------------------------------------------
      DATA TAUR, TAUI/ -.5, 0.866025403784439/ 
      DO K = 1, L1
         CR2 = CC(1,K,2) + CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1) + TAUR*CR2
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2 + DR3
            CI2 = DI2 + DI3
            CH(I-1,1,K) = CC(I-1,K,1) + CR2
            CH(I,1,K) = CC(I,K,1) + CI2
            TR2 = CC(I-1,K,1) + TAUR*CR2
            TI2 = CC(I,K,1) + TAUR*CI2
            TR3 = TAUI*(DI2 - DI3)
            TI3 = TAUI*(DR3 - DR2)
            CH(I-1,3,K) = TR2 + TR3
            CH(IC-1,2,K) = TR2 - TR3
            CH(I,3,K) = TI2 + TI3
            CH(IC,2,K) = TI3 - TI2
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADF3


      SUBROUTINE RADF4(IDO, L1, CC, CH, WA1, WA2, WA3)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,L1,4)
      REAL , INTENT(OUT) :: CH(IDO,4,L1)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
      REAL , INTENT(IN) :: WA3(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: HSQT2, TR1, TR2, CR2, CI2, CR3, CI3, CR4, CI4, TR4, TI1, 
     1   TI4, TI2, TI3, TR3
C-----------------------------------------------
      DATA HSQT2/ 0.7071067811865475/ 
      DO K = 1, L1
         TR1 = CC(1,K,2) + CC(1,K,4)
         TR2 = CC(1,K,1) + CC(1,K,3)
         CH(1,1,K) = TR1 + TR2
         CH(IDO,4,K) = TR2 - TR1
         CH(IDO,2,K) = CC(1,K,1) - CC(1,K,3)
         CH(1,3,K) = CC(1,K,4) - CC(1,K,2)
      END DO
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  CR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
                  CI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
                  CR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
                  CI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
                  CR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
                  CI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
                  TR1 = CR2 + CR4
                  TR4 = CR4 - CR2
                  TI1 = CI2 + CI4
                  TI4 = CI2 - CI4
                  TI2 = CC(I,K,1) + CI3
                  TI3 = CC(I,K,1) - CI3
                  TR2 = CC(I-1,K,1) + CR3
                  TR3 = CC(I-1,K,1) - CR3
                  CH(I-1,1,K) = TR1 + TR2
                  CH(IC-1,4,K) = TR2 - TR1
                  CH(I,1,K) = TI1 + TI2
                  CH(IC,4,K) = TI1 - TI2
                  CH(I-1,3,K) = TI4 + TR3
                  CH(IC-1,2,K) = TR3 - TI4
                  CH(I,3,K) = TR4 + TI3
                  CH(IC,2,K) = TR4 - TI3
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         DO K = 1, L1
            TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
            TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
            CH(IDO,1,K) = TR1 + CC(IDO,K,1)
            CH(IDO,3,K) = CC(IDO,K,1) - TR1
            CH(1,2,K) = TI1 - CC(IDO,K,3)
            CH(1,4,K) = TI1 + CC(IDO,K,3)
         END DO
      ENDIF
      RETURN 
      END SUBROUTINE RADF4


      SUBROUTINE RADF5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL , INTENT(IN) :: CC(IDO,L1,5)
      REAL , INTENT(OUT) :: CH(IDO,5,L1)
      REAL , INTENT(IN) :: WA1(*)
      REAL , INTENT(IN) :: WA2(*)
      REAL , INTENT(IN) :: WA3(*)
      REAL , INTENT(IN) :: WA4(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: TR11, TI11, TR12, TI12, CR2, CI5, CR3, CI4, DR2, DI2, DR3
     1   , DI3, DR4, DI4, DR5, DI5, CR5, CI2, CR4, CI3, TR2, TI2, TR3, 
     2   TI3, TR5, TI5, TR4, TI4
C-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, 0.951056516295154
     1   , -.809016994374947, 0.587785252292473/ 
      DO K = 1, L1
         CR2 = CC(1,K,5) + CC(1,K,2)
         CI5 = CC(1,K,5) - CC(1,K,2)
         CR3 = CC(1,K,4) + CC(1,K,3)
         CI4 = CC(1,K,4) - CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2 + CR3
         CH(IDO,2,K) = CC(1,K,1) + TR11*CR2 + TR12*CR3
         CH(1,3,K) = TI11*CI5 + TI12*CI4
         CH(IDO,4,K) = CC(1,K,1) + TR12*CR2 + TR11*CR3
         CH(1,5,K) = TI12*CI5 - TI11*CI4
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5) + WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5) - WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2 + DR5
            CI5 = DR5 - DR2
            CR5 = DI2 - DI5
            CI2 = DI2 + DI5
            CR3 = DR3 + DR4
            CI4 = DR4 - DR3
            CR4 = DI3 - DI4
            CI3 = DI3 + DI4
            CH(I-1,1,K) = CC(I-1,K,1) + CR2 + CR3
            CH(I,1,K) = CC(I,K,1) + CI2 + CI3
            TR2 = CC(I-1,K,1) + TR11*CR2 + TR12*CR3
            TI2 = CC(I,K,1) + TR11*CI2 + TR12*CI3
            TR3 = CC(I-1,K,1) + TR12*CR2 + TR11*CR3
            TI3 = CC(I,K,1) + TR12*CI2 + TR11*CI3
            TR5 = TI11*CR5 + TI12*CR4
            TI5 = TI11*CI5 + TI12*CI4
            TR4 = TI12*CR5 - TI11*CR4
            TI4 = TI12*CI5 - TI11*CI4
            CH(I-1,3,K) = TR2 + TR5
            CH(IC-1,2,K) = TR2 - TR5
            CH(I,3,K) = TI2 + TI5
            CH(IC,2,K) = TI5 - TI2
            CH(I-1,5,K) = TR3 + TR4
            CH(IC-1,4,K) = TR3 - TR4
            CH(I,5,K) = TI3 + TI4
            CH(IC,4,K) = TI4 - TI3
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADF5


      SUBROUTINE RADFG(IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL , INTENT(OUT) :: CC(IDO,IP,L1)
      REAL , INTENT(INOUT) :: C1(IDO,L1,IP)
      REAL , INTENT(INOUT) :: C2(IDL1,IP)
      REAL , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL , INTENT(IN) :: WA(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::IPPH,IPP2,IDP2,NBD,IK,J,K,IS,IDIJ,I,JC,L,LC,J2,IC
      REAL::TPI,DUM,ARG,DCP,DSP,AR1,AI1,AR1H,DC2,DS2,AR2,AI2,AR2H
C-----------------------------------------------
      TPI = 8.0*ATAN(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP + 1)/2
      IPP2 = IP + 2
      IDP2 = IDO + 2
      NBD = (IDO - 1)/2
      IF (IDO /= 1) THEN
         CH2(:,1) = C2(:,1)
         CH(1,:,2:IP) = C1(1,:,2:IP)
         IF (NBD <= L1) THEN
            IS = -IDO
            DO J = 2, IP
               IS = IS + IDO
               IDIJ = IS
               DO I = 3, IDO, 2
                  IDIJ = IDIJ + 2
                  CH(I-1,:,J)=WA(IDIJ-1)*C1(I-1,:,J)+WA(IDIJ)*C1(I,:,J)
                  CH(I,:,J)=WA(IDIJ-1)*C1(I,:,J)-WA(IDIJ)*C1(I-1,:,J)
               END DO
            END DO
         ELSE
            IS = -IDO
            DO J = 2, IP
               IS = IS + IDO
               DO K = 1, L1
                  IDIJ = IS
                  CH(2:IDO-1:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*C1(2:IDO-1
     1               :2,K,J) + WA(IDIJ+2:IDO-1+IDIJ:2)*C1(3:IDO:2,K,J)
                  CH(3:IDO:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*C1(3:IDO:2,K
     1               ,J) - WA(IDIJ+2:IDO-1+IDIJ:2)*C1(2:IDO-1:2,K,J)
               END DO
            END DO
         ENDIF
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               C1(2:IDO-1:2,:,J)=CH(2:IDO-1:2,:,J)+CH(2:IDO-1:2,:,JC)
               C1(2:IDO-1:2,:,JC) = CH(3:IDO:2,:,J) - CH(3:IDO:2,:,JC)
               C1(3:IDO:2,:,J) = CH(3:IDO:2,:,J) + CH(3:IDO:2,:,JC)
               C1(3:IDO:2,:,JC) = CH(2:IDO-1:2,:,JC) - CH(2:IDO-1:2,:,J)
            END DO
            GO TO 121
         ENDIF
         DO J = 2, IPPH
            JC = IPP2 - J
            C1(2:IDO-1:2,:,J) = CH(2:IDO-1:2,:,J) + CH(2:IDO-1:2,:,JC)
            C1(2:IDO-1:2,:,JC) = CH(3:IDO:2,:,J) - CH(3:IDO:2,:,JC)
            C1(3:IDO:2,:,J) = CH(3:IDO:2,:,J) + CH(3:IDO:2,:,JC)
            C1(3:IDO:2,:,JC) = CH(2:IDO-1:2,:,JC) - CH(2:IDO-1:2,:,J)
         END DO
         GO TO 121
      ENDIF
      C2(:,1) = CH2(:,1)
  121 CONTINUE
      DO J = 2, IPPH
         JC = IPP2 - J
         C1(1,:,J) = CH(1,:,J) + CH(1,:,JC)
         C1(1,:,JC) = CH(1,:,JC) - CH(1,:,J)
      END DO
C
      AR1 = 1.
      AI1 = 0.
      DO L = 2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         CH2(:,L) = C2(:,1) + AR1*C2(:,2)
         CH2(:,LC) = AI1*C2(:,IP)
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO J = 3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            CH2(:,L) = CH2(:,L) + AR2*C2(:,J)
            CH2(:,LC) = CH2(:,LC) + AI2*C2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + C2(:,J)
      END DO
C
      IF (IDO >= L1) THEN
         CC(:,1,:) = CH(:,:,1)
      ELSE
         CC(:,1,:) = CH(:,:,1)
      ENDIF
      CC(IDO,2:(IPPH-1)*2:2,:) = TRANSPOSE(CH(1,:,2:IPPH))
      CC(1,3:IPPH*2-1:2,:) = TRANSPOSE(CH(1,:,IPP2-2:IPP2-IPPH:(-1)))
      IF (IDO == 1) RETURN 
      IF (NBD >= L1) THEN
         CC(2:IDO-1:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:,
     1      2:IPPH)+CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO
     2      -1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(IDP2-4:IDP2-1-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = 
     1      CH(2:IDO-1:2,:,2:IPPH)-CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1))
     2      ,SHAPE = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(3:IDO:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(3:IDO:2,:,2:
     1      IPPH)+CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/
     2      2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(IDP2-3:IDP2-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH
     1      (3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1))-CH(3:IDO:2,:,2:IPPH),SHAPE
     2       = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         RETURN 
      ENDIF
      CC(2:IDO-1:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:,2:
     1   IPPH)+CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2
     2   ,IPPH-1,L1/),ORDER = (/1,3,2/))
      CC(IDP2-4:IDP2-1-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH(
     1   2:IDO-1:2,:,2:IPPH)-CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE
     2    = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
      CC(3:IDO:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(3:IDO:2,:,2:IPPH)
     1   +CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2,IPPH-1
     2   ,L1/),ORDER = (/1,3,2/))
      CC(IDP2-3:IDP2-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH(3:
     1   IDO:2,:,IPP2-2:IPP2-IPPH:(-1))-CH(3:IDO:2,:,2:IPPH),SHAPE = (/(
     2   IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
      RETURN 
      END SUBROUTINE RADFG
!     this function is define in the file comf.f
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June 2004 2004    fortran 90 updates
C-----------------------------------------------------------------------
c     END
C
C     file genbun.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR)
C
C
C DIMENSION OF           A(M),B(M),C(M),Y(IDIMY,N)
C ARGUMENTS
C
C LATEST REVISION        JUNE 2004
C
C PURPOSE                THE NAME OF THIS PACKAGE IS A MNEMONIC FOR THE
C                        GENERALIZED BUNEMAN ALGORITHM.
C
C                        IT SOLVES THE REAL LINEAR SYSTEM OF EQUATIONS
C
C                        A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C                        + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
C
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.
C
C                        INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
C                        I.E., X(0,J) = X(M,J) AND X(M+1,J) = X(1,J),
C                        AND X(I,0) MAY EQUAL 0, X(I,2), OR X(I,N),
C                        AND X(I,N+1) MAY EQUAL 0, X(I,N-1), OR X(I,1)
C                        DEPENDING ON AN INPUT PARAMETER.
C
C USAGE                  CALL GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
C                                     IERROR)
C
C ARGUMENTS
C
C ON INPUT               NPEROD
C
C                          INDICATES THE VALUES THAT X(I,0) AND
C                          X(I,N+1) ARE ASSUMED TO HAVE.
C
C                          = 0  IF X(I,0) = X(I,N) AND X(I,N+1) =
C                               X(I,1).
C                          = 1  IF X(I,0) = X(I,N+1) = 0  .
C                          = 2  IF X(I,0) = 0 AND X(I,N+1) = X(I,N-1).
C                          = 3  IF X(I,0) = X(I,2) AND X(I,N+1) =
C                               X(I,N-1).
C                          = 4  IF X(I,0) = X(I,2) AND X(I,N+1) = 0.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        MPEROD
C                          = 0 IF A(1) AND C(M) ARE NOT ZERO
C                          = 1 IF A(1) = C(M) = 0
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        A,B,C
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0
C                          THE ARRAY ELEMENTS MUST NOT DEPEND UPON
C                          THE INDEX I, BUT MUST BE CONSTANT.
C                          SPECIFICALLY, THE SUBROUTINE CHECKS THE
C                          FOLLOWING CONDITION .
C
C                            A(I) = C(1)
C                            C(I) = C(1)
C                            B(I) = B(1)
C
C                          FOR I=1,2,...,M.
C
C                        IDIMY
C                          THE ROW (OR FIRST) DIMENSION OF THE
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS
C                          IN THE PROGRAM CALLING GENBUN.
C                          THIS PARAMETER IS USED TO SPECIFY THE
C                          VARIABLE DIMENSION OF Y.
C                          IDIMY MUST BE AT LEAST M.
C
C                        Y
C                          A TWO-DIMENSIONAL COMPLEX ARRAY THAT
C                          SPECIFIES THE VALUES OF THE RIGHT SIDE
C                          OF THE LINEAR SYSTEM OF EQUATIONS GIVEN
C                          ABOVE.
C                          Y MUST BE DIMENSIONED AT LEAST M*N.
C
C
C  ON OUTPUT             Y
C
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG WHICH INDICATES INVALID
C                          INPUT PARAMETERS  EXCEPT FOR NUMBER
C                          ZERO, A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR.
C                          = 1  M .LE. 2  .
C                          = 2  N .LE. 2
C                          = 3  IDIMY .LT. M
C                          = 4  NPEROD .LT. 0 OR NPEROD .GT. 4
C                          = 5  MPEROD .LT. 0 OR MPEROD .GT. 1
C                          = 6  A(I) .NE. C(1) OR C(I) .NE. C(1) OR
C                               B(I) .NE. B(1) FOR
C                               SOME I=1,2,...,M.
C                          = 7  A(1) .NE. 0 OR C(M) .NE. 0 AND
C                                 MPEROD = 1
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C
C SPECIAL CONDITONS      NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED FILES         comf.f,gnbnaux.f,fish.f
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
C                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C ALGORITHM              THE LINEAR SYSTEM IS SOLVED BY A CYCLIC
C                        REDUCTION ALGORITHM DESCRIBED IN THE
C                        REFERENCE.
C
C PORTABILITY            FORTRAN 90 --
C                        THE MACHINE DEPENDENT CONSTANT PI IS
C                        DEFINED IN FUNCTION PIMACH.
C
C REFERENCES             SWEET, R., "A CYCLIC REDUCTION ALGORITHM FOR
C                        SOLVING BLOCK TRIDIAGONAL SYSTEMS OF ARBITRARY
C                        DIMENSIONS," SIAM J. ON NUMER. ANAL., 14 (1977)
C                        PP. 706-720.
C
C ACCURACY               THIS TEST WAS PERFORMED ON a platform with
c                        64 bit floating point arithmetic.
C                        A UNIFORM RANDOM NUMBER GENERATOR WAS USED
C                        TO CREATE A SOLUTION ARRAY X FOR THE SYSTEM
C                        GIVEN IN THE 'PURPOSE' DESCRIPTION ABOVE
C                        WITH
C                          A(I) = C(I) = -0.5*B(I) = 1, I=1,2,...,M
C
C                        AND, WHEN MPEROD = 1
C
C                          A(1) = C(M) = 0
C                          A(M) = C(1) = 2.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE
C                        GIVEN SYSTEM  AND, USING DOUBLE PRECISION
C                        A RIGHT SIDE Y WAS COMPUTED.
C                        USING THIS ARRAY Y, SUBROUTINE GENBUN
C                        WAS CALLED TO PRODUCE APPROXIMATE
C                        SOLUTION Z.  THEN RELATIVE ERROR
C                          E = MAX(ABS(Z(I,J)-X(I,J)))/
C                              MAX(ABS(X(I,J)))
C                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
C                        OVER I=1,2,...,M AND J=1,...,N.
C
C                        THE VALUE OF E IS GIVEN IN THE TABLE
C                        BELOW FOR SOME TYPICAL VALUES OF M AND N.
C
C                   M (=N)    MPEROD    NPEROD        E
C                   ------    ------    ------      ------
C
C                     31        0         0         6.E-14
C                     31        1         1         4.E-13
C                     31        1         3         3.E-13
C                     32        0         0         9.E-14
C                     32        1         1         3.E-13
C                     32        1         3         1.E-13
C                     33        0         0         9.E-14
C                     33        1         1         4.E-13
C                     33        1         3         1.E-13
C                     63        0         0         1.E-13
C                     63        1         1         1.E-12
C                     63        1         3         2.E-13
C                     64        0         0         1.E-13
C                     64        1         1         1.E-12
C                     64        1         3         6.E-13
C                     65        0         0         2.E-13
C                     65        1         1         1.E-12
C                     65        1         3         4.E-13
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GENBUN(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, IERROR)

      implicit none
      TYPE(fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL  :: A(*)
      REAL  :: B(*)
      REAL  :: C(*)
      REAL  :: Y(IDIMY,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
      IERROR = 0
!     check input arguments
      IF (M <= 2) then
	 ierror = 1
	 return
      end if
      IF (N <= 2) then
	 ierror = 2
	 return
      end if
      IF (IDIMY < M) then
	 ierror = 3
	 return
      end if
      IF (NPEROD<0 .OR. NPEROD>4) then
	 ierror = 4
	 return
      end if
      IF (MPEROD<0 .OR. MPEROD>1) then
	 ierror = 5
	 return
      end if
!     compute and allocate real work space for genbun
      CALL GEN_SPACE (N, M, IRWK)
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed (e.g., if n,m are too large)
      IF (IERROR == 20) RETURN 
      call genbunn(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE GENBUN

 
      SUBROUTINE GENBUNN(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER , INTENT(INOUT) :: IERROR
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL  :: Y(IDIMY,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, MP1, IWBA, IWBB, IWBC, IWB2, IWB3, IWW1, IWW2, IWW3
     1   , IWD, IWTCOS, IWP, K, J, MP, NP, IPSTOR, IREV, MH, MHM1, MODD
     2   , MHPI, MHMI, NBY2, MSKIP
      REAL :: ALL, A1

      SAVE ALL
C-----------------------------------------------
      iError = 0

      IF (MPEROD /= 1) THEN
         DO I = 2, M
            IF (A(I) /= C(1)) GO TO 103
            IF (C(I) /= C(1)) GO TO 103
            IF (B(I) /= B(1)) GO TO 103
         END DO
         GO TO 104
      ENDIF
      IF (A(1)/=0. .OR. C(M)/=0.) IERROR = 7
      GO TO 104
  103 CONTINUE
      IERROR = 6
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
      MP1 = M + 1
      IWBA = MP1
      IWBB = IWBA + M
      IWBC = IWBB + M
      IWB2 = IWBC + M
      IWB3 = IWB2 + M
      IWW1 = IWB3 + M
      IWW2 = IWW1 + M
      IWW3 = IWW2 + M
      IWD = IWW3 + M
      IWTCOS = IWD + M
      IWP = IWTCOS + 4*N
      W(IWBA:M-1+IWBA) = -A(:M)
      W(IWBC:M-1+IWBC) = -C(:M)
      W(IWBB:M-1+IWBB) = 2. - B(:M)
      Y(:M,:N) = -Y(:M,:N)
      MP = MPEROD + 1
      NP = NPEROD + 1
      GO TO (114,107) MP
  107 CONTINUE
      GO TO (108,109,110,111,123) NP
  108 CONTINUE
      CALL POISP2 (M, N, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W(IWB2)
     1   , W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), W(IWP)
     2   )
      GO TO 112
  109 CONTINUE
      CALL POISD2 (M, N, 1, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W(
     1   IWW1), W(IWD), W(IWTCOS), W(IWP))
      GO TO 112
  110 CONTINUE
      CALL POISN2 (M, N, 1, 2, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W
     1   (IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), 
     2   W(IWP))
      GO TO 112
  111 CONTINUE
      CALL POISN2 (M, N, 1, 1, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W
     1   (IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), 
     2   W(IWP))
  112 CONTINUE
      IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD == 4) GO TO 124
  113 CONTINUE
      GO TO (127,133) MP
  114 CONTINUE
      MH = (M + 1)/2
      MHM1 = MH - 1
      MODD = 1
      IF (MH*2 == M) MODD = 2
      DO J = 1, N
         W(:MHM1) = Y(MH-1:MH-MHM1:(-1),J) - Y(MH+1:MHM1+MH,J)
         W(MH+1:MHM1+MH) = Y(MH-1:MH-MHM1:(-1),J) + Y(MH+1:MHM1+MH,J)
         W(MH) = 2.*Y(MH,J)
         GO TO (117,116) MODD
  116    CONTINUE
         W(M) = 2.*Y(M,J)
  117    CONTINUE
         Y(:M,J) = W(:M)
      END DO
      K = IWBC + MHM1 - 1
      I = IWBA + MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      SELECT CASE (MODD) 
      CASE DEFAULT
         K = IWBB + MHM1 - 1
         W(K) = W(K) - W(I-1)
         W(IWBC-1) = W(IWBC-1) + W(IWBB-1)
      CASE (2) 
         W(IWBB-1) = W(K+1)
      END SELECT
  122 CONTINUE
      GO TO 107
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
  123 CONTINUE
      IREV = 1
      NBY2 = N/2
  124 CONTINUE
      DO J = 1, NBY2
         MSKIP = N + 1 - J
         DO I = 1, M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
         END DO
      END DO
      GO TO (110,113) IREV
  127 CONTINUE
      DO J = 1, N
         W(MH-1:MH-MHM1:(-1)) = 0.5*(Y(MH+1:MHM1+MH,J)+Y(:MHM1,J))
         W(MH+1:MHM1+MH) = 0.5*(Y(MH+1:MHM1+MH,J)-Y(:MHM1,J))
         W(MH) = 0.5*Y(MH,J)
         GO TO (130,129) MODD
  129    CONTINUE
         W(M) = 0.5*Y(M,J)
  130    CONTINUE
         Y(:M,J) = W(:M)
      END DO
  133 CONTINUE
      W(1) = IPSTOR + IWP - 1
      RETURN 
      END SUBROUTINE GENBUNN


      SUBROUTINE POISD2(MR,NR,ISTAG,BA,BB,BC,Q,IDIMQ,B,W,D,TCOS,P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MR
      INTEGER , INTENT(IN) :: NR
      INTEGER , INTENT(IN) :: ISTAG
      INTEGER , INTENT(IN) :: IDIMQ
      REAL  :: BA(*)
      REAL  :: BB(*)
      REAL  :: BC(*)
      REAL , INTENT(INOUT) :: Q(IDIMQ,1)
      REAL  :: B(*)
      REAL  :: W(*)
      REAL  :: D(*)
      REAL  :: TCOS(*)
      REAL , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: M, N, JSH, IP, IPSTOR, KR, IRREG, JSTSAV, I, LR, NUN, 
     1   JST, JSP, L, NODD, J, JM1, JP1, JM2, JP2, JM3, JP3, NODDPR, IP1
     2   , KRPI, IDEG, JDEG
      REAL :: ALL, FI, T

      SAVE ALL
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION FOR DIRICHLET BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A+I.
C
      M = MR
      N = NR
      JSH = 0
      FI = 1./FLOAT(ISTAG)
      IP = -M
      IPSTOR = 0
      SELECT CASE (ISTAG) 
      CASE DEFAULT
         KR = 0
         IRREG = 1
         IF (N > 1) GO TO 106
         TCOS(1) = 0.
      CASE (2) 
         KR = 1
         JSTSAV = 1
         IRREG = 2
         IF (N > 1) GO TO 106
         TCOS(1) = -1.
      END SELECT
  103 CONTINUE
      B(:M) = Q(:M,1)
      CALL TRIX (1, 0, M, BA, BB, BC, B, TCOS, D, W)
      Q(:M,1) = B(:M)
      GO TO 183
  106 CONTINUE
      LR = 0
      P(:M) = 0.
      NUN = N
      JST = 1
      JSP = N
C
C     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
C
  108 CONTINUE
      L = 2*JST
      NODD = 2 - 2*((NUN + 1)/2) + NUN
C
C     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
C
      SELECT CASE (NODD) 
      CASE DEFAULT
         JSP = JSP - L
      CASE (1) 
         JSP = JSP - JST
         IF (IRREG /= 1) JSP = JSP - L
      END SELECT
  111 CONTINUE
      CALL COSGEN (JST, 1, 0.5, 0.0, TCOS)
      IF (L <= JSP) THEN
         DO J = L, JSP, L
            JM1 = J - JSH
            JP1 = J + JSH
            JM2 = J - JST
            JP2 = J + JST
            JM3 = JM2 - JSH
            JP3 = JP2 + JSH
            IF (JST == 1) THEN
               B(:M) = 2.*Q(:M,J)
               Q(:M,J) = Q(:M,JM2) + Q(:M,JP2)
            ELSE
               DO I = 1, M
                  T = Q(I,J) - Q(I,JM1) - Q(I,JP1) + Q(I,JM2) + Q(I,JP2)
                  B(I) = T + Q(I,J) - Q(I,JM3) - Q(I,JP3)
                  Q(I,J) = T
               END DO
            ENDIF
            CALL TRIX (JST, 0, M, BA, BB, BC, B, TCOS, D, W)
            Q(:M,J) = Q(:M,J) + B(:M)
         END DO
      ENDIF
C
C     REDUCTION FOR LAST UNKNOWN
C
      SELECT CASE (NODD) 
      CASE DEFAULT
         GO TO (152,120) IRREG
C
C     ODD NUMBER OF UNKNOWNS
C
  120    CONTINUE
         JSP = JSP + L
         J = JSP
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         JM3 = JM2 - JSH
         GO TO (123,121) ISTAG
  121    CONTINUE
         IF (JST /= 1) GO TO 123
         B(:M) = Q(:M,J)
         Q(:M,J) = 0.
         GO TO 130
  123    CONTINUE
         SELECT CASE (NODDPR) 
         CASE DEFAULT
            B(:M) = 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3)) + P(IP+1:M+IP)
     1          + Q(:M,J)
         CASE (2) 
            B(:M) = 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3)) + Q(:M,JP2) - Q(
     1         :M,JP1) + Q(:M,J)
         END SELECT
  128    CONTINUE
         Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1))
  130    CONTINUE
         CALL TRIX (JST, 0, M, BA, BB, BC, B, TCOS, D, W)
         IP = IP + M
         IPSTOR = MAX0(IPSTOR,IP + M)
         P(IP+1:M+IP) = Q(:M,J) + B(:M)
         B(:M) = Q(:M,JP2) + P(IP+1:M+IP)
         IF (LR == 0) THEN
            DO I = 1, JST
               KRPI = KR + I
               TCOS(KRPI) = TCOS(I)
            END DO
         ELSE
            CALL COSGEN (LR, JSTSAV, 0., FI, TCOS(JST+1))
            CALL MERGE (TCOS, 0, JST, JST, LR, KR)
         ENDIF
         CALL COSGEN (KR, JSTSAV, 0.0, FI, TCOS)
         CALL TRIX (KR, KR, M, BA, BB, BC, B, TCOS, D, W)
         Q(:M,J) = Q(:M,JM2) + B(:M) + P(IP+1:M+IP)
         LR = KR
         KR = KR + L
C
C     EVEN NUMBER OF UNKNOWNS
C
      CASE (2) 
         JSP = JSP + L
         J = JSP
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         JM3 = JM2 - JSH
         SELECT CASE (IRREG) 
         CASE DEFAULT
            JSTSAV = JST
            IDEG = JST
            KR = L
         CASE (2) 
            CALL COSGEN (KR, JSTSAV, 0.0, FI, TCOS)
            CALL COSGEN (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
            IDEG = KR
            KR = KR + JST
         END SELECT
  139    CONTINUE
         IF (JST == 1) THEN
            IRREG = 2
            B(:M) = Q(:M,J)
            Q(:M,J) = Q(:M,JM2)
         ELSE
            B(:M) = Q(:M,J) + 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3))
            SELECT CASE (IRREG) 
            CASE DEFAULT
               Q(:M,J) = Q(:M,JM2) + 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1))
               IRREG = 2
            CASE (2) 
               SELECT CASE (NODDPR) 
               CASE DEFAULT
                  Q(:M,J) = Q(:M,JM2) + P(IP+1:M+IP)
                  IP = IP - M
               CASE (2) 
                  Q(:M,J) = Q(:M,JM2) + Q(:M,J) - Q(:M,JM1)
               END SELECT
            END SELECT
         ENDIF
  150    CONTINUE
         CALL TRIX (IDEG, LR, M, BA, BB, BC, B, TCOS, D, W)
         Q(:M,J) = Q(:M,J) + B(:M)
      END SELECT
  152 CONTINUE
      NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN >= 2) GO TO 108
C
C     START SOLUTION.
C
      J = JSP
      B(:M) = Q(:M,J)
      SELECT CASE (IRREG) 
      CASE DEFAULT
         CALL COSGEN (JST, 1, 0.5, 0.0, TCOS)
         IDEG = JST
      CASE (2) 
         KR = LR + JST
         CALL COSGEN (KR, JSTSAV, 0.0, FI, TCOS)
         CALL COSGEN (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
         IDEG = KR
      END SELECT
  156 CONTINUE
      CALL TRIX (IDEG, LR, M, BA, BB, BC, B, TCOS, D, W)
      JM1 = J - JSH
      JP1 = J + JSH
      SELECT CASE (IRREG) 
      CASE DEFAULT
         Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1)) + B(:M)
      CASE (2) 
         SELECT CASE (NODDPR) 
         CASE DEFAULT
            Q(:M,J) = P(IP+1:M+IP) + B(:M)
            IP = IP - M
         CASE (2) 
            Q(:M,J) = Q(:M,J) - Q(:M,JM1) + B(:M)
         END SELECT
      END SELECT
  164 CONTINUE
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN > N) GO TO 183
      DO J = JST, N, L
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         IF (J <= JST) THEN
            B(:M) = Q(:M,J) + Q(:M,JP2)
         ELSE
            IF (JP2 <= N) GO TO 168
            B(:M) = Q(:M,J) + Q(:M,JM2)
            IF (JST < JSTSAV) IRREG = 1
            GO TO (170,171) IRREG
  168       CONTINUE
            B(:M) = Q(:M,J) + Q(:M,JM2) + Q(:M,JP2)
         ENDIF
  170    CONTINUE
         CALL COSGEN (JST, 1, 0.5, 0.0, TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    CONTINUE
         IF (J + L > N) LR = LR - JST
         KR = JST + LR
         CALL COSGEN (KR, JSTSAV, 0.0, FI, TCOS)
         CALL COSGEN (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL TRIX (IDEG, JDEG, M, BA, BB, BC, B, TCOS, D, W)
         IF (JST <= 1) THEN
            Q(:M,J) = B(:M)
         ELSE
            IF (JP2 > N) GO TO 177
  175       CONTINUE
            Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1)) + B(:M)
            CYCLE 
  177       CONTINUE
            GO TO (175,178) IRREG
  178       CONTINUE
            IF (J + JSH <= N) THEN
               Q(:M,J) = B(:M) + P(IP+1:M+IP)
               IP = IP - M
            ELSE
               Q(:M,J) = B(:M) + Q(:M,J) - Q(:M,JM1)
            ENDIF
         ENDIF
      END DO
      L = L/2
      GO TO 164
  183 CONTINUE
      W(1) = IPSTOR
      RETURN 
      END SUBROUTINE POISD2


      SUBROUTINE POISN2(M, N, ISTAG, MIXBND, A, BB, C, Q, IDIMQ, B, B2, 
     1   B3, W, W2, W3, D, TCOS, P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: ISTAG
      INTEGER , INTENT(IN) :: MIXBND
      INTEGER , INTENT(IN) :: IDIMQ
      REAL  :: A(*)
      REAL  :: BB(*)
      REAL  :: C(*)
      REAL , INTENT(INOUT) :: Q(IDIMQ,*)
      REAL  :: B(*)
      REAL  :: B2(*)
      REAL  :: B3(*)
      REAL  :: W(*)
      REAL  :: W2(*)
      REAL  :: W3(*)
      REAL  :: D(*)
      REAL  :: TCOS(*)
      REAL , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER , DIMENSION(4) :: K
      INTEGER :: K1, K2, K3, K4, MR, IP, IPSTOR, I2R, JR, NR, NLAST, KR
     1   , LR, I, NROD, JSTART, JSTOP, I2RBY2, J, JP1, JP2, JP3, JM1, 
     2   JM2, JM3, NRODPR, II, I1, I2, JR2, NLASTP, JSTEP
      REAL :: ALL, FISTAG, FNUM, FDEN, FI, T

      SAVE ALL
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION WITH NEUMANN BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS A-I.
C     MIXBND = 1 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTH BOUNDARIES.
C     MIXBND = 2 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTTOM AND
C     DIRICHLET CONDITION AT TOP.  (FOR THIS CASE, MUST HAVE ISTAG = 1.)
C
      EQUIVALENCE (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
      FISTAG = 3 - ISTAG
      FNUM = 1./FLOAT(ISTAG)
      FDEN = 0.5*FLOAT(ISTAG - 1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103) ISTAG
  101 CONTINUE
      Q(:MR,N) = 0.5*Q(:MR,N)
      GO TO (103,104) MIXBND
  103 CONTINUE
      IF (N <= 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 == NR) NROD = 0
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         JSTART = 1
      CASE (2) 
         JSTART = JR
         NROD = 1 - NROD
      END SELECT
  107 CONTINUE
      JSTOP = NLAST - JR
      IF (NROD == 0) JSTOP = JSTOP - I2R
      CALL COSGEN (I2R, 1, 0.5, 0.0, TCOS)
      I2RBY2 = I2R/2
      IF (JSTOP < JSTART) THEN
         J = JR
      ELSE
         DO J = JSTART, JSTOP, JR
            JP1 = J + I2RBY2
            JP2 = J + I2R
            JP3 = JP2 + I2RBY2
            JM1 = J - I2RBY2
            JM2 = J - I2R
            JM3 = JM2 - I2RBY2
            IF (J == 1) THEN
               JM1 = JP1
               JM2 = JP2
               JM3 = JP3
            ENDIF
            IF (I2R == 1) THEN
               IF (J == 1) JM2 = JP2
               B(:MR) = 2.*Q(:MR,J)
               Q(:MR,J) = Q(:MR,JM2) + Q(:MR,JP2)
            ELSE
               DO I = 1, MR
                  FI = Q(I,J)
                  Q(I,J)=Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
                  B(I) = FI + Q(I,J) - Q(I,JM3) - Q(I,JP3)
               END DO
            ENDIF
            CALL TRIX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,J) = Q(:MR,J) + B(:MR)
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
         END DO
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
         J = JSTOP + JR
      ENDIF
      NLAST = J
      JM1 = J - I2RBY2
      JM2 = J - I2R
      JM3 = JM2 - I2RBY2
      IF (NROD /= 0) THEN
C
C     ODD NUMBER OF UNKNOWNS
C
         IF (I2R == 1) THEN
            B(:MR) = FISTAG*Q(:MR,J)
            Q(:MR,J) = Q(:MR,JM2)
         ELSE
            B(:MR) = Q(:MR,J) + 0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
            IF (NRODPR == 0) THEN
               Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP)
               IP = IP - MR
            ELSE
               Q(:MR,J) = Q(:MR,J) - Q(:MR,JM1) + Q(:MR,JM2)
            ENDIF
            IF (LR /= 0) THEN
               CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(KR+1))
            ELSE
               B(:MR) = FISTAG*B(:MR)
            ENDIF
         ENDIF
         CALL COSGEN (KR, 1, 0.5, FDEN, TCOS)
         CALL TRIX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         KR = KR + I2R
      ELSE
         JP1 = J + I2RBY2
         JP2 = J + I2R
         IF (I2R == 1) THEN
            B(:MR) = Q(:MR,J)
            CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            IP = 0
            IPSTOR = MR
            SELECT CASE (ISTAG) 
            CASE DEFAULT
               P(:MR) = B(:MR)
               B(:MR) = B(:MR) + Q(:MR,N)
               TCOS(1) = 1.
               TCOS(2) = 0.
               CALL TRIX (1, 1, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,JM2) + P(:MR) + B(:MR)
               GO TO 150
            CASE (1) 
               P(:MR) = B(:MR)
               Q(:MR,J) = Q(:MR,JM2) + 2.*Q(:MR,JP2) + 3.*B(:MR)
               GO TO 150
            END SELECT
         ENDIF
         B(:MR) = Q(:MR,J) + 0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
         IF (NRODPR == 0) THEN
            B(:MR) = B(:MR) + P(IP+1:MR+IP)
         ELSE
            B(:MR) = B(:MR) + Q(:MR,JP2) - Q(:MR,JP1)
         ENDIF
         CALL TRIX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
         IP = IP + MR
         IPSTOR = MAX0(IPSTOR,IP + MR)
         P(IP+1:MR+IP) = B(:MR) + 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         B(:MR) = P(IP+1:MR+IP) + Q(:MR,JP2)
         IF (LR /= 0) THEN
            CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(I2R+1))
            CALL MERGE (TCOS, 0, I2R, I2R, LR, KR)
         ELSE
            DO I = 1, I2R
               II = KR + I
               TCOS(II) = TCOS(I)
            END DO
         ENDIF
         CALL COSGEN (KR, 1, 0.5, FDEN, TCOS)
         IF (LR == 0) THEN
            GO TO (146,145) ISTAG
         ENDIF
  145    CONTINUE
         CALL TRIX (KR, KR, MR, A, BB, C, B, TCOS, D, W)
         GO TO 148
  146    CONTINUE
         B(:MR) = FISTAG*B(:MR)
  148    CONTINUE
         Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP) + B(:MR)
  150    CONTINUE
         LR = KR
         KR = KR + JR
      ENDIF
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         NR = (NLAST - 1)/JR + 1
         IF (NR <= 3) GO TO 155
      CASE (2) 
         NR = NLAST/JR
         IF (NR <= 1) GO TO 192
      END SELECT
  154 CONTINUE
      I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
      J = 1 + JR
      JM1 = J - I2R
      JP1 = J + I2R
      JM2 = NLAST - I2R
      IF (NR /= 2) THEN
         IF (LR /= 0) GO TO 170
         IF (N == 3) THEN
C
C     CASE N = 3.
C
            GO TO (156,168) ISTAG
  156       CONTINUE
            B(:MR) = Q(:MR,2)
            TCOS(1) = 0.
            CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,2) = B(:MR)
            B(:MR) = 4.*B(:MR) + Q(:MR,1) + 2.*Q(:MR,3)
            TCOS(1) = -2.
            TCOS(2) = 2.
            I1 = 2
            I2 = 0
            CALL TRIX (I1, I2, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,2) = Q(:MR,2) + B(:MR)
            B(:MR) = Q(:MR,1) + 2.*Q(:MR,2)
            TCOS(1) = 0.
            CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,1) = B(:MR)
            JR = 1
            I2R = 0
            GO TO 194
         ENDIF
C
C     CASE N = 2**P+1
C
         GO TO (162,170) ISTAG
  162    CONTINUE
         B(:MR) = Q(:MR,J) + 0.5*Q(:MR,1) - Q(:MR,JM1) + Q(:MR,NLAST) - 
     1      Q(:MR,JM2)
         CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1)) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,NLAST) + 4.*Q(:MR,J)
         JR2 = 2*JR
         CALL COSGEN (JR, 1, 0.0, 0.0, TCOS)
         TCOS(JR+1:JR*2) = -TCOS(JR:1:(-1))
         CALL TRIX (JR2, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,J)
         CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168    CONTINUE
         B(:MR) = Q(:MR,2)
         Q(:MR,2) = 0.
         B2(:MR) = Q(:MR,3)
         B3(:MR) = Q(:MR,1)
         JR = 1
         I2R = 0
         J = 2
         GO TO 177
  170    CONTINUE
         B(:MR) = 0.5*Q(:MR,1) - Q(:MR,JM1) + Q(:MR,J)
         IF (NROD == 0) THEN
            B(:MR) = B(:MR) + P(IP+1:MR+IP)
         ELSE
            B(:MR) = B(:MR) + Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
         DO I = 1, MR
            T = 0.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
            Q(I,J) = T
            B2(I) = Q(I,NLAST) + T
            B3(I) = Q(I,1) + 2.*T
         END DO
  177    CONTINUE
         K1 = KR + 2*JR - 1
         K2 = KR + JR
         TCOS(K1+1) = -2.
         K4 = K1 + 3 - ISTAG
         CALL COSGEN (K2 + ISTAG - 2, 1, 0.0, FNUM, TCOS(K4))
         K4 = K1 + K2 + 1
         CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K4))
         CALL MERGE (TCOS, K1, K2, K1 + K2, JR - 1, 0)
         K3 = K1 + K2 + LR
         CALL COSGEN (JR, 1, 0.5, 0.0, TCOS(K3+1))
         K4 = K3 + JR + 1
         CALL COSGEN (KR, 1, 0.5, FDEN, TCOS(K4))
         CALL MERGE (TCOS, K3, JR, K3 + JR, KR, K1)
         IF (LR /= 0) THEN
            CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(K4))
            CALL MERGE (TCOS, K3, JR, K3 + JR, LR, K3 - LR)
            CALL COSGEN (KR, 1, 0.5, FDEN, TCOS(K4))
         ENDIF
         K3 = KR
         K4 = KR
         CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
         B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
         TCOS(1) = 2.
         CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,J)
         CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         IF (JR == 1) THEN
            Q(:MR,1) = B(:MR)
            GO TO 194
         ENDIF
         Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 194
      ENDIF
      IF (N == 2) THEN
C
C     CASE  N = 2
C
         B(:MR) = Q(:MR,1)
         TCOS(1) = 0.
         CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = B(:MR)
         B(:MR) = 2.*(Q(:MR,2)+B(:MR))*FISTAG
         TCOS(1) = -FISTAG
         TCOS(2) = 2.
         CALL TRIX (2, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = Q(:MR,1) + B(:MR)
         JR = 1
         I2R = 0
         GO TO 194
      ENDIF
      B3(:MR) = 0.
      B(:MR) = Q(:MR,1) + 2.*P(IP+1:MR+IP)
      Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1)
      B2(:MR) = 2.*(Q(:MR,1)+Q(:MR,NLAST))
      K1 = KR + JR - 1
      TCOS(K1+1) = -2.
      K4 = K1 + 3 - ISTAG
      CALL COSGEN (KR + ISTAG - 2, 1, 0.0, FNUM, TCOS(K4))
      K4 = K1 + KR + 1
      CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K4))
      CALL MERGE (TCOS, K1, KR, K1 + KR, JR - 1, 0)
      CALL COSGEN (KR, 1, 0.5, FDEN, TCOS(K1+1))
      K2 = KR
      K4 = K1 + K2 + 1
      CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(K4))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
      B(:MR) = B(:MR) + B2(:MR)
      TCOS(1) = 2.
      CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
      Q(:MR,1) = Q(:MR,1) + B(:MR)
      GO TO 194
  192 CONTINUE
      B(:MR) = Q(:MR,NLAST)
      GO TO 196
  194 CONTINUE
      J = NLAST - JR
      B(:MR) = Q(:MR,NLAST) + Q(:MR,J)
  196 CONTINUE
      JM2 = NLAST - I2R
      IF (JR == 1) THEN
         Q(:MR,NLAST) = 0.
      ELSE
         IF (NROD == 0) THEN
            Q(:MR,NLAST) = P(IP+1:MR+IP)
            IP = IP - MR
         ELSE
            Q(:MR,NLAST) = Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
      ENDIF
      CALL COSGEN (KR, 1, 0.5, FDEN, TCOS)
      CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(KR+1))
      IF (LR == 0) THEN
         B(:MR) = FISTAG*B(:MR)
      ENDIF
      CALL TRIX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
      Q(:MR,NLAST) = Q(:MR,NLAST) + B(:MR)
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR == 0) GO TO 222
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         JSTART = 1 + JR
      CASE (2) 
         JSTART = JR
      END SELECT
  209 CONTINUE
      KR = KR - JR
      IF (NLAST + JR <= N) THEN
         KR = KR - JR
         NLAST = NLAST + JR
         JSTOP = NLAST - JSTEP
      ELSE
         JSTOP = NLAST - JR
      ENDIF
      LR = KR - JR
      CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
      DO J = JSTART, JSTOP, JSTEP
         JM2 = J - JR
         JP2 = J + JR
         IF (J == JR) THEN
            B(:MR) = Q(:MR,J) + Q(:MR,JP2)
         ELSE
            B(:MR) = Q(:MR,J) + Q(:MR,JM2) + Q(:MR,JP2)
         ENDIF
         IF (JR == 1) THEN
            Q(:MR,J) = 0.
         ELSE
            JM1 = J - I2R
            JP1 = J + I2R
            Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         ENDIF
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
      END DO
      NROD = 1
      IF (NLAST + I2R <= N) NROD = 0
      IF (NLASTP /= NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
      W(1) = IPSTOR
      RETURN 
      END SUBROUTINE POISN2


      SUBROUTINE POISP2(M,N,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,W3,D,TCOS,P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER  :: IDIMQ
      REAL  :: A(*)
      REAL  :: BB(*)
      REAL  :: C(*)
      REAL  :: Q(IDIMQ,1)
      REAL  :: B(*)
      REAL  :: B2(*)
      REAL  :: B3(*)
      REAL  :: W(*)
      REAL  :: W2(*)
      REAL  :: W3(*)
      REAL  :: D(*)
      REAL  :: TCOS(*)
      REAL  :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MR, NR, NRM1, J, NRMJ, NRPJ, I, IPSTOR, LH
      REAL :: ALL, S, T

      SAVE ALL
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON EQUATION WITH PERIODIC BOUNDARY
C     CONDITIONS.
C
      MR = M
      NR = (N + 1)/2
      NRM1 = NR - 1
      IF (2*NR == N) THEN
C
C     EVEN NUMBER OF UNKNOWNS
C
         DO J = 1, NRM1
            NRMJ = NR - J
            NRPJ = NR + J
            DO I = 1, MR
               S = Q(I,NRMJ) - Q(I,NRPJ)
               T = Q(I,NRMJ) + Q(I,NRPJ)
               Q(I,NRMJ) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 2.*Q(:MR,NR)
         Q(:MR,N) = 2.*Q(:MR,N)
         CALL POISD2 (MR, NRM1, 1, A, BB, C, Q, IDIMQ, B, W, D, TCOS, P)
         IPSTOR = W(1)
         CALL POISN2 (MR, NR + 1, 1, 1, A, BB, C, Q(1,NR), IDIMQ, B, B2
     1      , B3, W, W2, W3, D, TCOS, P)
         IPSTOR = MAX0(IPSTOR,INT(W(1)))
         DO J = 1, NRM1
            NRMJ = NR - J
            NRPJ = NR + J
            DO I = 1, MR
               S = 0.5*(Q(I,NRPJ)+Q(I,NRMJ))
               T = 0.5*(Q(I,NRPJ)-Q(I,NRMJ))
               Q(I,NRMJ) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 0.5*Q(:MR,NR)
         Q(:MR,N) = 0.5*Q(:MR,N)
      ELSE
         DO J = 1, NRM1
            NRPJ = N + 1 - J
            DO I = 1, MR
               S = Q(I,J) - Q(I,NRPJ)
               T = Q(I,J) + Q(I,NRPJ)
               Q(I,J) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 2.*Q(:MR,NR)
         LH = NRM1/2
         DO J = 1, LH
            NRMJ = NR - J
            DO I = 1, MR
               S = Q(I,J)
               Q(I,J) = Q(I,NRMJ)
               Q(I,NRMJ) = S
            END DO
         END DO
         CALL POISD2 (MR, NRM1, 2, A, BB, C, Q, IDIMQ, B, W, D, TCOS, P)
         IPSTOR = W(1)
         CALL POISN2 (MR, NR, 2, 1, A, BB, C, Q(1,NR), IDIMQ, B, B2, B3
     1      , W, W2, W3, D, TCOS, P)
         IPSTOR = MAX0(IPSTOR,INT(W(1)))
         DO J = 1, NRM1
            NRPJ = NR + J
            DO I = 1, MR
               S = 0.5*(Q(I,NRPJ)+Q(I,J))
               T = 0.5*(Q(I,NRPJ)-Q(I,J))
               Q(I,NRPJ) = T
               Q(I,J) = S
            END DO
         END DO
         Q(:MR,NR) = 0.5*Q(:MR,NR)
         DO J = 1, LH
            NRMJ = NR - J
            DO I = 1, MR
               S = Q(I,J)
               Q(I,J) = Q(I,NRMJ)
               Q(I,NRMJ) = S
            END DO
         END DO
      ENDIF
      W(1) = IPSTOR
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE POISP2
C
C     file gnbnaux.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C PACKAGE GNBNAUX
C
C LATEST REVISION        June 2004
C
C PURPOSE                TO PROVIDE AUXILIARY ROUTINES FOR FISHPACK
C                        ENTRIES GENBUN AND POISTG.
C
C USAGE                  THERE ARE NO USER ENTRIES IN THIS PACKAGE.
C                        THE ROUTINES IN THIS PACKAGE ARE NOT INTENDED
C                        TO BE CALLED BY USERS, BUT RATHER BY ROUTINES
C                        IN PACKAGES GENBUN AND POISTG.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
C                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
c                        Revised by John Adams in June 2004 incorporating
c                        Fortran 90 features
C
C PORTABILITY            FORTRAN 90
C ********************************************************************
      SUBROUTINE COSGEN(N, IJUMP, FNUM, FDEN, A)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IJUMP
      REAL , INTENT(IN) :: FNUM
      REAL , INTENT(IN) :: FDEN
      REAL , INTENT(OUT) :: A(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K3, K4, K, K1, K5, I, K2, NP1
      REAL :: PI, DUM, PIBYN, X, Y
C-----------------------------------------------
C
C
C     THIS SUBROUTINE COMPUTES REQUIRED COSINE VALUES IN ASCENDING
C     ORDER.  WHEN IJUMP .GT. 1 THE ROUTINE COMPUTES VALUES
C
C        2*COS(J*PI/L) , J=1,2,...,L AND J .NE. 0(MOD N/IJUMP+1)
C
C     WHERE L = IJUMP*(N/IJUMP+1).
C
C
C     WHEN IJUMP = 1 IT COMPUTES
C
C            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
C
C     WHERE
C        FNUM = 0.5, FDEN = 0.0,  FOR REGULAR REDUCTION VALUES
C        FNUM = 0.0, FDEN = 1.0, FOR B-R AND C-R WHEN ISTAG = 1
C        FNUM = 0.0, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
C        FNUM = 0.5, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
C                                IN POISN2 ONLY.
C
C
      PI = 4.0*ATAN(1.0)
      IF (N /= 0) THEN
         IF (IJUMP /= 1) THEN
            K3 = N/IJUMP + 1
            K4 = K3 - 1
            PIBYN = PI/FLOAT(N + IJUMP)
            DO K = 1, IJUMP
               K1 = (K - 1)*K3
               K5 = (K - 1)*K4
               DO I = 1, K4
                  X = K1 + I
                  K2 = K5 + I
                  A(K2) = -2.*COS(X*PIBYN)
               END DO
            END DO
         ELSE
            NP1 = N + 1
            Y = PI/(FLOAT(N) + FDEN)
            DO I = 1, N
               X = FLOAT(NP1 - I) - FNUM
               A(I) = 2.*COS(X*Y)
            END DO
         ENDIF
      ENDIF
      RETURN 
      END SUBROUTINE COSGEN


      SUBROUTINE MERGE(TCOS, I1, M1, I2, M2, I3)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: M1
      INTEGER , INTENT(IN) :: I2
      INTEGER , INTENT(IN) :: M2
      INTEGER , INTENT(IN) :: I3
      REAL , INTENT(INOUT) :: TCOS(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J11, J3, J1, J2, J, L, K, M
      REAL :: X, Y
C-----------------------------------------------
C
C     THIS SUBROUTINE MERGES TWO ASCENDING STRINGS OF NUMBERS IN THE
C     ARRAY TCOS.  THE FIRST STRING IS OF LENGTH M1 AND STARTS AT
C     TCOS(I1+1).  THE SECOND STRING IS OF LENGTH M2 AND STARTS AT
C     TCOS(I2+1).  THE MERGED STRING GOES INTO TCOS(I3+1).
C
C
      J1 = 1
      J2 = 1
      J = I3
      IF (M1 == 0) GO TO 107
      IF (M2 == 0) GO TO 104
  101 CONTINUE
      J11 = J1
      J3 = MAX(M1,J11)
      DO J1 = J11, J3
         J = J + 1
         L = J1 + I1
         X = TCOS(L)
         L = J2 + I2
         Y = TCOS(L)
         IF (X - Y > 0.) GO TO 103
         TCOS(J) = X
      END DO
      GO TO 106
  103 CONTINUE
      TCOS(J) = Y
      J2 = J2 + 1
      IF (J2 <= M2) GO TO 101
      IF (J1 > M1) GO TO 109
  104 CONTINUE
      K = J - J1 + 1
      DO J = J1, M1
         M = K + J
         L = J + I1
         TCOS(M) = TCOS(L)
      END DO
      GO TO 109
  106 CONTINUE
      IF (J2 > M2) GO TO 109
  107 CONTINUE
      K = J - J2 + 1
      DO J = J2, M2
         M = K + J
         L = J + I2
         TCOS(M) = TCOS(L)
      END DO
  109 CONTINUE
      RETURN 
      END SUBROUTINE MERGE


      SUBROUTINE TRIX(IDEGBR, IDEGCR, M, A, B, C, Y, TCOS, D, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDEGBR
      INTEGER , INTENT(IN) :: IDEGCR
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL , INTENT(INOUT) :: Y(*)
      REAL , INTENT(IN) :: TCOS(*)
      REAL , INTENT(INOUT) :: D(*)
      REAL , INTENT(INOUT) :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MM1, IFB, IFC, L, LINT, K, I, IP
      REAL :: X, XX, Z
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE A SYSTEM OF LINEAR EQUATIONS WHERE THE
C     COEFFICIENT MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
C     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ).
C
      MM1 = M - 1
      IFB = IDEGBR + 1
      IFC = IDEGCR + 1
      L = IFB/IFC
      LINT = 1
      DO K = 1, IDEGBR
         X = TCOS(K)
         IF (K == L) THEN
            I = IDEGBR + LINT
            XX = X - TCOS(I)
            W(:M) = Y(:M)
            Y(:M) = XX*Y(:M)
         ENDIF
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y(1) = Y(1)*Z
         DO I = 2, MM1
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y(I) = (Y(I)-A(I)*Y(I-1))*Z
         END DO
         Z = B(M) - X - A(M)*D(MM1)
         IF (Z == 0.) THEN
            Y(M) = 0.
         ELSE
            Y(M) = (Y(M)-A(M)*Y(MM1))/Z
         ENDIF
         DO IP = 1, MM1
            Y(M-IP) = Y(M-IP) - D(M-IP)*Y(M+1-IP)
         END DO
         IF (K /= L) CYCLE 
         Y(:M) = Y(:M) + W(:M)
         LINT = LINT + 1
         L = (LINT*IFB)/IFC
      END DO
      RETURN 
      END SUBROUTINE TRIX


      SUBROUTINE TRI3(M, A, B, C, K, Y1, Y2, Y3, TCOS, D, W1, W2, W3)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: K(4)
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL , INTENT(INOUT) :: Y1(*)
      REAL , INTENT(INOUT) :: Y2(*)
      REAL , INTENT(INOUT) :: Y3(*)
      REAL , INTENT(IN) :: TCOS(*)
      REAL , INTENT(INOUT) :: D(*)
      REAL , INTENT(INOUT) :: W1(*)
      REAL , INTENT(INOUT) :: W2(*)
      REAL , INTENT(INOUT) :: W3(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MM1, K1, K2, K3, K4, IF1, IF2, IF3, IF4, K2K3K4, L1, L2
     1   , L3, LINT1, LINT2, LINT3, KINT1, KINT2, KINT3, N, I, IP
      REAL :: X, Z, XX
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE THREE LINEAR SYSTEMS WHOSE COMMON COEFFICIENT
C     MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
C
C                  TRIDIAGONAL (...,A(I),B(I),C(I),...)
C
      MM1 = M - 1
      K1 = K(1)
      K2 = K(2)
      K3 = K(3)
      K4 = K(4)
      IF1 = K1 + 1
      IF2 = K2 + 1
      IF3 = K3 + 1
      IF4 = K4 + 1
      K2K3K4 = K2 + K3 + K4
      IF (K2K3K4 /= 0) THEN
         L1 = IF1/IF2
         L2 = IF1/IF3
         L3 = IF1/IF4
         LINT1 = 1
         LINT2 = 1
         LINT3 = 1
         KINT1 = K1
         KINT2 = KINT1 + K2
         KINT3 = KINT2 + K3
      ENDIF
      DO N = 1, K1
         X = TCOS(N)
         IF (K2K3K4 /= 0) THEN
            IF (N == L1) THEN
               W1(:M) = Y1(:M)
            ENDIF
            IF (N == L2) THEN
               W2(:M) = Y2(:M)
            ENDIF
            IF (N == L3) THEN
               W3(:M) = Y3(:M)
            ENDIF
         ENDIF
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y1(1) = Y1(1)*Z
         Y2(1) = Y2(1)*Z
         Y3(1) = Y3(1)*Z
         DO I = 2, M
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y1(I) = (Y1(I)-A(I)*Y1(I-1))*Z
            Y2(I) = (Y2(I)-A(I)*Y2(I-1))*Z
            Y3(I) = (Y3(I)-A(I)*Y3(I-1))*Z
         END DO
         DO IP = 1, MM1
            Y1(M-IP) = Y1(M-IP) - D(M-IP)*Y1(M+1-IP)
            Y2(M-IP) = Y2(M-IP) - D(M-IP)*Y2(M+1-IP)
            Y3(M-IP) = Y3(M-IP) - D(M-IP)*Y3(M+1-IP)
         END DO
         IF (K2K3K4 == 0) CYCLE 
         IF (N == L1) THEN
            I = LINT1 + KINT1
            XX = X - TCOS(I)
            Y1(:M) = XX*Y1(:M) + W1(:M)
            LINT1 = LINT1 + 1
            L1 = (LINT1*IF1)/IF2
         ENDIF
         IF (N == L2) THEN
            I = LINT2 + KINT2
            XX = X - TCOS(I)
            Y2(:M) = XX*Y2(:M) + W2(:M)
            LINT2 = LINT2 + 1
            L2 = (LINT2*IF1)/IF3
         ENDIF
         IF (N /= L3) CYCLE 
         I = LINT3 + KINT3
         XX = X - TCOS(I)
         Y3(:M) = XX*Y3(:M) + W3(:M)
         LINT3 = LINT3 + 1
         L3 = (LINT3*IF1)/IF4
      END DO
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C OCTOBER   1980    CHANGED SEVERAL DIVIDES OF FLOATING INTEGERS
C                   TO INTEGER DIVIDES TO ACCOMODATE CRAY-1 ARITHMETIC.
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
C-----------------------------------------------------------------------
      END SUBROUTINE TRI3
C
C     file hstcrt.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                 SOLVES THE STANDARD FIVE-POINT FINITE
C                         DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C                         EQUATION
C                           (D/DX)(DU/DX) + (D/DY)(DU/DY) + LAMBDA*U
C                           = F(X,Y)
C                         ON A STAGGERED GRID IN CARTESIAN COORDINATES.
C
C USAGE                   CALL HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D
C                                      N,NBDCND,BDC,BDD,ELMBDA,
C                                      F,IDIMF,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT
C
C                        A,B
C                          THE RANGE OF X, I.E. A .LE. X .LE. B.
C                          A MUST BE LESS THAN B.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE
C                          INTERVAL (A,B).  THE GRID POINTS
C                          IN THE X-DIRECTION ARE GIVEN BY
C                          X(I) = A + (I-0.5)DX FOR I=1,2,...,M
C                          WHERE DX =(B-A)/M.  M MUST BE GREATER
C                          THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT X = A AND X = B.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN X,
C                               U(M+I,J) = U(I,J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND X = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO X
C                               IS SPECIFIED AT X = B.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED
C                               AT X = A  AND X = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED
C                               AT X = A  AND THE SOLUTION IS
C                               SPECIFIED AT X = B.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
C                          THAT SPECIFIES THE BOUNDARY VALUES
C                          (IF ANY) OF THE SOLUTION AT X = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,Y(J)) ,         J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DX)U(A,Y(J)) ,   J=1,2,...,N.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
C                          THAT SPECIFIES THE BOUNDARY VALUES
C                          OF THE SOLUTION AT X = B.
C
C                          WHEN MBDCND = 1 OR 4
C                            BDB(J) = U(B,Y(J)) ,        J=1,2,...,N.
C
C                          WHEN MBDCND = 2 OR 3
C                            BDB(J) = (D/DX)U(B,Y(J)) ,  J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF Y, I.E. C .LE. Y .LE. D.
C                          C MUST BE LESS THAN D.
C
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE Y-DIRECTION
C                          ARE GIVEN BY Y(J) = C + (J-0.5)DY,
C                          J=1,2,...,N, WHERE DY = (D-C)/N.
C                          N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Y = C   AND Y = D.
C
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
C                               U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT Y = C
C                               AND Y = D.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT Y = C
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = D.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND Y = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND THE SOLUTION IS SPECIFIED
C                               AT Y = D.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Y = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(X(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DY)U(X(I),C),   I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Y = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(X(I),D) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DY)U(X(I),D) ,  I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
C                          A SOLUTION MAY NOT EXIST. HOWEVER,
C                          HSTCRT WILL  ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE
C                          HELMHOLTZ EQUATION.  FOR I=1,2,...,M
C                          AND J=1,2,...,N
C
C                            F(I,J) = F(X(I),Y(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTCRT.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (X(I),Y(J)) FOR  I=1,2,...,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  HSTCRT
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION.  THIS SOLUTION PLUS ANY
C                          CONSTANT IS ALSO A SOLUTION; HENCE, THE
C                          SOLUTION IS NOT UNIQUE.  THE VALUE OF
C                          PERTRB SHOULD BE SMALL COMPARED TO THE
C                          RIGHT SIDE F.  OTHERWISE, A SOLUTION IS
C                          OBTAINED TO AN ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT TO NUMBERS 0 AND  6,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .GE. B
C
C                          =  2  MBDCND .LT. 0 OR MBDCND .GT. 4
C
C                          =  3  C .GE. D
C
C                          =  4  N .LE. 2
C
C                         =  5  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                         =  6  LAMBDA .GT. 0
C
C                         =  7  IDIMF .LT. M
C
C                         =  8  M .LE. 2
C
C                         SINCE THIS IS THE ONLY MEANS OF INDICATING
C                         A POSSIBLY INCORRECT CALL TO HSTCRT, THE
C                         USER SHOULD TEST IERROR AFTER THE CALL.
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
C                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR
C                        AND CALLS EITHER POISTG OR GENBUN WHICH SOLVES
C                        THE LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A
C                        LOSS OF NO MORE THAN FOUR SIGNIFICANT DIGITS
C                        FOR N AND M AS LARGE AS 64.  MORE DETAILED
C                        INFORMATION ABOUT ACCURACY CAN BE FOUND IN
C                        THE DOCUMENTATION FOR PACKAGE POISTG WHICH
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET,"A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
C                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
C                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
C                        PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTCRT(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
C
C     CHECK FOR INVALID PARAMETERS.
C
      IERROR = 0
      IF (A >= B) IERROR = 1
      IF (MBDCND<0 .OR. MBDCND>4) IERROR = 2
      IF (C >= D) IERROR = 3
      IF (N <= 2) IERROR = 4
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 5
      IF (IDIMF < M) IERROR = 7
      IF (M <= 2) IERROR = 8
      IF (IERROR /= 0) RETURN 
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hstcrtt(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     +             elmbda,f,idimf,pertrb,ierror,w%rew)
      RETURN 
      END SUBROUTINE HSTCRT


 
      SUBROUTINE HSTCRTT(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER  :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NPEROD, MPEROD, NP, MP, ID2, ID3, ID4, I, J, IERR1
      REAL::DELTAX,TWDELX,DELXSQ,DELTAY,TWDELY,DELYSQ,TWDYSQ,S,ST2
C-----------------------------------------------
 
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND > 0) MPEROD = 1
      DELTAX = (B - A)/FLOAT(M)
      TWDELX = 1./DELTAX
      DELXSQ = 2./DELTAX**2
      DELTAY = (D - C)/FLOAT(N)
      TWDELY = 1./DELTAY
      DELYSQ = DELTAY**2
      TWDYSQ = 2./DELYSQ
      NP = NBDCND + 1
      MP = MBDCND + 1
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = M
      ID3 = ID2 + M
      ID4 = ID3 + M
      S = (DELTAY/DELTAX)**2
      ST2 = 2.*S
      W(:M) = S
      W(ID2+1:M+ID2) = (-ST2) + ELMBDA*DELYSQ
      W(ID3+1:M+ID3) = S
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (111,102,102,104,104) MP
  102 CONTINUE
      F(1,:N) = F(1,:N) - BDA(:N)*DELXSQ
      W(ID2+1) = W(ID2+1) - W(1)
      GO TO 106
  104 CONTINUE
      F(1,:N) = F(1,:N) + BDA(:N)*TWDELX
      W(ID2+1) = W(ID2+1) + W(1)
  106 CONTINUE
      GO TO (111,107,109,109,107) MP
  107 CONTINUE
      F(M,:N) = F(M,:N) - BDB(:N)*DELXSQ
      W(ID3) = W(ID3) - W(1)
      GO TO 111
  109 CONTINUE
      F(M,:N) = F(M,:N) - BDB(:N)*TWDELX
      W(ID3) = W(ID3) + W(1)
  111 CONTINUE
      GO TO (121,112,112,114,114) NP
  112 CONTINUE
      F(:M,1) = F(:M,1) - BDC(:M)*TWDYSQ
      GO TO 116
  114 CONTINUE
      F(:M,1) = F(:M,1) + BDC(:M)*TWDELY
  116 CONTINUE
      GO TO (121,117,119,119,117) NP
  117 CONTINUE
      F(:M,N) = F(:M,N) - BDD(:M)*TWDYSQ
      GO TO 121
  119 CONTINUE
      F(:M,N) = F(:M,N) - BDD(:M)*TWDELY
  121 CONTINUE
      F(:M,:N) = F(:M,:N)*DELYSQ
      IF (MPEROD /= 0) THEN
         W(1) = 0.
         W(ID4) = 0.
      ENDIF
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 6
         ELSE
            GO TO (127,133,133,127,133) MP
  127       CONTINUE
            GO TO (128,133,133,128,133) NP
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
  128       CONTINUE
            S = 0.
            DO J = 1, N
               S = S + SUM(F(:M,J))
            END DO
            PERTRB = S/FLOAT(M*N)
            F(:M,:N) = F(:M,:N) - PERTRB
            PERTRB = PERTRB/DELYSQ
C
C     SOLVE THE EQUATION.
C
         ENDIF
      ENDIF
  133 CONTINUE
      IERR1 = 0
      IF (NPEROD /= 0) THEN
         CALL POISTGG (NPEROD, N, MPEROD, M, W(1), W(ID2+1), W(ID3+1), 
     1      IDIMF, F, IERR1, W(ID4+1))
      ELSE
         CALL GENBUNN (NPEROD, N, MPEROD, M, W(1), W(ID2+1), W(ID3+1), 
     1      IDIMF, F, IERR1, W(ID4+1))
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HSTCRTT
C
C     file hstcsp.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HSTCSP (INTL,A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,
C    +                   BDD,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED
C                        GRID TO THE MODIFIED HELMHOLTZ EQUATION IN
C                        SPHERICAL COORDINATES ASSUMING AXISYMMETRY
C                        (NO DEPENDENCE ON LONGITUDE).
C
C                        THE EQUATION IS
C
C                           (1/R**2)(D/DR)(R**2(DU/DR)) +
C                           1/(R**2*SIN(THETA))(D/DTHETA)
C                           (SIN(THETA)(DU/DTHETA)) +
C                           (LAMBDA/(R*SIN(THETA))**2)U  =  F(THETA,R)
C
C                        WHERE THETA IS COLATITUDE AND R IS THE
C                        RADIAL COORDINATE. THIS TWO-DIMENSIONAL
C                        MODIFIED HELMHOLTZ EQUATION RESULTS FROM
C                        THE FOURIER TRANSFORM OF THE THREE-
C                        DIMENSIONAL POISSON EQUATION.
C
C
C USAGE                  CALL HSTCSP (INTL,A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR,W)
C
C ARGUMENTS
C  ON INPUT              INTL
C
C                          = 0  ON INITIAL ENTRY TO HSTCSP OR IF ANY
C                               OF THE ARGUMENTS C, D, N, OR NBDCND
C                               ARE CHANGED FROM A PREVIOUS CALL
C
C                          = 1  IF C, D, N, AND NBDCND ARE ALL
C                               UNCHANGED FROM PREVIOUS CALL TO HSTCSP
C
C                          NOTE:
C                          A CALL WITH INTL = 0 TAKES APPROXIMATELY
C                          1.5 TIMES AS MUCH TIME AS A CALL WITH
C                          INTL = 1.  ONCE A CALL WITH INTL = 0
C                          HAS BEEN MADE THEN SUBSEQUENT SOLUTIONS
C                          CORRESPONDING TO DIFFERENT F, BDA, BDB,
C                          BDC, AND BDD CAN BE OBTAINED FASTER WITH
C                          INTL = 1 SINCE INITIALIZATION IS NOT
C                          REPEATED.
C
C                        A,B
C                          THE RANGE OF THETA (COLATITUDE),
C                          I.E. A .LE. THETA .LE. B.  A
C                          MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.  A AND B ARE IN RADIANS.
C                          A = 0 CORRESPONDS TO THE NORTH POLE AND
C                          B = PI CORRESPONDS TO THE SOUTH POLE.
C
C                          * * *  IMPORTANT  * * *
C
C                          IF B IS EQUAL TO PI, THEN B MUST BE
C                          COMPUTED USING THE STATEMENT
C                              B = PIMACH(DUM)
C                          THIS INSURES THAT B IN THE USER'S PROGRAM
C                          IS EQUAL TO PI IN THIS PROGRAM, PERMITTING
C                          SEVERAL TESTS OF THE INPUT PARAMETERS THAT
C                          OTHERWISE WOULD NOT BE POSSIBLE.
C
C                          * * * * * * * * * * * *
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE THETA-
C                          DIRECTION ARE GIVEN BY
C                            THETA(I) = A + (I-0.5)DTHETA
C                          FOR I=1,2,...,M WHERE DTHETA =(B-A)/M.
C                          M MUST BE GREATER THAN 4.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = A AND THETA = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THETA = B.
C                               (SEE NOTES 1, 2 BELOW)
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = B
C                               (SEE NOTES 1, 2 BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = A (SEE NOTES 1, 2 BELOW)
C                               AND THETA = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED AT
C                               THETA = A (SEE NOTES 1, 2 BELOW) AND
C                               THE SOLUTION IS SPECIFIED AT THETA = B.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT THETA = B.
C                               (SEE NOTE 2 BELOW)
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = B
C                               (SEE NOTE 2 BELOW).
C
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THE SOLUTION IS
C                               UNSPECIFIED AT THETA = B = PI.
C
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED AT
C                               THETA = A (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS UNSPECIFIED AT
C                               THETA = B = PI.
C
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                                THETA = A = 0 AND THETA = B = PI.
C
C                          NOTE 1:
C                          IF A = 0, DO NOT USE MBDCND = 1,2,3,4,7
C                          OR 8, BUT INSTEAD USE MBDCND = 5, 6, OR 9.
C
C                          NOTE 2:
C                          IF B = PI, DO NOT USE MBDCND = 1,2,3,4,5,
C                          OR 6, BUT INSTEAD USE MBDCND = 7, 8, OR 9.
C
C                          NOTE 3:
C                          WHEN A = 0  AND/OR B = PI THE ONLY
C                          MEANINGFUL BOUNDARY CONDITION IS
C                          DU/DTHETA = 0.   SEE D. GREENSPAN,
C                          'NUMERICAL ANALYSIS OF ELLIPTIC
C                           BOUNDARY VALUE PROBLEMS,'
C                          HARPER AND ROW, 1965, CHAPTER 5.)
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
C                          THE SOLUTION AT THETA = A.
C
C                          WHEN  MBDCND = 1, 2, OR 7,
C                            BDA(J) = U(A,R(J)),   J=1,2,...,N.
C
C                          WHEN MBDCND = 3, 4, OR 8,
C                            BDA(J) = (D/DTHETA)U(A,R(J)), J=1,2,...,N.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS A
C                          DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = B.
C
C                          WHEN MBDCND = 1, 4, OR 5,
C                            BDB(J) = U(B,R(J)),     J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DTHETA)U(B,R(J)), J=1,2,...,N.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF R , I.E. C .LE. R .LE. D.
C                          C MUST BE LESS THAN D AND NON-NEGATIVE.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE R-DIRECTION
C                          ARE GIVEN BY R(J) = C + (J-0.5)DR,
C                          J=1,2,...,N, WHERE DR = (D-C)/N.
C                          N MUST BE GREATER THAN 4.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = C AND R = D.
C
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = C AND R = D.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS
C                               SPECIFIED AT R = D. (SEE NOTE 1 BELOW)
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = C AND R = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS
C                               SPECIFIED AT R = C AND THE SOLUTION
C                               IS SPECIFIED AT R = D.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = C = 0 (SEE NOTE 2 BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT R = D.
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = C = 0 (SEE NOTE 2 BELOW)
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = D.
C
C                          NOTE 1:
C                          IF C = 0 AND MBDCND = 3,6,8 OR 9, THE
C                          SYSTEM OF EQUATIONS TO BE SOLVED IS
C                          SINGULAR.  THE UNIQUE SOLUTION IS
C                          DETERMINED BY EXTRAPOLATION TO THE
C                          SPECIFICATION OF U(THETA(1),C).
C                          BUT IN THESE CASES THE RIGHT SIDE OF THE
C                          SYSTEM WILL BE PERTURBED BY THE CONSTANT
C                          PERTRB.
C
C                          NOTE 2:
C                          NBDCND = 5 OR 6 CANNOT BE USED WITH
C                          MBDCND =1, 2, 4, 5, OR 7
C                          (THE FORMER INDICATES THAT THE SOLUTION IS
C                          UNSPECIFIED AT R = 0; THE LATTER INDICATES
C                          SOLUTION IS SPECIFIED).
C                          USE INSTEAD NBDCND = 1 OR 2.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = C.  WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(THETA(I),C),    I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DR)U(THETA(I),C), I=1,2,...,M.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = D.  WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(THETA(I),D) ,    I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DR)U(THETA(I),D), I=1,2,...,M.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE MODIFIED
C                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
C                          THAN 0, A SOLUTION MAY NOT EXIST.
C                          HOWEVER, HSTCSP WILL ATTEMPT TO FIND A
C                          SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE MODIFIED
C                          HELMHOLTZ EQUATION.  FOR I=1,2,...,M AND
C                          J=1,2,...,N
C
C                                F(I,J) = F(THETA(I),R(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTCSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C                        W
c                          A fortran 90 derived TYPE (fishworkspace) variable
c                          that must be declared by the user.  The first
c                          two declarative statements in the user program
c                          calling HSTCSP must be:
c
c                               USE fish
c                               TYPE (fishworkspace) :: W
c
c                          The first statement makes the fishpack module
c                          defined in the file "fish.f" available to the
c                          user program calling HSTCSP.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in BLKTRI to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (IFLG=0) call to HSTCSP.  These
c                          must be preserved on non-initial calls (INTL=1)
c                          to HSTCSP.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling HSTCSP should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          HSTCSP.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
C
C                                                                       
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),R(J)) FOR I=1,2,..,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTCSP THEN COMPUTES THIS
C                          SOLUTION, WHICH IS A LEAST SQUARES SOLUTION
C                          TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION; HENCE, THE SOLUTION IS NOT
C                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
C                          SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS. EXCEPT FOR NUMBERS 0 AND 10,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .LT. 0 OR B .GT. PI
C
C                          =  2  A .GE. B
C
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 9
C
C                          =  4  C .LT. 0
C
C                          =  5  C .GE. D
C
C                          =  6  NBDCND .LT. 1 OR NBDCND .GT. 6
C
C                          =  7  N .LT. 5
C
C                          =  8  NBDCND = 5 OR 6 AND
C                                MBDCND = 1, 2, 4, 5, OR 7
C
C                          =  9  C .GT. 0 AND NBDCND .GE. 5
C
C                          = 10  ELMBDA .GT. 0
C
C                          = 11  IDIMF .LT. M
C
C                          = 12  M .LT. 5
C
C                          = 13  A = 0 AND MBDCND =1,2,3,4,7 OR 8
C
C                          = 14  B = PI AND MBDCND .LE. 6
C
C                          = 15  A .GT. 0 AND MBDCND = 5, 6, OR 9
C
C                          = 16  B .LT. PI AND MBDCND .GE. 7
C
C                          = 17  LAMBDA .NE. 0 AND NBDCND .GE. 5
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTCSP,
C                          THE USER SHOULD TEST IERROR AFTER THE CALL.
C
C                        = 20 If the dynamic allocation of real and
C                             complex work space in the derived type
C                             (fishworkspace) variable W fails (e.g.,
c                             if N,M are too large for the platform used)
C                                                                       
C                        W
c                             The derived type (fishworkspace) variable W
c                             contains real and complex values that must not
C                             be destroyed if HSTCSP is called again with
C                             IFLG=1.
C                                                                       
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       fish.f,blktri.f,comf.f
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980. Revised by John Adams in June
C                        2004 using Fortan 90 dynamically allocated work
c                        space and derived data types to eliminate mixed
c                        mode conflicts in the earlier versions.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
C                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR
C                        AND CALLS BLKTRI WHICH SOLVES THE LINEAR
C                        SYSTEM OF EQUATIONS.
C
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT IS
C                        ROUGHLY PROPORTIONAL TO M*N*LOG2(N).  THE
C                        TIMING ALSO DEPENDS ON INPUT PARAMETER INTL.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE BLKTRI WHICH IS THE ROUTINE
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             P.N. SWARZTRAUBER, "A DIRECT METHOD FOR
C                        THE DISCRETE SOLUTION OF SEPARABLE ELLIPTIC
C                        EQUATIONS",
C                        SIAM J. NUMER. ANAL. 11(1974), PP. 1136-1150.
C
C                        U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
C                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
C                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
C                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
C                        PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTCSP(INTL, A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND
     1   , BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: INTL
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IW1, IWBM, IWCM, IWAN, IWBN, IWCN, IWSNTH, IWRSQ, K, L
     1   , NP, IRWK, ICWK, IERR1
      REAL :: PI, DUM

      SAVE IW1, IWBM, IWCM, IWAN, IWBN, IWCN, IWSNTH, IWRSQ
C-----------------------------------------------
c     USE fish
c     TYPE (fishworkspace) :: w
      PI = 4.0*ATAN(1.0)
C
C     CHECK FOR INVALID INPUT PARAMETERS
C
      IERROR = 0
      IF (A<0. .OR. B>PI) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<1 .OR. MBDCND>9) IERROR = 3
      IF (C < 0.) IERROR = 4
      IF (C >= D) IERROR = 5
      IF (NBDCND<1 .OR. NBDCND>6) IERROR = 6
      IF (N < 5) IERROR = 7
      IF ((NBDCND==5 .OR. NBDCND==6) .AND. (MBDCND==1 .OR. MBDCND==2
     1    .OR. MBDCND==4 .OR. MBDCND==5 .OR. MBDCND==7)) IERROR = 8
      IF (C>0. .AND. NBDCND>=5) IERROR = 9
      IF (IDIMF < M) IERROR = 11
      IF (M < 5) IERROR = 12
      IF(A==0..AND.MBDCND/=5.AND.MBDCND/=6.AND.MBDCND/=9)IERROR=13
      IF (B==PI .AND. MBDCND<=6) IERROR = 14
      IF(A>0..AND.(MBDCND==5.OR.MBDCND==6.OR.MBDCND==9))IERROR=15
      IF (B<PI .AND. MBDCND>=7) IERROR = 16
      IF (ELMBDA/=0. .AND. NBDCND>=5) IERROR = 17
      IF (IERROR == 0) THEN
         IF (INTL == 0) THEN
!     allocate required work space
            K = M + 1
            L = N + 1
            NP = NBDCND
!          compute blktri requirements in irwk,icwk
            CALL BLK_SPACE (N, M, IRWK, ICWK)
!     set work space indices
            IW1 = IRWK + 1
            IWBM = IW1 + M
            IWCM = IWBM + M
            IWAN = IWCM + M
            IWBN = IWAN + N
            IWCN = IWBN + N
            IWSNTH = IWCN + N
            IWRSQ = IWSNTH + M
!     allocate hstcsp required work spac
            IRWK = IWRSQ + N
            ICWK = ICWK + 3*K
            CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
            IF (IERROR == 20) RETURN 
         ENDIF
         IERR1 = 0
      CALL HSTCS1 (INTL,A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +ELMBDA,F,IDIMF,PERTRB,IERR1,w%rew(iw1),w%rew(IWBM),w%rew(IWCM),
     +w%rew(IWAN),w%rew(IWBN),w%rew(IWCN),w%rew(IWSNTH),w%rew(IWRSQ),
     +w%rew,w%cxw)
         IERROR = IERR1
      ENDIF
      RETURN 
      END SUBROUTINE HSTCSP


      SUBROUTINE HSTCS1(INTL, A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND
     1   , BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERR1, AM, BM, CM, AN, BN
     2   , CN, SNTH, RSQ, W, WC)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: INTL
      INTEGER  :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER  :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERR1
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL , INTENT(INOUT) :: SNTH(*)
      REAL , INTENT(INOUT) :: RSQ(*)
      REAL  :: W(*)
      COMPLEX  :: WC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J, ISW, NB
      REAL :: DTH, DTHSQ, DR, X, Y, A2, A1, A3
C-----------------------------------------------
      DTH = (B - A)/FLOAT(M)
      DTHSQ = DTH*DTH
      DO I = 1, M
         SNTH(I) = SIN(A + (FLOAT(I) - 0.5)*DTH)
      END DO
      DR = (D - C)/FLOAT(N)
      DO J = 1, N
         RSQ(J) = (C + (FLOAT(J) - 0.5)*DR)**2
      END DO
C
C     MULTIPLY RIGHT SIDE BY R(J)**2
C
      DO J = 1, N
         X = RSQ(J)
         F(:M,J) = X*F(:M,J)
      END DO
C
C      DEFINE COEFFICIENTS AM,BM,CM
C
      X = 1./(2.*COS(DTH/2.))
      AM(2:M) = (SNTH(:M-1)+SNTH(2:M))*X
      CM(:M-1) = AM(2:M)
      AM(1) = SIN(A)
      CM(M) = SIN(B)
      DO I = 1, M
         X = 1./SNTH(I)
         Y = X/DTHSQ
         AM(I) = AM(I)*Y
         CM(I) = CM(I)*Y
         BM(I) = ELMBDA*X*X - AM(I) - CM(I)
      END DO
C
C     DEFINE COEFFICIENTS AN,BN,CN
C
      X = C/DR
      DO J = 1, N
         AN(J) = (X + FLOAT(J - 1))**2
         CN(J) = (X + FLOAT(J))**2
         BN(J) = -(AN(J)+CN(J))
      END DO
      ISW = 1
      NB = NBDCND
      IF (C==0. .AND. NB==2) NB = 6
C
C     ENTER DATA ON THETA BOUNDARIES
C
      GO TO (108,108,110,110,112,112,108,110,112) MBDCND
  108 CONTINUE
      BM(1) = BM(1) - AM(1)
      X = 2.*AM(1)
      F(1,:N) = F(1,:N) - X*BDA(:N)
      GO TO 112
  110 CONTINUE
      BM(1) = BM(1) + AM(1)
      X = DTH*AM(1)
      F(1,:N) = F(1,:N) + X*BDA(:N)
  112 CONTINUE
      GO TO (113,115,115,113,113,115,117,117,117) MBDCND
  113 CONTINUE
      BM(M) = BM(M) - CM(M)
      X = 2.*CM(M)
      F(M,:N) = F(M,:N) - X*BDB(:N)
      GO TO 117
  115 CONTINUE
      BM(M) = BM(M) + CM(M)
      X = DTH*CM(M)
      F(M,:N) = F(M,:N) - X*BDB(:N)
  117 CONTINUE
      GO TO (118,118,120,120,122,122) NB
  118 CONTINUE
      BN(1) = BN(1) - AN(1)
      X = 2.*AN(1)
      F(:M,1) = F(:M,1) - X*BDC(:M)
      GO TO 122
  120 CONTINUE
      BN(1) = BN(1) + AN(1)
      X = DR*AN(1)
      F(:M,1) = F(:M,1) + X*BDC(:M)
  122 CONTINUE
      GO TO (123,125,125,123,123,125) NB
  123 CONTINUE
      BN(N) = BN(N) - CN(N)
      X = 2.*CN(N)
      F(:M,N) = F(:M,N) - X*BDD(:M)
      GO TO 127
  125 CONTINUE
      BN(N) = BN(N) + CN(N)
      X = DR*CN(N)
      F(:M,N) = F(:M,N) - X*BDD(:M)
  127 CONTINUE
      PERTRB = 0.
      GO TO (137,137,128,137,137,128,137,128,128) MBDCND
  128 CONTINUE
      GO TO (137,137,129,137,137,129) NB
  129 CONTINUE
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERR1 = 10
         ELSE
            ISW = 2
            DO I = 1, M
               X = 0.
               X = SUM(F(I,:N))
               PERTRB = PERTRB + X*SNTH(I)
            END DO
            X = 0.
            X = SUM(RSQ(:N))
            PERTRB = 2.*(PERTRB*SIN(DTH/2.))/(X*(COS(A) - COS(B)))
            DO J = 1, N
               X = RSQ(J)*PERTRB
               F(:M,J) = F(:M,J) - X
            END DO
         ENDIF
      ENDIF
  137 CONTINUE
      A2 = SUM(F(:M,1))
      A2 = A2/RSQ(1)
C
C     INITIALIZE BLKTRI
C
      IERR1 = 0
      IF (INTL == 0) CALL BLKTRII (0, 1, N, AN, BN, CN, 1, M, AM, BM, CM
     1   , IDIMF, F, IERR1, W, WC)
      CALL BLKTRII(1,1,N,AN,BN,CN,1,M,AM,BM,CM,IDIMF,F,IERR1,W,WC)
      IF (.NOT.(ISW/=2 .OR. C/=0. .OR. NBDCND/=2)) THEN
         A3 = 0.
         A1 = DOT_PRODUCT(SNTH(:M),F(:M,1))
         A3 = SUM(SNTH(:M))
         A1 = A1 + RSQ(1)*A2/2.
         IF(MBDCND==3)A1=A1+(SIN(B)*BDB(1)-SIN(A)*BDA(1))/(2.*(B-A))
         A1 = A1/A3
         A1 = BDC(1) - A1
         F(:M,:N) = F(:M,:N) + A1
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HSTCS1
C
C     file hstcyl.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED
C                        GRID TO THE MODIFIED HELMHOLTZ EQUATION
C                        IN CYLINDRICAL COORDINATES. THIS EQUATION
C
C                          (1/R)(D/DR)(R(DU/DR)) + (D/DZ)(DU/DZ)
C
C                          + LAMBDA*(1/R**2)*U = F(R,Z)
C
C                        IS A TWO-DIMENSIONAL MODIFIED HELMHOLTZ
C                        EQUATION RESULTING FROM THE FOURIER TRANSFORM
C                        OF A THREE-DIMENSIONAL POISSON EQUATION.
C
C USAGE                  CALL HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               A,B
C
C                          THE RANGE OF R, I.E. A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          BE NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE R-DIRECTION
C                          R-DIRECTION ARE GIVEN BY
C                          R(I) = A + (I-0.5)DR FOR I=1,2,...,M
C                          WHERE DR =(B-A)/M.
C                          M MUST BE GREATER THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
C                               (SEE NOTE BELOW) AND R = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
C                               (SEE NOTE BELOW) AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO R IS
C                               SPECIFIED AT R = B.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND R = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT R = B.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT R = B.
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS SPECIFIED
C                               AT R = B.
C
C                          NOTE:
C                          IF A = 0, DO NOT USE MBDCND = 1,2,3, OR 4,
C                          BUT INSTEAD USE MBDCND = 5 OR 6.
C                          THE RESULTING APPROXIMATION GIVES THE ONLY
C                          MEANINGFUL BOUNDARY CONDITION,
C                          I.E. DU/DR = 0.
C                          (SEE D. GREENSPAN, 'INTRODUCTORY NUMERICAL
C                          ANALYSIS OF ELLIPTIC BOUNDARY VALUE
C                          PROBLEMS,' HARPER AND ROW, 1965, CHAPTER 5.)
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY)
C                          OF THE SOLUTION AT R = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,Z(J)) ,       J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,Z(J)) ,   J=1,2,...,N.
C
C                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
C                          VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = B.
C
C                          WHEN MBDCND = 1,4,OR 5,
C                            BDB(J) = U(B,Z(J)) ,        J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,Z(J)) ,   J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF Z, I.E. C .LE. Z .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE Z-DIRECTION
C                          ARE GIVEN BY Z(J) = C + (J-0.5)DZ,
C                          J=1,2,...,N, WHERE DZ = (D-C)/N.
C                          N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Z = C  AND Z = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
C                               U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT Z = C
C                               AND Z = D.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT Z = C
C                               AND THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = D.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = C
C                               AND Z = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = C AND
C                               THE SOLUTION IS SPECIFIED AT Z = D.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Z = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(R(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DZ)U(R(I),C),    I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Z = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(R(I),D) ,       I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DZ)U(R(I),D) ,   I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE MODIFIED
C                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
C                          THAN 0, A SOLUTION MAY NOT EXIST.
C                          HOWEVER, HSTCYL WILL ATTEMPT TO FIND A
C                          SOLUTION.  LAMBDA MUST BE ZERO WHEN
C                          MBDCND = 5 OR 6.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE
C                          MODIFIED HELMHOLTZ EQUATION.
C                          FOR I=1,2,...,M   AND J=1,2,...,N
C                            F(I,J) = F(R(I),Z(J)) .
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTCYL.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M.
C
C ON OUTPUT
C
C                        F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),Z(J)) FOR  I=1,2,...,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTCYL THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION; HENCE, THE SOLUTION IS NOT
C                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
C                          SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 11,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .LT. 0
C
C                          =  2  A .GE. B
C
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6
C
C                          =  4  C .GE. D
C
C                          =  5  N .LE. 2
C
C                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                          =  7  A = 0 AND MBDCND = 1,2,3, OR 4
C
C                          =  8  A .GT. 0 AND MBDCND .GE. 5
C
C                          =  9  M .LE. 2
C
C                          = 10  IDIMF .LT. M
C
C                          = 11  LAMBDA .GT. 0
C
C                          = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTCYL, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
C                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR AND
C                        CALLS EITHER POISTG OR GENBUN WHICH SOLVES THE
C                        LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS RESULTS IN A LOSS
C                        OF NO MORE THAN FOUR SIGNIFICANT DIGITS
C                        FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE POISTG WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
C                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
C                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
C                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
C                        PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTCYL(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
C
      IERROR = 0
      IF (A < 0.) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 2) IERROR = 5
      IF (NBDCND<0 .OR. NBDCND>=5) IERROR = 6
      IF (A==0. .AND. MBDCND/=5 .AND. MBDCND/=6) IERROR = 7
      IF (A>0. .AND. MBDCND>=5) IERROR = 8
      IF (IDIMF < M) IERROR = 10
      IF (M <= 2) IERROR = 9
      IF (A==0. .AND. MBDCND>=5 .AND. ELMBDA/=0.) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     allocate real work space
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hstcyll(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +             ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HSTCYL


 
      SUBROUTINE HSTCYLL(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER  :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NP, IWB, IWC, IWR, I, J, K, LP, IERR1
      REAL :: DELTAR, DLRSQ, DELTHT, DLTHSQ, A1
C-----------------------------------------------
      DELTAR = (B - A)/FLOAT(M)
      DLRSQ = DELTAR**2
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB + M
      IWR = IWC + M
      DO I = 1, M
         J = IWR + I
         W(J) = A + (FLOAT(I) - 0.5)*DELTAR
         W(I) = (A + FLOAT(I - 1)*DELTAR)/(DLRSQ*W(J))
         K = IWC + I
         W(K) = (A + FLOAT(I)*DELTAR)/(DLRSQ*W(J))
         K = IWB + I
         W(K) = ELMBDA/W(J)**2 - 2./DLRSQ
      END DO
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (102,102,104,104,106,106) MBDCND
  102 CONTINUE
      A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1) - W(1)
      F(1,:N) = F(1,:N) - A1*BDA(:N)
      GO TO 106
  104 CONTINUE
      A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1) + W(1)
      F(1,:N) = F(1,:N) + A1*BDA(:N)
  106 CONTINUE
      GO TO (107,109,109,107,107,109) MBDCND
  107 CONTINUE
      W(IWC) = W(IWC) - W(IWR)
      A1 = 2.*W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
      GO TO 111
  109 CONTINUE
      W(IWC) = W(IWC) + W(IWR)
      A1 = DELTAR*W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  111 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (121,112,112,114,114) NP
  112 CONTINUE
      F(:M,1) = F(:M,1) - A1*BDC(:M)
      GO TO 116
  114 CONTINUE
      A1 = 1./DELTHT
      F(:M,1) = F(:M,1) + A1*BDC(:M)
  116 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (121,117,119,119,117) NP
  117 CONTINUE
      F(:M,N) = F(:M,N) - A1*BDD(:M)
      GO TO 121
  119 CONTINUE
      A1 = 1./DELTHT
      F(:M,N) = F(:M,N) - A1*BDD(:M)
  121 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 11
         ELSE
            GO TO (130,130,124,130,130,124) MBDCND
  124       CONTINUE
            GO TO (125,130,130,125,130) NP
  125       CONTINUE
            DO I = 1, M
               A1 = 0.
               A1 = SUM(F(I,:N))
               J = IWR + I
               PERTRB = PERTRB + A1*W(J)
            END DO
            PERTRB = PERTRB/(FLOAT(M*N)*0.5*(A + B))
            F(:M,:N) = F(:M,:N) - PERTRB
         ENDIF
      ENDIF
  130 CONTINUE
      W(:M) = W(:M)*DLTHSQ
      W(IWC+1:M+IWC) = W(IWC+1:M+IWC)*DLTHSQ
      W(IWB+1:M+IWB) = W(IWB+1:M+IWB)*DLTHSQ
      F(:M,:N) = F(:M,:N)*DLTHSQ
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     SOLVE THE SYSTEM OF EQUATIONS.
C
      IERR1 = 0
      IF (NBDCND /= 0) THEN
         CALL POISTGG (LP, N, 1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ELSE
         CALL GENBUNN (LP, N, 1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HSTCYLL
C
C     file hstplr.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HSTPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED
C                        GRID TO THE HELMHOLTZ EQUATION IN POLAR
C                        COORDINATES.  THE EQUATION IS
C
C                           (1/R)(D/DR)(R(DU/DR)) +
C                           (1/R**2)(D/DTHETA)(DU/DTHETA) +
C                           LAMBDA*U = F(R,THETA)
C
C USAGE                  CALL HSTPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,
C                                     IDIMF,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               A,B
C
C                          THE RANGE OF R, I.E. A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE R-DIRECTION
C                          ARE GIVEN BY R(I) = A + (I-0.5)DR FOR
C                          I=1,2,...,M WHERE DR =(B-A)/M.
C                          M MUST BE GREATER THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
C                               AND R = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT R = B.
C                               (SEE NOTE 1 BELOW)
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE 2 BELOW) AND R = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               SPECIFIED AT R = A (SEE NOTE 2 BELOW)
C                               AND THE SOLUTION IS SPECIFIED AT R = B.
C
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT R = B.
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS SPECIFIED
C                               AT R = B.
C
C                          NOTE 1:
C                          IF A = 0, MBDCND = 2, AND NBDCND = 0 OR 3,
C                          THE SYSTEM OF EQUATIONS TO BE SOLVED IS
C                          SINGULAR.  THE UNIQUE SOLUTION IS
C                          IS DETERMINED BY EXTRAPOLATION TO THE
C                          SPECIFICATION OF U(0,THETA(1)).
C                          BUT IN THIS CASE THE RIGHT SIDE OF THE
C                          SYSTEM WILL BE PERTURBED BY THE CONSTANT
C                          PERTRB.
C
C                          NOTE 2:
C                          IF A = 0, DO NOT USE MBDCND = 3 OR 4,
C                          BUT INSTEAD USE MBDCND = 1,2,5, OR 6.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
C                          THE SOLUTION AT R = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,THETA(J)) ,     J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,THETA(J)) ,
C                            J=1,2,...,N.
C
C                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
C                          VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = B.
C
C                          WHEN MBDCND = 1,4, OR 5,
C                            BDB(J) = U(B,THETA(J)) ,     J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,THETA(J)) ,
C                            J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF THETA, I.E. C .LE. THETA .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE THETA-
C                          DIRECTION ARE GIVEN BY THETA(J) = C +
C                          (J-0.5)DT,   J=1,2,...,N, WHERE
C                          DT = (D-C)/N.  N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = C  AND THETA = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN THETA,
C                               I.E. U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THETA = D
C                               (SEE NOTE BELOW).
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = C AND THETA = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = C AND THE SOLUTION IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C
C                          NOTE:
C                          WHEN NBDCND = 1, 2, OR 4, DO NOT USE
C                          MBDCND = 5 OR 6 (THE FORMER INDICATES THAT
C                          THE SOLUTION IS SPECIFIED AT R =  0; THE
C                          LATTER INDICATES THE SOLUTION IS UNSPECIFIED
C                          AT R = 0).  USE INSTEAD MBDCND = 1 OR 2.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(R(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DTHETA)U(R(I),C),
C                            I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(R(I),D) ,         I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) =(D/DTHETA)U(R(I),D), I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
C                          A SOLUTION MAY NOT EXIST.  HOWEVER, HSTPLR
C                          WILL ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION.
C
C                          FOR I=1,2,...,M AND J=1,2,...,N
C                            F(I,J) = F(R(I),THETA(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTPLR.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C
C ON OUTPUT
C
C                        F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),THETA(J)) FOR I=1,2,...,M,
C                          J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTPLR THEN COMPUTES THIS
C                          SOLUTION, WHICH IS A LEAST SQUARES SOLUTION
C                          TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION; HENCE, THE SOLUTION IS NOT
C                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
C                          SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 11,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .LT. 0
C
C                          =  2  A .GE. B
C
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6
C
C                          =  4  C .GE. D
C
C                          =  5  N .LE. 2
C
C                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                          =  7  A = 0 AND MBDCND = 3 OR 4
C
C                          =  8  A .GT. 0 AND MBDCND .GE. 5
C
C                          =  9  MBDCND .GE. 5 AND NBDCND .NE. 0 OR 3
C
C                          = 10  IDIMF .LT. M
C
C                          = 11  LAMBDA .GT. 0
C
C                          = 12  M .LE. 2
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTPLR, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED FILES         fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-
C                        DIFFERENCE EQUATIONS, INCORPORATES BOUNDARY
C                        DATA, ADJUSTS THE RIGHT SIDE WHEN THE SYSTEM
C                        IS SINGULAR AND CALLS EITHER POISTG OR GENBUN
C                        WHICH SOLVES THE LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        ROUTINE POISTG WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
C                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
C                        GRID OF ARBITRARY SIZE," J. COMP. PHYS.
C                        20(1976), PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTPLR(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
      IERROR = 0
      IF (A < 0.) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 2) IERROR = 5
      IF (NBDCND<0 .OR. NBDCND>=5) IERROR = 6
      IF (A==0. .AND. (MBDCND==3 .OR. MBDCND==4)) IERROR = 7
      IF (A>0. .AND. MBDCND>=5) IERROR = 8
      IF (MBDCND>=5 .AND. NBDCND/=0 .AND. NBDCND/=3) IERROR = 9
      IF (IDIMF < M) IERROR = 10
      IF (M <= 2) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hstplrr(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +                   ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
      RETURN 
      END SUBROUTINE HSTPLR


      SUBROUTINE HSTPLRR(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER  :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NP, ISW, MB, IWB, IWC, IWR, I, J, K, LP, IERR1
      REAL :: DELTAR, DLRSQ, DELTHT, DLTHSQ, A1, A2
C-----------------------------------------------
      DELTAR = (B - A)/FLOAT(M)
      DLRSQ = DELTAR**2
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
      ISW = 1
      MB = MBDCND
      IF (A==0. .AND. MBDCND==2) MB = 6
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB + M
      IWR = IWC + M
      DO I = 1, M
         J = IWR + I
         W(J) = A + (FLOAT(I) - 0.5)*DELTAR
         W(I) = (A + FLOAT(I - 1)*DELTAR)/DLRSQ
         K = IWC + I
         W(K) = (A + FLOAT(I)*DELTAR)/DLRSQ
         K = IWB + I
         W(K) = (ELMBDA - 2./DLRSQ)*W(J)
      END DO
      DO I = 1, M
         J = IWR + I
         A1 = W(J)
         F(I,:N) = A1*F(I,:N)
      END DO
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (104,104,106,106,108,108) MB
  104 CONTINUE
      A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1) - W(1)
      F(1,:N) = F(1,:N) - A1*BDA(:N)
      GO TO 108
  106 CONTINUE
      A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1) + W(1)
      F(1,:N) = F(1,:N) + A1*BDA(:N)
  108 CONTINUE
      GO TO (109,111,111,109,109,111) MB
  109 CONTINUE
      A1 = 2.*W(IWR)
      W(IWC) = W(IWC) - W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
      GO TO 113
  111 CONTINUE
      A1 = DELTAR*W(IWR)
      W(IWC) = W(IWC) + W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  113 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (123,114,114,116,116) NP
  114 CONTINUE
      F(:M,1) = F(:M,1) - A1*BDC(:M)/W(IWR+1:M+IWR)
      GO TO 118
  116 CONTINUE
      A1 = 1./DELTHT
      F(:M,1) = F(:M,1) + A1*BDC(:M)/W(IWR+1:M+IWR)
  118 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (123,119,121,121,119) NP
  119 CONTINUE
      F(:M,N) = F(:M,N) - A1*BDD(:M)/W(IWR+1:M+IWR)
      GO TO 123
  121 CONTINUE
      A1 = 1./DELTHT
      F(:M,N) = F(:M,N) - A1*BDD(:M)/W(IWR+1:M+IWR)
  123 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 11
         ELSE
            GO TO (133,133,126,133,133,126) MB
  126       CONTINUE
            GO TO (127,133,133,127,133) NP
  127       CONTINUE
            ISW = 2
            DO J = 1, N
               PERTRB = PERTRB + SUM(F(:M,J))
            END DO
            PERTRB = PERTRB/(FLOAT(M*N)*0.5*(A + B))
            DO I = 1, M
               J = IWR + I
               A1 = PERTRB*W(J)
               F(I,:N) = F(I,:N) - A1
            END DO
            A2 = SUM(F(1,:N))
            A2 = A2/W(IWR+1)
         ENDIF
      ENDIF
  133 CONTINUE
      DO I = 1, M
         J = IWR + I
         A1 = DLTHSQ*W(J)
         W(I) = A1*W(I)
         J = IWC + I
         W(J) = A1*W(J)
         J = IWB + I
         W(J) = A1*W(J)
         F(I,:N) = A1*F(I,:N)
      END DO
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IERR1 = 0
      IF (LP /= 0) THEN
         CALL POISTGG (LP, N, 1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ELSE
         CALL GENBUNN (LP, N, 1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ENDIF
      IF (.NOT.(A/=0. .OR. MBDCND/=2 .OR. ISW/=2)) THEN
         A1 = SUM(F(1,:N))
         A1 = (A1 - DLRSQ*A2/16.)/FLOAT(N)
         IF (NBDCND == 3) A1 = A1 + (BDD(1)-BDC(1))/(D - C)
         A1 = BDA(1) - A1
         F(:M,:N) = F(:M,:N) + A1
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HSTPLRR
C
C     file hstssp.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HSTSSP (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED GRID
C                        TO THE HELMHOLTZ EQUATION IN SPHERICAL
C                        COORDINATES AND ON THE SURFACE OF THE UNIT
C                        SPHERE (RADIUS OF 1).  THE EQUATION IS
C
C                          (1/SIN(THETA))(D/DTHETA)(SIN(THETA)
C                          (DU/DTHETA)) + (1/SIN(THETA)**2)
C                          (D/DPHI)(DU/DPHI) + LAMBDA*U = F(THETA,PHI)
C
C                        WHERE THETA IS COLATITUDE AND PHI IS
C                        LONGITUDE.
C
C USAGE                  CALL HSTSSP (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR)
C
C
C ARGUMENTS
C ON INPUT
C
C                        A,B
C                          THE RANGE OF THETA (COLATITUDE),
C                          I.E. A .LE. THETA .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.  A AND B ARE IN RADIANS.
C                          A = 0 CORRESPONDS TO THE NORTH POLE AND
C                          B = PI CORRESPONDS TO THE SOUTH POLE.
C
C
C                            * * *  IMPORTANT  * * *
C
C                          IF B IS EQUAL TO PI, THEN B MUST BE
C                          COMPUTED USING THE STATEMENT
C                            B = PIMACH(DUM)
C
C                          THIS INSURES THAT B IN THE USER"S PROGRAM
C                          IS EQUAL TO PI IN THIS PROGRAM WHICH
C                          PERMITS SEVERAL TESTS OF THE INPUT
C                          PARAMETERS THAT OTHERWISE WOULD NOT BE
C                          POSSIBLE.
C
C                            * * * * * * * * * * * *
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE THETA
C                          DIRECTION ARE GIVEN BY
C                          THETA(I) = A + (I-0.5)DTHETA
C                          FOR I=1,2,...,M WHERE DTHETA =(B-A)/M.
C                          M MUST BE GREATER THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = A AND THETA = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THETA = B.
C                               (SEE NOTE 3 BELOW)
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = B
C                               (SEE NOTES 2 AND 3 BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                             WITH RESPECT TO THETA IS SPECIFIED
C                             AT THETA = A
C                             (SEE NOTES 1, 2 BELOW) AND THETA = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = A
C                               (SEE NOTES 1 AND 2 BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT THETA = B.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT THETA = B.
C                               (SEE NOTE 3 BELOW)
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO THETA
C                               IS SPECIFIED AT THETA = B
C                               (SEE NOTE 2 BELOW).
C
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THE SOLUTION IS
C                               UNSPECIFIED AT THETA = B = PI.
C                               (SEE NOTE 3 BELOW)
C
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED AT
C                               THETA = A (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS UNSPECIFIED AT
C                               THETA = B = PI.
C
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THETA = B = PI.
C
C                          NOTE 1:
C                          IF A = 0, DO NOT USE MBDCND = 3, 4, OR 8,
C                          BUT INSTEAD USE MBDCND = 5, 6, OR 9.
C
C                          NOTE 2:
C                          IF B = PI, DO NOT USE MBDCND = 2, 3, OR 6,
C                          BUT INSTEAD USE MBDCND = 7, 8, OR 9.
C
C                          NOTE 3:
C                          WHEN THE SOLUTION IS SPECIFIED AT
C                          THETA = 0 AND/OR THETA = PI AND THE OTHER
C                          BOUNDARY CONDITIONS ARE COMBINATIONS
C                          OF UNSPECIFIED, NORMAL DERIVATIVE, OR
C                          PERIODICITY A SINGULAR SYSTEM RESULTS.
C                          THE UNIQUE SOLUTION IS DETERMINED BY
C                          EXTRAPOLATION TO THE SPECIFICATION OF THE
C                          SOLUTION AT EITHER THETA = 0 OR THETA = PI.
C                          BUT IN THESE CASES THE RIGHT SIDE OF THE
C                          SYSTEM  WILL BE PERTURBED BY THE CONSTANT
C                          PERTRB.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
C                          THE SOLUTION AT THETA = A.
C
C                          WHEN MBDCND = 1, 2, OR 7,
C                            BDA(J) = U(A,PHI(J)) ,      J=1,2,...,N.
C
C                          WHEN MBDCND = 3, 4, OR 8,
C                            BDA(J) = (D/DTHETA)U(A,PHI(J)) ,
C                            J=1,2,...,N.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE,
C                          BDA IS A DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = B.
C
C                          WHEN MBDCND = 1,4, OR 5,
C                            BDB(J) = U(B,PHI(J)) ,       J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DTHETA)U(B,PHI(J)) ,
C                            J=1,2,...,N.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF PHI (LONGITUDE),
C                          I.E. C .LE. PHI .LE. D.
C                          C MUST BE LESS THAN D.  IF D-C = 2*PI,
C                          PERIODIC BOUNDARY CONDITIONS ARE USUALLY
C                          USUALLY PRESCRIBED.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE PHI-DIRECTION
C                          ARE GIVEN BY PHI(J) = C + (J-0.5)DPHI,
C                          J=1,2,...,N, WHERE DPHI = (D-C)/N.
C                          N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT PHI = C  AND PHI = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN PHI,
C                               I.E.  U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = C AND PHI = D
C                               (SEE NOTE BELOW).
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO PHI IS
C                               SPECIFIED AT PHI = D
C                               (SEE NOTE BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = C AND PHI = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = C AND THE SOLUTION IS
C                               SPECIFIED AT PHI = D
C                               (SEE NOTE BELOW).
C
C                          NOTE:
C                          WHEN NBDCND = 1, 2, OR 4, DO NOT USE
C                          MBDCND = 5, 6, 7, 8, OR 9
C                          (THE FORMER INDICATES THAT THE SOLUTION
C                          IS SPECIFIED AT A POLE; THE LATTER
C                          INDICATES THE SOLUTION IS UNSPECIFIED).
C                          USE INSTEAD MBDCND = 1 OR 2.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT PHI = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(THETA(I),C) ,     I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DPHI)U(THETA(I),C),
C                            I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT PHI = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(THETA(I),D) ,     I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DPHI)U(THETA(I),D) ,
C                            I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
C                          A SOLUTION MAY NOT EXIST.  HOWEVER,
C                          HSTSSP WILL ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE
C                          HELMHOLTZ EQUATION.
C                          FOR I=1,2,...,M AND J=1,2,...,N
C
C                            F(I,J) = F(THETA(I),PHI(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTSSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),PHI(J)) FOR
C                          I=1,2,...,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTSSP THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION; HENCE, THE SOLUTION IS NOT
C                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
C                          SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 14,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .LT. 0 OR B .GT. PI
C
C                          =  2  A .GE. B
C
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 9
C
C                          =  4  C .GE. D
C
C                          =  5  N .LE. 2
C
C                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                          =  7  A .GT. 0 AND MBDCND = 5, 6, OR 9
C
C                          =  8  A = 0 AND MBDCND = 3, 4, OR 8
C
C                          =  9  B .LT. PI AND MBDCND .GE. 7
C
C                          = 10  B = PI AND MBDCND = 2,3, OR 6
C
C                          = 11  MBDCND .GE. 5 AND NDBCND = 1, 2, OR 4
C
C                          = 12  IDIMF .LT. M
C
C                          = 13  M .LE. 2
C
C                          = 14  LAMBDA .GT. 0
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTSSP, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED FILES         fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90.
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-
C                        DIFFERENCE EQUATIONS, INCORPORATES BOUNDARY
C                        DATA, ADJUSTS THE RIGHT SIDE WHEN THE SYSTEM
C                        IS SINGULAR AND CALLS EITHER POISTG OR GENBUN
C                        WHICH SOLVES THE LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        ROUTINE POISTG WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET,"A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
C                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
C                        GRID OF ARBITRARY SIZE," J. COMP. PHYS.
C                        20(1976), PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTSSP(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
      REAL :: PI, DUM
C-----------------------------------------------
C
      IERROR = 0
      PI = 4.0*ATAN(1.0)
      IF (A<0. .OR. B>PI) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>9) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 2) IERROR = 5
      IF (NBDCND<0 .OR. NBDCND>=5) IERROR = 6
      IF(A>0..AND.(MBDCND==5.OR.MBDCND==6.OR.MBDCND==9))IERROR=7
      IF(A==0..AND.(MBDCND==3.OR.MBDCND==4.OR.MBDCND==8))IERROR=8
      IF (B<PI .AND. MBDCND>=7) IERROR = 9
      IF(B==PI.AND.(MBDCND==2.OR.MBDCND==3.OR.MBDCND==6))IERROR=10
      IF (MBDCND>=5 .AND. (NBDCND==1 .OR. NBDCND==2 .OR. NBDCND==4)) 
     1   IERROR = 11
      IF (IDIMF < M) IERROR = 12
      IF (M <= 2) IERROR = 13
      IF (IERROR /= 0) RETURN 
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hstsspp(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +                   ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HSTSSP


 
      SUBROUTINE HSTSSPP(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER  :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,1)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::NP,ISW,JSW,MB,IWB,IWC,IWR,IWS,I,J,MM1,K,LP,IERR1,I1
      REAL :: DELTAR, DLRSQ, DELTHT, DLTHSQ, PI, A1, A2, A3
C-----------------------------------------------
 
      DELTAR = (B - A)/FLOAT(M)
      DLRSQ = DELTAR**2
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
      ISW = 1
      JSW = 1
      MB = MBDCND
      IF (ELMBDA == 0.) THEN
         GO TO (101,102,105,103,101,105,101,105,105) MBDCND
  101    CONTINUE
         IF (A/=0. .OR. B/=PI) GO TO 105
         MB = 9
         GO TO 104
  102    CONTINUE
         IF (A /= 0.) GO TO 105
         MB = 6
         GO TO 104
  103    CONTINUE
         IF (B /= PI) GO TO 105
         MB = 8
  104    CONTINUE
         JSW = 2
      ENDIF
  105 CONTINUE
      IWB = M
      IWC = IWB + M
      IWR = IWC + M
      IWS = IWR + M
      DO I = 1, M
         J = IWR + I
         W(J) = SIN(A + (FLOAT(I) - 0.5)*DELTAR)
         W(I) = SIN(A + FLOAT(I - 1)*DELTAR)/DLRSQ
      END DO
      MM1 = M - 1
      W(IWC+1:MM1+IWC) = W(2:MM1+1)
      W(IWB+1:MM1+IWB) = ELMBDA*W(IWR+1:MM1+IWR) - (W(:MM1)+W(2:MM1+1))
      W(IWR) = SIN(B)/DLRSQ
      W(IWC) = ELMBDA*W(IWS) - (W(M)+W(IWR))
      DO I = 1, M
         J = IWR + I
         A1 = W(J)
         F(I,:N) = A1*F(I,:N)
      END DO
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
      GO TO (110,110,112,112,114,114,110,112,114) MB
  110 CONTINUE
      A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1) - W(1)
      F(1,:N) = F(1,:N) - A1*BDA(:N)
      GO TO 114
  112 CONTINUE
      A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1) + W(1)
      F(1,:N) = F(1,:N) + A1*BDA(:N)
  114 CONTINUE
      GO TO (115,117,117,115,115,117,119,119,119) MB
  115 CONTINUE
      A1 = 2.*W(IWR)
      W(IWC) = W(IWC) - W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
      GO TO 119
  117 CONTINUE
      A1 = DELTAR*W(IWR)
      W(IWC) = W(IWC) + W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
C
C     ENTER BOUNDARY DATA FOR PHI-BOUNDARIES.
C
  119 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (129,120,120,122,122) NP
  120 CONTINUE
      F(:M,1) = F(:M,1) - A1*BDC(:M)/W(IWR+1:M+IWR)
      GO TO 124
  122 CONTINUE
      A1 = 1./DELTHT
      F(:M,1) = F(:M,1) + A1*BDC(:M)/W(IWR+1:M+IWR)
  124 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (129,125,127,127,125) NP
  125 CONTINUE
      F(:M,N) = F(:M,N) - A1*BDD(:M)/W(IWR+1:M+IWR)
      GO TO 129
  127 CONTINUE
      A1 = 1./DELTHT
      F(:M,N) = F(:M,N) - A1*BDD(:M)/W(IWR+1:M+IWR)
  129 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 14
         ELSE
            GO TO (139,139,132,139,139,132,139,132,132) MB
  132       CONTINUE
            GO TO (133,139,139,133,139) NP
  133       CONTINUE
            ISW = 2
            DO J = 1, N
               PERTRB = PERTRB + SUM(F(:M,J))
            END DO
            A1 = FLOAT(N)*(COS(A) - COS(B))/(2.*SIN(0.5*DELTAR))
            PERTRB = PERTRB/A1
            DO I = 1, M
               J = IWR + I
               A1 = PERTRB*W(J)
               F(I,:N) = F(I,:N) - A1
            END DO
            A2 = SUM(F(1,:N))
            A3 = SUM(F(M,:N))
            A2 = A2/W(IWR+1)
            A3 = A3/W(IWS)
         ENDIF
      ENDIF
  139 CONTINUE
      DO I = 1, M
         J = IWR + I
         A1 = DLTHSQ*W(J)
         W(I) = A1*W(I)
         J = IWC + I
         W(J) = A1*W(J)
         J = IWB + I
         W(J) = A1*W(J)
         F(I,:N) = A1*F(I,:N)
      END DO
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     CALL POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IERR1 = 0
      I1 = 1
      IF (NBDCND /= 0) THEN
         CALL POISTGG (LP, N, I1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ELSE
         CALL GENBUNN (LP, N, I1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ENDIF
      IF (ISW==2 .AND. JSW==2) THEN
         IF (MB == 8) THEN
            A1 = SUM(F(M,:N))
            A1 = (A1 - DLRSQ*A3/16.)/FLOAT(N)
            IF (NBDCND == 3) A1 = A1 + (BDD(M)-BDC(M))/(D - C)
            A1 = BDB(1) - A1
         ELSE
            A1 = SUM(F(1,:N))
            A1 = (A1 - DLRSQ*A2/16.)/FLOAT(N)
            IF (NBDCND == 3) A1 = A1 + (BDD(1)-BDC(1))/(D - C)
            A1 = BDA(1) - A1
         ENDIF
         F(:M,:N) = F(:M,:N) + A1
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HSTSSPP
C
C     file hw3crt.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,
C    +                   BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF,
C    +                   MDIMF,F,PERTRB,IERROR)
C
C
C DIMENSION OF           BDXS(MDIMF,N+1),    BDXF(MDIMF,N+1),
C ARGUMENTS              BDYS(LDIMF,N+1),    BDYF(LDIMF,N+1),
C                        BDZS(LDIMF,M+1),    BDZF(LDIMF,M+1),
C                        F(LDIMF,MDIMF,N+1)
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C                        EQUATION IN CARTESIAN COORDINATES.  THIS
C                        EQUATION IS
C
C                          (D/DX)(DU/DX) + (D/DY)(DU/DY) +
C                          (D/DZ)(DU/DZ) + LAMBDA*U = F(X,Y,Z) .
C
C USAGE                  CALL HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,
C                                     MBDCND,BDYS,BDYF,ZS,ZF,N,NBDCND,
C                                     BDZS,BDZF,ELMBDA,LDIMF,MDIMF,F,
C                                     PERTRB,IERROR)
C
C ARGUMENTS
C
C ON INPUT               XS,XF
C
C                          THE RANGE OF X, I.E. XS .LE. X .LE. XF .
C                          XS MUST BE LESS THAN XF.
C
C                        L
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (XS,XF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE L+1 GRID POINTS
C                          IN THE X-DIRECTION GIVEN BY
C                          X(I) = XS+(I-1)DX FOR I=1,2,...,L+1,
C                          WHERE DX = (XF-XS)/L IS THE PANEL WIDTH.
C                          L MUST BE AT LEAST 5.
C
C                        LBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT X = XS AND X = XF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN X,
C                               I.E. U(L+I,J,K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               X = XS AND X = XF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               X = XS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO X IS
C                               SPECIFIED AT X = XF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               X = XS AND X = XF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               X = XS AND THE SOLUTION IS SPECIFIED
C                               AT X=XF.
C
C                        BDXS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE DERIVATIVE OF THE SOLUTION
C                          WITH RESPECT TO X AT X = XS.
C
C                          WHEN LBDCND = 3 OR 4,
C
C                            BDXS(J,K) = (D/DX)U(XS,Y(J),Z(K)),
C                            J=1,2,...,M+1,      K=1,2,...,N+1.
C
C                          WHEN LBDCND HAS ANY OTHER VALUE, BDXS
C                          IS A DUMMY VARIABLE. BDXS MUST BE
C                          DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                        BDXF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE DERIVATIVE OF THE SOLUTION
C                          WITH RESPECT TO X AT X = XF.
C
C                          WHEN LBDCND = 2 OR 3,
C
C                            BDXF(J,K) = (D/DX)U(XF,Y(J),Z(K)),
C                            J=1,2,...,M+1,      K=1,2,...,N+1.
C
C                          WHEN LBDCND HAS ANY OTHER VALUE, BDXF IS
C                          A DUMMY VARIABLE.  BDXF MUST BE
C                          DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                        YS,YF
C                          THE RANGE OF Y, I.E. YS .LE. Y .LE. YF.
C                          YS MUST BE LESS THAN YF.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (YS,YF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS IN
C                          THE Y-DIRECTION GIVEN BY Y(J) = YS+(J-1)DY
C                          FOR J=1,2,...,M+1,
C                          WHERE DY = (YF-YS)/M IS THE PANEL WIDTH.
C                          M MUST BE AT LEAST 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Y = YS AND Y = YF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
C                               U(I,M+J,K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Y = YS AND Y = YF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Y = YS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO Y IS
C                               SPECIFIED AT Y = YF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = YS AND Y = YF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               AT Y = YS AND THE SOLUTION IS
C                               SPECIFIED AT Y=YF.
C
C                        BDYS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Y AT Y = YS.
C
C                          WHEN MBDCND = 3 OR 4,
C
C                            BDYS(I,K) = (D/DY)U(X(I),YS,Z(K)),
C                            I=1,2,...,L+1,      K=1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDYS
C                          IS A DUMMY VARIABLE. BDYS MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(N+1).
C
C                        BDYF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Y AT Y = YF.
C
C                          WHEN MBDCND = 2 OR 3,
C
C                            BDYF(I,K) = (D/DY)U(X(I),YF,Z(K)),
C                            I=1,2,...,L+1,      K=1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDYF
C                          IS A DUMMY VARIABLE. BDYF MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(N+1).
C
C                        ZS,ZF
C                          THE RANGE OF Z, I.E. ZS .LE. Z .LE. ZF.
C                          ZS MUST BE LESS THAN ZF.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (ZS,ZF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE Z-DIRECTION GIVEN BY
C                          Z(K) = ZS+(K-1)DZ FOR K=1,2,...,N+1,
C                          WHERE DZ = (ZF-ZS)/N IS THE PANEL WIDTH.
C                          N MUST BE AT LEAST 5.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Z = ZS AND Z = ZF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
C                               U(I,J,N+K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Z = ZS AND Z = ZF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Z = ZS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO Z IS
C                               SPECIFIED AT Z = ZF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS SPECIFIED AT
C                               Z = ZS AND Z = ZF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS SPECIFIED AT
C                               Z = ZS AND THE SOLUTION IS SPECIFIED
C                               AT Z=ZF.
C
C                        BDZS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Z AT Z = ZS.
C
C                          WHEN NBDCND = 3 OR 4,
C
C                            BDZS(I,J) = (D/DZ)U(X(I),Y(J),ZS),
C                            I=1,2,...,L+1,      J=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDZS
C                          IS A DUMMY VARIABLE. BDZS MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(M+1).
C
C                        BDZF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Z AT Z = ZF.
C
C                          WHEN NBDCND = 2 OR 3,
C
C                            BDZF(I,J) = (D/DZ)U(X(I),Y(J),ZF),
C                            I=1,2,...,L+1,      J=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDZF
C                          IS A DUMMY VARIABLE. BDZF MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(M+1).
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION. IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HW3CRT WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        LDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE
C                          ARRAYS F,BDYS,BDYF,BDZS,AND BDZF AS IT
C                          APPEARS IN THE PROGRAM CALLING HW3CRT.
C                          THIS PARAMETER IS USED TO SPECIFY THE
C                          VARIABLE DIMENSION OF THESE ARRAYS.
C                          LDIMF MUST BE AT LEAST L+1.
C
C                        MDIMF
C                          THE COLUMN (OR SECOND) DIMENSION OF THE
C                          ARRAY F AND THE ROW (OR FIRST) DIMENSION
C                          OF THE ARRAYS BDXS AND BDXF AS IT APPEARS
C                          IN THE PROGRAM CALLING HW3CRT.  THIS
C                          PARAMETER IS USED TO SPECIFY THE VARIABLE
C                          DIMENSION OF THESE ARRAYS.
C                          MDIMF MUST BE AT LEAST M+1.
C
C                        F
C                          A THREE-DIMENSIONAL ARRAY OF DIMENSION AT
C                          AT LEAST (L+1)*(M+1)*(N+1), SPECIFYING THE
C                          VALUES OF THE RIGHT SIDE OF THE HELMHOLZ
C                          EQUATION AND BOUNDARY VALUES (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I=2,3,...,L,  J=2,3,...,M,
C                          AND K=2,3,...,N
C                          F(I,J,K) = F(X(I),Y(J),Z(K)).
C
C                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
C                          FOR J=1,2,...,M+1,  K=1,2,...,N+1,
C                          AND I=1,2,...,L+1
C
C                          LBDCND      F(1,J,K)         F(L+1,J,K)
C                          ------   ---------------   ---------------
C
C                            0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K))
C                            1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
C                            2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))
C                            3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))
C                            4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
C
C                          MBDCND      F(I,1,K)         F(I,M+1,K)
C                          ------   ---------------   ---------------
C
C                            0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K))
C                            1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K))
C                            2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))
C                            3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))
C                            4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K))
C
C                          NBDCND      F(I,J,1)         F(I,J,N+1)
C                          ------   ---------------   ---------------
C
C                            0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS)
C                            1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
C                            2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)
C                            3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)
C                            4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F ON A BOUNDARY,
C                          THEN THE SOLUTION MUST BE SPECIFIED.
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J,K) OF THE
C                          FINITE DIFFERENCE APPROXIMATION FOR THE
C                          GRID POINT (X(I),Y(J),Z(K)) FOR
C                          I=1,2,...,L+1, J=1,2,...,M+1,
C                          AND K=1,2,...,N+1.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  PWSCRT
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION.  THIS SOLUTION IS NOT
C                          UNIQUE AND IS UNNORMALIZED.  THE VALUE OF
C                          PERTRB SHOULD BE SMALL COMPARED TO THE
C                          THE RIGHT SIDE F.  OTHERWISE, A SOLUTION
C                          IS OBTAINED TO AN ESSENTIALLY DIFFERENT
C                          PROBLEM.  THIS COMPARISON SHOULD ALWAYS
C                          BE MADE TO INSURE THAT A MEANINGFUL
C                          SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 12,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C                          =  1  XS .GE. XF
C                          =  2  L .LT. 5
C                          =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4
C                          =  4  YS .GE. YF
C                          =  5  M .LT. 5
C                          =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4
C                          =  7  ZS .GE. ZF
C                          =  8  N .LT. 5
C                          =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4
C                          = 10  LDIMF .LT. L+1
C                          = 11  MDIMF .LT. M+1
C                          = 12  LAMBDA .GT. 0
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HW3CRT, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED Files         fish.f,pois3d.f,fftpack.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS AND
C                        THEN CALLS POIS3D TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          L*M*N*(LOG2(L)+LOG2(M)+5),
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS LBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR L, M AND N AS LARGE AS 32.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        ROUTINE POIS3D WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             NONE
C***********************************************************************
      SUBROUTINE HW3CRT(XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, MBDCND
     1   , BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA, LDIMF, 
     2   MDIMF, F, PERTRB, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: L
      INTEGER  :: LBDCND
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL  :: XS
      REAL  :: XF
      REAL  :: YS
      REAL  :: YF
      REAL  :: ZS
      REAL  :: ZF
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDXS(MDIMF,*)
      REAL  :: BDXF(MDIMF,*)
      REAL  :: BDYS(LDIMF,*)
      REAL  :: BDYF(LDIMF,*)
      REAL  :: BDZS(LDIMF,*)
      REAL  :: BDZF(LDIMF,*)
      REAL  :: F(LDIMF,MDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
C
C     CHECK FOR INVALID INPUT.
C
      IERROR = 0
      IF (XF <= XS) IERROR = 1
      IF (L < 5) IERROR = 2
      IF (LBDCND<0 .OR. LBDCND>4) IERROR = 3
      IF (YF <= YS) IERROR = 4
      IF (M < 5) IERROR = 5
      IF (MBDCND<0 .OR. MBDCND>4) IERROR = 6
      IF (ZF <= ZS) IERROR = 7
      IF (N < 5) IERROR = 8
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 9
      IF (LDIMF < L + 1) IERROR = 10
      IF (MDIMF < M + 1) IERROR = 11
c     IF (IERROR .NE. 0) GO TO 188
      IF (IERROR /= 0) RETURN 
!     allocate required work space length (generous estimate)
      IRWK=30+L+M+5*N+MAX0(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hw3crtt(XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,
     +             BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF,
     +             MDIMF,F,PERTRB,IERROR,w%rew,w%ifac,w%ifac2)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HW3CRT


 
      SUBROUTINE HW3CRTT(XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, 
     1   MBDCND, BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA, 
     2   LDIMF, MDIMF, F, PERTRB, IERROR, W, IFACX, IFACY)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: L
      INTEGER  :: LBDCND
      INTEGER , INTENT(IN) :: M
      INTEGER  :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: XS
      REAL , INTENT(IN) :: XF
      REAL , INTENT(IN) :: YS
      REAL , INTENT(IN) :: YF
      REAL , INTENT(IN) :: ZS
      REAL , INTENT(IN) :: ZF
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDXS(MDIMF,*)
      REAL , INTENT(IN) :: BDXF(MDIMF,*)
      REAL , INTENT(IN) :: BDYS(LDIMF,*)
      REAL , INTENT(IN) :: BDYF(LDIMF,*)
      REAL , INTENT(IN) :: BDZS(LDIMF,*)
      REAL , INTENT(IN) :: BDZF(LDIMF,*)
      REAL  :: F(LDIMF,MDIMF,*)
      REAL  :: W(*)
      INTEGER:: IFACX(*)
      INTEGER:: IFACY(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MSTART, MSTOP, MP1, MP, MUNK, NP, NP1, NSTART, NSTOP, 
     1   NUNK, LP1, LP, LSTART, LSTOP, J, K, LUNK, I, IWB, IWC, IWW, 
     2   MSTPM1, LSTPM1, NSTPM1, NPEROD, IR
      REAL::DY,TWBYDY,C2,DZ,TWBYDZ,C3,DX,C1,TWBYDX,XLP,YLP,ZLP,S1,S2,S
C-----------------------------------------------
 
      DY = (YF - YS)/FLOAT(M)
      TWBYDY = 2./DY
      C2 = 1./DY**2
      MSTART = 1
      MSTOP = M
      MP1 = M + 1
      MP = MBDCND + 1
      GO TO (104,101,101,102,102) MP
  101 CONTINUE
      MSTART = 2
  102 CONTINUE
      GO TO (104,104,103,103,104) MP
  103 CONTINUE
      MSTOP = MP1
  104 CONTINUE
      MUNK = MSTOP - MSTART + 1
      DZ = (ZF - ZS)/FLOAT(N)
      TWBYDZ = 2./DZ
      NP = NBDCND + 1
      C3 = 1./DZ**2
      NP1 = N + 1
      NSTART = 1
      NSTOP = N
      GO TO (108,105,105,106,106) NP
  105 CONTINUE
      NSTART = 2
  106 CONTINUE
      GO TO (108,108,107,107,108) NP
  107 CONTINUE
      NSTOP = NP1
  108 CONTINUE
      NUNK = NSTOP - NSTART + 1
      LP1 = L + 1
      DX = (XF - XS)/FLOAT(L)
      C1 = 1./DX**2
      TWBYDX = 2./DX
      LP = LBDCND + 1
      LSTART = 1
      LSTOP = L
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (122,109,109,112,112) LP
  109 CONTINUE
      LSTART = 2
      F(2,MSTART:MSTOP,NSTART:NSTOP) = F(2,MSTART:MSTOP,NSTART:NSTOP) - 
     1   C1*F(1,MSTART:MSTOP,NSTART:NSTOP)
      GO TO 115
  112 CONTINUE
      F(1,MSTART:MSTOP,NSTART:NSTOP) = F(1,MSTART:MSTOP,NSTART:NSTOP) + 
     1   TWBYDX*BDXS(MSTART:MSTOP,NSTART:NSTOP)
  115 CONTINUE
      GO TO (122,116,119,119,116) LP
  116 CONTINUE
      F(L,MSTART:MSTOP,NSTART:NSTOP) = F(L,MSTART:MSTOP,NSTART:NSTOP) - 
     1   C1*F(LP1,MSTART:MSTOP,NSTART:NSTOP)
      GO TO 122
  119 CONTINUE
      LSTOP = LP1
      F(LP1,MSTART:MSTOP,NSTART:NSTOP) = F(LP1,MSTART:MSTOP,NSTART:NSTOP
     1   ) - TWBYDX*BDXF(MSTART:MSTOP,NSTART:NSTOP)
  122 CONTINUE
      LUNK = LSTOP - LSTART + 1
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (136,123,123,126,126) MP
  123 CONTINUE
      F(LSTART:LSTOP,2,NSTART:NSTOP) = F(LSTART:LSTOP,2,NSTART:NSTOP) - 
     1   C2*F(LSTART:LSTOP,1,NSTART:NSTOP)
      GO TO 129
  126 CONTINUE
      F(LSTART:LSTOP,1,NSTART:NSTOP) = F(LSTART:LSTOP,1,NSTART:NSTOP) + 
     1   TWBYDY*BDYS(LSTART:LSTOP,NSTART:NSTOP)
  129 CONTINUE
      GO TO (136,130,133,133,130) MP
  130 CONTINUE
      F(LSTART:LSTOP,M,NSTART:NSTOP) = F(LSTART:LSTOP,M,NSTART:NSTOP) - 
     1   C2*F(LSTART:LSTOP,MP1,NSTART:NSTOP)
      GO TO 136
  133 CONTINUE
      F(LSTART:LSTOP,MP1,NSTART:NSTOP) = F(LSTART:LSTOP,MP1,NSTART:NSTOP
     1   ) - TWBYDY*BDYF(LSTART:LSTOP,NSTART:NSTOP)
  136 CONTINUE
      GO TO (150,137,137,140,140) NP
  137 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,2) = F(LSTART:LSTOP,MSTART:MSTOP,2) - 
     1   C3*F(LSTART:LSTOP,MSTART:MSTOP,1)
      GO TO 143
  140 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,1) = F(LSTART:LSTOP,MSTART:MSTOP,1) + 
     1   TWBYDZ*BDZS(LSTART:LSTOP,MSTART:MSTOP)
  143 CONTINUE
      GO TO (150,144,147,147,144) NP
  144 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,N) = F(LSTART:LSTOP,MSTART:MSTOP,N) - 
     1   C3*F(LSTART:LSTOP,MSTART:MSTOP,NP1)
      GO TO 150
  147 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,NP1) = F(LSTART:LSTOP,MSTART:MSTOP,NP1
     1   ) - TWBYDZ*BDZF(LSTART:LSTOP,MSTART:MSTOP)
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
  150 CONTINUE
      IWB = NUNK + 1
      IWC = IWB + NUNK
      IWW = IWC + NUNK
      W(:NUNK) = C3
      W(IWC:NUNK-1+IWC) = C3
      W(IWB:NUNK-1+IWB) = (-2.*C3) + ELMBDA
      GO TO (155,155,153,152,152) NP
  152 CONTINUE
      W(IWC) = 2.*C3
  153 CONTINUE
      GO TO (155,155,154,154,155) NP
  154 CONTINUE
      W(IWB-1) = 2.*C3
  155 CONTINUE
      PERTRB = 0.
C
C     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.
C
      GO TO (156,172,172,156,172) LP
  156 CONTINUE
      GO TO (157,172,172,157,172) MP
  157 CONTINUE
      GO TO (158,172,172,158,172) NP
  158 CONTINUE
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 12
         ELSE
            MSTPM1 = MSTOP - 1
            LSTPM1 = LSTOP - 1
            NSTPM1 = NSTOP - 1
            XLP = (2 + LP)/3
            YLP = (2 + MP)/3
            ZLP = (2 + NP)/3
            S1 = 0.
            DO K = 2, NSTPM1
               DO J = 2, MSTPM1
                  S1 = S1 + SUM(F(2:LSTPM1,J,K))
                  S1 = S1 + (F(1,J,K)+F(LSTOP,J,K))/XLP
               END DO
               S2 = SUM(F(2:LSTPM1,1,K)+F(2:LSTPM1,MSTOP,K))
               S2 = (S2 + (F(1,1,K)+F(1,MSTOP,K)+F(LSTOP,1,K)+F(LSTOP,
     1            MSTOP,K))/XLP)/YLP
               S1 = S1 + S2
            END DO
            S = (F(1,1,1)+F(LSTOP,1,1)+F(1,1,NSTOP)+F(LSTOP,1,NSTOP)+F(1
     1         ,MSTOP,1)+F(LSTOP,MSTOP,1)+F(1,MSTOP,NSTOP)+F(LSTOP,MSTOP
     2         ,NSTOP))/(XLP*YLP)
            DO J = 2, MSTPM1
               S = S + SUM(F(2:LSTPM1,J,1)+F(2:LSTPM1,J,NSTOP))
            END DO
            S2 = 0.
            S2 = SUM(F(2:LSTPM1,1,1)+F(2:LSTPM1,1,NSTOP)+F(2:LSTPM1,
     1         MSTOP,1)+F(2:LSTPM1,MSTOP,NSTOP))
            S = S2/YLP + S
            S2 = 0.
            S2 = SUM(F(1,2:MSTPM1,1)+F(1,2:MSTPM1,NSTOP)+F(LSTOP,2:
     1         MSTPM1,1)+F(LSTOP,2:MSTPM1,NSTOP))
            S = S2/XLP + S
            PERTRB = (S/ZLP + S1)/((FLOAT(LUNK + 1) - XLP)*(FLOAT(MUNK
     1          + 1) - YLP)*(FLOAT(NUNK + 1) - ZLP))
            F(:LUNK,:MUNK,:NUNK) = F(:LUNK,:MUNK,:NUNK) - PERTRB
         ENDIF
      ENDIF
  172 CONTINUE
      NPEROD = 0
      IF (NBDCND /= 0) THEN
         NPEROD = 1
         W(1) = 0.
         W(IWW-1) = 0.
      ENDIF
      CALL POIS3DD (LBDCND, LUNK, C1, MBDCND, MUNK, C2, NPEROD, NUNK, W
     1   , W(IWB), W(IWC), LDIMF, MDIMF, F(LSTART,MSTART,NSTART), IR,
     2   W(IWW), IFACX, IFACY)
C
C     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
C
      IF (LP == 1) THEN
         IF (MP == 1) THEN
            F(1,MP1,NSTART:NSTOP) = F(1,1,NSTART:NSTOP)
            MSTOP = MP1
         ENDIF
         IF (NP == 1) THEN
            F(1,MSTART:MSTOP,NP1) = F(1,MSTART:MSTOP,1)
            NSTOP = NP1
         ENDIF
         F(LP1,MSTART:MSTOP,NSTART:NSTOP) = F(1,MSTART:MSTOP,NSTART:
     1      NSTOP)
      ENDIF
      IF (MP == 1) THEN
         IF (NP == 1) THEN
            F(LSTART:LSTOP,1,NP1) = F(LSTART:LSTOP,1,1)
            NSTOP = NP1
         ENDIF
         F(LSTART:LSTOP,MP1,NSTART:NSTOP) = F(LSTART:LSTOP,1,NSTART:
     1      NSTOP)
      ENDIF
      IF (NP == 1) THEN
         F(LSTART:LSTOP,MSTART:MSTOP,NP1) = F(LSTART:LSTOP,MSTART:MSTOP,
     1      1)
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HW3CRTT
C
C     file hwscrt.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HWSCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N),      BDB(N),   BDC(M),BDD(M),
C ARGUMENTS              F(IDIMF,N)
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C                        EQUATION IN CARTESIAN COORDINATES.  THIS
C                        EQUATION IS
C
C                          (D/DX)(DU/DX) + (D/DY)(DU/DY)
C                          + LAMBDA*U = F(X,Y).
C
C USAGE                  CALL HWSCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               A,B
C
C                          THE RANGE OF X, I.E., A .LE. X .LE. B.
C                          A MUST BE LESS THAN B.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (A,B) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS
C                          IN THE X-DIRECTION GIVEN BY
C                          X(I) = A+(I-1)DX FOR I = 1,2,...,M+1,
C                          WHERE DX = (B-A)/M IS THE PANEL WIDTH.
C                          M MUST BE GREATER THAN 3.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT X = A AND X = B.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN X,
C                               I.E., U(I,J) = U(M+I,J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND X = B.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO X IS
C                               SPECIFIED AT X = B.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               AT X = A AND X = B.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               X = A AND THE SOLUTION IS SPECIFIED
C                               AT X = B.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO X AT X = A.
C
C                          WHEN MBDCND = 3 OR 4,
C
C                            BDA(J) = (D/DX)U(A,Y(J)), J = 1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO X AT X = B.
C
C                          WHEN MBDCND = 2 OR 3,
C
C                            BDB(J) = (D/DX)U(B,Y(J)), J = 1,2,...,N+1
C
C                          WHEN MBDCND HAS ANY OTHER VALUE BDB IS A
C                          DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF Y, I.E., C .LE. Y .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE
C                          Y-DIRECTION GIVEN BY Y(J) = C+(J-1)DY
C                          FOR J = 1,2,...,N+1, WHERE
C                          DY = (D-C)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 3.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS AT
C                          Y = C AND Y = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Y,
C                               I.E., U(I,J) = U(I,N+J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Y = C AND Y = D.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Y = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO Y IS
C                               SPECIFIED AT Y = D.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND Y = D.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND THE SOLUTION IS SPECIFIED
C                               AT Y = D.
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO Y AT Y = C.
C
C                          WHEN NBDCND = 3 OR 4,
C
C                            BDC(I) = (D/DY)U(X(I),C), I = 1,2,...,M+1
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO Y AT Y = D.
C
C                          WHEN NBDCND = 2 OR 3,
C
C                            BDD(I) = (D/DY)U(X(I),D), I = 1,2,...,M+1
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSCRT WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES OF THE
C                          RIGHT SIDE OF THE HELMHOLTZ  EQUATION AND
C                          BOUNDARY VALUES (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(X(I),Y(J)).
C
C                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
C                          FOR J=1,2,...,N+1,  I=1,2,...,M+1,
C
C                          MBDCND     F(1,J)        F(M+1,J)
C                          ------     ---------     --------
C
C                            0        F(A,Y(J))     F(A,Y(J))
C                            1        U(A,Y(J))     U(B,Y(J))
C                            2        U(A,Y(J))     F(B,Y(J))
C                            3        F(A,Y(J))     F(B,Y(J))
C                            4        F(A,Y(J))     U(B,Y(J))
C
C
C                          NBDCND     F(I,1)        F(I,N+1)
C                          ------     ---------     --------
C
C                            0        F(X(I),C)     F(X(I),C)
C                            1        U(X(I),C)     U(X(I),D)
C                            2        U(X(I),C)     F(X(I),D)
C                            3        F(X(I),C)     F(X(I),D)
C                            4        F(X(I),C)     U(X(I),D)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION U
C                          AND THE RIGHT SIDE F AT A CORNER THEN THE
C                          SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSCRT.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1  .
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (X(I),Y(J)), I = 1,2,...,M+1,
C                          J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  HWSCRT
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION.  THIS SOLUTION PLUS ANY
C                          CONSTANT IS ALSO A SOLUTION.  HENCE, THE
C                          SOLUTION IS NOT UNIQUE.  THE VALUE OF
C                          PERTRB SHOULD BE SMALL COMPARED TO THE
C                          RIGHT SIDE F.  OTHERWISE, A SOLUTION IS
C                          OBTAINED TO AN ESSENTIALLY DIFFERENT
C                          PROBLEM. THIS COMPARISON SHOULD ALWAYS
C                          BE MADE TO INSURE THAT A MEANINGFUL
C                          SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 6,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR
C                          = 1  A .GE. B
C                          = 2  MBDCND .LT. 0 OR MBDCND .GT. 4
C                          = 3  C .GE. D
C                          = 4  N .LE. 3
C                          = 5  NBDCND .LT. 0 OR NBDCND .GT. 4
C                          = 6  LAMBDA .GT. 0
C                          = 7  IDIMF .LT. M+1
C                          = 8  M .LE. 3
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HWSCRT, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,genbun.f,gnbnaux.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS GENBUN TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE  M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          M*N*(LOG2(N)
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS NBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
C                        OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR N
C                        AND M AS LARGE AS 64.  MORE DETAILS ABOUT
C                        ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE GENBUN WHICH IS THE ROUTINE THAT
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS"
C                          NCAR TN/IA-109, JULY, 1975, 138 PP.
C***********************************************************************
      SUBROUTINE HWSCRT(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
C
C     CHECK FOR INVALID PARAMETERS.
C
      IERROR = 0
      IF (A >= B) IERROR = 1
      IF (MBDCND<0 .OR. MBDCND>4) IERROR = 2
      IF (C >= D) IERROR = 3
      IF (N <= 3) IERROR = 4
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 5
      IF (IDIMF < M + 1) IERROR = 7
      IF (M <= 3) IERROR = 8
      IF (IERROR /= 0) RETURN 
!     allocate required work space length (generous estimate)
      IRWK=4*(N+1)+(13+INT(ALOG(FLOAT(N+1)/ALOG(2.0))))*(M+1)
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hwscrtt(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +             ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HWSCRT


 
      SUBROUTINE HWSCRTT(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NPEROD, MPEROD, NP, NP1, MP, MP1, NSTART, NSTOP, NSKIP
     1   , NUNK, MSTART, MSTOP, MSKIP, J, MUNK, I, ID2, ID3, ID4, MSP1, 
     2   MSTM1, NSP1, NSTM1, IERR1
      REAL::DELTAX,TWDELX,DELXSQ,DELTAY,TWDELY,DELYSQ,S,ST2,A1,A2,S1
C-----------------------------------------------
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND > 0) MPEROD = 1
      DELTAX = (B - A)/FLOAT(M)
      TWDELX = 2./DELTAX
      DELXSQ = 1./DELTAX**2
      DELTAY = (D - C)/FLOAT(N)
      TWDELY = 2./DELTAY
      DELYSQ = 1./DELTAY**2
      NP = NBDCND + 1
      NP1 = N + 1
      MP = MBDCND + 1
      MP1 = M + 1
      NSTART = 1
      NSTOP = N
      NSKIP = 1
      GO TO (104,101,102,103,104) NP
  101 CONTINUE
      NSTART = 2
      GO TO 104
  102 CONTINUE
      NSTART = 2
  103 CONTINUE
      NSTOP = NP1
      NSKIP = 2
  104 CONTINUE
      NUNK = NSTOP - NSTART + 1
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      MSTART = 1
      MSTOP = M
      MSKIP = 1
      GO TO (117,105,106,109,110) MP
  105 CONTINUE
      MSTART = 2
      GO TO 107
  106 CONTINUE
      MSTART = 2
      MSTOP = MP1
      MSKIP = 2
  107 CONTINUE
      F(2,NSTART:NSTOP) = F(2,NSTART:NSTOP) - F(1,NSTART:NSTOP)*DELXSQ
      GO TO 112
  109 CONTINUE
      MSTOP = MP1
      MSKIP = 2
  110 CONTINUE
      F(1,NSTART:NSTOP) = F(1,NSTART:NSTOP) + BDA(NSTART:NSTOP)*TWDELX
  112 CONTINUE
      SELECT CASE (MSKIP) 
      CASE DEFAULT
         F(M,NSTART:NSTOP) = F(M,NSTART:NSTOP) - F(MP1,NSTART:NSTOP)*
     1      DELXSQ
      CASE (2) 
         F(MP1,NSTART:NSTOP) = F(MP1,NSTART:NSTOP) - BDB(NSTART:NSTOP)*
     1      TWDELX
      END SELECT
  117 CONTINUE
      MUNK = MSTOP - MSTART + 1
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (127,118,118,120,120) NP
  118 CONTINUE
      F(MSTART:MSTOP,2) = F(MSTART:MSTOP,2) - F(MSTART:MSTOP,1)*DELYSQ
      GO TO 122
  120 CONTINUE
      F(MSTART:MSTOP,1) = F(MSTART:MSTOP,1) + BDC(MSTART:MSTOP)*TWDELY
  122 CONTINUE
      SELECT CASE (NSKIP) 
      CASE DEFAULT
         F(MSTART:MSTOP,N) = F(MSTART:MSTOP,N) - F(MSTART:MSTOP,NP1)*
     1      DELYSQ
      CASE (2) 
         F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,NP1) - BDD(MSTART:MSTOP)*
     1      TWDELY
      END SELECT
C
C    MULTIPLY RIGHT SIDE BY DELTAY**2.
C
  127 CONTINUE
      DELYSQ = DELTAY*DELTAY
      F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:NSTOP)*DELYSQ
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2 + MUNK
      ID4 = ID3 + MUNK
      S = DELYSQ*DELXSQ
      ST2 = 2.*S
      W(:MUNK) = S
      W(ID2+1:MUNK+ID2) = (-ST2) + ELMBDA*DELYSQ
      W(ID3+1:MUNK+ID3) = S
      IF (MP /= 1) THEN
         W(1) = 0.
         W(ID4) = 0.
      ENDIF
      GO TO (135,135,132,133,134) MP
  132 CONTINUE
      W(ID2) = ST2
      GO TO 135
  133 CONTINUE
      W(ID2) = ST2
  134 CONTINUE
      W(ID3+1) = ST2
  135 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 6
         ELSE
            IF ((NBDCND==0 .OR. NBDCND==3) .AND. (MBDCND==0 .OR. MBDCND
     1         ==3)) THEN
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
               A1 = 1.
               A2 = 1.
               IF (NBDCND == 3) A2 = 2.
               IF (MBDCND == 3) A1 = 2.
               S1 = 0.
               MSP1 = MSTART + 1
               MSTM1 = MSTOP - 1
               NSP1 = NSTART + 1
               NSTM1 = NSTOP - 1
               DO J = NSP1, NSTM1
                  S = 0.
                  S = SUM(F(MSP1:MSTM1,J))
                  S1 = S1 + S*A1 + F(MSTART,J) + F(MSTOP,J)
               END DO
               S1 = A2*S1
               S = 0.
               S = SUM(F(MSP1:MSTM1,NSTART)+F(MSP1:MSTM1,NSTOP))
               S1 = S1 + S*A1 + F(MSTART,NSTART) + F(MSTART,NSTOP) + F(
     1            MSTOP,NSTART) + F(MSTOP,NSTOP)
               S = (2. + FLOAT(NUNK - 2)*A2)*(2. + FLOAT(MUNK - 2)*A1)
               PERTRB = S1/S
               F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:
     1            NSTOP) - PERTRB
               PERTRB = PERTRB/DELYSQ
C
C     SOLVE THE EQUATION.
C
            ENDIF
         ENDIF
      ENDIF
      IERR1 = 0
      CALL GENBUNN (NPEROD, NUNK, MPEROD, MUNK, W(1), W(ID2+1), W(ID3+1)
     1   , IDIMF, F(MSTART,NSTART), IERR1, W(ID4+1))
C
C     FILL IN IDENTICAL VALUES WHEN HAVE PERIODIC BOUNDARY CONDITIONS.
C
      IF (NBDCND == 0) THEN
         F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,1)
      ENDIF
      IF (MBDCND == 0) THEN
         F(MP1,NSTART:NSTOP) = F(1,NSTART:NSTOP)
         IF (NBDCND == 0) F(MP1,NP1) = F(1,NP1)
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HWSCRTT
C
C     file hwscsp.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,
C    +                   BDRS,BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C
C DIMENSION OF           BDTS(N+1),     BDTF(N+1), BDRS(M+1), BDRF(M+1),
C ARGUMENTS              F(IDIMF,N+1)
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION
C                        TO THE MODIFIED HELMHOLTZ EQUATION IN
C                        SPHERICAL COORDINATES ASSUMING AXISYMMETRY
C                        (NO DEPENDENCE ON LONGITUDE).  THE EQUATION
C                        IS
C
C                          (1/R**2)(D/DR)((R**2)(D/DR)U) +
C
C                          (1/(R**2)SIN(THETA))(D/DTHETA)
C
C                          (SIN(THETA)(D/DTHETA)U) +
C
C                          (LAMBDA/(RSIN(THETA))**2)U = F(THETA,R).
C
C                        THIS TWO DIMENSIONAL MODIFIED HELMHOLTZ
C                        EQUATION RESULTS FROM THE FOURIER TRANSFORM
C                        OF THE THREE DIMENSIONAL POISSON EQUATION.
C
C USAGE                  CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,
C                                     RS,RF,N,NBDCND,BDRS,BDRF,ELMBDA,
C                                     F,IDIMF,PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               INTL
C                          = 0  ON INITIAL ENTRY TO HWSCSP OR IF ANY
C                               OF THE ARGUMENTS RS, RF, N, NBDCND
C                               ARE CHANGED FROM A PREVIOUS CALL.
C                          = 1  IF RS, RF, N, NBDCND ARE ALL UNCHANGED
C                               FROM PREVIOUS CALL TO HWSCSP.
C
C                          NOTE:
C                          A CALL WITH INTL=0 TAKES APPROXIMATELY
C                          1.5 TIMES AS MUCH TIME AS A CALL WITH
C                          INTL = 1  .  ONCE A CALL WITH INTL = 0
C                          HAS BEEN MADE THEN SUBSEQUENT SOLUTIONS
C                          CORRESPONDING TO DIFFERENT F, BDTS, BDTF,
C                          BDRS, BDRF CAN BE OBTAINED FASTER WITH
C                          INTL = 1 SINCE INITIALIZATION IS NOT
C                          REPEATED.
C
C                        TS,TF
C                          THE RANGE OF THETA (COLATITUDE), I.E.,
C                          TS .LE. THETA .LE. TF. TS MUST BE LESS
C                          THAN TF.  TS AND TF ARE IN RADIANS. A TS OF
C                          ZERO CORRESPONDS TO THE NORTH POLE AND A
C                          TF OF PI CORRESPONDS TO THE SOUTH POLE.
C
C                          **** IMPORTANT ****
C
C                          IF TF IS EQUAL TO PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          TF = PIMACH(DUM). THIS INSURES THAT TF
C                          IN THE USER'S PROGRAM IS EQUAL TO PI IN
C                          THIS PROGRAM WHICH PERMITS SEVERAL TESTS
C                          OF THE  INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (TS,TF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS
C                          IN THE THETA-DIRECTION GIVEN BY
C                          THETA(K) = (I-1)DTHETA+TS FOR
C                          I = 1,2,...,M+1, WHERE DTHETA = (TF-TS)/M
C                          IS THE PANEL WIDTH.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT THETA = TS AND  THETA = TF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THETA = TF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS AND THETA = TF
C                               (SEE NOTES 1,2 BELOW).
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW) AND
C                               SOLUTION IS SPECIFIED AT THETA = TF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE SOLUTION IS
C                                SPECIFIED AT THETA = TF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO THETA
C                               IS SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE SOLUTION IS
C                                UNSPECIFIED AT THETA = TF = PI.
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TF = PI.
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THETA = TF = PI.
C
C                          NOTE 1:
C                          IF TS = 0, DO NOT USE MBDCND = 3,4, OR 8,
C                          BUT INSTEAD USE MBDCND = 5,6, OR 9  .
C
C                          NOTE 2:
C                          IF TF = PI, DO NOT USE MBDCND = 2,3, OR 6,
C                          BUT INSTEAD USE MBDCND = 7,8, OR 9  .
C
C                        BDTS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TS.  WHEN MBDCND = 3,4, OR 8,
C
C                            BDTS(J) = (D/DTHETA)U(TS,R(J)),
C                            J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTS IS
C                          A DUMMY VARIABLE.
C
C                        BDTF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TF.  WHEN MBDCND = 2,3, OR 6,
C
C                          BDTF(J) = (D/DTHETA)U(TF,R(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTF IS
C                          A DUMMY VARIABLE.
C
C                        RS,RF
C                          THE RANGE OF R, I.E., RS .LE. R .LT. RF.
C                          RS MUST BE LESS THAN RF.  RS MUST BE
C                          NON-NEGATIVE.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (RS,RF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS IN THE
C                          R-DIRECTION GIVEN BY R(J) = (J-1)DR+RS
C                          FOR J = 1,2,...,N+1, WHERE DR = (RF-RS)/N
C                          IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 2
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT R = RS AND R = RF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = RS AND R = RF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = RS AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO R
C                               IS SPECIFIED AT R = RF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = RS AND R = RF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               RS AND THE SOLUTION IS SPECIFIED AT
C                               R = RF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = RS = 0 (SEE NOTE BELOW)  AND THE
C                               SOLUTION IS SPECIFIED AT R = RF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = RS = 0 (SEE NOTE BELOW) AND THE
C                               DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO R IS SPECIFIED AT R = RF.
C
C                          NOTE:
C                          NBDCND = 5 OR 6 CANNOT BE USED WITH
C                          MBDCND = 1,2,4,5, OR 7.  THE FORMER
C                          INDICATES THAT THE SOLUTION IS UNSPECIFIED
C                          AT R = 0, THE LATTER INDICATES THAT THE
C                          SOLUTION IS SPECIFIED).
C                          USE INSTEAD   NBDCND = 1 OR 2  .
C
C                        BDRS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = RS.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDRS(I) = (D/DR)U(THETA(I),RS),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDRS IS
C                          A DUMMY VARIABLE.
C
C                        BDRF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUES OF THE
C                          DERIVATIVE OF THE SOLUTION WITH RESPECT
C                          TO R AT R = RF.
C
C                          WHEN NBDCND = 2,3, OR 6,
C                            BDRF(I) = (D/DR)U(THETA(I),RF),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDRF IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSCSP WILL
C                          ATTEMPT TO FIND A SOLUTION.  IF NBDCND = 5
C                          OR 6 OR  MBDCND = 5,6,7,8, OR 9, ELMBDA
C                          MUST BE ZERO.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES OF THE
C                          RIGHT SIDE OF THE HELMHOLTZ EQUATION AND
C                          BOUNDARY VALUES (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(THETA(I),R(J)).
C
C                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
C                          FOR J=1,2,...,N+1,  I=1,2,...,M+1,
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   ----------        ----------
C
C                            1      U(TS,R(J))        U(TF,R(J))
C                            2      U(TS,R(J))        F(TF,R(J))
C                            3      F(TS,R(J))        F(TF,R(J))
C                            4      F(TS,R(J))        U(TF,R(J))
C                            5      F(0,R(J))         U(TF,R(J))
C                            6      F(0,R(J))         F(TF,R(J))
C                            7      U(TS,R(J))        F(PI,R(J))
C                            8      F(TS,R(J))        F(PI,R(J))
C                            9      F(0,R(J))         F(PI,R(J))
C
C                            NBDCND   F(I,1)            F(I,N+1)
C                            ------   --------------    --------------
C
C                              1      U(THETA(I),RS)    U(THETA(I),RF)
C                              2      U(THETA(I),RS)    F(THETA(I),RF)
C                              3      F(THETA(I),RS)    F(THETA(I),RF)
C                              4      F(THETA(I),RS)    U(THETA(I),RF)
C                              5      F(TS,0)           U(THETA(I),RF)
C                              6      F(TS,0)           F(THETA(I),RF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F AT A CORNER THEN
C                          THE SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSCSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1  .
C
C                        W
c                          A fortran 90 derived TYPE (fishworkspace) variable
c                          that must be declared by the user.  The first
c                          two declarative statements in the user program
c                          calling SEPELI must be:
c
c                               USE fish
c                               TYPE (fishworkspace) :: W
c
c                          The first statement makes the fishpack module
c                          defined in the file "fish.f" available to the
c                          user program calling HWSCSP.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in HWSCSP to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (INTL=0) call to HWSCSP.  These
c                          must be preserved on non-initial calls (INTL=1)
c                          to HWSCSP.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling HWSCSP should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          HWSCSP.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
c
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),R(J)),  I = 1,2,...,M+1,
C                                            J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  HWSCSP
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION. THIS SOLUTION IS NOT UNIQUE
C                          AND IS UNNORMALIZED. THE VALUE OF PERTRB
C                          SHOULD BE SMALL COMPARED TO THE RIGHT SIDE
C                          F. OTHERWISE , A SOLUTION IS OBTAINED TO
C                          AN ESSENTIALLY DIFFERENT PROBLEM. THIS
C                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
C                          THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 10,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          = 1  TS.LT.0. OR TF.GT.PI
C                          = 2  TS.GE.TF
C                          = 3  M.LT.5
C                          = 4  MBDCND.LT.1 OR MBDCND.GT.9
C                          = 5  RS.LT.0
C                          = 6  RS.GE.RF
C                          = 7  N.LT.5
C                          = 8  NBDCND.LT.1 OR NBDCND.GT.6
C                          = 9  ELMBDA.GT.0
C                          = 10 IDIMF.LT.M+1
C                          = 11 ELMBDA.NE.0 AND MBDCND.GE.5
C                          = 12 ELMBDA.NE.0 AND NBDCND EQUALS 5 OR 6
C                          = 13 MBDCND EQUALS 5,6 OR 9 AND TS.NE.0
C                          = 14 MBDCND.GE.7 AND TF.NE.PI
C                          = 15 TS.EQ.0 AND MBDCND EQUALS 3,4 OR 8
C                          = 16 TF.EQ.PI AND MBDCND EQUALS 2,3 OR 6
C                          = 17 NBDCND.GE.5 AND RS.NE.0
C                          = 18 NBDCND.GE.5 AND MBDCND EQUALS 1,2,4,5 OR
C                          = 20 If the dynamic allocation of real and
C                               complex work space in the derived type
C                               (fishworkspace) variable W fails (e.g.,
c                               if N,M are too large for the platform used)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSLIBY INCORRECT CALL TO HWSCSP, THE
C                          USER SHOULD TEST IERROR AFTER A CALL.
C
C                        W
c                          The derived type (fishworkspace) variable W
c                          contains real and complex values that must not
C                          be destroyed if HWSCSP is called again with
C                          INTL=1.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,blktri.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980. Revised by John
c                        Adams in June 2004 using Fortran 90 dynamically
C                        allocated work space and derived datat types
c                        to eliminate mixed mode conflicts in the earlier
c                        versions.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS BLKTRI TO SOLVE THE SYSTEM.
C
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS"
C                          NCAR TN/IA-109, JULY, 1975, 138 PP.
C***********************************************************************
      SUBROUTINE HWSCSP(INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, 
     1   NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: INTL
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: TS
      REAL  :: TF
      REAL  :: RS
      REAL  :: RF
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDTS(*)
      REAL  :: BDTF(*)
      REAL  :: BDRS(*)
      REAL  :: BDRF(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, NCK, L, K, NP
     1   , IRWK, ICWK, NP1, MP1
      REAL :: PI, DUM

      SAVE I1, I2, I3, I4, I5, I6, I7, I8, I9, I10
C-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      IERROR = 0
      IF (TS<0. .OR. TF>PI) IERROR = 1
      IF (TS >= TF) IERROR = 2
      IF (M < 5) IERROR = 3
      IF (MBDCND<1 .OR. MBDCND>9) IERROR = 4
      IF (RS < 0.) IERROR = 5
      IF (RS >= RF) IERROR = 6
      IF (N < 5) IERROR = 7
      IF (NBDCND<1 .OR. NBDCND>6) IERROR = 8
      IF (ELMBDA > 0.) IERROR = 9
      IF (IDIMF < M + 1) IERROR = 10
      IF (ELMBDA/=0. .AND. MBDCND>=5) IERROR = 11
      IF (ELMBDA/=0. .AND. (NBDCND==5 .OR. NBDCND==6)) IERROR = 12
      IF((MBDCND==5.OR.MBDCND==6.OR.MBDCND==9).AND.TS/=0.)IERROR=13
      IF (MBDCND>=7 .AND. TF/=PI) IERROR = 14
      IF(TS==0..AND.(MBDCND==4.OR.MBDCND==8.OR.MBDCND==3))IERROR=15
      IF(TF==PI.AND.(MBDCND==2.OR.MBDCND==3.OR.MBDCND==6))IERROR=16
      IF (NBDCND>=5 .AND. RS/=0.) IERROR = 17
      IF (NBDCND>=5 .AND. (MBDCND==1 .OR. MBDCND==2 .OR. MBDCND==5 .OR. 
     1   MBDCND==7)) IERROR = 18
      IF (IERROR/=0 .AND. IERROR/=9) RETURN 
      NCK = N
      GO TO (101,103,102,103,101,103) NBDCND
  101 CONTINUE
      NCK = NCK - 1
      GO TO 103
  102 CONTINUE
      NCK = NCK + 1
  103 CONTINUE
      L = 2
      K = 1
      L = L + L
      K = K + 1
      DO WHILE(NCK - L > 0)
         L = L + L
         K = K + 1
      END DO
      L = L + L
 
      IF (INTL == 0) THEN
!          compute blktri work space lengths
         NP = NBDCND
         CALL BLK_SPACE (N, M, IRWK, ICWK)
         NP1 = N + 1
         MP1 = M + 1
         I1 = (K - 2)*L + K + MAX0(2*N,6*M) + 13
         I2 = I1 + NP1
         I3 = I2 + NP1
         I4 = I3 + NP1
         I5 = I4 + NP1
         I6 = I5 + NP1
         I7 = I6 + MP1
         I8 = I7 + MP1
         I9 = I8 + MP1
         I10 = I9 + MP1
!          set real and complex work space requirements
         IRWK = I10 + MP1
         ICWK = ICWK + 3*M
!          allocate work space
         CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!          return if allocation fails
         IF (IERROR == 20) RETURN 
      ENDIF
      CALL HWSCS1 (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     +BDRF,ELMBDA,F,IDIMF,PERTRB,w%rew,w%cxw,w%rew(i1),w%rew(i2),
     +w%rew(i3),w%rew(i4),w%rew(i5),w%rew(i6),w%rew(i7),w%rew(i8),
     +w%rew(i9),w%rew(i10),ierror)
      RETURN 
      END SUBROUTINE HWSCSP


      SUBROUTINE HWSCS1(INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, 
     1   NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, W, WC, S, AN, BN
     2   , CN, R, AM, BM, CM, SINT, BMH, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: INTL
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL , INTENT(IN) :: TS
      REAL , INTENT(IN) :: TF
      REAL , INTENT(IN) :: RS
      REAL , INTENT(IN) :: RF
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDTS(*)
      REAL , INTENT(IN) :: BDTF(*)
      REAL , INTENT(IN) :: BDRS(*)
      REAL , INTENT(IN) :: BDRF(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
      REAL , INTENT(INOUT) :: S(*)
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL , INTENT(INOUT) :: R(*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL , INTENT(INOUT) :: SINT(*)
      REAL , INTENT(INOUT) :: BMH(*)
      COMPLEX  :: WC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MP1, I, NP1, J, MP, NP, ITS, ITF, ITSP, ITFM, ICTR, JRS
     1   , L, JRF, JRSP, JRFM, MUNK, NUNK, ISING, IFLG
      REAL :: PI, DUM, EPS, DTH, TDT, HDTH, SDTS, THETA, T1, DR, HDR, 
     1   TDR, DR2, CZR, AT, CT, WTS, WTF, AR, WTNM, YPS, CR, WRS, WRF, 
     2   WRZ, SUM, R2, HNE, YHLD, RS2, RF2, RSQ, XP, YPH, XPS
!!!      REAL :: EPMACH
C-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      EPS = EPMACH(DUM)
      MP1 = M + 1
      DTH = (TF - TS)/FLOAT(M)
      TDT = DTH + DTH
      HDTH = DTH/2.
      SDTS = 1./(DTH*DTH)
      DO I = 1, MP1
         THETA = TS + FLOAT(I - 1)*DTH
         SINT(I) = SIN(THETA)
         IF (SINT(I) == 0.) CYCLE 
         T1 = SDTS/SINT(I)
         AM(I) = T1*SIN(THETA - HDTH)
         CM(I) = T1*SIN(THETA + HDTH)
         BM(I) = -(AM(I)+CM(I))
      END DO
      NP1 = N + 1
      DR = (RF - RS)/FLOAT(N)
      HDR = DR/2.
      TDR = DR + DR
      DR2 = DR*DR
      CZR = 6.*DTH/(DR2*(COS(TS) - COS(TF)))
      DO J = 1, NP1
         R(J) = RS + FLOAT(J - 1)*DR
         AN(J) = (R(J)-HDR)**2/DR2
         CN(J) = (R(J)+HDR)**2/DR2
         BN(J) = -(AN(J)+CN(J))
      END DO
      MP = 1
      NP = 1
C
C BOUNDARY CONDITION AT PHI=PS
C
      GO TO (104,104,105,105,106,106,104,105,106) MBDCND
  104 CONTINUE
      AT = AM(2)
      ITS = 2
      GO TO 107
  105 CONTINUE
      AT = AM(1)
      ITS = 1
      CM(1) = CM(1) + AM(1)
      GO TO 107
  106 CONTINUE
      ITS = 1
      BM(1) = -4.*SDTS
      CM(1) = -BM(1)
C
C BOUNDARY CONDITION AT PHI=PF
C
  107 CONTINUE
      GO TO (108,109,109,108,108,109,110,110,110) MBDCND
  108 CONTINUE
      CT = CM(M)
      ITF = M
      GO TO 111
  109 CONTINUE
      CT = CM(M+1)
      AM(M+1) = AM(M+1) + CM(M+1)
      ITF = M + 1
      GO TO 111
  110 CONTINUE
      ITF = M + 1
      AM(M+1) = 4.*SDTS
      BM(M+1) = -AM(M+1)
  111 CONTINUE
      WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
      WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
      ITSP = ITS + 1
      ITFM = ITF - 1
C
C BOUNDARY CONDITION AT R=RS
C
      ICTR = 0
      SELECT CASE (NBDCND) 
      CASE DEFAULT
         AR = AN(2)
         JRS = 2
      CASE (3:4) 
         AR = AN(1)
         JRS = 1
         CN(1) = CN(1) + AN(1)
      CASE (5:6) 
         JRS = 2
         ICTR = 1
         S(N) = AN(N)/BN(N)
         DO J = 3, N
            L = N - J + 2
            S(L) = AN(L)/(BN(L)-CN(L)*S(L+1))
         END DO
         S(2) = -S(2)
         DO J = 3, N
            S(J) = -S(J)*S(J-1)
         END DO
         WTNM = WTS + WTF
         DO I = ITSP, ITFM
            WTNM = WTNM + SINT(I)
         END DO
         YPS = CZR*WTNM*(S(2)-1.)
      END SELECT
C
C BOUNDARY CONDITION AT R=RF
C
  118 CONTINUE
      GO TO (119,120,120,119,119,120) NBDCND
  119 CONTINUE
      CR = CN(N)
      JRF = N
      GO TO 121
  120 CONTINUE
      CR = CN(N+1)
      AN(N+1) = AN(N+1) + CN(N+1)
      JRF = N + 1
  121 CONTINUE
      WRS = AN(JRS+1)*R(JRS)**2/CN(JRS)
      WRF = CN(JRF-1)*R(JRF)**2/AN(JRF)
      WRZ = AN(JRS)/CZR
      JRSP = JRS + 1
      JRFM = JRF - 1
      MUNK = ITF - ITS + 1
      NUNK = JRF - JRS + 1
      BMH(ITS:ITF) = BM(ITS:ITF)
      ISING = 0
      GO TO (132,132,123,132,132,123) NBDCND
  123 CONTINUE
      GO TO (132,132,124,132,132,124,132,124,124) MBDCND
  124 CONTINUE
      IF (ELMBDA >= 0.) THEN
         ISING = 1
         SUM = WTS*WRS + WTS*WRF + WTF*WRS + WTF*WRF
         IF (ICTR /= 0) THEN
            SUM = SUM + WRZ
         ENDIF
         DO J = JRSP, JRFM
            R2 = R(J)**2
            DO I = ITSP, ITFM
               SUM = SUM + R2*SINT(I)
            END DO
         END DO
         DO J = JRSP, JRFM
            SUM = SUM + (WTS + WTF)*R(J)**2
         END DO
         DO I = ITSP, ITFM
            SUM = SUM + (WRS + WRF)*SINT(I)
         END DO
         HNE = SUM
      ENDIF
  132 CONTINUE
      GO TO (133,133,133,133,134,134,133,133,134) MBDCND
  133 CONTINUE
      BM(ITS) = BMH(ITS) + ELMBDA/SINT(ITS)**2
  134 CONTINUE
      GO TO (135,135,135,135,135,135,136,136,136) MBDCND
  135 CONTINUE
      BM(ITF) = BMH(ITF) + ELMBDA/SINT(ITF)**2
  136 CONTINUE
      BM(ITSP:ITFM) = BMH(ITSP:ITFM) + ELMBDA/SINT(ITSP:ITFM)**2
      GO TO (138,138,140,140,142,142,138,140,142) MBDCND
  138 CONTINUE
      F(2,JRS:JRF) = F(2,JRS:JRF) - AT*F(1,JRS:JRF)/R(JRS:JRF)**2
      GO TO 142
  140 CONTINUE
      F(1,JRS:JRF) = F(1,JRS:JRF) + TDT*BDTS(JRS:JRF)*AT/R(JRS:JRF)**2
  142 CONTINUE
      GO TO (143,145,145,143,143,145,147,147,147) MBDCND
  143 CONTINUE
      F(M,JRS:JRF) = F(M,JRS:JRF) - CT*F(M+1,JRS:JRF)/R(JRS:JRF)**2
      GO TO 147
  145 CONTINUE
      F(M+1,JRS:JRF)=F(M+1,JRS:JRF)-TDT*BDTF(JRS:JRF)*CT/R(JRS:JRF)**2
  147 CONTINUE
      SELECT CASE (NBDCND) 
      CASE DEFAULT
         IF (MBDCND - 3 /= 0) GO TO 155
         YHLD = F(ITS,1) - CZR/TDT*(SIN(TF)*BDTF(2)-SIN(TS)*BDTS(2))
         F(:MP1,1) = YHLD
      CASE (1:2) 
         RS2 = (RS + DR)**2
         F(ITS:ITF,2) = F(ITS:ITF,2) - AR*F(ITS:ITF,1)/RS2
      CASE (3:4) 
         F(ITS:ITF,1) = F(ITS:ITF,1) + TDR*BDRS(ITS:ITF)*AR/RS**2
      END SELECT
  155 CONTINUE
      GO TO (156,158,158,156,156,158) NBDCND
  156 CONTINUE
      RF2 = (RF - DR)**2
      F(ITS:ITF,N) = F(ITS:ITF,N) - CR*F(ITS:ITF,N+1)/RF2
      GO TO 160
  158 CONTINUE
      F(ITS:ITF,N+1) = F(ITS:ITF,N+1) - TDR*BDRF(ITS:ITF)*CR/RF**2
  160 CONTINUE
      PERTRB = 0.
      IF (ISING /= 0) THEN
         SUM = WTS*WRS*F(ITS,JRS) + WTS*WRF*F(ITS,JRF) + WTF*WRS*F(ITF,
     1      JRS) + WTF*WRF*F(ITF,JRF)
         IF (ICTR /= 0) THEN
            SUM = SUM + WRZ*F(ITS,1)
         ENDIF
         DO J = JRSP, JRFM
            R2 = R(J)**2
            DO I = ITSP, ITFM
               SUM = SUM + R2*SINT(I)*F(I,J)
            END DO
         END DO
         SUM = SUM + DOT_PRODUCT(R(JRSP:JRFM)**2,WTS*F(ITS,JRSP:JRFM)+
     1      WTF*F(ITF,JRSP:JRFM))
         SUM = SUM + DOT_PRODUCT(SINT(ITSP:ITFM),WRS*F(ITSP:ITFM,JRS)+
     1      WRF*F(ITSP:ITFM,JRF))
         PERTRB = SUM/HNE
         F(:MP1,:NP1) = F(:MP1,:NP1) - PERTRB
      ENDIF
      DO J = JRS, JRF
         RSQ = R(J)**2
         F(ITS:ITF,J) = RSQ*F(ITS:ITF,J)
      END DO
      IFLG = INTL
      CALL BLKTRII (IFLG, NP, NUNK, AN(JRS), BN(JRS), CN(JRS), MP, MUNK
     1   , AM(ITS), BM(ITS), CM(ITS), IDIMF, F(ITS,JRS), IERROR, W, WC)
      IFLG = IFLG + 1
      DO WHILE(IFLG - 1 == 0)
         CALL BLKTRII (IFLG, NP, NUNK, AN(JRS), BN(JRS), CN(JRS), MP, 
     1      MUNK, AM(ITS), BM(ITS), CM(ITS), IDIMF, F(ITS,JRS), IERROR, 
     2      W, WC)
         IFLG = IFLG + 1
      END DO
      IF (NBDCND == 0) THEN
         F(:MP1,JRF+1) = F(:MP1,JRS)
      ENDIF
      IF (MBDCND == 0) THEN
         F(ITF+1,:NP1) = F(ITS,:NP1)
      ENDIF
      XP = 0.
      IF (ICTR /= 0) THEN
         IF (ISING == 0) THEN
            SUM = WTS*F(ITS,2) + WTF*F(ITF,2)
            SUM = SUM + DOT_PRODUCT(SINT(ITSP:ITFM),F(ITSP:ITFM,2))
            YPH = CZR*SUM
            XP = (F(ITS,1)-YPH)/YPS
            DO J = JRS, JRF
               XPS = XP*S(J)
               F(ITS:ITF,J) = F(ITS:ITF,J) + XPS
            END DO
         ENDIF
         F(:MP1,1) = XP
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE HWSCS1
C
C     file hwscyl.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HWSCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N+1)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION
C                        TO THE HELMHOLTZ EQUATION IN CYLINDRICAL
C                        COORDINATES.  THIS MODIFIED HELMHOLTZ EQUATION
C
C                          (1/R)(D/DR)(R(DU/DR)) + (D/DZ)(DU/DZ)
C
C                          + (LAMBDA/R**2)U = F(R,Z)
C
C                        RESULTS FROM THE FOURIER TRANSFORM OF THE
C                        THREE-DIMENSIONAL POISSON EQUATION.
C
C USAGE                  CALL HWSCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               A,B
C                          THE RANGE OF R, I.E., A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (A,B) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE
C                          R-DIRECTION GIVEN BY R(I) = A+(I-1)DR,
C                          FOR I = 1,2,...,M+1, WHERE DR = (B-A)/M
C                          IS THE PANEL WIDTH. M MUST BE GREATER
C                          THAN 3.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = A AND R = B.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = A AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS
C                               SPECIFIED AT R = B.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND R = B.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT R = B.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT R = B.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS SPECIFIED
C                               AT R = B.
C
C                          IF A = 0, DO NOT USE MBDCND = 3 OR 4,
C                          BUT INSTEAD USE MBDCND = 1,2,5, OR 6  .
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = A.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,Z(J)), J = 1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO R AT R = B.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,Z(J)), J = 1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF Z, I.E., C .LE. Z .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE
C                          Z-DIRECTION GIVEN BY Z(J) = C+(J-1)DZ,
C                          FOR J = 1,2,...,N+1,
C                          WHERE DZ = (D-C)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 3.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Z = C AND Z = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Z,
C                               I.E., U(I,1) = U(I,N+1).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Z = C AND Z = D.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Z = C AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO Z IS
C                               SPECIFIED AT Z = D.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS
C                               SPECIFIED AT Z = C AND Z = D.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS SPECIFIED AT
C                               Z = C AND THE SOLUTION IS SPECIFIED
C                               AT Z = D.
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO Z AT Z = C.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DZ)U(R(I),C), I = 1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO Z AT Z = D.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DZ)U(R(I),D), I = 1,2,...,M+1
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSCYL WILL
C                          ATTEMPT TO FIND A SOLUTION.  LAMBDA MUST
C                          BE ZERO WHEN MBDCND = 5 OR 6  .
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES
C                          OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION AND BOUNDARY DATA (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(R(I),Z(J)).
C
C                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
C                          FOR J = 1,2,...,N+1 AND I = 1,2,...,M+1
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   ---------         ---------
C
C                            1      U(A,Z(J))         U(B,Z(J))
C                            2      U(A,Z(J))         F(B,Z(J))
C                            3      F(A,Z(J))         F(B,Z(J))
C                            4      F(A,Z(J))         U(B,Z(J))
C                            5      F(0,Z(J))         U(B,Z(J))
C                            6      F(0,Z(J))         F(B,Z(J))
C
C                          NBDCND   F(I,1)            F(I,N+1)
C                          ------   ---------         ---------
C
C                            0      F(R(I),C)         F(R(I),C)
C                            1      U(R(I),C)         U(R(I),D)
C                            2      U(R(I),C)         F(R(I),D)
C                            3      F(R(I),C)         F(R(I),D)
C                            4      F(R(I),C)         U(R(I),D)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F AT A CORNER THEN
C                          THE SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSCYL.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1  .
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),Z(J)), I =1,2,...,M+1, J =1,2,...,N+1.
C
C                        PERTRB
C                          IF ONE SPECIFIES A COMBINATION OF PERIODIC,
C                          DERIVATIVE, AND UNSPECIFIED BOUNDARY
C                          CONDITIONS FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HWSCYL THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION.  HENCE, THE SOLUTION IS NOT
C                          UNIQUE. THE VALUE OF PERTRB SHOULD BE
C                          SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM.  THIS
C                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
C                          THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG WHICH INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 11,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR.
C                          =  1  A .LT. 0  .
C                          =  2  A .GE. B.
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6  .
C                          =  4  C .GE. D.
C                          =  5  N .LE. 3
C                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4  .
C                          =  7  A = 0, MBDCND = 3 OR 4  .
C                          =  8  A .GT. 0, MBDCND .GE. 5  .
C                          =  9  A = 0, LAMBDA .NE. 0, MBDCND .GE. 5  .
C                          = 10  IDIMF .LT. M+1  .
C                          = 11  LAMBDA .GT. 0  .
C                          = 12  M .LE. 3
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HWSCYL, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,genbun.f,gnbnaux.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS GENBUN TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE  M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          M*N*(LOG2(N)
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS NBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
C                        OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR N
C                        AND M AS LARGE AS 64.  MORE DETAILS ABOUT
C                        ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE GENBUN WHICH IS THE ROUTINE THAT
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS"
C                          NCAR TN/IA-109, JULY, 1975, 138 PP.
C***********************************************************************
      SUBROUTINE HWSCYL(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
C
C     CHECK FOR INVALID PARAMETERS.
C
      IERROR = 0
      IF (A < 0.) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 3) IERROR = 5
      IF (NBDCND<=(-1) .OR. NBDCND>=5) IERROR = 6
      IF (A==0. .AND. (MBDCND==3 .OR. MBDCND==4)) IERROR = 7
      IF (A>0. .AND. MBDCND>=5) IERROR = 8
      IF (A==0. .AND. ELMBDA/=0. .AND. MBDCND>=5) IERROR = 9
      IF (IDIMF < M + 1) IERROR = 10
      IF (M <= 3) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     allocate real work space
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hwscyll(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +             ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HWSCYL


 
      SUBROUTINE HWSCYLL(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MP1, NP1, NP, MSTART, MSTOP, MUNK, NSTART, NSTOP, NUNK
     1   , ID2, ID3, ID4, ID5, ID6, ISTART, IJ, I, J, K, L, NSP1, NSTM1
     2   , IERR1, I1
      REAL::DELTAR,DLRBY2,DLRSQ,DELTHT,DLTHSQ,A1,R,A2,S,S1,S2
C-----------------------------------------------
      MP1 = M + 1
      DELTAR = (B - A)/FLOAT(M)
      DLRBY2 = DELTAR/2.
      DLRSQ = DELTAR**2
      NP1 = N + 1
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
C
C     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I,J).
C
      MSTART = 2
      MSTOP = M
      GO TO (104,103,102,101,101,102) MBDCND
  101 CONTINUE
      MSTART = 1
      GO TO 104
  102 CONTINUE
      MSTART = 1
  103 CONTINUE
      MSTOP = MP1
  104 CONTINUE
      MUNK = MSTOP - MSTART + 1
      NSTART = 1
      NSTOP = N
      GO TO (108,105,106,107,108) NP
  105 CONTINUE
      NSTART = 2
      GO TO 108
  106 CONTINUE
      NSTART = 2
  107 CONTINUE
      NSTOP = NP1
  108 CONTINUE
      NUNK = NSTOP - NSTART + 1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2 + MUNK
      ID4 = ID3 + MUNK
      ID5 = ID4 + MUNK
      ID6 = ID5 + MUNK
      ISTART = 1
      A1 = 2./DLRSQ
      IJ = 0
      IF (MBDCND==3 .OR. MBDCND==4) IJ = 1
      IF (MBDCND > 4) THEN
         W(1) = 0.
         W(ID2+1) = -2.*A1
         W(ID3+1) = 2.*A1
         ISTART = 2
         IJ = 1
      ENDIF
      DO I = ISTART, MUNK
         R = A + FLOAT(I - IJ)*DELTAR
         J = ID5 + I
         W(J) = R
         J = ID6 + I
         W(J) = 1./R**2
         W(I) = (R - DLRBY2)/(R*DLRSQ)
         J = ID3 + I
         W(J) = (R + DLRBY2)/(R*DLRSQ)
         K = ID6 + I
         J = ID2 + I
         W(J) = (-A1) + ELMBDA*W(K)
      END DO
      GO TO (114,111,112,113,114,112) MBDCND
  111 CONTINUE
      W(ID2) = A1
      GO TO 114
  112 CONTINUE
      W(ID2) = A1
  113 CONTINUE
      W(ID3+1) = A1*FLOAT(ISTART)
  114 CONTINUE
      GO TO (115,115,117,117,119,119) MBDCND
  115 CONTINUE
      A1 = W(1)
      F(2,NSTART:NSTOP) = F(2,NSTART:NSTOP) - A1*F(1,NSTART:NSTOP)
      GO TO 119
  117 CONTINUE
      A1 = 2.*DELTAR*W(1)
      F(1,NSTART:NSTOP) = F(1,NSTART:NSTOP) + A1*BDA(NSTART:NSTOP)
  119 CONTINUE
      GO TO (120,122,122,120,120,122) MBDCND
  120 CONTINUE
      A1 = W(ID4)
      F(M,NSTART:NSTOP) = F(M,NSTART:NSTOP) - A1*F(MP1,NSTART:NSTOP)
      GO TO 124
  122 CONTINUE
      A1 = 2.*DELTAR*W(ID4)
      F(MP1,NSTART:NSTOP) = F(MP1,NSTART:NSTOP) - A1*BDB(NSTART:NSTOP)
C
C     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
C
  124 CONTINUE
      A1 = 1./DLTHSQ
      L = ID5 - MSTART + 1
      GO TO (134,125,125,127,127) NP
  125 CONTINUE
      F(MSTART:MSTOP,2) = F(MSTART:MSTOP,2) - A1*F(MSTART:MSTOP,1)
      GO TO 129
  127 CONTINUE
      A1 = 2./DELTHT
      F(MSTART:MSTOP,1) = F(MSTART:MSTOP,1) + A1*BDC(MSTART:MSTOP)
  129 CONTINUE
      A1 = 1./DLTHSQ
      GO TO (134,130,132,132,130) NP
  130 CONTINUE
      F(MSTART:MSTOP,N) = F(MSTART:MSTOP,N) - A1*F(MSTART:MSTOP,NP1)
      GO TO 134
  132 CONTINUE
      A1 = 2./DELTHT
      F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,NP1) - A1*BDD(MSTART:MSTOP)
  134 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 11
         ELSE
            W(ID5+1) = 0.5*(W(ID5+2)-DLRBY2)
            GO TO (146,146,138,146,146,137) MBDCND
  137       CONTINUE
            W(ID5+1) = 0.5*W(ID5+1)
  138       CONTINUE
            GO TO (140,146,146,139,146) NP
  139       CONTINUE
            A2 = 2.
            GO TO 141
  140       CONTINUE
            A2 = 1.
  141       CONTINUE
            K = ID5 + MUNK
            W(K) = 0.5*(W(K-1)+DLRBY2)
            S = 0.
            DO I = MSTART, MSTOP
               S1 = 0.
               NSP1 = NSTART + 1
               NSTM1 = NSTOP - 1
               S1 = SUM(F(I,NSP1:NSTM1))
               K = I + L
               S = S + (A2*S1 + F(I,NSTART)+F(I,NSTOP))*W(K)
            END DO
            S2 = FLOAT(M)*A + (0.75 + FLOAT((M - 1)*(M + 1)))*DLRBY2
            IF (MBDCND == 3) S2 = S2 + 0.25*DLRBY2
            S1 = (2. + A2*FLOAT(NUNK - 2))*S2
            PERTRB = S/S1
            F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:NSTOP)
     1          - PERTRB
         ENDIF
      ENDIF
  146 CONTINUE
      W(:MSTOP-MSTART+1) = W(:MSTOP-MSTART+1)*DLTHSQ
      W(ID2+1:MSTOP-MSTART+1+ID2) = W(ID2+1:MSTOP-MSTART+1+ID2)*DLTHSQ
      W(ID3+1:MSTOP-MSTART+1+ID3) = W(ID3+1:MSTOP-MSTART+1+ID3)*DLTHSQ
      F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:NSTOP)*DLTHSQ
      W(1) = 0.
      W(ID4) = 0.
C
C     SOLVE THE SYSTEM OF EQUATIONS.
C
      IERR1 = 0
      I1 = 1
      CALL GENBUNN (NBDCND, NUNK, I1, MUNK, W(1), W(ID2+1), W(ID3+1), 
     1   IDIMF, F(MSTART,NSTART), IERR1, W(ID4+1))
      IF (NBDCND == 0) THEN
         F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,1)
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE HWSCYLL
C
C     file hwsplr.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HWSPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N+1)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION TO
C                        THE HELMHOLTZ EQUATION IN POLAR COORDINATES.
C                        THE EQUATION IS
C
C                            (1/R)(D/DR)(R(DU/DR)) +
C                            (1/R**2)(D/DTHETA)(DU/DTHETA) +
C                            LAMBDA*U = F(R,THETA).
C
C USAGE                  CALL HWSPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               A,B
C                          THE RANGE OF R, I.E., A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (A,B) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE
C                          R-DIRECTION GIVEN BY R(I) = A+(I-1)DR,
C                          FOR I = 1,2,...,M+1,
C                          WHERE DR = (B-A)/M IS THE PANEL WIDTH.
C                          M MUST BE GREATER THAN 3.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = A AND R = B.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = A AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO R IS
C                               SPECIFIED AT R = B.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND R = B.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT R = B.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT R = B.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS SPECIFIED
C                               AT R = B.
C
C                          NOTE:
C                          IF A = 0, DO NOT USE MBDCND = 3 OR 4, BUT
C                          INSTEAD USE MBDCND = 1,2,5, OR 6  .
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = A.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,THETA(J)),
C                            J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = B.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,THETA(J)),
C                            J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF THETA, I.E., C .LE.
C                          THETA .LE. D.  C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE
C                          THETA-DIRECTION GIVEN BY
C                          THETA(J) = C+(J-1)DTHETA FOR
C                          J = 1,2,...,N+1, WHERE
C                          DTHETA = (D-C)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 3.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = C AND AT THETA = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN THETA,
C                               I.E., U(I,J) = U(I,N+J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THETA = D
C                               (SEE NOTE BELOW).
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = C AND THE SOLUTION IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C
C                          NOTE:
C                          WHEN NBDCND = 1,2, OR 4, DO NOT USE
C                          MBDCND = 5 OR 6
C                          (THE FORMER INDICATES THAT THE SOLUTION
C                          IS SPECIFIED AT R = 0, THE LATTER INDICATES
C                          THE SOLUTION IS UNSPECIFIED AT R = 0).
C                          USE INSTEAD MBDCND = 1 OR 2  .
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = C.  WHEN NBDCND = 3 OR 4,
C
C                            BDC(I) = (D/DTHETA)U(R(I),C),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = D.  WHEN NBDCND = 2 OR 3,
C
C                            BDD(I) = (D/DTHETA)U(R(I),D),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .LT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSPLR WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES
C                          OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION AND BOUNDARY DATA (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(R(I),THETA(J)).
C
C                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
C                          FOR J = 1,2,...,N+1 AND I = 1,2,...,M+1
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   -------------     -------------
C
C                            1      U(A,THETA(J))     U(B,THETA(J))
C                            2      U(A,THETA(J))     F(B,THETA(J))
C                            3      F(A,THETA(J))     F(B,THETA(J))
C                            4      F(A,THETA(J))     U(B,THETA(J))
C                            5      F(0,0)            U(B,THETA(J))
C                            6      F(0,0)            F(B,THETA(J))
C
C                          NBDCND   F(I,1)            F(I,N+1)
C                          ------   ---------         ---------
C
C                            0      F(R(I),C)         F(R(I),C)
C                            1      U(R(I),C)         U(R(I),D)
C                            2      U(R(I),C)         F(R(I),D)
C                            3      F(R(I),C)         F(R(I),D)
C                            4      F(R(I),C)         U(R(I),D)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F AT A CORNER THEN
C                          THEN THE SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSPLR.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1.
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),THETA(J)),
C                          I = 1,2,...,M+1, J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HWSPLR THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION.  HENCE, THE SOLUTION IS NOT
C                          UNIQUE.  PERTRB SHOULD BE SMALL COMPARED
C                          TO THE RIGHT SIDE. OTHERWISE, A SOLUTION
C                          IS OBTAINED TO AN ESSENTIALLY DIFFERENT
C                          PROBLEM.  THIS COMPARISON SHOULD ALWAYS
C                          BE MADE TO INSURE THAT A MEANINGFUL
C                          SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 11,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR.
C                          =  1  A .LT. 0  .
C                          =  2  A .GE. B.
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6  .
C                          =  4  C .GE. D.
C                          =  5  N .LE. 3
C                          =  6  NBDCND .LT. 0 OR .GT. 4  .
C                          =  7  A = 0, MBDCND = 3 OR 4  .
C                          =  8  A .GT. 0, MBDCND .GE. 5  .
C                          =  9  MBDCND .GE. 5, NBDCND .NE. 0
C                                AND NBDCND .NE. 3  .
C                          = 10  IDIMF .LT. M+1  .
C                          = 11  LAMBDA .GT. 0  .
C                          = 12  M .LE. 3
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HWSPLR, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,genbun.f,gnbnaux.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS GENBUN TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE  M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          M*N*(LOG2(N)
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS NBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
C                        OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR N
C                        AND M AS LARGE AS 64.  MORE DETAILS ABOUT
C                        ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE GENBUN WHICH IS THE ROUTINE THAT
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS"
C                          NCAR TN/IA-109, JULY, 1975, 138 PP.
C***********************************************************************
      SUBROUTINE HWSPLR(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)
!     Insert first required declarative statements for dynamically
!     allocated work space

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
C
C     CHECK FOR INVALID PARAMETERS.
C
      IERROR = 0
      IF (A < 0.) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 3) IERROR = 5
      IF (NBDCND<=(-1) .OR. NBDCND>=5) IERROR = 6
      IF (A==0. .AND. (MBDCND==3 .OR. MBDCND==4)) IERROR = 7
      IF (A>0. .AND. MBDCND>=5) IERROR = 8
      IF (MBDCND>=5 .AND. NBDCND/=0 .AND. NBDCND/=3) IERROR = 9
      IF (IDIMF < M + 1) IERROR = 10
      IF (M <= 3) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     compute and allocate required work space
      IRWK=4*(N+1)+(M+1)*(13+INT(ALOG(FLOAT(N+1))/ALOG(2.0)))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hwsplrr(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +                   ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HWSPLR


 
      SUBROUTINE HWSPLRR(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MP1, NP1, NP, MSTART, MSTOP, MUNK, NSTART, NSTOP, NUNK
     1   , ID2, ID3, ID4, ID5, ID6, IJ, I, J, L, LP, K, I1, IERR1, IP
      REAL :: ALL, DELTAR, DLRBY2, DLRSQ, DELTHT, DLTHSQ, A1, R, S2, A2
     1   , S, S1, YPOLE

      SAVE ALL
C-----------------------------------------------
      MP1 = M + 1
      DELTAR = (B - A)/FLOAT(M)
      DLRBY2 = DELTAR/2.
      DLRSQ = DELTAR**2
      NP1 = N + 1
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
C
C     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I,J).
C
      MSTART = 2
      MSTOP = MP1
      GO TO (101,105,102,103,104,105) MBDCND
  101 CONTINUE
      MSTOP = M
      GO TO 105
  102 CONTINUE
      MSTART = 1
      GO TO 105
  103 CONTINUE
      MSTART = 1
  104 CONTINUE
      MSTOP = M
  105 CONTINUE
      MUNK = MSTOP - MSTART + 1
      NSTART = 1
      NSTOP = N
      GO TO (109,106,107,108,109) NP
  106 CONTINUE
      NSTART = 2
      GO TO 109
  107 CONTINUE
      NSTART = 2
  108 CONTINUE
      NSTOP = NP1
  109 CONTINUE
      NUNK = NSTOP - NSTART + 1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2 + MUNK
      ID4 = ID3 + MUNK
      ID5 = ID4 + MUNK
      ID6 = ID5 + MUNK
      A1 = 2./DLRSQ
      IJ = 0
      IF (MBDCND==3 .OR. MBDCND==4) IJ = 1
      DO I = 1, MUNK
         R = A + FLOAT(I - IJ)*DELTAR
         J = ID5 + I
         W(J) = R
         J = ID6 + I
         W(J) = 1./R**2
         W(I) = (R - DLRBY2)/(R*DLRSQ)
         J = ID3 + I
         W(J) = (R + DLRBY2)/(R*DLRSQ)
         J = ID2 + I
         W(J) = (-A1) + ELMBDA
      END DO
      GO TO (114,111,112,113,114,111) MBDCND
  111 CONTINUE
      W(ID2) = A1
      GO TO 114
  112 CONTINUE
      W(ID2) = A1
  113 CONTINUE
      W(ID3+1) = A1
  114 CONTINUE
      GO TO (115,115,117,117,119,119) MBDCND
  115 CONTINUE
      A1 = W(1)
      F(2,NSTART:NSTOP) = F(2,NSTART:NSTOP) - A1*F(1,NSTART:NSTOP)
      GO TO 119
  117 CONTINUE
      A1 = 2.*DELTAR*W(1)
      F(1,NSTART:NSTOP) = F(1,NSTART:NSTOP) + A1*BDA(NSTART:NSTOP)
  119 CONTINUE
      GO TO (120,122,122,120,120,122) MBDCND
  120 CONTINUE
      A1 = W(ID4)
      F(M,NSTART:NSTOP) = F(M,NSTART:NSTOP) - A1*F(MP1,NSTART:NSTOP)
      GO TO 124
  122 CONTINUE
      A1 = 2.*DELTAR*W(ID4)
      F(MP1,NSTART:NSTOP) = F(MP1,NSTART:NSTOP) - A1*BDB(NSTART:NSTOP)
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  124 CONTINUE
      A1 = 1./DLTHSQ
      L = ID5 - MSTART + 1
      LP = ID6 - MSTART + 1
      GO TO (134,125,125,127,127) NP
  125 CONTINUE
      F(MSTART:MSTOP,2) = F(MSTART:MSTOP,2) - A1*W(MSTART+LP:MSTOP+LP)*F
     1   (MSTART:MSTOP,1)
      GO TO 129
  127 CONTINUE
      A1 = 2./DELTHT
      F(MSTART:MSTOP,1) = F(MSTART:MSTOP,1) + A1*W(MSTART+LP:MSTOP+LP)*
     1   BDC(MSTART:MSTOP)
  129 CONTINUE
      A1 = 1./DLTHSQ
      GO TO (134,130,132,132,130) NP
  130 CONTINUE
      F(MSTART:MSTOP,N) = F(MSTART:MSTOP,N) - A1*W(MSTART+LP:MSTOP+LP)*F
     1   (MSTART:MSTOP,NP1)
      GO TO 134
  132 CONTINUE
      A1 = 2./DELTHT
      F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,NP1) - A1*W(MSTART+LP:MSTOP+
     1   LP)*BDD(MSTART:MSTOP)
  134 CONTINUE
      IF (MBDCND>=5 .AND. NBDCND==3) F(1,1) = F(1,1) - (BDD(2)-BDC(2))*
     1   4./(FLOAT(N)*DELTHT*DLRSQ)
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 11
         ELSE
            IF (NBDCND==0 .OR. NBDCND==3) THEN
               S2 = 0.
               GO TO (144,144,137,144,144,138) MBDCND
  137          CONTINUE
               W(ID5+1) = 0.5*(W(ID5+2)-DLRBY2)
               S2 = 0.25*DELTAR
  138          CONTINUE
               A2 = 2.
               IF (NBDCND == 0) A2 = 1.
               J = ID5 + MUNK
               W(J) = 0.5*(W(J-1)+DLRBY2)
               S = 0.
               DO I = MSTART, MSTOP
                  S1 = 0.
                  IJ = NSTART + 1
                  K = NSTOP - 1
                  S1 = SUM(F(I,IJ:K))
                  J = I + L
                  S = S + (A2*S1 + F(I,NSTART)+F(I,NSTOP))*W(J)
               END DO
               S2=FLOAT(M)*A+DELTAR*(FLOAT((M-1)*(M+1))*0.5+0.25)+S2
               S1 = (2. + A2*FLOAT(NUNK - 2))*S2
               IF (MBDCND /= 3) THEN
                  S2 = FLOAT(N)*A2*DELTAR/8.
                  S = S + F(1,1)*S2
                  S1 = S1 + S2
               ENDIF
               PERTRB = S/S1
               F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:
     1            NSTOP) - PERTRB
            ENDIF
         ENDIF
      ENDIF
  144 CONTINUE
      DO I = MSTART, MSTOP
         K = I - MSTART + 1
         J = I + LP
         A1 = DLTHSQ/W(J)
         W(K) = A1*W(K)
         J = ID2 + K
         W(J) = A1*W(J)
         J = ID3 + K
         W(J) = A1*W(J)
         F(I,NSTART:NSTOP) = A1*F(I,NSTART:NSTOP)
      END DO
      W(1) = 0.
      W(ID4) = 0.
C
C     SOLVE THE SYSTEM OF EQUATIONS.
C
      I1 = 1
      IERR1 = 0
      CALL GENBUNN (NBDCND, NUNK, I1, MUNK, W(1), W(ID2+1), W(ID3+1), 
     1   IDIMF, F(MSTART,NSTART), IERR1, W(ID4+1))
      GO TO (157,157,157,157,148,147) MBDCND
C
C     ADJUST THE SOLUTION AS NECESSARY FOR THE PROBLEMS WHERE A = 0.
C
  147 CONTINUE
      IF (ELMBDA /= 0.) GO TO 148
      YPOLE = 0.
      GO TO 155
  148 CONTINUE
      J = ID5 + MUNK
      W(J) = W(ID2)/W(ID3)
      DO IP = 3, MUNK
         I = MUNK - IP + 2
         J = ID5 + I
         LP = ID2 + I
         K = ID3 + I
         W(J) = W(I)/(W(LP)-W(K)*W(J+1))
      END DO
      W(ID5+1) = -0.5*DLTHSQ/(W(ID2+1)-W(ID3+1)*W(ID5+2))
      DO I = 2, MUNK
         J = ID5 + I
         W(J) = -W(J)*W(J-1)
      END DO
      S = 0.
      S = SUM(F(2,NSTART:NSTOP))
      A2 = NUNK
      IF (NBDCND /= 0) THEN
         S = S - 0.5*(F(2,NSTART)+F(2,NSTOP))
         A2 = A2 - 1.
      ENDIF
      YPOLE = (0.25*DLRSQ*F(1,1)-S/A2)/(W(ID5+1)-1.+ELMBDA*DLRSQ*0.25)
      DO I = MSTART, MSTOP
         K = L + I
         F(I,NSTART:NSTOP) = F(I,NSTART:NSTOP) + YPOLE*W(K)
      END DO
  155 CONTINUE
      F(1,:NP1) = YPOLE
  157 CONTINUE
      IF (NBDCND == 0) THEN
         F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,1)
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE HWSPLRR
C
C     file hwsssp.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE HWSSSP (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,
C    +                   BDPF,ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDTS(N+1),    BDTF(N+1), BDPS(M+1), BDPF(M+1),
C ARGUMENTS              F(IDIMF,N+1)
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION TO
C                        THE HELMHOLTZ EQUATION IN SPHERICAL
C                        COORDINATES AND ON THE SURFACE OF THE UNIT
C                        SPHERE (RADIUS OF 1).  THE EQUATION IS
C
C                          (1/SIN(THETA))(D/DTHETA)(SIN(THETA)
C                          (DU/DTHETA)) + (1/SIN(THETA)**2)(D/DPHI)
C                          (DU/DPHI)  + LAMBDA*U = F(THETA,PHI)
C
C                        WHERE THETA IS COLATITUDE AND PHI IS
C                        LONGITUDE.
C
C USAGE                  CALL HWSSSP (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,
C                                     N,NBDCND,BDPS,BDPF,ELMBDA,F,
C                                     IDIMF,PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               TS,TF
C
C                          THE RANGE OF THETA (COLATITUDE), I.E.,
C                          TS .LE. THETA .LE. TF. TS MUST BE LESS
C                          THAN TF.  TS AND TF ARE IN RADIANS.
C                          A TS OF ZERO CORRESPONDS TO THE NORTH
C                          POLE AND A TF OF PI CORRESPONDS TO
C                          THE SOUTH POLE.
C
C                          * * * IMPORTANT * * *
C
C                          IF TF IS EQUAL TO PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          TF = PIMACH(DUM). THIS INSURES THAT TF
C                          IN THE USER'S PROGRAM IS EQUAL TO PI IN
C                          THIS PROGRAM WHICH PERMITS SEVERAL TESTS
C                          OF THE INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (TS,TF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS IN THE
C                          THETA-DIRECTION GIVEN BY
C                          THETA(I) = (I-1)DTHETA+TS FOR
C                          I = 1,2,...,M+1, WHERE
C                          DTHETA = (TF-TS)/M IS THE PANEL WIDTH.
C                          M MUST BE GREATER THAN 5
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT THETA = TS AND THETA = TF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THETA = TF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               SPECIFIED AT THETA = TS AND
C                               THETA = TF (SEE NOTES 1,2 BELOW).
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS SPECIFIED AT
C                               THETA = TF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE SOLUTION
C                               IS SPECIFIED AT THETA = TF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO THETA
C                               IS SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE SOLUTION IS
C                               IS UNSPECIFIED AT THETA = TF = PI.
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW) AND
C                               THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TF = PI.
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THETA = TF = PI.
C
C                          NOTES:
C                          IF TS = 0, DO NOT USE MBDCND = 3,4, OR 8,
C                          BUT INSTEAD USE MBDCND = 5,6, OR 9  .
C
C                          IF TF = PI, DO NOT USE MBDCND = 2,3, OR 6,
C                          BUT INSTEAD USE MBDCND = 7,8, OR 9  .
C
C                        BDTS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TS.  WHEN MBDCND = 3,4, OR 8,
C
C                          BDTS(J) = (D/DTHETA)U(TS,PHI(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTS IS
C                          A DUMMY VARIABLE.
C
C                        BDTF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TF.  WHEN MBDCND = 2,3, OR 6,
C
C                          BDTF(J) = (D/DTHETA)U(TF,PHI(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTF IS
C                          A DUMMY VARIABLE.
C
C                        PS,PF
C                          THE RANGE OF PHI (LONGITUDE), I.E.,
C                          PS .LE. PHI .LE. PF.  PS MUST BE LESS
C                          THAN PF.  PS AND PF ARE IN RADIANS.
C                          IF PS = 0 AND PF = 2*PI, PERIODIC
C                          BOUNDARY CONDITIONS ARE USUALLY PRESCRIBED.
C
C                          * * * IMPORTANT * * *
C
C                          IF PF IS EQUAL TO 2*PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          PF = 2.*PIMACH(DUM). THIS INSURES THAT
C                          PF IN THE USERS PROGRAM IS EQUAL TO
C                          2*PI IN THIS PROGRAM WHICH PERMITS TESTS
C                          OF THE INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (PS,PF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE PHI-DIRECTION GIVEN BY
C                          PHI(J) = (J-1)DPHI+PS  FOR
C                          J = 1,2,...,N+1, WHERE
C                          DPHI = (PF-PS)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 4
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT PHI = PS AND PHI = PF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN PHI,
C                               I.U., U(I,J) = U(I,N+J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = PS AND PHI = PF
C                               (SEE NOTE BELOW).
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = PS (SEE NOTE BELOW)
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = PF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = PS AND PHI = PF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PS AND THE SOLUTION IS SPECIFIED
C                               AT PHI = PF
C
C                          NOTE:
C                          NBDCND = 1,2, OR 4 CANNOT BE USED WITH
C                          MBDCND = 5,6,7,8, OR 9.  THE FORMER INDICATES
C                          THAT THE SOLUTION IS SPECIFIED AT A POLE, THE
C                          LATTER INDICATES THAT THE SOLUTION IS NOT
C                          SPECIFIED.  USE INSTEAD  MBDCND = 1 OR 2.
C
C                        BDPS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO PHI AT
C                          PHI = PS.  WHEN NBDCND = 3 OR 4,
C
C                            BDPS(I) = (D/DPHI)U(THETA(I),PS),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDPS IS
C                          A DUMMY VARIABLE.
C
C                        BDPF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO PHI AT
C                          PHI = PF.  WHEN NBDCND = 2 OR 3,
C
C                            BDPF(I) = (D/DPHI)U(THETA(I),PF),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDPF IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSSSP WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUE OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION AND BOUNDARY VALUES (IF ANY).
C                          F MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(THETA(I),PHI(J)).
C
C                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
C                          FOR J = 1,2,...,N+1 AND I = 1,2,...,M+1
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   ------------      ------------
C
C                            1      U(TS,PHI(J))      U(TF,PHI(J))
C                            2      U(TS,PHI(J))      F(TF,PHI(J))
C                            3      F(TS,PHI(J))      F(TF,PHI(J))
C                            4      F(TS,PHI(J))      U(TF,PHI(J))
C                            5      F(0,PS)           U(TF,PHI(J))
C                            6      F(0,PS)           F(TF,PHI(J))
C                            7      U(TS,PHI(J))      F(PI,PS)
C                            8      F(TS,PHI(J))      F(PI,PS)
C                            9      F(0,PS)           F(PI,PS)
C
C                          NBDCND   F(I,1)            F(I,N+1)
C                          ------   --------------    --------------
C
C                            0      F(THETA(I),PS)    F(THETA(I),PS)
C                            1      U(THETA(I),PS)    U(THETA(I),PF)
C                            2      U(THETA(I),PS)    F(THETA(I),PF)
C                            3      F(THETA(I),PS)    F(THETA(I),PF)
C                            4      F(THETA(I),PS)    U(THETA(I),PF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION U
C                          AND THE RIGHT SIDE F AT A CORNER THEN THE
C                          SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSSSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F. IDIMF MUST BE
C                          AT LEAST M+1  .
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),PHI(J)),  I = 1,2,...,M+1  AND
C                          J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF ONE SPECIFIES A COMBINATION OF PERIODIC,
C                          DERIVATIVE OR UNSPECIFIED BOUNDARY
C                          CONDITIONS FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HWSSSP THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION IS NOT UNIQUE AND IS
C                          UNNORMALIZED. THE VALUE OF PERTRB SHOULD
C                          BE SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE , A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM. THIS
C                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
C                          THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 8,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR
C                          = 1  TS.LT.0 OR TF.GT.PI
C                          = 2  TS.GE.TF
C                          = 3  MBDCND.LT.1 OR MBDCND.GT.9
C                          = 4  PS.LT.0 OR PS.GT.PI+PI
C                          = 5  PS.GE.PF
C                          = 6  N.LT.5
C                          = 7  M.LT.5
C                          = 8  NBDCND.LT.0 OR NBDCND.GT.4
C                          = 9  ELMBDA.GT.0
C                          = 10 IDIMF.LT.M+1
C                          = 11 NBDCND EQUALS 1,2 OR 4 AND MBDCND.GE.5
C                          = 12 TS.EQ.0 AND MBDCND EQUALS 3,4 OR 8
C                          = 13 TF.EQ.PI AND MBDCND EQUALS 2,3 OR 6
C                          = 14 MBDCND EQUALS 5,6 OR 9 AND TS.NE.0
C                          = 15 MBDCND.GE.7 AND TF.NE.PI
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,genbun.f,gnbnaux.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS GENBUN TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE  M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          M*N*(LOG2(N)
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS NBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
C                        OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR N
C                        AND M AS LARGE AS 64.  MORE DETAILS ABOUT
C                        ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE GENBUN WHICH IS THE ROUTINE THAT
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             P. N. SWARZTRAUBER, "THE DIRECT SOLUTION OF
C                        THE DISCRETE POISSON EQUATION ON THE SURFACE OF
C                        A SPHERE", S.I.A.M. J. NUMER. ANAL.,15(1974),
C                        PP 212-215.
C
C                        SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS", NCAR TN/IA-109, JULY,
C                        1975, 138 PP.
C***********************************************************************
      SUBROUTINE HWSSSP(TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND
     1   , BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: TS
      REAL  :: TF
      REAL  :: PS
      REAL  :: PF
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDTS(*)
      REAL  :: BDTF(*)
      REAL  :: BDPS(*)
      REAL  :: BDPF(*)
      REAL  :: F(IDIMF,1)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NBR, IRWK, ICWK
      REAL :: PI, DUM, TPI
C-----------------------------------------------
C
      NBR = NBDCND + 1
      PI = 4.0*ATAN(1.0)
      TPI = 2.*PI
      IERROR = 0
      IF (TS<0. .OR. TF>PI) IERROR = 1
      IF (TS >= TF) IERROR = 2
      IF (MBDCND<1 .OR. MBDCND>9) IERROR = 3
      IF (PS<0. .OR. PF>TPI) IERROR = 4
      IF (PS >= PF) IERROR = 5
      IF (N < 5) IERROR = 6
      IF (M < 5) IERROR = 7
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 8
      IF (ELMBDA > 0.) IERROR = 9
      IF (IDIMF < M + 1) IERROR = 10
      IF ((NBDCND==1 .OR. NBDCND==2 .OR. NBDCND==4) .AND. MBDCND>=5) 
     1   IERROR = 11
      IF(TS==0..AND.(MBDCND==3.OR.MBDCND==4.OR.MBDCND==8))IERROR=12
      IF(TF==PI.AND.(MBDCND==2.OR.MBDCND==3.OR.MBDCND==6))IERROR=13
      IF((MBDCND==5.OR.MBDCND==6.OR.MBDCND==9).AND.TS/=0.)IERROR=14
      IF (MBDCND>=7 .AND. TF/=PI) IERROR = 15
      IF (IERROR/=0 .AND. IERROR/=9) RETURN 
!     allocate generous work space estimate
      IRWK=4*(N+1)+(16+INT(ALOG(FLOAT(N+1))/ALOG(2.0)))*(M+1)
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hwssspp(TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,
     +             BDPF,ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HWSSSP


 
      SUBROUTINE HWSSSPP(TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, 
     1   NBDCND, BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: TS
      REAL  :: TF
      REAL  :: PS
      REAL  :: PF
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDTS(*)
      REAL  :: BDTF(*)
      REAL  :: BDPS(*)
      REAL  :: BDPF(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
      CALL HWSSS1 (TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND, 
     1   BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, W, W(M+2), W(2*M+3), W(3*
     2   M+4), W(4*M+5), W(5*M+6), W(6*M+7))
      RETURN 
      END SUBROUTINE HWSSSPP


 
      SUBROUTINE HWSSS1(TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND
     1   , BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, AM, BM, CM, SN, SS, 
     2   SINT, D)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      REAL , INTENT(IN) :: TS
      REAL , INTENT(IN) :: TF
      REAL , INTENT(IN) :: PS
      REAL , INTENT(IN) :: PF
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDTS(*)
      REAL , INTENT(IN) :: BDTF(*)
      REAL , INTENT(IN) :: BDPS(*)
      REAL , INTENT(IN) :: BDPF(*)
      REAL  :: F(IDIMF,*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL , INTENT(INOUT) :: SN(*)
      REAL , INTENT(INOUT) :: SS(*)
      REAL , INTENT(INOUT) :: SINT(*)
      REAL  :: D(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MP1, NP1, I, INP, ISP, MBR, ITS, ITF, ITSP, ITFM, MUNK
     1   , IID, II, NBR, JPS, JPF, JPSP, JPFM, NUNK, ISING, J, IERROR
      REAL :: PI, DUM, TPI, HPI, FN, FM, DTH, HDTH, TDT, DPHI, TDP, 
     1   DPHI2, EDP2, DTH2, CP, WP, FIM1, THETA, T1, AT, CT, WTS, WTF, 
     2   WPS, WPF, FJJ, CF, SUM, SUM1, HNE, YHLD, SUM2, DFN, DNN, DSN, 
     3   CNP, HLD, DFS, DSS, DNS, CSP, RTN, RTS, DEN
C-----------------------------------------------
C
      PI = 4.0*ATAN(1.0)
      TPI = PI + PI
      HPI = PI/2.
      MP1 = M + 1
      NP1 = N + 1
      FN = N
      FM = M
      DTH = (TF - TS)/FM
      HDTH = DTH/2.
      TDT = DTH + DTH
      DPHI = (PF - PS)/FN
      TDP = DPHI + DPHI
      DPHI2 = DPHI*DPHI
      EDP2 = ELMBDA*DPHI2
      DTH2 = DTH*DTH
      CP = 4./(FN*DTH2)
      WP = FN*SIN(HDTH)/4.
      DO I = 1, MP1
         FIM1 = I - 1
         THETA = FIM1*DTH + TS
         SINT(I) = SIN(THETA)
         IF (SINT(I) == 0.) CYCLE 
         T1 = 1./(DTH2*SINT(I))
         AM(I) = T1*SIN(THETA - HDTH)
         CM(I) = T1*SIN(THETA + HDTH)
         BM(I) = (-AM(I)) - CM(I) + ELMBDA
      END DO
      INP = 0
      ISP = 0
C
C BOUNDARY CONDITION AT THETA=TS
C
      MBR = MBDCND + 1
      GO TO (103,104,104,105,105,106,106,104,105,106) MBR
  103 CONTINUE
      ITS = 1
      GO TO 107
  104 CONTINUE
      AT = AM(2)
      ITS = 2
      GO TO 107
  105 CONTINUE
      AT = AM(1)
      ITS = 1
      CM(1) = AM(1) + CM(1)
      GO TO 107
  106 CONTINUE
      AT = AM(2)
      INP = 1
      ITS = 2
C
C BOUNDARY CONDITION THETA=TF
C
  107 CONTINUE
      GO TO (108,109,110,110,109,109,110,111,111,111) MBR
  108 CONTINUE
      ITF = M
      GO TO 112
  109 CONTINUE
      CT = CM(M)
      ITF = M
      GO TO 112
  110 CONTINUE
      CT = CM(M+1)
      AM(M+1) = AM(M+1) + CM(M+1)
      ITF = M + 1
      GO TO 112
  111 CONTINUE
      ITF = M
      ISP = 1
      CT = CM(M)
C
C COMPUTE HOMOGENEOUS SOLUTION WITH SOLUTION AT POLE EQUAL TO ONE
C
  112 CONTINUE
      ITSP = ITS + 1
      ITFM = ITF - 1
      WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
      WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
      MUNK = ITF - ITS + 1
      IF (ISP > 0) THEN
         D(ITS) = CM(ITS)/BM(ITS)
         DO I = ITSP, M
            D(I) = CM(I)/(BM(I)-AM(I)*D(I-1))
         END DO
         SS(M) = -D(M)
         IID = M - ITS
         DO II = 1, IID
            I = M - II
            SS(I) = -D(I)*SS(I+1)
         END DO
         SS(M+1) = 1.
      ENDIF
      IF (INP > 0) THEN
         SN(1) = 1.
         D(ITF) = AM(ITF)/BM(ITF)
         IID = ITF - 2
         DO II = 1, IID
            I = ITF - II
            D(I) = AM(I)/(BM(I)-CM(I)*D(I+1))
         END DO
         SN(2) = -D(2)
         DO I = 3, ITF
            SN(I) = -D(I)*SN(I-1)
         END DO
      ENDIF
C
C BOUNDARY CONDITIONS AT PHI=PS
C
      NBR = NBDCND + 1
      WPS = 1.
      WPF = 1.
      SELECT CASE (NBR) 
      CASE DEFAULT
         JPS = 1
      CASE (2:3) 
         JPS = 2
      CASE (4:5) 
         JPS = 1
         WPS = 0.5
      END SELECT
C
C BOUNDARY CONDITION AT PHI=PF
C
  124 CONTINUE
      GO TO (125,126,127,127,126) NBR
  125 CONTINUE
      JPF = N
      GO TO 128
  126 CONTINUE
      JPF = N
      GO TO 128
  127 CONTINUE
      WPF = 0.5
      JPF = N + 1
  128 CONTINUE
      JPSP = JPS + 1
      JPFM = JPF - 1
      NUNK = JPF - JPS + 1
      FJJ = JPFM - JPSP + 1
C
C SCALE COEFFICIENTS FOR SUBROUTINE GENBUN
C
      DO I = ITS, ITF
         CF = DPHI2*SINT(I)*SINT(I)
         AM(I) = CF*AM(I)
         BM(I) = CF*BM(I)
         CM(I) = CF*CM(I)
      END DO
      AM(ITS) = 0.
      CM(ITF) = 0.
      ISING = 0
      GO TO (130,138,138,130,138,138,130,138,130,130) MBR
  130 CONTINUE
      GO TO (131,138,138,131,138) NBR
  131 CONTINUE
      IF (ELMBDA >= 0.) THEN
         ISING = 1
         SUM = WTS*WPS + WTS*WPF + WTF*WPS + WTF*WPF
         IF (INP > 0) THEN
            SUM = SUM + WP
         ENDIF
         IF (ISP > 0) THEN
            SUM = SUM + WP
         ENDIF
         SUM1 = 0.
         DO I = ITSP, ITFM
            SUM1 = SUM1 + SINT(I)
         END DO
         SUM = SUM + FJJ*(SUM1 + WTS + WTF)
         SUM = SUM + (WPS + WPF)*SUM1
         HNE = SUM
      ENDIF
  138 CONTINUE
      GO TO (146,142,142,144,144,139,139,142,144,139) MBR
  139 CONTINUE
      IF (NBDCND - 3 /= 0) GO TO 146
      YHLD = F(1,JPS) - 4./(FN*DPHI*DTH2)*(BDPF(2)-BDPS(2))
      F(1,:NP1) = YHLD
      GO TO 146
  142 CONTINUE
      F(2,JPS:JPF) = F(2,JPS:JPF) - AT*F(1,JPS:JPF)
      GO TO 146
  144 CONTINUE
      F(1,JPS:JPF) = F(1,JPS:JPF) + TDT*BDTS(JPS:JPF)*AT
  146 CONTINUE
      GO TO (154,150,152,152,150,150,152,147,147,147) MBR
  147 CONTINUE
      IF (NBDCND - 3 /= 0) GO TO 154
      YHLD = F(M+1,JPS) - 4./(FN*DPHI*DTH2)*(BDPF(M)-BDPS(M))
      F(M+1,:NP1) = YHLD
      GO TO 154
  150 CONTINUE
      F(M,JPS:JPF) = F(M,JPS:JPF) - CT*F(M+1,JPS:JPF)
      GO TO 154
  152 CONTINUE
      F(M+1,JPS:JPF) = F(M+1,JPS:JPF) - TDT*BDTF(JPS:JPF)*CT
  154 CONTINUE
      GO TO (159,155,155,157,157) NBR
  155 CONTINUE
      F(ITS:ITF,2) = F(ITS:ITF,2) - F(ITS:ITF,1)/(DPHI2*SINT(ITS:ITF)*
     1   SINT(ITS:ITF))
      GO TO 159
  157 CONTINUE
      F(ITS:ITF,1) = F(ITS:ITF,1) + TDP*BDPS(ITS:ITF)/(DPHI2*SINT(ITS:
     1   ITF)*SINT(ITS:ITF))
  159 CONTINUE
      GO TO (164,160,162,162,160) NBR
  160 CONTINUE
      F(ITS:ITF,N) = F(ITS:ITF,N) - F(ITS:ITF,N+1)/(DPHI2*SINT(ITS:ITF)*
     1   SINT(ITS:ITF))
      GO TO 164
  162 CONTINUE
      F(ITS:ITF,N+1) = F(ITS:ITF,N+1) - TDP*BDPF(ITS:ITF)/(DPHI2*SINT(
     1   ITS:ITF)*SINT(ITS:ITF))
  164 CONTINUE
      PERTRB = 0.
      IF (ISING /= 0) THEN
         SUM = WTS*WPS*F(ITS,JPS) + WTS*WPF*F(ITS,JPF) + WTF*WPS*F(ITF,
     1      JPS) + WTF*WPF*F(ITF,JPF)
         IF (INP > 0) THEN
            SUM = SUM + WP*F(1,JPS)
         ENDIF
         IF (ISP > 0) THEN
            SUM = SUM + WP*F(M+1,JPS)
         ENDIF
         DO I = ITSP, ITFM
            SUM1 = 0.
            DO J = JPSP, JPFM
               SUM1 = SUM1 + F(I,J)
            END DO
            SUM = SUM + SINT(I)*SUM1
         END DO
         SUM1 = 0.
         SUM2 = 0.
         DO J = JPSP, JPFM
            SUM1 = SUM1 + F(ITS,J)
            SUM2 = SUM2 + F(ITF,J)
         END DO
         SUM = SUM + WTS*SUM1 + WTF*SUM2
         SUM1 = 0.
         SUM2 = 0.
         SUM1 = DOT_PRODUCT(SINT(ITSP:ITFM),F(ITSP:ITFM,JPS))
         SUM2 = DOT_PRODUCT(SINT(ITSP:ITFM),F(ITSP:ITFM,JPF))
         SUM = SUM + WPS*SUM1 + WPF*SUM2
         PERTRB = SUM/HNE
         F(:MP1,:NP1) = F(:MP1,:NP1) - PERTRB
      ENDIF
C
C SCALE RIGHT SIDE FOR SUBROUTINE GENBUN
C
      DO I = ITS, ITF
         CF = DPHI2*SINT(I)*SINT(I)
         F(I,JPS:JPF) = CF*F(I,JPS:JPF)
      END DO
      CALL GENBUNN (NBDCND, NUNK, 1, MUNK, AM(ITS), BM(ITS), CM(ITS), 
     1   IDIMF, F(ITS,JPS), IERROR, D)
      IF (ISING <= 0) GO TO 186
      IF (INP > 0) THEN
         IF (ISP > 0) GO TO 186
         F(1,:NP1) = 0.
         GO TO 209
      ENDIF
      IF (ISP <= 0) GO TO 186
      F(M+1,:NP1) = 0.
      GO TO 209
  186 CONTINUE
      IF (INP > 0) THEN
         SUM = WPS*F(ITS,JPS) + WPF*F(ITS,JPF)
         DO J = JPSP, JPFM
            SUM = SUM + F(ITS,J)
         END DO
         DFN = CP*SUM
         DNN = CP*((WPS + WPF + FJJ)*(SN(2)-1.)) + ELMBDA
         DSN = CP*(WPS + WPF + FJJ)*SN(M)
         IF (ISP > 0) GO TO 194
         CNP = (F(1,1)-DFN)/DNN
         DO I = ITS, ITF
            HLD = CNP*SN(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
         END DO
         F(1,:NP1) = CNP
         GO TO 209
      ENDIF
      IF (ISP <= 0) GO TO 209
  194 CONTINUE
      SUM = WPS*F(ITF,JPS) + WPF*F(ITF,JPF)
      DO J = JPSP, JPFM
         SUM = SUM + F(ITF,J)
      END DO
      DFS = CP*SUM
      DSS = CP*((WPS + WPF + FJJ)*(SS(M)-1.)) + ELMBDA
      DNS = CP*(WPS + WPF + FJJ)*SS(2)
      IF (INP <= 0) THEN
         CSP = (F(M+1,1)-DFS)/DSS
         DO I = ITS, ITF
            HLD = CSP*SS(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
         END DO
         F(M+1,:NP1) = CSP
      ELSE
         RTN = F(1,1) - DFN
         RTS = F(M+1,1) - DFS
         IF (ISING > 0) THEN
            CSP = 0.
            CNP = RTN/DNN
         ELSE
            IF (ABS(DNN) - ABS(DSN) > 0.) THEN
               DEN = DSS - DNS*DSN/DNN
               RTS = RTS - RTN*DSN/DNN
               CSP = RTS/DEN
               CNP = (RTN - CSP*DNS)/DNN
            ELSE
               DEN = DNS - DSS*DNN/DSN
               RTN = RTN - RTS*DNN/DSN
               CSP = RTN/DEN
               CNP = (RTS - DSS*CSP)/DSN
            ENDIF
         ENDIF
         DO I = ITS, ITF
            HLD = CNP*SN(I) + CSP*SS(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
         END DO
         F(1,:NP1) = CNP
         F(M+1,:NP1) = CSP
      ENDIF
  209 CONTINUE
      IF (NBDCND == 0) THEN
         F(:MP1,JPF+1) = F(:MP1,JPS)
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE HWSSS1
C
C     file pois3d.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,
C    +                   MDIMF,F,IERROR)
C
C
C DIMENSION OF           A(N), B(N), C(N), F(LDIMF,MDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
C                        FOR UNKNOWN X VALUES, WHERE I=1,2,...,L,
C                        J=1,2,...,M, AND K=1,2,...,N
C
C                        C1*(X(I-1,J,K) -2.*X(I,J,K) +X(I+1,J,K)) +
C                        C2*(X(I,J-1,K) -2.*X(I,J,K) +X(I,J+1,K)) +
C                        A(K)*X(I,J,K-1) +B(K)*X(I,J,K)+ C(K)*X(I,J,K+1)
C                        = F(I,J,K)
C
C                        THE INDICES K-1 AND K+1 ARE EVALUATED MODULO N,
C                        I.E. X(I,J,0)=X(I,J,N) AND X(I,J,N+1)=X(I,J,1).
C                        THE UNKNOWNS
C                        X(0,J,K), X(L+1,J,K), X(I,0,K), AND X(I,M+1,K)
C                        ARE ASSUMED TO TAKE ON CERTAIN PRESCRIBED
C                        VALUES DESCRIBED BELOW.
C
C USAGE                  CALL POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,
C                        N,A,B,C,LDIMF,MDIMF,F,IERROR)
C
C ARGUMENTS
C
C ON INPUT
C                        LPEROD
C                          INDICATES THE VALUES THAT X(0,J,K) AND
C                          X(L+1,J,K) ARE ASSUMED TO HAVE.
C                          = 0  X(0,J,K)=X(L,J,K), X(L+1,J,K)=X(1,J,K)
C                          = 1  X(0,J,K) = 0,      X(L+1,J,K) = 0
C                          = 2  X(0,J,K)=0,        X(L+1,J,K)=X(L-1,J,K)
C                          = 3  X(0,J,K)=X(2,J,K), X(L+1,J,K)=X(L-1,J,K)
C                          = 4  X(0,J,K)=X(2,J,K), X(L+1,J,K) = 0.
C
C                        L
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                          L MUST BE AT LEAST 3.
C
C                        C1
C                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
C                          OF EQUATIONS TO BE SOLVED.
C
C                        MPEROD
C                          INDICATES THE VALUES THAT X(I,0,K) AND
C                          X(I,M+1,K) ARE ASSUMED TO HAVE.
C                          = 0  X(I,0,K)=X(I,M,K), X(I,M+1,K)=X(I,1,K)
C                          = 1  X(I,0,K)=0,        X(I,M+1,K)=0
C                          = 2  X(I,0,K)=0,        X(I,M+1,K)=X(I,M-1,K)
C                          = 3  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=X(I,M-1,K)
C                          = 4  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=0
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          M MUST BE AT LEAST 3.
C
C                        C2
C                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
C                          OF EQUATIONS TO BE SOLVED.
C
C                        NPEROD
C                          = 0  IF A(1) AND C(N) ARE NOT ZERO.
C                          = 1  IF A(1) = C(N) = 0.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE K-DIRECTION.
C                          N MUST BE AT LEAST 3.
C
C                        A, B, C
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.
C
C                          IF NPEROD = 0 THE ARRAY ELEMENTS MUST NOT
C                          DEPEND UPON INDEX K, BUT MUST BE CONSTANT.
C                          SPECIFICALLY,THE SUBROUTINE CHECKS THE
C                          FOLLOWING CONDITION
C                            A(K) = C(1)
C                            C(K) = C(1)
C                            B(K) = B(1)
C                          FOR K=1,2,...,N.
C
C                        LDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE THREE-
C                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
C                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
C                          USED TO SPECIFY THE VARIABLE DIMENSION
C                          OF F.  LDIMF MUST BE AT LEAST L.
C
C                        MDIMF
C                          THE COLUMN (OR SECOND) DIMENSION OF THE THREE
C                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
C                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
C                          USED TO SPECIFY THE VARIABLE DIMENSION
C                          OF F.  MDIMF MUST BE AT LEAST M.
C
C                        F
C                          A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
C                          OF EQUATIONS GIVEN ABOVE.  F MUST BE
C                          DIMENSIONED AT LEAST L X M X N.
C
C ON OUTPUT
C
C                        F
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
C                          SOLUTION IS NOT ATTEMPTED.
C                          = 0  NO ERROR
C                          = 1  IF LPEROD .LT. 0 OR .GT. 4
C                          = 2  IF L .LT. 3
C                          = 3  IF MPEROD .LT. 0 OR .GT. 4
C                          = 4  IF M .LT. 3
C                          = 5  IF NPEROD .LT. 0 OR .GT. 1
C                          = 6  IF N .LT. 3
C                          = 7  IF LDIMF .LT. L
C                          = 8  IF MDIMF .LT. M
C                          = 9  IF A(K) .NE. C(1) OR C(K) .NE. C(1)
C                               OR B(I) .NE.B(1) FOR SOME K=1,2,...,N.
C                          = 10 IF NPEROD = 1 AND A(1) .NE. 0
C                               OR C(N) .NE. 0
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
C                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
C                          SHOULD TEST IERROR AFTER THE CALL.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,comf.f,fftpack.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY, 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE SOLVES THREE-DIMENSIONAL BLOCK
C                        TRIDIAGONAL LINEAR SYSTEMS ARISING FROM FINITE
C                        DIFFERENCE APPROXIMATIONS TO THREE-DIMENSIONAL
C                        POISSON EQUATIONS USING THE FFT PACKAGE
C                        FFTPACK WRITTEN BY PAUL SWARZTRAUBER.
C
C TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          L*M*N*(LOG2(L)+LOG2(M)+5)
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS LPEROD
C                        AND MPEROD.
C
C ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
C                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
C                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
C                        IN THE 'PURPOSE' SECTION WITH
C                          A(K) = C(K) = -0.5*B(K) = 1,  K=1,2,...,N
C                        AND, WHEN NPEROD = 1
C                          A(1) = C(N) = 0
C                          A(N) = C(1) = 2.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
C                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT
C                        SIDE Y WAS COMPUTED.  USING THIS ARRAY Y
C                        SUBROUTINE POIS3D WAS CALLED TO PRODUCE AN
C                        APPROXIMATE SOLUTION Z.  RELATIVE ERROR
C
C                        E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K
C
C                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
C                        OVER I=1,2,...,L, J=1,2,...,M AND K=1,2,...,N.
C                        VALUES OF E ARE GIVEN IN THE TABLE BELOW FOR
C                        SOME TYPICAL VALUES OF L,M AND N.
C
C                        L(=M=N)   LPEROD    MPEROD       E
C                        ------    ------    ------     ------
C
C                          16        0         0        1.E-13
C                          15        1         1        4.E-13
C                          17        3         3        2.E-13
C                          32        0         0        2.E-13
C                          31        1         1        2.E-12
C                          33        3         3        7.E-13
C
C REFERENCES              NONE
C ********************************************************************
      SUBROUTINE POIS3D(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, C
     1   , LDIMF, MDIMF, F, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: LPEROD
      INTEGER  :: L
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL  :: C1
      REAL  :: C2
      REAL  :: A(*)
      REAL  :: B(*)
      REAL  :: C(*)
      REAL  :: F(LDIMF,MDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: LP, MP, NP, K, IRWK, ICWK
      REAL, DIMENSION(6) :: SAVE
C-----------------------------------------------
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
C
C     CHECK FOR INVALID INPUT.
C
      IERROR = 0
      IF (LP<1 .OR. LP>5) IERROR = 1
      IF (L < 3) IERROR = 2
      IF (MP<1 .OR. MP>5) IERROR = 3
      IF (M < 3) IERROR = 4
      IF (NP<1 .OR. NP>2) IERROR = 5
      IF (N < 3) IERROR = 6
      IF (LDIMF < L) IERROR = 7
      IF (MDIMF < M) IERROR = 8
      IF (NP == 1) THEN
         DO K = 1, N
            IF (A(K) /= C(1)) GO TO 102
            IF (C(K) /= C(1)) GO TO 102
            IF (B(K) /= B(1)) GO TO 102
         END DO
         GO TO 104
  102    CONTINUE
         IERROR = 9
      ENDIF
      IF (NPEROD==1 .AND. (A(1)/=0. .OR. C(N)/=0.)) IERROR = 10
c 104 IF (IERROR .NE. 0) GO TO 122
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
!     allocate required work space length (generous estimate)
      IRWK=30+L+M+2*N+MAX0(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call pois3dd(LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,
     +             MDIMF,F,IERROR,w%rew,w%ifac,w%ifac2)
!     release work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE POIS3D


 
      SUBROUTINE POIS3DD(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, 
     1   C, LDIMF, MDIMF, F, IERROR, W, IFACX, IFACY)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: LPEROD
      INTEGER  :: L
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL  :: C1
      REAL  :: C2
      REAL  :: A(*)
      REAL  :: B(*)
      REAL  :: C(*)
      REAL  :: F(LDIMF,MDIMF,*)
      REAL  :: W(*)
      INTEGER:: IFACX(*)
      INTEGER:: IFACY(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: LP, MP, NP, IWYRT, IWT, IWD, IWBB, IWX, IWY, NH, NHM1, 
     1   NODD, I, J, K, NHPK, NHMK
      REAL, DIMENSION(6) :: SAVE
C-----------------------------------------------
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
      IWYRT = L + 1
      IWT = IWYRT + M
      IWD = IWT + MAX0(L,M,N) + 1
      IWBB = IWD + N
      IWX = IWBB + N
      IWY = IWX + 7*((L + 1)/2) + 15
      GO TO (105,114) NP
C
C     REORDER UNKNOWNS WHEN NPEROD = 0.
C
  105 CONTINUE
      NH = (N + 1)/2
      NHM1 = NH - 1
      NODD = 1
      IF (2*NH == N) NODD = 2
      DO I = 1, L
         DO J = 1, M
            DO K = 1, NHM1
               W(K) = F(I,J,NH-K) - F(I,J,K+NH)
               W(K+NH) = F(I,J,NH-K) + F(I,J,K+NH)
            END DO
            W(NH) = 2.*F(I,J,NH)
            GO TO (108,107) NODD
  107       CONTINUE
            W(N) = 2.*F(I,J,N)
  108       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      SAVE(1) = C(NHM1)
      SAVE(2) = A(NH)
      SAVE(3) = C(NH)
      SAVE(4) = B(NHM1)
      SAVE(5) = B(N)
      SAVE(6) = A(N)
      C(NHM1) = 0.
      A(NH) = 0.
      C(NH) = 2.*C(NH)
      SELECT CASE (NODD) 
      CASE DEFAULT
         B(NHM1) = B(NHM1) - A(NH-1)
         B(N) = B(N) + A(N)
      CASE (2) 
         A(N) = C(NH)
      END SELECT
  114 CONTINUE

      call POS3D1(LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, W, 
     1       W(IWYRT), W(IWT), W(IWD), W(IWX), IFACX, W(IWY), IFACY,
     2       C1, C2, W(IWBB))

      GO TO (115,122) NP
  115 CONTINUE
      DO I = 1, L
         DO J = 1, M
            W(NH-1:NH-NHM1:(-1))=0.5*(F(I,J,NH+1:NHM1+NH)+F(I,J,:NHM1))
            W(NH+1:NHM1+NH) = 0.5*(F(I,J,NH+1:NHM1+NH)-F(I,J,:NHM1))
            W(NH) = 0.5*F(I,J,NH)
            GO TO (118,117) NODD
  117       CONTINUE
            W(N) = 0.5*F(I,J,N)
  118       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      C(NHM1) = SAVE(1)
      A(NH) = SAVE(2)
      C(NH) = SAVE(3)
      B(NHM1) = SAVE(4)
      B(N) = SAVE(5)
      A(N) = SAVE(6)
  122 CONTINUE
      RETURN 
      END SUBROUTINE POIS3DD


      SUBROUTINE POS3D1(LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, XRT, 
     1   YRT, T, D, WX, IFACX, WY, IFACY, C1, C2, BB)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: LP
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: MP
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: LDIMF
      INTEGER , INTENT(IN) :: MDIMF
      REAL , INTENT(IN) :: C1
      REAL , INTENT(IN) :: C2
      REAL  :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL  :: C(*)
      REAL , INTENT(INOUT) :: F(LDIMF,MDIMF,1)
      REAL , INTENT(INOUT) :: XRT(*)
      REAL , INTENT(INOUT) :: YRT(*)
      REAL  :: T(*)
      REAL  :: D(*)
      REAL  :: WX(*)
      INTEGER:: IFACX(*)
      REAL  :: WY(*)
      INTEGER:: IFACY(*)
      REAL  :: BB(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: LR, MR, NR, LRDEL, I, MRDEL, J, IFWRD, IS, K
      REAL :: PI, DUM, SCALX, DX, DI, SCALY, DY, DJ
C-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      LR = L
      MR = M
      NR = N
C
C     GENERATE TRANSFORM ROOTS
C
      LRDEL = ((LP - 1)*(LP - 3)*(LP - 5))/3
      SCALX = LR + LRDEL
      DX = PI/(2.*SCALX)
      GO TO (108,103,101,102,101) LP
  101 CONTINUE
      DI = 0.5
      SCALX = 2.*SCALX
      GO TO 104
  102 CONTINUE
      DI = 1.0
      GO TO 104
  103 CONTINUE
      DI = 0.0
  104 CONTINUE
      DO I = 1, LR
         XRT(I) = -4.*C1*SIN((FLOAT(I) - DI)*DX)**2
      END DO
      SCALX = 2.*SCALX
      GO TO (112,106,110,107,111) LP
  106 CONTINUE
      CALL SINTI (LR, WX, IFACX)
      GO TO 112
  107 CONTINUE
      CALL COSTI (LR, WX, IFACX)
      GO TO 112
  108 CONTINUE
      XRT(1) = 0.
      XRT(LR) = -4.*C1
      DO I = 3, LR, 2
         XRT(I-1) = -4.*C1*SIN(FLOAT(I - 1)*DX)**2
         XRT(I) = XRT(I-1)
      END DO
      CALL RFFTI (LR, WX, IFACX)
      GO TO 112
  110 CONTINUE
      CALL SINQI (LR, WX, IFACX)
      GO TO 112
  111 CONTINUE
      CALL COSQI (LR, WX, IFACX)
  112 CONTINUE
      MRDEL = ((MP - 1)*(MP - 3)*(MP - 5))/3
      SCALY = MR + MRDEL
      DY = PI/(2.*SCALY)
      GO TO (120,115,113,114,113) MP
  113 CONTINUE
      DJ = 0.5
      SCALY = 2.*SCALY
      GO TO 116
  114 CONTINUE
      DJ = 1.0
      GO TO 116
  115 CONTINUE
      DJ = 0.0
  116 CONTINUE
      DO J = 1, MR
         YRT(J) = -4.*C2*SIN((FLOAT(J) - DJ)*DY)**2
      END DO
      SCALY = 2.*SCALY
      GO TO (124,118,122,119,123) MP
  118 CONTINUE
      CALL SINTI (MR, WY, IFACY)
      GO TO 124
  119 CONTINUE
      CALL COSTI (MR, WY, IFACY)
      GO TO 124
  120 CONTINUE
      YRT(1) = 0.
      YRT(MR) = -4.*C2
      DO J = 3, MR, 2
         YRT(J-1) = -4.*C2*SIN(FLOAT(J - 1)*DY)**2
         YRT(J) = YRT(J-1)
      END DO
      CALL RFFTI (MR, WY, IFACY)
      GO TO 124
  122 CONTINUE
      CALL SINQI (MR, WY, IFACY)
      GO TO 124
  123 CONTINUE
      CALL COSQI (MR, WY, IFACY)
  124 CONTINUE
      IFWRD = 1
      IS = 1
  125 CONTINUE
C
C     TRANSFORM X
C
      DO J=1,MR
	 DO K=1,NR
	    DO I=1,LR
               T(I) = F(I,J,K)
	    END DO
            GO TO (127,130,131,134,135),LP
  127       GO TO (128,129),IFWRD
  128       CALL RFFTF (LR,T,WX, IFACX)
            GO TO 138
  129       CALL RFFTB (LR,T,WX, IFACX)
            GO TO 138
  130       CALL SINT (LR,T,WX,IFACX)
            GO TO 138
  131       GO TO (132,133),IFWRD
  132       CALL SINQF (LR,T,WX, IFACX)
            GO TO 138
  133       CALL SINQB (LR,T,WX, IFACX)
            GO TO 138
  134       CALL COST (LR,T,WX, IFACX)
            GO TO 138
  135       GO TO (136,137),IFWRD
  136       CALL COSQF (LR,T,WX, IFACX)
            GO TO 138
  137       CALL COSQB (LR,T,WX, IFACX)
  138       CONTINUE
	    DO I=1,LR
               F(I,J,K) = T(I)
	    END DO
	 END DO
      END DO
      GO TO (142,164) IFWRD
C
C     TRANSFORM Y
C
  142 CONTINUE
      DO I=1,LR
	 DO K=1,NR
	    DO J=1,MR
               T(J) = F(I,J,K)
	    END DO
            GO TO (144,147,148,151,152),MP
  144       GO TO (145,146),IFWRD
  145       CALL RFFTF (MR,T,WY, IFACY)
            GO TO 155
  146       CALL RFFTB (MR,T,WY, IFACY)
            GO TO 155
  147       CALL SINT (MR,T,WY,IFACY)
            GO TO 155
  148       GO TO (149,150),IFWRD
  149       CALL SINQF (MR,T,WY, IFACY)
            GO TO 155
  150       CALL SINQB (MR,T,WY, IFACY)
            GO TO 155
  151       CALL COST (MR,T,WY, IFACY)
            GO TO 155
  152       GO TO (153,154),IFWRD
  153       CALL COSQF (MR,T,WY, IFACY)
            GO TO 155
  154       CALL COSQB (MR,T,WY, IFACY)
  155       CONTINUE
	    DO J=1,MR
               F(I,J,K) = T(J)
	    END DO
	 END DO
      END DO
      GO TO (159,125) IFWRD
  159 CONTINUE
      DO I = 1, LR
         DO J = 1, MR
            BB(:NR) = B(:NR) + XRT(I) + YRT(J)
            T(:NR) = F(I,J,:NR)
            CALL TRID (NR, A, BB, C, T, D)
            F(I,J,:NR) = T(:NR)
         END DO
      END DO
      IFWRD = 2
      IS = -1
      GO TO 142
  164 CONTINUE
      F(:LR,:MR,:NR) = F(:LR,:MR,:NR)/(SCALX*SCALY)
      RETURN 
      END SUBROUTINE POS3D1


      SUBROUTINE TRID(MR, A, B, C, Y, D)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MR
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL , INTENT(INOUT) :: Y(*)
      REAL , INTENT(INOUT) :: D(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: M, MM1, I, IP
      REAL :: Z
C-----------------------------------------------
      M = MR
      MM1 = M - 1
      Z = 1./B(1)
      D(1) = C(1)*Z
      Y(1) = Y(1)*Z
      DO I = 2, MM1
         Z = 1./(B(I)-A(I)*D(I-1))
         D(I) = C(I)*Z
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z
      END DO
      Z = B(M) - A(M)*D(MM1)
      IF (Z == 0.) THEN
         Y(M) = 0.
      ELSE
         Y(M) = (Y(M)-A(M)*Y(MM1))/Z
      ENDIF
      DO IP = 1, MM1
         I = M - IP
         Y(I) = Y(I) - D(I)*Y(I+1)
      END DO
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE TRID
C
C     file poistg.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE POISTG (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR)
C
C
C DIMENSION OF           A(M),  B(M),  C(M),  Y(IDIMY,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
C                        FOR UNKNOWN X VALUES, WHERE I=1,2,...,M
C                        AND J=1,2,...,N
C
C                        A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C                        + X(I,J-1) - 2.*X(I,J) + X(I,J+1)
C                        = Y(I,J)
C
C                        THE INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
C                        I.E. X(0,J) = X(M,J) AND X(M+1,J) = X(1,J), AND
C                        X(I,0) MAY BE EQUAL TO X(I,1) OR -X(I,1), AND
C                        X(I,N+1) MAY BE EQUAL TO X(I,N) OR -X(I,N),
C                        DEPENDING ON AN INPUT PARAMETER.
C
C USAGE                  CALL POISTG (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
C                                     IERROR)
C
C ARGUMENTS
C
C ON INPUT
C
C                        NPEROD
C                          INDICATES VALUES WHICH X(I,0) AND X(I,N+1)
C                          ARE ASSUMED TO HAVE.
C                          = 1 IF X(I,0) = -X(I,1) AND X(I,N+1) = -X(I,N
C                          = 2 IF X(I,0) = -X(I,1) AND X(I,N+1) =  X(I,N
C                          = 3 IF X(I,0) =  X(I,1) AND X(I,N+1) =  X(I,N
C                          = 4 IF X(I,0) =  X(I,1) AND X(I,N+1) = -X(I,N
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        MPEROD
C                          = 0 IF A(1) AND C(M) ARE NOT ZERO
C                          = 1 IF A(1) = C(M) = 0
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                          M MUST BE GREATER THAN 2.
C
C                        A,B,C
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0 THE
C                          ARRAY ELEMENTS MUST NOT DEPEND ON INDEX I,
C                          BUT MUST BE CONSTANT.  SPECIFICALLY, THE
C                          SUBROUTINE CHECKS THE FOLLOWING CONDITION
C                            A(I) = C(1)
C                            B(I) = B(1)
C                            C(I) = C(1)
C                          FOR I = 1, 2, ..., M.
C
C                        IDIMY
C                          THE ROW (OR FIRST) DIMENSION OF THE TWO-
C                          DIMENSIONAL ARRAY Y AS IT APPEARS IN THE
C                          PROGRAM CALLING POISTG.  THIS PARAMETER IS
C                          USED TO SPECIFY THE VARIABLE DIMENSION OF Y.
C                          IDIMY MUST BE AT LEAST M.
C
C                        Y
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
C                          OF EQUATIONS GIVEN ABOVE.
C                          Y MUST BE DIMENSIONED AT LEAST M X N.
C
C ON OUTPUT
C
C                        Y
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
C                          SOLUTION IS NOT ATTEMPTED.
C                          = 0  NO ERROR
C                          = 1  IF M .LE. 2
C                          = 2  IF N .LE. 2
C                          = 3  IDIMY .LT. M
C                          = 4  IF NPEROD .LT. 1 OR NPEROD .GT. 4
C                          = 5  IF MPEROD .LT. 0 OR MPEROD .GT. 1
C                          = 6  IF MPEROD = 0 AND A(I) .NE. C(1)
C                               OR B(I) .NE. B(1) OR C(I) .NE. C(1)
C                               FOR SOME I = 1, 2, ..., M.
C                          = 7  IF MPEROD .EQ. 1 .AND.
C                               (A(1).NE.0 .OR. C(M).NE.0)
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
C                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
C                          SHOULD TEST IERROR AFTER THE CALL.
C
C
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,gnbnaux.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY, 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE IS AN IMPLEMENTATION OF THE
C                        ALGORITHM PRESENTED IN THE REFERENCE BELOW.
C
C TIMING                 FOR LARGE M AND N, THE EXECUTION TIME IS
C                        ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
C                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
C                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
C                        IN THE 'PURPOSE' SECTION ABOVE, WITH
C                          A(I) = C(I) = -0.5*B(I) = 1,    I=1,2,...,M
C                        AND, WHEN MPEROD = 1
C                          A(1) = C(M) = 0
C                          B(1) = B(M) =-1.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
C                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT SID
C                        Y WAS COMPUTED.  USING THIS ARRAY Y SUBROUTINE
C                        POISTG WAS CALLED TO PRODUCE AN APPROXIMATE
C                        SOLUTION Z.  THEN THE RELATIVE ERROR, DEFINED A
C                          E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J)))
C                        WHERE THE TWO MAXIMA ARE TAKEN OVER I=1,2,...,M
C                        AND J=1,2,...,N, WAS COMPUTED.  VALUES OF E ARE
C                        GIVEN IN THE TABLE BELOW FOR SOME TYPICAL VALUE
C                        OF M AND N.
C
C                        M (=N)    MPEROD    NPEROD      E
C                        ------    ------    ------    ------
C
C                          31        0-1       1-4     9.E-13
C                          31        1         1       4.E-13
C                          31        1         3       3.E-13
C                          32        0-1       1-4     3.E-12
C                          32        1         1       3.E-13
C                          32        1         3       1.E-13
C                          33        0-1       1-4     1.E-12
C                          33        1         1       4.E-13
C                          33        1         3       1.E-13
C                          63        0-1       1-4     3.E-12
C                          63        1         1       1.E-12
C                          63        1         3       2.E-13
C                          64        0-1       1-4     4.E-12
C                          64        1         1       1.E-12
C                          64        1         3       6.E-13
C                          65        0-1       1-4     2.E-13
C                          65        1         1       1.E-11
C                          65        1         3       4.E-13
C
C REFERENCES             SCHUMANN, U. AND R. SWEET,"A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON"S EQUATION WITH
C                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
C                        GRID OF ARBITRARY SIZE," J. COMP. PHYS.
C                        20(1976), PP. 171-182.
C *********************************************************************
      SUBROUTINE POISTG(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, IERROR)

      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL  :: A(*)
      REAL  :: B(*)
      REAL  :: C(*)
      REAL  :: Y(IDIMY,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, IRWK, ICWK
C-----------------------------------------------
      IERROR = 0
      IF (M <= 2) IERROR = 1
      IF (N <= 2) IERROR = 2
      IF (IDIMY < M) IERROR = 3
      IF (NPEROD<1 .OR. NPEROD>4) IERROR = 4
      IF (MPEROD<0 .OR. MPEROD>1) IERROR = 5
      IF (MPEROD /= 1) THEN
         DO I = 1, M
            IF (A(I) /= C(1)) GO TO 102
            IF (C(I) /= C(1)) GO TO 102
            IF (B(I) /= B(1)) GO TO 102
         END DO
         GO TO 104
  102    CONTINUE
         IERROR = 6
         RETURN 
      ENDIF
      IF (A(1)/=0. .OR. C(M)/=0.) IERROR = 7
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
!     compute and allocate real work space for poistg
      IRWK = M*(9 + INT(ALOG(FLOAT(N))/ALOG(2.0))) + 4*N
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed (e.g., if n,m are too large)
      IF (IERROR == 20) RETURN 
      call  poistgg(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,w%rew)
!     release work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE POISTG


 
      SUBROUTINE POISTGG(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL  :: Y(IDIMY,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IWBA, IWBB, IWBC, IWB2, IWB3, IWW1, IWW2, IWW3, IWD, 
     1   IWTCOS, IWP, I, K, J, NP, MP, IPSTOR, IREV, MH, MHM1, MODD, 
     2   MHPI, MHMI, NBY2, MSKIP
      REAL :: A1
C-----------------------------------------------
      IWBA = M + 1
      IWBB = IWBA + M
      IWBC = IWBB + M
      IWB2 = IWBC + M
      IWB3 = IWB2 + M
      IWW1 = IWB3 + M
      IWW2 = IWW1 + M
      IWW3 = IWW2 + M
      IWD = IWW3 + M
      IWTCOS = IWD + M
      IWP = IWTCOS + 4*N
      DO I = 1, M
         K = IWBA + I - 1
         W(K) = -A(I)
         K = IWBC + I - 1
         W(K) = -C(I)
         K = IWBB + I - 1
         W(K) = 2. - B(I)
         Y(I,:N) = -Y(I,:N)
      END DO
      NP = NPEROD
      MP = MPEROD + 1
      GO TO (110,107) MP
  107 CONTINUE
      GO TO (108,108,108,119) NPEROD
  108 CONTINUE
      CALL POSTG2 (NP, N, M, W(IWBA), W(IWBB), W(IWBC), IDIMY, Y, W, W(
     1   IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), W
     2   (IWP))
      IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD == 4) GO TO 120
  109 CONTINUE
      GO TO (123,129) MP
  110 CONTINUE
      MH = (M + 1)/2
      MHM1 = MH - 1
      MODD = 1
      IF (MH*2 == M) MODD = 2
      DO J = 1, N
         DO I = 1, MHM1
            W(I) = Y(MH-I,J) - Y(I+MH,J)
            W(I+MH) = Y(MH-I,J) + Y(I+MH,J)
         END DO
         W(MH) = 2.*Y(MH,J)
         GO TO (113,112) MODD
  112    CONTINUE
         W(M) = 2.*Y(M,J)
  113    CONTINUE
         Y(:M,J) = W(:M)
      END DO
      K = IWBC + MHM1 - 1
      I = IWBA + MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      SELECT CASE (MODD) 
      CASE DEFAULT
         K = IWBB + MHM1 - 1
         W(K) = W(K) - W(I-1)
         W(IWBC-1) = W(IWBC-1) + W(IWBB-1)
      CASE (2) 
         W(IWBB-1) = W(K+1)
      END SELECT
  118 CONTINUE
      GO TO 107
  119 CONTINUE
      IREV = 1
      NBY2 = N/2
      NP = 2
  120 CONTINUE
      DO J = 1, NBY2
         MSKIP = N + 1 - J
         DO I = 1, M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
         END DO
      END DO
      GO TO (108,109) IREV
  123 CONTINUE
      DO J = 1, N
         W(MH-1:MH-MHM1:(-1)) = 0.5*(Y(MH+1:MHM1+MH,J)+Y(:MHM1,J))
         W(MH+1:MHM1+MH) = 0.5*(Y(MH+1:MHM1+MH,J)-Y(:MHM1,J))
         W(MH) = 0.5*Y(MH,J)
         GO TO (126,125) MODD
  125    CONTINUE
         W(M) = 0.5*Y(M,J)
  126    CONTINUE
         Y(:M,J) = W(:M)
      END DO
  129 CONTINUE
      W(1) = IPSTOR + IWP - 1
      RETURN 
      END SUBROUTINE POISTGG


      SUBROUTINE POSTG2(NPEROD, N, M, A, BB, C, IDIMQ, Q, B, B2, B3, W, 
     1   W2, W3, D, TCOS, P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: IDIMQ
      REAL  :: A(*)
      REAL  :: BB(*)
      REAL  :: C(*)
      REAL , INTENT(INOUT) :: Q(IDIMQ,*)
      REAL  :: B(*)
      REAL  :: B2(*)
      REAL  :: B3(*)
      REAL  :: W(*)
      REAL  :: W2(*)
      REAL  :: W3(*)
      REAL  :: D(*)
      REAL  :: TCOS(*)
      REAL , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER , DIMENSION(4) :: K
      INTEGER :: K1, K2, K3, K4, NP, MR, IP, IPSTOR, I2R, JR, NR, NLAST
     1   , KR, LR, NROD, JSTART, JSTOP, I2RBY2, J, IJUMP, JP1, JP2, JP3
     2   , JM1, JM2, JM3, I, NRODPR, II, NLASTP, JSTEP
      REAL :: FNUM, FNUM2, FI, T
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION ON A STAGGERED GRID.
C
C
      EQUIVALENCE (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
      NP = NPEROD
      FNUM = 0.5*FLOAT(NP/3)
      FNUM2 = 0.5*FLOAT(NP/2)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      IF (NR > 3) THEN
  101    CONTINUE
         JR = 2*I2R
         NROD = 1
         IF ((NR/2)*2 == NR) NROD = 0
         JSTART = 1
         JSTOP = NLAST - JR
         IF (NROD == 0) JSTOP = JSTOP - I2R
         I2RBY2 = I2R/2
         IF (JSTOP < JSTART) THEN
            J = JR
         ELSE
            IJUMP = 1
            DO J = JSTART, JSTOP, JR
               JP1 = J + I2RBY2
               JP2 = J + I2R
               JP3 = JP2 + I2RBY2
               JM1 = J - I2RBY2
               JM2 = J - I2R
               JM3 = JM2 - I2RBY2
               IF (J == 1) THEN
                  CALL COSGEN (I2R, 1, FNUM, 0.5, TCOS)
                  IF (I2R == 1) THEN
                     B(:MR) = Q(:MR,1)
                     Q(:MR,1) = Q(:MR,2)
                     GO TO 112
                  ENDIF
                  B(:MR) = Q(:MR,1) + 0.5*(Q(:MR,JP2)-Q(:MR,JP1)-Q(:MR,
     1               JP3))
                  Q(:MR,1) = Q(:MR,JP2) + Q(:MR,1) - Q(:MR,JP1)
                  GO TO 112
               ENDIF
               GO TO (107,108) IJUMP
  107          CONTINUE
               IJUMP = 2
               CALL COSGEN (I2R, 1, 0.5, 0.0, TCOS)
  108          CONTINUE
               IF (I2R == 1) THEN
                  B(:MR) = 2.*Q(:MR,J)
                  Q(:MR,J) = Q(:MR,JM2) + Q(:MR,JP2)
               ELSE
                  DO I = 1, MR
                     FI = Q(I,J)
                     Q(I,J)=Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
                     B(I) = FI + Q(I,J) - Q(I,JM3) - Q(I,JP3)
                  END DO
               ENDIF
  112          CONTINUE
               CALL TRIX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,J) + B(:MR)
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
            END DO
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
            J = JSTOP + JR
         ENDIF
         NLAST = J
         JM1 = J - I2RBY2
         JM2 = J - I2R
         JM3 = JM2 - I2RBY2
         IF (NROD /= 0) THEN
C
C     ODD NUMBER OF UNKNOWNS
C
            IF (I2R == 1) THEN
               B(:MR) = Q(:MR,J)
               Q(:MR,J) = Q(:MR,JM2)
            ELSE
               B(:MR)=Q(:MR,J)+0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
               IF (NRODPR == 0) THEN
                  Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP)
                  IP = IP - MR
               ELSE
                  Q(:MR,J) = Q(:MR,J) - Q(:MR,JM1) + Q(:MR,JM2)
               ENDIF
               IF (LR /= 0) CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(KR+1))
            ENDIF
            CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS)
            CALL TRIX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,J) = Q(:MR,J) + B(:MR)
            KR = KR + I2R
         ELSE
            JP1 = J + I2RBY2
            JP2 = J + I2R
            IF (I2R == 1) THEN
               B(:MR) = Q(:MR,J)
               TCOS(1) = 0.
               CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
               IP = 0
               IPSTOR = MR
               P(:MR) = B(:MR)
               B(:MR) = B(:MR) + Q(:MR,N)
               TCOS(1) = (-1.) + 2.*FLOAT(NP/2)
               TCOS(2) = 0.
               CALL TRIX (1, 1, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,JM2) + P(:MR) + B(:MR)
            ELSE
               B(:MR)=Q(:MR,J)+0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
               IF (NRODPR == 0) THEN
                  B(:MR) = B(:MR) + P(IP+1:MR+IP)
               ELSE
                  B(:MR) = B(:MR) + Q(:MR,JP2) - Q(:MR,JP1)
               ENDIF
               CALL COSGEN (I2R, 1, 0.5, 0.0, TCOS)
               CALL TRIX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
               IP = IP + MR
               IPSTOR = MAX0(IPSTOR,IP + MR)
               P(IP+1:MR+IP) = B(:MR) + 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,
     1            JP1))
               B(:MR) = P(IP+1:MR+IP) + Q(:MR,JP2)
               IF (LR /= 0) THEN
                  CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(I2R+1))
                  CALL MERGE (TCOS, 0, I2R, I2R, LR, KR)
               ELSE
                  DO I = 1, I2R
                     II = KR + I
                     TCOS(II) = TCOS(I)
                  END DO
               ENDIF
               CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS)
               CALL TRIX (KR, KR, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP) + B(:MR)
            ENDIF
            LR = KR
            KR = KR + JR
         ENDIF
         NR = (NLAST - 1)/JR + 1
         IF (NR <= 3) GO TO 142
         I2R = JR
         NRODPR = NROD
         GO TO 101
      ENDIF
  142 CONTINUE
      J = 1 + JR
      JM1 = J - I2R
      JP1 = J + I2R
      JM2 = NLAST - I2R
      IF (NR /= 2) THEN
         IF (LR == 0) THEN
            IF (N == 3) THEN
C
C     CASE N = 3.
C
               GO TO (143,148,143) NP
  143          CONTINUE
               B(:MR) = Q(:MR,2)
               B2(:MR) = Q(:MR,1) + Q(:MR,3)
               B3(:MR) = 0.
               SELECT CASE (NP) 
               CASE DEFAULT
                  TCOS(1) = -1.
                  TCOS(2) = 1.
                  K1 = 1
               CASE (1:2) 
                  TCOS(1) = -2.
                  TCOS(2) = 1.
                  TCOS(3) = -1.
                  K1 = 2
               END SELECT
  147          CONTINUE
               K2 = 1
               K3 = 0
               K4 = 0
               GO TO 150
  148          CONTINUE
               B(:MR) = Q(:MR,2)
               B2(:MR) = Q(:MR,3)
               B3(:MR) = Q(:MR,1)
               CALL COSGEN (3, 1, 0.5, 0.0, TCOS)
               TCOS(4) = -1.
               TCOS(5) = 1.
               TCOS(6) = -1.
               TCOS(7) = 1.
               K1 = 3
               K2 = 2
               K3 = 1
               K4 = 1
  150          CONTINUE
               CALL TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
               B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
               GO TO (153,153,152) NP
  152          CONTINUE
               TCOS(1) = 2.
               CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
  153          CONTINUE
               Q(:MR,2) = B(:MR)
               B(:MR) = Q(:MR,1) + B(:MR)
               TCOS(1) = (-1.) + 4.*FNUM
               CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,1) = B(:MR)
               JR = 1
               I2R = 0
               GO TO 188
            ENDIF
C
C     CASE N = 2**P+1
C
            B(:MR)=Q(:MR,J)+Q(:MR,1)-Q(:MR,JM1)+Q(:MR,NLAST)-Q(:MR,JM2)
            GO TO (158,160,158) NP
  158       CONTINUE
            B2(:MR) = Q(:MR,1) + Q(:MR,NLAST) + Q(:MR,J) - Q(:MR,JM1) - 
     1         Q(:MR,JP1)
            B3(:MR) = 0.
            K1 = NLAST - 1
            K2 = NLAST + JR - 1
            CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(NLAST))
            TCOS(K2) = 2.*FLOAT(NP - 2)
            CALL COSGEN (JR, 1, 0.5 - FNUM, 0.5, TCOS(K2+1))
            K3 = (3 - NP)/2
            CALL MERGE (TCOS, K1, JR - K3, K2 - K3, JR + K3, 0)
            K1 = K1 - 1 + K3
            CALL COSGEN (JR, 1, FNUM, 0.5, TCOS(K1+1))
            K2 = JR
            K3 = 0
            K4 = 0
            GO TO 162
  160       CONTINUE
            DO I = 1, MR
               FI = (Q(I,J)-Q(I,JM1)-Q(I,JP1))/2.
               B2(I) = Q(I,1) + FI
               B3(I) = Q(I,NLAST) + FI
            END DO
            K1 = NLAST + JR - 1
            K2 = K1 + JR - 1
            CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K1+1))
            CALL COSGEN (NLAST, 1, 0.5, 0.0, TCOS(K2+1))
            CALL MERGE (TCOS, K1, JR - 1, K2, NLAST, 0)
            K3 = K1 + NLAST - 1
            K4 = K3 + JR
            CALL COSGEN (JR, 1, 0.5, 0.5, TCOS(K3+1))
            CALL COSGEN (JR, 1, 0.0, 0.5, TCOS(K4+1))
            CALL MERGE (TCOS, K3, JR, K4, JR, K1)
            K2 = NLAST - 1
            K3 = JR
            K4 = JR
  162       CONTINUE
            CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
            B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
            IF (NP == 3) THEN
               TCOS(1) = 2.
               CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            ENDIF
            Q(:MR,J) = B(:MR) + 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
            B(:MR) = Q(:MR,J) + Q(:MR,1)
            CALL COSGEN (JR, 1, FNUM, 0.5, TCOS)
            CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,1) = Q(:MR,1) - Q(:MR,JM1) + B(:MR)
            GO TO 188
         ENDIF
C
C     CASE OF GENERAL N WITH NR = 3 .
C
         B(:MR) = Q(:MR,1) - Q(:MR,JM1) + Q(:MR,J)
         IF (NROD == 0) THEN
            B(:MR) = B(:MR) + P(IP+1:MR+IP)
         ELSE
            B(:MR) = B(:MR) + Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
         DO I = 1, MR
            T = 0.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
            Q(I,J) = T
            B2(I) = Q(I,NLAST) + T
            B3(I) = Q(I,1) + T
         END DO
         K1 = KR + 2*JR
         CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K1+1))
         K2 = K1 + JR
         TCOS(K2) = 2.*FLOAT(NP - 2)
         K4 = (NP - 1)*(3 - NP)
         K3 = K2 + 1 - K4
         CALL COSGEN(KR+JR+K4,1,FLOAT(K4)/2.,1.-FLOAT(K4),TCOS(K3))
         K4 = 1 - NP/3
         CALL MERGE (TCOS, K1, JR - K4, K2 - K4, KR + JR + K4, 0)
         IF (NP == 3) K1 = K1 - 1
         K2 = KR + JR
         K4 = K1 + K2
         CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS(K4+1))
         K3 = K4 + KR
         CALL COSGEN (JR, 1, FNUM, 0.5, TCOS(K3+1))
         CALL MERGE (TCOS, K4, KR, K3, JR, K1)
         K4 = K3 + JR
         CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(K4+1))
         CALL MERGE (TCOS, K3, JR, K4, LR, K1 + K2)
         CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS(K3+1))
         K3 = KR
         K4 = KR
         CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
         B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
         IF (NP == 3) THEN
            TCOS(1) = 2.
            CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         ENDIF
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + Q(:MR,J)
         CALL COSGEN (JR, 1, FNUM, 0.5, TCOS)
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         IF (JR == 1) THEN
            Q(:MR,1) = B(:MR)
            GO TO 188
         ENDIF
         Q(:MR,1) = Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 188
      ENDIF
      B3(:MR) = 0.
      B(:MR) = Q(:MR,1) + P(IP+1:MR+IP)
      Q(:MR,1) = Q(:MR,1) - Q(:MR,JM1)
      B2(:MR) = Q(:MR,1) + Q(:MR,NLAST)
      K1 = KR + JR
      K2 = K1 + JR
      CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K1+1))
      GO TO (182,183,182) NP
  182 CONTINUE
      TCOS(K2) = 2.*FLOAT(NP - 2)
      CALL COSGEN (KR, 1, 0.0, 1.0, TCOS(K2+1))
      GO TO 184
  183 CONTINUE
      CALL COSGEN (KR + 1, 1, 0.5, 0.0, TCOS(K2))
  184 CONTINUE
      K4 = 1 - NP/3
      CALL MERGE (TCOS, K1, JR - K4, K2 - K4, KR + K4, 0)
      IF (NP == 3) K1 = K1 - 1
      K2 = KR
      CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS(K1+1))
      K4 = K1 + KR
      CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(K4+1))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
      B(:MR) = B(:MR) + B2(:MR)
      IF (NP == 3) THEN
         TCOS(1) = 2.
         CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
      ENDIF
      Q(:MR,1) = Q(:MR,1) + B(:MR)
  188 CONTINUE
      J = NLAST - JR
      B(:MR) = Q(:MR,NLAST) + Q(:MR,J)
      JM2 = NLAST - I2R
      IF (JR == 1) THEN
         Q(:MR,NLAST) = 0.
      ELSE
         IF (NROD == 0) THEN
            Q(:MR,NLAST) = P(IP+1:MR+IP)
            IP = IP - MR
         ELSE
            Q(:MR,NLAST) = Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
      ENDIF
      CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS)
      CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(KR+1))
      CALL TRIX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
      Q(:MR,NLAST) = Q(:MR,NLAST) + B(:MR)
      NLASTP = NLAST
  197 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR == 0) GO TO 210
      JSTART = 1 + JR
      KR = KR - JR
      IF (NLAST + JR <= N) THEN
         KR = KR - JR
         NLAST = NLAST + JR
         JSTOP = NLAST - JSTEP
      ELSE
         JSTOP = NLAST - JR
      ENDIF
      LR = KR - JR
      CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
      DO J = JSTART, JSTOP, JSTEP
         JM2 = J - JR
         JP2 = J + JR
         IF (J == JR) THEN
            B(:MR) = Q(:MR,J) + Q(:MR,JP2)
         ELSE
            B(:MR) = Q(:MR,J) + Q(:MR,JM2) + Q(:MR,JP2)
         ENDIF
         IF (JR == 1) THEN
            Q(:MR,J) = 0.
         ELSE
            JM1 = J - I2R
            JP1 = J + I2R
            Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         ENDIF
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
      END DO
      NROD = 1
      IF (NLAST + I2R <= N) NROD = 0
      IF (NLASTP /= NLAST) GO TO 188
      GO TO 197
  210 CONTINUE
      W(1) = IPSTOR
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE POSTG2
C
C     file sepaux.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C PACKAGE SEPAUX         CONTAINS NO USER ENTRY POINTS.
C
C LATEST REVISION        June 2004
C
C PURPOSE                THIS PACKAGE CONTAINS AUXILIARY ROUTINES FOR
C                        THE FISHPACK SOLVERS SEPELI AND SEPX4.
C
C USAGE                  SINCE THIS PACKAGE CONTAINS NO USER ENTRIES,
C                        NO USAGE INSTRUCTIONS OR ARGUMENT DESCRIPTIONS
C                        ARE GIVEN HERE.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       NONE
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                DEVELOPED IN THE LATE 1970'S BY JOHN C. ADAMS
C                        OF NCAR'S SCIENTTIFIC COMPUTING DIVISION.
c                        Revised in June 2004 incorporating fortran 90
c                        features
C
C PORTABILITY            FORTRAN 90
C **********************************************************************
      SUBROUTINE SEPORT(USOL, IDMN, ZN, ZM, PERTRB)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDMN
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(INOUT) :: USOL(IDMN,1)
      REAL , INTENT(IN) :: ZN(*)
      REAL , INTENT(IN) :: ZM(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ISTR, IFNL, JSTR, JFNL, I, II, J, JJ
      REAL :: UTE, ETE
C-----------------------------------------------
C
C     THIS SUBROUTINE ORTHOGANALIZES THE ARRAY USOL WITH RESPECT TO
C     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM
C
      ISTR = IS
      IFNL = MS
      JSTR = JS
      JFNL = NS
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO I = IS, MS
         II = I - IS + 1
         ETE = ETE + SUM(ZM(II)*ZN(:NS-JS+1))
         UTE = UTE + SUM(USOL(I,JS:NS)*ZM(II)*ZN(:NS-JS+1))
      END DO
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      USOL(ISTR:IFNL,JSTR:JFNL) = USOL(ISTR:IFNL,JSTR:JFNL) - PERTRB
      RETURN 
      END SUBROUTINE SEPORT


      SUBROUTINE SEPMIN(USOL, IDMN, ZN, ZM, PERTB)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDMN
      REAL  :: PERTB
      REAL , INTENT(INOUT) :: USOL(IDMN,1)
      REAL , INTENT(IN) :: ZN(*)
      REAL , INTENT(IN) :: ZM(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ISTR, IFNL, JSTR, JFNL, I, II, J, JJ
      REAL :: UTE, ETE, PERTRB
C-----------------------------------------------
C
C     THIS SUBROUTINE ORHTOGONALIZES THE ARRAY USOL WITH RESPECT TO
C     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM
C
C
C     ENTRY AT SEPMIN OCCURRS WHEN THE FINAL SOLUTION IS
C     TO BE MINIMIZED WITH RESPECT TO THE WEIGHTED
C     LEAST SQUARES NORM
C
      ISTR = 1
      IFNL = K
      JSTR = 1
      JFNL = L
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO I = IS, MS
         II = I - IS + 1
         ETE = ETE + SUM(ZM(II)*ZN(:NS-JS+1))
         UTE = UTE + SUM(USOL(I,JS:NS)*ZM(II)*ZN(:NS-JS+1))
      END DO
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      USOL(ISTR:IFNL,JSTR:JFNL) = USOL(ISTR:IFNL,JSTR:JFNL) - PERTRB
      RETURN 
      END SUBROUTINE SEPMIN


      SUBROUTINE SEPTRI(N, A, B, C, D, U, Z)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: A(N)
      REAL , INTENT(IN) :: B(N)
      REAL , INTENT(IN) :: C(N)
      REAL , INTENT(INOUT) :: D(N)
      REAL , INTENT(INOUT) :: U(N)
      REAL , INTENT(INOUT) :: Z(N)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NM2, J, NM1, K
      REAL :: BN, V, DEN, AN
C-----------------------------------------------
C
C     THIS SUBROUTINE SOLVES FOR A NON-ZERO EIGENVECTOR CORRESPONDING
C     TO THE ZERO EIGENVALUE OF THE TRANSPOSE OF THE RANK
C     DEFICIENT ONE MATRIX WITH SUBDIAGONAL A, DIAGONAL B, AND
C     SUPERDIAGONAL C , WITH A(1) IN THE (1,N) POSITION, WITH
C     C(N) IN THE (N,1) POSITION, AND ALL OTHER ELEMENTS ZERO.
C
      BN = B(N)
      D(1) = A(2)/B(1)
      V = A(1)
      U(1) = C(N)/B(1)
      NM2 = N - 2
      DO J = 2, NM2
         DEN = B(J) - C(J-1)*D(J-1)
         D(J) = A(J+1)/DEN
         U(J) = -C(J-1)*U(J-1)/DEN
         BN = BN - V*U(J-1)
         V = -V*D(J-1)
      END DO
      DEN = B(N-1) - C(N-2)*D(N-2)
      D(N-1) = (A(N)-C(N-2)*U(N-2))/DEN
      AN = C(N-1) - V*D(N-2)
      BN = BN - V*U(N-2)
      DEN = BN - AN*D(N-1)
C
C     SET LAST COMPONENT EQUAL TO ONE
C
      Z(N) = 1.0
      Z(N-1) = -D(N-1)
      NM1 = N - 1
      DO J = 2, NM1
         K = N - J
         Z(K) = (-D(K)*Z(K+1)) - U(K)*Z(N)
      END DO
      RETURN 
      END SUBROUTINE SEPTRI


      SUBROUTINE SEPDX(U, IDMN, I, J, UXXX, UXXXX)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDMN
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: J
      REAL , INTENT(OUT) :: UXXX
      REAL , INTENT(OUT) :: UXXXX
      REAL , INTENT(IN) :: U(IDMN,1)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE
C     APPROXIMATIONS TO THE THIRD AND FOURTH X
C     PARTIAL DERIVATIVES OF U AT THE (I,J) MESH POINT
C
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A
C
      IF (I == 1) THEN
         IF (KSWX /= 1) THEN
            UXXX = ((-5.0*U(1,J))+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)-
     1         3.0*U(5,J))/TDLX3
            UXXXX = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+11.0
     1         *U(5,J)-2.0*U(6,J))/DLX4
            RETURN 
         ELSE
C
C     PERIODIC AT X=A
C
            UXXX = ((-U(K-2,J))+2.0*U(K-1,J)-2.0*U(2,J)+U(3,J))/TDLX3
            UXXXX = (U(K-2,J)-4.0*U(K-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J))
     1         /DLX4
            RETURN 
         ENDIF
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX
C
      ELSE IF (I == 2) THEN
         IF (KSWX /= 1) THEN
            UXXX = ((-3.0*U(1,J))+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-U(5
     1         ,J))/TDLX3
            UXXXX = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*U
     1         (5,J)-U(6,J))/DLX4
            RETURN 
         ELSE
C
C     PERIODIC AT X=A+DLX
C
            UXXX = ((-U(K-1,J))+2.0*U(1,J)-2.0*U(3,J)+U(4,J))/TDLX3
            UXXXX = (U(K-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J))/
     1         DLX4
            RETURN 
         ENDIF
      ELSE IF (I>2 .AND. I<K-1) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
         UXXX = ((-U(I-2,J))+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J))/TDLX3
         UXXXX = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J)
     1      )/DLX4
         RETURN 
      ELSE IF (I == K - 1) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX
C
         IF (KSWX /= 1) THEN
            UXXX = (U(K-4,J)-6.0*U(K-3,J)+12.0*U(K-2,J)-10.0*U(K-1,J)+
     1         3.0*U(K,J))/TDLX3
            UXXXX = ((-U(K-5,J))+6.0*U(K-4,J)-14.0*U(K-3,J)+16.0*U(K-2,J
     1         )-9.0*U(K-1,J)+2.0*U(K,J))/DLX4
            RETURN 
         ELSE
C
C     PERIODIC AT X=B-DLX
C
            UXXX = ((-U(K-3,J))+2.0*U(K-2,J)-2.0*U(1,J)+U(2,J))/TDLX3
            UXXXX = (U(K-3,J)-4.0*U(K-2,J)+6.0*U(K-1,J)-4.0*U(1,J)+U(2,J
     1         ))/DLX4
            RETURN 
         ENDIF
      ELSE IF (I == K) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B
C
         UXXX = -(3.0*U(K-4,J)-14.0*U(K-3,J)+24.0*U(K-2,J)-18.0*U(K-1,J)
     1      +5.0*U(K,J))/TDLX3
         UXXXX = ((-2.0*U(K-5,J))+11.0*U(K-4,J)-24.0*U(K-3,J)+26.0*U(K-2
     1      ,J)-14.0*U(K-1,J)+3.0*U(K,J))/DLX4
         RETURN 
      ENDIF
      RETURN 
      END SUBROUTINE SEPDX


      SUBROUTINE SEPDY(U, IDMN, I, J, UYYY, UYYYY)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDMN
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: J
      REAL , INTENT(OUT) :: UYYY
      REAL , INTENT(OUT) :: UYYYY
      REAL , INTENT(IN) :: U(IDMN,6)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE
C     APPROXIMATIONS TO THE THIRD AND FOURTH Y
C     PARTIAL DERIVATIVES OF U AT THE (I,J) MESH POINT
C
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
C
      IF (J == 1) THEN
         IF (KSWY /= 1) THEN
            UYYY = ((-5.0*U(I,1))+18.0*U(I,2)-24.0*U(I,3)+14.0*U(I,4)-
     1         3.0*U(I,5))/TDLY3
            UYYYY = (3.0*U(I,1)-14.0*U(I,2)+26.0*U(I,3)-24.0*U(I,4)+11.0
     1         *U(I,5)-2.0*U(I,6))/DLY4
            RETURN 
         ELSE
C
C     PERIODIC AT X=A
C
            UYYY = ((-U(I,L-2))+2.0*U(I,L-1)-2.0*U(I,2)+U(I,3))/TDLY3
            UYYYY = (U(I,L-2)-4.0*U(I,L-1)+6.0*U(I,1)-4.0*U(I,2)+U(I,3))
     1         /DLY4
            RETURN 
         ENDIF
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY
C
      ELSE IF (J == 2) THEN
         IF (KSWY /= 1) THEN
            UYYY = ((-3.0*U(I,1))+10.0*U(I,2)-12.0*U(I,3)+6.0*U(I,4)-U(I
     1         ,5))/TDLY3
            UYYYY = (2.0*U(I,1)-9.0*U(I,2)+16.0*U(I,3)-14.0*U(I,4)+6.0*U
     1         (I,5)-U(I,6))/DLY4
            RETURN 
         ELSE
C
C     PERIODIC AT Y=C+DLY
C
            UYYY = ((-U(I,L-1))+2.0*U(I,1)-2.0*U(I,3)+U(I,4))/TDLY3
            UYYYY = (U(I,L-1)-4.0*U(I,1)+6.0*U(I,2)-4.0*U(I,3)+U(I,4))/
     1         DLY4
            RETURN 
         ENDIF
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
      ELSE IF (J>2 .AND. J<L-1) THEN
         UYYY = ((-U(I,J-2))+2.0*U(I,J-1)-2.0*U(I,J+1)+U(I,J+2))/TDLY3
         UYYYY = (U(I,J-2)-4.0*U(I,J-1)+6.0*U(I,J)-4.0*U(I,J+1)+U(I,J+2)
     1      )/DLY4
         RETURN 
      ELSE IF (J == L - 1) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY
C
         IF (KSWY /= 1) THEN
            UYYY = (U(I,L-4)-6.0*U(I,L-3)+12.0*U(I,L-2)-10.0*U(I,L-1)+
     1         3.0*U(I,L))/TDLY3
            UYYYY = ((-U(I,L-5))+6.0*U(I,L-4)-14.0*U(I,L-3)+16.0*U(I,L-2
     1         )-9.0*U(I,L-1)+2.0*U(I,L))/DLY4
            RETURN 
         ELSE
C
C     PERIODIC AT Y=D-DLY
C
            UYYY = ((-U(I,L-3))+2.0*U(I,L-2)-2.0*U(I,1)+U(I,2))/TDLY3
            UYYYY = (U(I,L-3)-4.0*U(I,L-2)+6.0*U(I,L-1)-4.0*U(I,1)+U(I,2
     1         ))/DLY4
            RETURN 
         ENDIF
      ELSE IF (J == L) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D
C
         UYYY = -(3.0*U(I,L-4)-14.0*U(I,L-3)+24.0*U(I,L-2)-18.0*U(I,L-1)
     1      +5.0*U(I,L))/TDLY3
         UYYYY = ((-2.0*U(I,L-5))+11.0*U(I,L-4)-24.0*U(I,L-3)+26.0*U(I,L
     1      -2)-14.0*U(I,L-1)+3.0*U(I,L))/DLY4
         RETURN 
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    version 5.0, fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE SEPDY
C
C     file sepeli.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE SEPELI (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,
C    +                   D,N,NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,GRHS,
C    +                   USOL,IDMN,W,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C ARGUMENTS              USOL(IDMN,N+1),GRHS(IDMN,N+1),
C
C LATEST REVISION        JUNE 2004
C
C PURPOSE                SEPELI SOLVES FOR EITHER THE SECOND-ORDER
C                        FINITE DIFFERENCE APPROXIMATION OR A
C                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
C                        ELLIPTIC EQUATION
C
C                                 2    2
C                          AF(X)*D U/DX + BF(X)*DU/DX  + CF(X)*U +
C                                 2    2
C                          DF(Y)*D U/DY  + EF(Y)*DU/DY + FF(Y)*U
C
C                          = G(X,Y)
C
C                        ON A RECTANGLE (X GREATER THAN OR EQUAL TO A
C                        AND LESS THAN OR EQUAL TO B; Y GREATER THAN
C                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
C                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
C                        CONDITIONS IS ALLOWED.
C
C                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
C                        IN THE X-DIRECTION:
C                        (0) PERIODIC, U(X+B-A,Y)=U(X,Y) FOR ALL
C                            Y,X (1) U(A,Y), U(B,Y) ARE SPECIFIED FOR
C                            ALL Y
C                        (2) U(A,Y), DU(B,Y)/DX+BETA*U(B,Y) ARE
C                            SPECIFIED FOR ALL Y
C                        (3) DU(A,Y)/DX+ALPHA*U(A,Y),DU(B,Y)/DX+
C                            BETA*U(B,Y) ARE SPECIFIED FOR ALL Y
C                        (4) DU(A,Y)/DX+ALPHA*U(A,Y),U(B,Y) ARE
C                            SPECIFIED FOR ALL Y
C
C                        IN THE Y-DIRECTION:
C                        (0) PERIODIC, U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                        (1) U(X,C),U(X,D) ARE SPECIFIED FOR ALL X
C                        (2) U(X,C),DU(X,D)/DY+XNU*U(X,D) ARE
C                            SPECIFIED FOR ALL X
C                        (3) DU(X,C)/DY+GAMA*U(X,C),DU(X,D)/DY+
C                            XNU*U(X,D) ARE SPECIFIED FOR ALL X
C                        (4) DU(X,C)/DY+GAMA*U(X,C),U(X,D) ARE
C                            SPECIFIED FOR ALL X
C
C USAGE                  CALL SEPELI (INTL,IORDER,A,B,M,MBDCND,BDA,
C                                     ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,
C                                     GAMA,BDD,XNU,COFX,COFY,GRHS,USOL,
C                                     IDMN,W,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               INTL
C                          = 0 ON INITIAL ENTRY TO SEPELI OR IF ANY
C                              OF THE ARGUMENTS C,D, N, NBDCND, COFY
C                              ARE CHANGED FROM A PREVIOUS CALL
C                          = 1 IF C, D, N, NBDCND, COFY ARE UNCHANGED
C                              FROM THE PREVIOUS CALL.
C
C                        IORDER
C                          = 2 IF A SECOND-ORDER APPROXIMATION
C                              IS SOUGHT
C                          = 4 IF A FOURTH-ORDER APPROXIMATION
C                              IS SOUGHT
C
C                        A,B
C                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
C                          I.E., X IS GREATER THAN OR EQUAL TO A
C                          AND LESS THAN OR EQUAL TO B.  A MUST BE
C                          LESS THAN B.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL [A,B] IS SUBDIVIDED. HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE X-
C                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
C                          FOR I=1,2,...,M+1 WHERE DLX=(B-A)/M IS
C                          THE PANEL WIDTH.  M MUST BE LESS THAN
C                          IDMN AND GREATER THAN 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT X=A AND X=B
C
C                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
C                              U(X+B-A,Y)=U(X,Y)  FOR ALL Y,X
C                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND X=B, I.E., U(A,Y) AND U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A AND
C                              THE BOUNDARY CONDITION IS MIXED AT X=B,
C                              I.E., U(A,Y) AND DU(B,Y)/DX+BETA*U(B,Y)
C                              ARE SPECIFIED FOR ALL Y
C                          = 3 IF THE BOUNDARY CONDITIONS AT X=A AND
C                              X=B ARE MIXED, I.E.,
C                              DU(A,Y)/DX+ALPHA*U(A,Y) AND
C                              DU(B,Y)/DX+BETA*U(B,Y) ARE SPECIFIED
C                              FOR ALL Y
C                          = 4 IF THE BOUNDARY CONDITION AT X=A IS
C                              MIXED AND THE SOLUTION IS SPECIFIED
C                              AT X=B, I.E., DU(A,Y)/DX+ALPHA*U(A,Y)
C                              AND U(B,Y) ARE SPECIFIED FOR ALL Y
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF
C                          DU(A,Y)/DX+ ALPHA*U(A,Y) AT X=A, WHEN
C                          MBDCND=3 OR 4.
C                          BDA(J) = DU(A,YJ)/DX+ALPHA*U(A,YJ),
C                          J=1,2,...,N+1. WHEN MBDCND HAS ANY OTHER
C                          OTHER VALUE, BDA IS A DUMMY PARAMETER.
C
C                        ALPHA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT X=A
C                          (SEE ARGUMENT BDA).  IF MBDCND IS NOT
C                          EQUAL TO 3 OR 4 THEN ALPHA IS A DUMMY
C                          PARAMETER.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF
C                          DU(B,Y)/DX+ BETA*U(B,Y) AT X=B.
C                          WHEN MBDCND=2 OR 3
C                          BDB(J) = DU(B,YJ)/DX+BETA*U(B,YJ),
C                          J=1,2,...,N+1. WHEN MBDCND HAS ANY OTHER
C                          OTHER VALUE, BDB IS A DUMMY PARAMETER.
C
C                        BETA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          X=B (SEE ARGUMENT BDB).  IF MBDCND IS
C                          NOT EQUAL TO 2 OR 3 THEN BETA IS A DUMMY
C                          PARAMETER.
C
C                        C,D
C                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
C                          I.E., Y IS GREATER THAN OR EQUAL TO C
C                          AND LESS THAN OR EQUAL TO D.  C MUST BE
C                          LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL [C,D] IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE Y-DIRECTION GIVEN BY
C                          YJ=C+(J-1)*DLY FOR J=1,2,...,N+1 WHERE
C                          DLY=(D-C)/N IS THE PANEL WIDTH.
C                          IN ADDITION, N MUST BE GREATER THAN 4.
C
C                        NBDCND
C                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
C                          AT Y=C AND Y=D
C
C                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
C                              I.E., U(X,Y+D-C)=U(X,Y)  FOR ALL X,Y
C                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND Y = D, I.E., U(X,C) AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND THE BOUNDARY CONDITION IS MIXED
C                              AT Y=D, I.E., U(X,C) AND
C                              DU(X,D)/DY+XNU*U(X,D) ARE SPECIFIED
C                              FOR ALL X
C                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
C                              AT Y=C AND Y=D, I.E.,
C                              DU(X,D)/DY+GAMA*U(X,C) AND
C                              DU(X,D)/DY+XNU*U(X,D) ARE SPECIFIED
C                              FOR ALL X
C                          = 4 IF THE BOUNDARY CONDITION IS MIXED
C                              AT Y=C AND THE SOLUTION IS SPECIFIED
C                              AT Y=D, I.E. DU(X,C)/DY+GAMA*U(X,C)
C                              AND U(X,D) ARE SPECIFIED FOR ALL X
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUE OF
C                          DU(X,C)/DY+GAMA*U(X,C) AT Y=C.
C                          WHEN NBDCND=3 OR 4 BDC(I) = DU(XI,C)/DY +
C                          GAMA*U(XI,C), I=1,2,...,M+1.
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC
C                          IS A DUMMY PARAMETER.
C
C                        GAMA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          Y=C (SEE ARGUMENT BDC).  IF NBDCND IS
C                          NOT EQUAL TO 3 OR 4 THEN GAMA IS A DUMMY
C                          PARAMETER.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUE OF
C                          DU(X,D)/DY + XNU*U(X,D) AT Y=C.
C                          WHEN NBDCND=2 OR 3 BDD(I) = DU(XI,D)/DY +
C                          XNU*U(XI,D), I=1,2,...,M+1.
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD
C                          IS A DUMMY PARAMETER.
C
C                        XNU
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          Y=D (SEE ARGUMENT BDD).  IF NBDCND IS
C                          NOT EQUAL TO 2 OR 3 THEN XNU IS A
C                          DUMMY PARAMETER.
C
C                        COFX
C                          A USER-SUPPLIED SUBPROGRAM WITH
C                          PARAMETERS X, AFUN, BFUN, CFUN WHICH
C                          RETURNS THE VALUES OF THE X-DEPENDENT
C                          COEFFICIENTS AF(X), BF(X), CF(X) IN THE
C                          ELLIPTIC EQUATION AT X.
C
C                        COFY
C                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
C                          Y, DFUN, EFUN, FFUN WHICH RETURNS THE
C                          VALUES OF THE Y-DEPENDENT COEFFICIENTS
C                          DF(Y), EF(Y), FF(Y) IN THE ELLIPTIC
C                          EQUATION AT Y.
C
C                          NOTE:  COFX AND COFY MUST BE DECLARED
C                          EXTERNAL IN THE CALLING ROUTINE.
C                          THE VALUES RETURNED IN AFUN AND DFUN
C                          MUST SATISFY AFUN*DFUN GREATER THAN 0
C                          FOR A LESS THAN X LESS THAN B, C LESS
C                          THAN Y LESS THAN D (SEE IERROR=10).
C                          THE COEFFICIENTS PROVIDED MAY LEAD TO A
C                          MATRIX EQUATION WHICH IS NOT DIAGONALLY
C                          DOMINANT IN WHICH CASE SOLUTION MAY FAIL
C                          (SEE IERROR=4).
C
C                        GRHS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT-HAND SIDE OF THE
C                          ELLIPTIC EQUATION, I.E.,
C                          GRHS(I,J)=G(XI,YI), FOR I=2,...,M,
C                          J=2,...,N.  AT THE BOUNDARIES, GRHS IS
C                          DEFINED BY
C
C                          MBDCND   GRHS(1,J)   GRHS(M+1,J)
C                          ------   ---------   -----------
C                            0      G(A,YJ)     G(B,YJ)
C                            1         *           *
C                            2         *        G(B,YJ)  J=1,2,...,N+1
C                            3      G(A,YJ)     G(B,YJ)
C                            4      G(A,YJ)        *
C
C                          NBDCND   GRHS(I,1)   GRHS(I,N+1)
C                          ------   ---------   -----------
C                            0      G(XI,C)     G(XI,D)
C                            1         *           *
C                            2         *        G(XI,D)  I=1,2,...,M+1
C                            3      G(XI,C)     G(XI,D)
C                            4      G(XI,C)        *
C
C                          WHERE * MEANS THESE QUANTITIES ARE NOT USED.
C                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
C                        USOL
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE SOLUTION ALONG THE BOUNDARIES.
C                          AT THE BOUNDARIES, USOL IS DEFINED BY
C
C                          MBDCND   USOL(1,J)   USOL(M+1,J)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(A,YJ)     U(B,YJ)
C                            2      U(A,YJ)        *     J=1,2,...,N+1
C                            3         *           *
C                            4         *        U(B,YJ)
C
C                          NBDCND   USOL(I,1)   USOL(I,N+1)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(XI,C)     U(XI,D)
C                            2      U(XI,C)        *     I=1,2,...,M+1
C                            3         *           *
C                            4         *        U(XI,D)
C
C                          WHERE * MEANS THE QUANTITIES ARE NOT USED
C                          IN THE SOLUTION.
C
C                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
C                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
C                          CASE THE TABLES SPECIFYING THE BOUNDARIES
C                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
C                          BOUNDARIES UNIQUELY EXCEPT AT THE CORNERS.
C                          IF THE TABLES CALL FOR BOTH G(X,Y) AND
C                          U(X,Y) AT A CORNER THEN THE SOLUTION MUST
C                          BE CHOSEN.  FOR EXAMPLE, IF MBDCND=2 AND
C                          NBDCND=4, THEN U(A,C), U(A,D), U(B,D) MUST
C                          BE CHOSEN AT THE CORNERS IN ADDITION
C                          TO G(B,C).
C
C                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
C                          GRHS, MUST BE DISTINCT.
C
C                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
C                        IDMN
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAYS
C                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
C                          CALLING SEPELI.  THIS PARAMETER IS USED
C                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
C                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
C                          GREATER THAN OR EQUAL TO M+1.
C
C                        W
c                          A fortran 90 derived TYPE (fishworkspace) variable
c                          that must be declared by the user.  The first
c                          two declarative statements in the user program
c                          calling SEPELI must be:
c
c                               USE fish
c                               TYPE (fishworkspace) :: W
c
c                          The first statement makes the fishpack module
c                          defined in the file "fish.f" available to the
c                          user program calling SEPELI.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in SEPELI to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (INTL=0) call to SEPELI.  These
c                          must be preserved on non-initial calls (INTL=1)
c                          to SEPELI.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling SEPELI should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          SEPELI.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
c
C ON OUTPUT              USOL
C                          CONTAINS THE APPROXIMATE SOLUTION TO THE
C                          ELLIPTIC EQUATION.
C                          USOL(I,J) IS THE APPROXIMATION TO U(XI,YJ)
C                          FOR I=1,2...,M+1 AND J=1,2,...,N+1.
C                          THE APPROXIMATION HAS ERROR
C                          O(DLX**2+DLY**2) IF CALLED WITH IORDER=2
C                          AND O(DLX**4+DLY**4) IF CALLED WITH
C                          IORDER=4.
C
C                        W
c                          The derived type (fishworkspace) variable W
c                          contains real and complex values that must not
C                          be destroyed if SEPELI is called again with
C                          INTL=1.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS
C                          (I.E., ALPHA=BETA=0 IF MBDCND=3;
C                          GAMA=XNU=0 IF NBDCND=3) IS SPECIFIED
C                          AND IF THE COEFFICIENTS OF U(X,Y) IN THE
C                          SEPARABLE ELLIPTIC EQUATION ARE ZERO
C                          (I.E., CF(X)=0 FOR X GREATER THAN OR EQUAL
C                          TO A AND LESS THAN OR EQUAL TO B;
C                          FF(Y)=0 FOR Y GREATER THAN OR EQUAL TO C
C                          AND LESS THAN OR EQUAL TO D) THEN A
C                          SOLUTION MAY NOT EXIST.  PERTRB IS A
C                          CONSTANT CALCULATED AND SUBTRACTED FROM
C                          THE RIGHT-HAND SIDE OF THE MATRIX EQUATIONS
C                          GENERATED BY SEPELI WHICH INSURES THAT A
C                          SOLUTION EXISTS. SEPELI THEN COMPUTES THIS
C                          SOLUTION WHICH IS A WEIGHTED MINIMAL LEAST
C                          SQUARES SOLUTION TO THE ORIGINAL PROBLEM.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS OR FAILURE TO FIND A SOLUTION
C                          = 0 NO ERROR
C                          = 1 IF A GREATER THAN B OR C GREATER THAN D
C                          = 2 IF MBDCND LESS THAN 0 OR MBDCND GREATER
C                              THAN 4
C                          = 3 IF NBDCND LESS THAN 0 OR NBDCND GREATER
C                              THAN 4
C                          = 4 IF ATTEMPT TO FIND A SOLUTION FAILS.
C                              (THE LINEAR SYSTEM GENERATED IS NOT
C                              DIAGONALLY DOMINANT.)
C                          = 5 IF IDMN IS TOO SMALL
C                              (SEE DISCUSSION OF IDMN)
C                          = 6 IF M IS TOO SMALL OR TOO LARGE
C                              (SEE DISCUSSION OF M)
C                          = 7 IF N IS TOO SMALL (SEE DISCUSSION OF N)
C                          = 8 IF IORDER IS NOT 2 OR 4
C                          = 9 IF INTL IS NOT 0 OR 1
C                          = 10 IF AFUN*DFUN LESS THAN OR EQUAL TO 0
C                               FOR SOME INTERIOR MESH POINT (XI,YJ)
C                          = 20 If the dynamic allocation of real and
C                               complex work space in the derived type
C                               (fishworkspace) variable W fails (e.g.,
c                               if N,M are too large for the platform used)
c
C                          NOTE (CONCERNING IERROR=4):  FOR THE
C                          COEFFICIENTS INPUT THROUGH COFX, COFY,
C                          THE DISCRETIZATION MAY LEAD TO A BLOCK
C                          TRIDIAGONAL LINEAR SYSTEM WHICH IS NOT
C                          DIAGONALLY DOMINANT (FOR EXAMPLE, THIS
C                          HAPPENS IF CFUN=0 AND BFUN/(2.*DLX) GREATER
C                          THAN AFUN/DLX**2).  IN THIS CASE SOLUTION
C                          MAY FAIL.  THIS CANNOT HAPPEN IN THE LIMIT
C                          AS DLX, DLY APPROACH ZERO.  HENCE, THE
C                          CONDITION MAY BE REMEDIED BY TAKING LARGER
C                          VALUES FOR M OR N.
C
C SPECIAL CONDITIONS     SEE COFX, COFY ARGUMENT DESCRIPTIONS ABOVE.
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED FILES         blktri.f,comf.f,sepaux.f,fish.f
C
C LANGUAGE               Fortran 90
C
C HISTORY                DEVELOPED AT NCAR DURING 1975-76 BY
C                        JOHN C. ADAMS OF THE SCIENTIFIC COMPUTING
C                        DIVISION.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980. Revised in June
C                        2004 using Fortan 90 dynamically allocated work
c                        space and derived data types to eliminate mixed
c                        mode conflicts in the earlier versions. All
c                        statement labels, arithmetic if statements and
c                        computed GO TO statements have been removed from
c                        the current version of SEPELI.
C
C ALGORITHM              SEPELI AUTOMATICALLY DISCRETIZES THE
C                        SEPARABLE ELLIPTIC EQUATION WHICH IS THEN
C                        SOLVED BY A GENERALIZED CYCLIC REDUCTION
C                        ALGORITHM IN THE SUBROUTINE, BLKTRI.  THE
C                        FOURTH-ORDER SOLUTION IS OBTAINED USING
C                        'DEFERRED CORRECTIONS' WHICH IS DESCRIBED
C                        AND REFERENCED IN SECTIONS, REFERENCES AND
C                        METHOD.
C
C TIMING                 THE OPERATIONAL COUNT IS PROPORTIONAL TO
C                        M*N*LOG2(N).
C
C ACCURACY               THE FOLLOWING ACCURACY RESULTS WERE OBTAINED
C                        using 64 bit floating point arithmetic.  Note
C                        THAT THE FOURTH-ORDER accuracy is not realized
C                        UNTIL THE MESH IS sufficiently refined.
C
C                                     SECOND-ORDER  FOURTH-ORDER
C                            M    N     ERROR         ERROR
C
C                             6    6    6.8E-1        1.2E0
C                            14   14    1.4E-1        1.8E-1
C                            30   30    3.2E-2        9.7E-3
C                            62   62    7.5E-3        3.0E-4
C                           126  126    1.8E-3        3.5E-6
C
C
C REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
C                        BOUNDARY-VALUE PROBLEMS, BLAISDEL (1968),
C                        WALTHAM, MASS.
C
C                        SWARZTRAUBER, P., AND R. SWEET (1975):
C                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
C                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
C                        EQUATIONS.  NCAR TECHNICAL NOTE
C                        NCAR-TN/IA-109, PP. 135-137.
C***********************************************************************
 
      SUBROUTINE SEPELI(INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, 
     1   BETA, C, D, N, NBDCND, BDC, GAMA, BDD, XNU, COFX, COFY, GRHS, 
     2   USOL, IDMN, W, PERTRB, IERROR)

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: INTL
      INTEGER  :: IORDER
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDMN
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: ALPHA
      REAL  :: BETA
      REAL  :: C
      REAL  :: D
      REAL  :: GAMA
      REAL  :: XNU
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: GRHS(IDMN,*)
      REAL  :: USOL(IDMN,*)

      TYPE (fishworkspace) :: W
      EXTERNAL :: COFX,COFY
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,K,L,NP,IRWK,ICWK

      SAVE I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12
C-----------------------------------------------
!     save local variable work space pointers for noninitial call
!     check input arguments
      CALL CHKPRM (INTL, IORDER, A, B, M, MBDCND, C, D, N, NBDCND,
     1   COFX, COFY, IDMN, IERROR)
      IF (IERROR /= 0) RETURN 
      IF (INTL == 0) THEN
!     allocate space and set work space indices on initial call only
         K = M + 1
         L = N + 1
!          compute required blktri work space lengths
         NP = NBDCND
         CALL BLK_SPACE (N, M, IRWK, ICWK)
C
C     SET WORK SPACE INDICES
C
         I1 = IRWK + 1
         I2 = I1 + L
         I3 = I2 + L
         I4 = I3 + L
         I5 = I4 + L
         I6 = I5 + L
         I7 = I6 + L
         I8 = I7 + K
         I9 = I8 + K
         I10 = I9 + K
         I11 = I10 + K
         I12 = I11 + K
!          set sepeli work space requirements
         IRWK = I12 + K
         ICWK = ICWK + 3*(M + 1)
!          allocate required real and complex work space
         CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!          return if allocation failure
         IF (IERROR == 20) RETURN 
      ENDIF
      IERROR = 0
!     compute second or fourth order solution
      CALL SPELIP (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     +NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,w%rew(i1),w%rew(i2),w%rew(i3),
     +w%rew(i4),w%rew(i5),w%rew(i6),w%rew(i7),w%rew(i8),w%rew(i9),
     +w%rew(i10),w%rew(i11),w%rew(i12),grhs,usol,idmn,w%rew,w%cxw,
     +pertrb,ierror)
      RETURN 
      END SUBROUTINE SEPELI


      SUBROUTINE SPELIP(INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, 
     1   BETA, C, D, N, NBDCND, BDC, GAMA, BDD, XNU, COFX, COFY, AN, BN
     2   , CN, DN, UN, ZN, AM, BM, CM, DM, UM, ZM, GRHS, USOL, IDMN, W, 
     3   WC, PERTRB, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: INTL
      INTEGER , INTENT(IN) :: IORDER
      INTEGER , INTENT(IN) :: M
      INTEGER  :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDMN
      INTEGER  :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL  :: ALPHA
      REAL  :: BETA
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL  :: GAMA
      REAL  :: XNU
      REAL  :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL  :: DN(*)
      REAL  :: UN(*)
      REAL  :: ZN(*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL  :: DM(*)
      REAL  :: UM(*)
      REAL  :: ZM(*)
      REAL  :: GRHS(IDMN,*)
      REAL  :: USOL(IDMN,*)
      REAL  :: W(*)
      COMPLEX  :: WC(*)
      EXTERNAL :: COFX, COFY
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J, I1, MP, NP
      REAL :: XI, AI, BI, CI, AXI, BXI, CXI, YJ, DJ, EJ, FJ, DYJ, EYJ, 
     1   FYJ, AX1, CXM, DY1, FYN, PRTRB
      LOGICAL :: SINGLR
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C     SPELIP SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
C     AND COMPUTES A SECOND ORDER SOLUTION IN USOL.  A RETURN JUMP TO
C     SEPELI OCCURRS IF IORDER=2.  IF IORDER=4 A FOURTH ORDER
C     SOLUTION IS GENERATED IN USOL.
C
C
C     SET PARAMETERS INTERNALLY
C
      KSWX = MBDCND + 1
      KSWY = NBDCND + 1
      K = M + 1
      L = N + 1
      AIT = A
      BIT = B
      CIT = C
      DIT = D
C
C     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
C     AND NON-SPECIFIED BOUNDARIES.
C
      USOL(2:M,2:N) = GRHS(2:M,2:N)
      IF (KSWX/=2 .AND. KSWX/=3) THEN
         USOL(1,2:N) = GRHS(1,2:N)
      ENDIF
      IF (KSWX/=2 .AND. KSWX/=5) THEN
         USOL(K,2:N) = GRHS(K,2:N)
      ENDIF
      IF (KSWY/=2 .AND. KSWY/=3) THEN
         USOL(2:M,1) = GRHS(2:M,1)
      ENDIF
      IF (KSWY/=2 .AND. KSWY/=5) THEN
         USOL(2:M,L) = GRHS(2:M,L)
      ENDIF
      IF (KSWX/=2 .AND. KSWX/=3 .AND. KSWY/=2 .AND. KSWY/=3) USOL(1,1)
     1    = GRHS(1,1)
      IF (KSWX/=2 .AND. KSWX/=5 .AND. KSWY/=2 .AND. KSWY/=3) USOL(K,1)
     1    = GRHS(K,1)
      IF (KSWX/=2 .AND. KSWX/=3 .AND. KSWY/=2 .AND. KSWY/=5) USOL(1,L)
     1    = GRHS(1,L)
      IF (KSWX/=2 .AND. KSWX/=5 .AND. KSWY/=2 .AND. KSWY/=5) USOL(K,L)
     1    = GRHS(K,L)
      I1 = 1
C
C     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
C
      MP = 1
      NP = 1
      IF (KSWX == 1) MP = 0
      IF (KSWY == 1) NP = 0
C
C     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
C     IN NINT,MINT
C
      DLX = (BIT - AIT)/FLOAT(M)
      MIT = K - 1
      IF (KSWX == 2) MIT = K - 2
      IF (KSWX == 4) MIT = K
      DLY = (DIT - CIT)/FLOAT(N)
      NIT = L - 1
      IF (KSWY == 2) NIT = L - 2
      IF (KSWY == 4) NIT = L
      TDLX3 = 2.0*DLX**3
      DLX4 = DLX**4
      TDLY3 = 2.0*DLY**3
      DLY4 = DLY**4
C
C     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
C
      IS = 1
      JS = 1
      IF (KSWX==2 .OR. KSWX==3) IS = 2
      IF (KSWY==2 .OR. KSWY==3) JS = 2
      NS = NIT + JS - 1
      MS = MIT + IS - 1
C
C     SET X - DIRECTION
C
      DO I = 1, MIT
         XI = AIT + FLOAT(IS + I - 2)*DLX
         CALL COFX (XI, AI, BI, CI)
         AXI = (AI/DLX - 0.5*BI)/DLX
         BXI = (-2.*AI/DLX**2) + CI
         CXI = (AI/DLX + 0.5*BI)/DLX
         AM(I) = AXI
         BM(I) = BXI
         CM(I) = CXI
      END DO
C
C     SET Y DIRECTION
C
      DO J = 1, NIT
         YJ = CIT + FLOAT(JS + J - 2)*DLY
         CALL COFY (YJ, DJ, EJ, FJ)
         DYJ = (DJ/DLY - 0.5*EJ)/DLY
         EYJ = (-2.*DJ/DLY**2) + FJ
         FYJ = (DJ/DLY + 0.5*EJ)/DLY
         AN(J) = DYJ
         BN(J) = EYJ
         CN(J) = FYJ
      END DO
C
C     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
C
      AX1 = AM(1)
      CXM = CM(MIT)
      SELECT CASE (KSWX) 
      CASE (2) 
C
C     DIRICHLET-DIRICHLET IN X DIRECTION
C
         AM(1) = 0.0
         CM(MIT) = 0.0
      CASE (5) 
C
C     MIXED-DIRICHLET IN X DIRECTION
C
         AM(1) = 0.0
         BM(1) = BM(1) + 2.*ALPHA*DLX*AX1
         CM(1) = CM(1) + AX1
         CM(MIT) = 0.0
      CASE (3) 
C
C     DIRICHLET-MIXED IN X DIRECTION
C
         AM(1) = 0.0
         AM(MIT) = AM(MIT) + CXM
         BM(MIT) = BM(MIT) - 2.*BETA*DLX*CXM
         CM(MIT) = 0.0
C
C     MIXED - MIXED IN X DIRECTION
C
      CASE (4) 
         AM(1) = 0.0
         BM(1) = BM(1) + 2.*DLX*ALPHA*AX1
         CM(1) = CM(1) + AX1
         AM(MIT) = AM(MIT) + CXM
         BM(MIT) = BM(MIT) - 2.*DLX*BETA*CXM
         CM(MIT) = 0.0
      END SELECT
C
C     ADJUST IN Y DIRECTION UNLESS PERIODIC
C
      DY1 = AN(1)
      FYN = CN(NIT)
      SELECT CASE (KSWY) 
      CASE (2) 
C
C     DIRICHLET-DIRICHLET IN Y DIRECTION
C
         AN(1) = 0.0
         CN(NIT) = 0.0
      CASE (5) 
C
C     MIXED-DIRICHLET IN Y DIRECTION
C
         AN(1) = 0.0
         BN(1) = BN(1) + 2.*DLY*GAMA*DY1
         CN(1) = CN(1) + DY1
         CN(NIT) = 0.0
      CASE (3) 
C
C     DIRICHLET-MIXED IN Y DIRECTION
C
         AN(1) = 0.0
         AN(NIT) = AN(NIT) + FYN
         BN(NIT) = BN(NIT) - 2.*DLY*XNU*FYN
         CN(NIT) = 0.0
      CASE (4) 
C
C     MIXED - MIXED DIRECTION IN Y DIRECTION
C
         AN(1) = 0.0
         BN(1) = BN(1) + 2.*DLY*GAMA*DY1
         CN(1) = CN(1) + DY1
         AN(NIT) = AN(NIT) + FYN
         BN(NIT) = BN(NIT) - 2.0*DLY*XNU*FYN
         CN(NIT) = 0.0
      END SELECT
      IF (KSWX /= 1) THEN
C
C     ADJUST USOL ALONG X EDGE
C
         IF (KSWX==2 .OR. KSWX==3) THEN
            IF (KSWX==2 .OR. KSWX==5) THEN
               USOL(IS,JS:NS) = USOL(IS,JS:NS) - AX1*USOL(1,JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - CXM*USOL(K,JS:NS)
            ELSE
               USOL(IS,JS:NS) = USOL(IS,JS:NS) - AX1*USOL(1,JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - 2.0*DLX*CXM*BDB(JS:NS)
            ENDIF
         ELSE
            IF (KSWX==2 .OR. KSWX==5) THEN
               USOL(IS,JS:NS) = USOL(IS,JS:NS) + 2.0*DLX*AX1*BDA(JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - CXM*USOL(K,JS:NS)
            ELSE
               USOL(IS,JS:NS) = USOL(IS,JS:NS) + 2.0*DLX*AX1*BDA(JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - 2.0*DLX*CXM*BDB(JS:NS)
            ENDIF
         ENDIF
      ENDIF
      IF (KSWY /= 1) THEN
C
C     ADJUST USOL ALONG Y EDGE
C
         IF (KSWY==2 .OR. KSWY==3) THEN
            IF (KSWY==2 .OR. KSWY==5) THEN
               USOL(IS:MS,JS) = USOL(IS:MS,JS) - DY1*USOL(IS:MS,1)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - FYN*USOL(IS:MS,L)
            ELSE
               USOL(IS:MS,JS) = USOL(IS:MS,JS) - DY1*USOL(IS:MS,1)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - 2.0*DLY*FYN*BDD(IS:MS)
            ENDIF
         ELSE
            IF (KSWY==2 .OR. KSWY==5) THEN
               USOL(IS:MS,JS) = USOL(IS:MS,JS) + 2.0*DLY*DY1*BDC(IS:MS)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - FYN*USOL(IS:MS,L)
            ELSE
               USOL(IS:MS,JS) = USOL(IS:MS,JS) + 2.0*DLY*DY1*BDC(IS:MS)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - 2.0*DLY*FYN*BDD(IS:MS)
            ENDIF
         ENDIF
      ENDIF
C
C     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
C
      IF (IORDER == 4) THEN
         GRHS(IS,JS:NS) = USOL(IS,JS:NS)
         GRHS(MS,JS:NS) = USOL(MS,JS:NS)
         GRHS(IS:MS,JS) = USOL(IS:MS,JS)
         GRHS(IS:MS,NS) = USOL(IS:MS,NS)
      ENDIF
      PERTRB = 0.0
C
C     CHECK IF OPERATOR IS SINGULAR
C
      CALL CHKSNG(MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,COFX,COFY,SINGLR)
C
C     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
C     IF SINGULAR
C
      IF (SINGLR) CALL SEPTRI (MIT, AM, BM, CM, DM, UM, ZM)
      IF (SINGLR) CALL SEPTRI (NIT, AN, BN, CN, DN, UN, ZN)
C
C     MAKE INITIALIZATION CALL TO blktrii
C
      IF (INTL == 0) THEN
         CALL BLKTRII (INTL, NP, NIT, AN, BN, CN, MP, MIT, AM, BM, CM, 
     1      IDMN, USOL(IS,JS), IERROR, W, WC)
         IF (IERROR /= 0) RETURN 
      ENDIF
C
C     ADJUST RIGHT HAND SIDE IF NECESSARY
C
      IF (SINGLR) CALL SEPORT (USOL, IDMN, ZN, ZM, PERTRB)
C
C     COMPUTE SOLUTION
C
      CALL BLKTRII (I1, NP, NIT, AN, BN, CN, MP, MIT, AM, BM, CM, IDMN, 
     1   USOL(IS,JS), IERROR, W, WC)
      IF (IERROR /= 0) RETURN 
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX == 1) THEN
         USOL(K,:L) = USOL(1,:L)
      ENDIF
      IF (KSWY == 1) THEN
         USOL(:K,L) = USOL(:K,1)
      ENDIF
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C
      IF (SINGLR) CALL SEPMIN (USOL, IDMN, ZN, ZM, PRTRB)
C
C     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
C     NOT FLAGGED
C
      IF (IORDER == 2) RETURN 
C
C     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
C
      CALL DEFER (COFX, COFY, IDMN, USOL, GRHS)
      IF (SINGLR) CALL SEPORT (USOL, IDMN, ZN, ZM, PERTRB)
C
C     COMPUTE fourth order SOLUTION
C
      CALL BLKTRII (I1, NP, NIT, AN, BN, CN, MP, MIT, AM, BM, CM, IDMN, 
     1   USOL(IS,JS), IERROR, W, WC)
      IF (IERROR /= 0) RETURN 
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX == 1) THEN
         USOL(K,:L) = USOL(1,:L)
      ENDIF
      IF (KSWY == 1) THEN
         USOL(:K,L) = USOL(:K,1)
      ENDIF
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C
      IF (SINGLR) CALL SEPMIN (USOL, IDMN, ZN, ZM, PRTRB)
      RETURN 
      END SUBROUTINE SPELIP


 
      SUBROUTINE CHKPRM(INTL, IORDER, A, B, M, MBDCND, C, D, N, NBDCND, 
     1   COFX, COFY, IDMN, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: INTL
      INTEGER , INTENT(IN) :: IORDER
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER , INTENT(IN) :: IDMN
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      TYPE (fishworkspace) :: w
      EXTERNAL :: COFX,COFY
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J
      REAL :: DLX, DLY, XI, AI, BI, CI, YJ, DJ, EJ, FJ
C-----------------------------------------------
C
C     THIS PROGRAM CHECKS THE INPUT arguments FOR ERRORS
C
C
C     CHECK DEFINITION OF SOLUTION REGION
C
      IF (A>=B .OR. C>=D) THEN
         IERROR = 1
         RETURN 
      ENDIF
c
c     check boundary condition arguments
c
      IF (MBDCND<0 .OR. MBDCND>4) THEN
         IERROR = 2
         RETURN 
      ENDIF
      IF (NBDCND<0 .OR. NBDCND>4) THEN
         IERROR = 3
         RETURN 
      ENDIF
C
C     CHECK FIRST DIMENSION IN CALLING ROUTINE
C
      IF (IDMN < 7) THEN
         IERROR = 5
         RETURN 
      ENDIF
C
C     CHECK M,N
C
      IF (M>IDMN - 1 .OR. M<6) THEN
         IERROR = 6
         RETURN 
      ENDIF
      IF (N < 5) THEN
         IERROR = 7
         RETURN 
      ENDIF
C
C     CHECK IORDER
C
      IF (IORDER/=2 .AND. IORDER/=4) THEN
         IERROR = 8
         RETURN 
      ENDIF
C
C     CHECK INTL
C
      IF (INTL/=0 .AND. INTL/=1) THEN
         IERROR = 9
         RETURN 
      ENDIF
C
C     CHECK THAT EQUATION IS ELLIPTIC (only on initial call)
C
      IF (INTL == 0) THEN
         DLX = (B - A)/FLOAT(M)
         DLY = (D - C)/FLOAT(N)
         DO I = 2, M
            XI = A + FLOAT(I - 1)*DLX
            CALL COFX (XI, AI, BI, CI)
            DO J = 2, N
               YJ = C + FLOAT(J - 1)*DLY
               CALL COFY (YJ, DJ, EJ, FJ)
               IF (AI*DJ > 0.0) CYCLE 
               IERROR = 10
               RETURN 
            END DO
         END DO
      ENDIF
C
C     NO ERROR FOUND
C
      IERROR = 0
      RETURN 
      END SUBROUTINE CHKPRM


      SUBROUTINE CHKSNG(MBDCND, NBDCND, ALPHA, BETA, GAMA, XNU, 
     1   COFX, COFY, SINGLR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: NBDCND
      REAL , INTENT(IN) :: ALPHA
      REAL , INTENT(IN) :: BETA
      REAL , INTENT(IN) :: GAMA
      REAL , INTENT(IN) :: XNU
      LOGICAL , INTENT(OUT) :: SINGLR
      EXTERNAL :: COFX, COFY
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J
      REAL :: XI, AI, BI, CI, YJ, DJ, EJ, FJ
C-----------------------------------------------
C
C     THIS SUBROUTINE CHECKS IF THE PDE   SEPELI
C     MUST SOLVE IS A SINGULAR OPERATOR
C
      SINGLR = .FALSE.
C
C     CHECK IF THE BOUNDARY CONDITIONS ARE
C     ENTIRELY PERIODIC AND/OR MIXED
C
      IF(MBDCND/=0.AND.MBDCND/=3.OR.NBDCND/=0.AND.NBDCND/=3)RETURN
C
C     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
C
      IF (MBDCND == 3) THEN
         IF (ALPHA/=0.0 .OR. BETA/=0.0) RETURN 
      ENDIF
 
      IF (NBDCND == 3) THEN
         IF (GAMA/=0.0 .OR. XNU/=0.0) RETURN 
      ENDIF
C
C     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
C     ARE ZERO
C
      DO I = IS, MS
         XI = AIT + FLOAT(I - 1)*DLX
         CALL COFX (XI, AI, BI, CI)
         IF (CI == 0.0) CYCLE 
         RETURN 
      END DO
      DO J = JS, NS
         YJ = CIT + FLOAT(J - 1)*DLY
         CALL COFY (YJ, DJ, EJ, FJ)
         IF (FJ == 0.0) CYCLE 
         RETURN 
      END DO
C
C     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
C
      SINGLR = .TRUE.
      RETURN 
      END SUBROUTINE CHKSNG


      SUBROUTINE DEFER(COFX, COFY, IDMN, USOL, GRHS)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IDMN
      REAL  :: USOL(IDMN,*)
      REAL , INTENT(INOUT) :: GRHS(IDMN,*)
      EXTERNAL:: COFX, COFY
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, I
      REAL::YJ,DJ,EJ,FJ,XI,AI,BI,CI,UXXX,UXXXX,UYYY,UYYYY,TX,TY
C-----------------------------------------------
C
C     THIS SUBROUTINE FIRST APPROXIMATES THE TRUNCATION ERROR GIVEN BY
C     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY WHERE
C     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 ON THE INTERIOR AND
C     AT THE BOUNDARIES IF PERIODIC(HERE UXXX,UXXXX ARE THE THIRD
C     AND FOURTH PARTIAL DERIVATIVES OF U WITH RESPECT TO X).
C     TX IS OF THE FORM AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
C     AT X=A OR X=B IF THE BOUNDARY CONDITION THERE IS MIXED.
C     TX=0.0 ALONG SPECIFIED BOUNDARIES.  TY HAS SYMMETRIC FORM
C     IN Y WITH X,AFUN(X),BFUN(X) REPLACED BY Y,DFUN(Y),EFUN(Y).
C     THE SECOND ORDER SOLUTION IN USOL IS USED TO APPROXIMATE
C     (VIA SECOND ORDER FINITE DIFFERENCING) THE TRUNCATION ERROR
C     AND THE RESULT IS ADDED TO THE RIGHT HAND SIDE IN GRHS
C     AND THEN TRANSFERRED TO USOL TO BE USED AS A NEW RIGHT
C     HAND SIDE WHEN CALLING BLKTRI FOR A FOURTH ORDER SOLUTION.
C
C
C     COMPUTE TRUNCATION ERROR APPROXIMATION OVER THE ENTIRE MESH
C
      DO J = JS, NS
         YJ = CIT + FLOAT(J - 1)*DLY
         CALL COFY (YJ, DJ, EJ, FJ)
         DO I = IS, MS
            XI = AIT + FLOAT(I - 1)*DLX
            CALL COFX (XI, AI, BI, CI)
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
C
            CALL SEPDX (USOL, IDMN, I, J, UXXX, UXXXX)
            CALL SEPDY (USOL, IDMN, I, J, UYYY, UYYYY)
            TX = AI*UXXXX/12.0 + BI*UXXX/6.0
            TY = DJ*UYYYY/12.0 + EJ*UYYY/6.0
C
C     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
C
            IF (KSWX/=1 .AND. (I==1 .OR. I==K)) TX = AI/3.0*(UXXXX/4.0
     1          + UXXX/DLX)
            IF (KSWY/=1 .AND. (J==1 .OR. J==L)) TY = DJ/3.0*(UYYYY/4.0
     1          + UYYY/DLY)
            GRHS(I,J) = GRHS(I,J) + DLX**2*TX + DLY**2*TY
         END DO
      END DO
C
C     RESET THE RIGHT HAND SIDE IN USOL
C
      USOL(IS:MS,JS:NS) = GRHS(IS:MS,JS:NS)
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    version 5.0, fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE DEFER
C
C     file sepx4.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE SEPX4 (IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
C    +                  NBDCND,BDC,BDD,COFX,GRHS,USOL,IDMN,PERTRB,
C    +                  IERROR)
C
C
C
C DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C ARGUMENTS              USOL(IDMN,N+1),     GRHS(IDMN,N+1),
C
C
C LATEST REVISION        June 2004
C
C PURPOSE                SEPX4 SOLVES FOR EITHER THE SECOND-ORDER
C                        FINITE DIFFERENCE APPROXIMATION OR A
C                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
C                        ELLIPTIC EQUATION
C
C                          AF(X)*UXX+BF(X)*UX+CF(X)*U+UYY = G(X,Y)
C
C                        ON A RECTANGLE (X GREATER THAN OR EQUAL TO
C                        A AND LESS THAN OR EQUAL TO B, Y GREATER THAN
C                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
C                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
C                        CONDITIONS IS ALLOWED.  IF BOUNDARY
C                        CONDITIONS IN THE X DIRECTION ARE PERIODIC
C                        (SEE MBDCND=0 BELOW) THEN THE COEFFICIENTS
C                        MUST SATISFY
C
C                          AF(X)=C1,BF(X)=0,CF(X)=C2 FOR ALL X.
C
C                        HERE C1,C2 ARE CONSTANTS, C1.GT.0.
C
C                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
C                        IN THE X-DIRECTION:
C                          (0) PERIODIC, U(X+B-A,Y)=U(X,Y) FOR
C                              ALL Y,X
C                          (1) U(A,Y), U(B,Y) ARE SPECIFIED FOR ALL Y
C                          (2) U(A,Y), DU(B,Y)/DX+BETA*U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          (3) DU(A,Y)/DX+ALPHA*U(A,Y),DU(B,Y)/DX+
C                              BETA*U(B,Y) ARE SPECIFIED FOR ALL Y
C                          (4) DU(A,Y)/DX+ALPHA*U(A,Y),U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C
C                        IN THE Y-DIRECTION:
C                          (0) PERIODIC, U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                          (1) U(X,C),U(X,D) ARE SPECIFIED FOR ALL X
C                          (2) U(X,C),DU(X,D)/DY ARE SPECIFIED FOR
C                              ALL X
C                          (3) DU(X,C)/DY,DU(X,D)/DY ARE SPECIFIED FOR
C                              ALL X
C                          (4) DU(X,C)/DY,U(X,D) ARE SPECIFIED FOR
C                              ALL X
C
C USAGE                  CALL SEPX4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,
C                                   BETA,C,D,N,NBDCND,BDC,BDD,COFX,
C                                   GRHS,USOL,IDMN,W,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               IORDER
C                          = 2 IF A SECOND-ORDER APPROXIMATION IS
C                              SOUGHT
C                          = 4 IF A FOURTH-ORDER APPROXIMATION IS
C                              SOUGHT
C
c *** caution ***          GRHS SHOULD BE RESET IF SEPX4 WAS FIRST CALLED
C                          WITH IORDER=2 AND WILL BE CALLED AGAIN WITH
C                          IORDER=4.  VALUES IN GRHS ARE DESTROYED BY THE
C                          IORDER=2 CALL.
C
C
C                        A,B
C                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
C                          I.E., X IS GREATER THAN OR EQUAL TO A
C                          AND LESS THAN OR EQUAL TO B.  A MUST BE
C                          LESS THAN B.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (A,B) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE X-
C                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
C                          FOR I=1,2,...,M+1 WHERE DLX=(B-A)/M IS
C                          THE PANEL WIDTH.  M MUST BE LESS THAN
C                          IDMN AND GREATER THAN 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT X=A AND X=B
C                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
C                              U(X+B-A,Y)=U(X,Y)  FOR ALL Y,X
C                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND X=B, I.E., U(A,Y) AND U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND THE BOUNDARY CONDITION IS MIXED AT
C                              X=B, I.E., U(A,Y) AND
C                              DU(B,Y)/DX+BETA*U(B,Y) ARE SPECIFIED
C                              FOR ALL Y
C                          = 3 IF THE BOUNDARY CONDITIONS AT X=A AND
C                              X=B ARE MIXED, I.E.,
C                              DU(A,Y)/DX+ALPHA*U(A,Y) AND
C                              DU(B,Y)/DX+BETA*U(B,Y) ARE SPECIFIED
C                              FOR ALL Y
C                          = 4 IF THE BOUNDARY CONDITION AT X=A IS
C                              MIXED AND THE SOLUTION IS SPECIFIED
C                              AT X=B, I.E., DU(A,Y)/DX+ALPHA*U(A,Y)
C                              AND U(B,Y) ARE SPECIFIED FOR ALL Y
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF
C                          DU(A,Y)/DX+ ALPHA*U(A,Y) AT X=A, WHEN
C                          MBDCND=3 OR 4.
C                          BDA(J) = DU(A,YJ)/DX+ALPHA*U(A,YJ),
C                          J=1,2,...,N+1
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY PARAMETER.
C
C                        ALPHA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN CASE
C                          OF A MIXED BOUNDARY CONDITION AT X=A
C                          (SEE ARGUMENT BDA).  IF MBDCND IS NOT EQUAL
C                          TO EITHER 3 OR 4, THEN ALPHA IS A DUMMY
C                          PARAMETER.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF
C                          DU(B,Y)/DX+ BETA*U(B,Y) AT X=B.
C                          WHEN MBDCND=2 OR 3
C                          BDB(J) = DU(B,YJ)/DX+BETA*U(B,YJ),
C                          J=1,2,...,N+1
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY PARAMETER.
C
C                        BETA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT X=B
C                          (SEE ARGUMENT BDB).  IF MBDCND IS NOT EQUAL
C                          TO 2 OR 3, THEN BETA IS A DUMMY PARAMETER.
C
C                        C,D
C                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
C                          I.E., Y IS GREATER THAN OR EQUAL TO C AND
C                          LESS THAN OR EQUAL TO D.  C MUST BE LESS
C                          THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE Y-
C                          DIRECTION GIVEN BY YJ=C+(J-1)*DLY FOR
C                          J=1,2,...,N+1 WHERE DLY=(D-C)/N IS THE
C                          PANEL WIDTH.  IN ADDITION, N MUST BE
C                          GREATER THAN 4.
C
C                        NBDCND
C                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
C                          AT Y=C AND Y=D
C                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
C                              I.E., U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND Y = D, I.E., U(X,C)  AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND THE BOUNDARY CONDITION IS MIXED
C                              AT Y=D, I.E., DU(X,C)/DY AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
C                              AT Y=CAND Y=D I.E.,
C                              DU(X,D)/DY AND DU(X,D)/DY ARE
C                              SPECIFIED FOR ALL X
C                          = 4 IF THE BOUNDARY CONDITION IS MIXED
C                              AT Y=C AND THE SOLUTION IS SPECIFIED
C                              AT Y=D, I.E. DU(X,C)/DY+GAMA*U(X,C)
C                              AND U(X,D) ARE SPECIFIED FOR ALL X
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUE DU(X,C)/DY AT Y=C.
C
C                          WHEN NBDCND=3 OR 4
C                            BDC(I) = DU(XI,C)/DY I=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY PARAMETER.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIED THE VALUE OF DU(X,D)/DY AT Y=D.
C
C                          WHEN NBDCND=2 OR 3
C                            BDD(I)=DU(XI,D)/DY I=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY PARAMETER.
C
C                        COFX
C                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
C                          X, AFUN, BFUN, CFUN WHICH RETURNS THE
C                          VALUES OF THE X-DEPENDENT COEFFICIENTS
C                          AF(X), BF(X), CF(X) IN THE ELLIPTIC
C                          EQUATION AT X.  IF BOUNDARY CONDITIONS IN
C                          THE X DIRECTION ARE PERIODIC THEN THE
C                          COEFFICIENTS MUST SATISFY AF(X)=C1,BF(X)=0,
C                          CF(X)=C2 FOR ALL X.  HERE C1.GT.0
C                          AND C2 ARE CONSTANTS.
C
C                          NOTE THAT COFX MUST BE DECLARED EXTERNAL
C                          IN THE CALLING ROUTINE.
C
C                        GRHS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT-HAND SIDE OF THE
C                          ELLIPTIC EQUATION, I.E.,GRHS(I,J)=G(XI,YI),
C                          FOR I=2,...,M, J=2,...,N.  AT THE
C                          BOUNDARIES, GRHS IS DEFINED BY
C
C                          MBDCND   GRHS(1,J)   GRHS(M+1,J)
C                          ------   ---------   -----------
C                            0      G(A,YJ)     G(B,YJ)
C                            1         *           *
C                            2         *        G(B,YJ)  J=1,2,...,N+1
C                            3      G(A,YJ)     G(B,YJ)
C                            4      G(A,YJ)        *
C
C                          NBDCND   GRHS(I,1)   GRHS(I,N+1)
C                          ------   ---------   -----------
C                            0      G(XI,C)     G(XI,D)
C                            1         *           *
C                            2         *        G(XI,D)  I=1,2,...,M+1
C                            3      G(XI,C)     G(XI,D)
C                            4      G(XI,C)        *
C
C                          WHERE * MEANS THESE QUANTITES ARE NOT USED.
C                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
c *** caution              GRHS SHOULD BE RESET IF SEPX4 WAS FIRST CALLED
C                          WITH IORDER=2 AND WILL BE CALLED AGAIN WITH
C                          IORDER=4.  VALUES IN GRHS ARE DESTROYED BY THE
C                          IORDER=2 CALL.
C
C                        USOL
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE SOLUTION ALONG THE BOUNDARIES.
C                          AT THE BOUNDARIES, USOL IS DEFINED BY
C
C                          MBDCND   USOL(1,J)   USOL(M+1,J)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(A,YJ)     U(B,YJ)
C                            2      U(A,YJ)        *     J=1,2,...,N+1
C                            3         *           *
C                            4         *        U(B,YJ)
C
C                          NBDCND   USOL(I,1)   USOL(I,N+1)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(XI,C)     U(XI,D)
C                            2      U(XI,C)        *     I=1,2,...,M+1
C                            3         *           *
C                            4         *        U(XI,D)
C
C                          WHERE * MEANS THE QUANTITES ARE NOT USED
C                          IN THE SOLUTION.
C
C                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
C                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
C                          CASE THE TABLES SPECIFYING THE BOUNDARIES
C                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
C                          BOUNDARIES UNIQUELY EXCEPT AT THE CORNERS.
C                          IF THE TABLES CALL FOR BOTH G(X,Y) AND
C                          U(X,Y) AT A CORNER THEN THE SOLUTION MUST
C                          BE CHOSEN.
C                          FOR EXAMPLE, IF MBDCND=2 AND NBDCND=4,
C                          THEN U(A,C), U(A,D),U(B,D) MUST BE CHOSEN
C                          AT THE CORNERS IN ADDITION TO G(B,C).
C
C                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
C                          GRHS, MUST BE DISTINCT.
C
C                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
C                        IDMN
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAYS
C                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
C                          CALLING SEPELI.  THIS PARAMETER IS USED
C                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
C                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
C                          GREATER THAN OR EQUAL TO M+1.
C
C
C ON OUTPUT              USOL
C                          CONTAINS THE APPROXIMATE SOLUTION TO THE
C                          ELLIPTIC EQUATION. USOL(I,J) IS THE
C                          APPROXIMATION TO U(XI,YJ) FOR I=1,2...,M+1
C                          AND J=1,2,...,N+1.  THE APPROXIMATION HAS
C                          ERROR O(DLX**2+DLY**2) IF CALLED WITH
C                          IORDER=2 AND O(DLX**4+DLY**4) IF CALLED
C                          WITH IORDER=4.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS (I.E., ALPHA=BETA=0 IF
C                          MBDCND=3) IS SPECIFIED AND IF CF(X)=0 FOR
C                          ALL X THEN A SOLUTION TO THE DISCRETIZED
C                          MATRIX EQUATION MAY NOT EXIST
C                          (REFLECTING THE NON-UNIQUENESS OF SOLUTIONS
C                          TO THE PDE).
C                          PERTRB IS A CONSTANT CALCULATED AND
C                          SUBTRACTED FROM THE RIGHT HAND SIDE OF THE
C                          MATRIX EQUATION INSURING THE EXISTENCE OF A
C                          SOLUTION.  SEPX4 COMPUTES THIS SOLUTION
C                          WHICH IS A WEIGHTED MINIMAL LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL PROBLEM.  IF
C                          SINGULARITY IS NOT DETECTED PERTRB=0.0 IS
C                          RETURNED BY SEPX4.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS OR FAILURE TO FIND A SOLUTION
C
C                          =  0 NO ERROR
C                          =  1 IF A GREATER THAN B OR C GREATER
C                               THAN D
C                          =  2 IF MBDCND LESS THAN 0 OR MBDCND
C                               GREATER THAN 4
C                          =  3 IF NBDCND LESS THAN 0 OR NBDCND
C                               GREATER THAN 4
C                          =  4 IF ATTEMPT TO FIND A SOLUTION FAILS.
C                               (THE LINEAR SYSTEM GENERATED IS NOT
C                               DIAGONALLY DOMINANT.)
C                          =  5 IF IDMN IS TOO SMALL (SEE DISCUSSION
C                               OF IDMN)
C                          =  6 IF M IS TOO SMALL OR TOO LARGE
C                               (SEE DISCUSSION OF M)
C                          =  7 IF N IS TOO SMALL (SEE DISCUSSION OF N)
C                          =  8 IF IORDER IS NOT 2 OR 4
C                          =  9 IF INTL IS NOT 0 OR 1
C                          = 10 IF AFUN IS LESS THAN OR EQUAL TO ZERO
C                               FOR SOME INTERIOR MESH POINT XI SOME
C                               INTERIOR MESH POINT (XI,YJ)
C                          = 12 IF MBDCND=0 AND AF(X)=CF(X)=CONSTANT
C                               OR BF(X)=0 FOR ALL X IS NOT TRUE.
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C REQUIRED files         fish.f,comf.f,genbun.f,gnbnaux.f,sepaux.f
C
C
C PRECISION              SINGLE
C
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                SEPX4 WAS DEVELOPED AT NCAR BY JOHN C.
C                        ADAMS OF THE SCIENTIFIC COMPUTING DIVISION
C                        IN OCTOBER 1978.  THE BASIS OF THIS CODE IS
C                        NCAR ROUTINE SEPELI.  BOTH PACKAGES WERE
C                        RELEASED ON NCAR'S PUBLIC LIBRARIES IN
C                        JANUARY 1980. SEPX4 was modified in June 2004
c                        incorporating fortran 90 dynamical storage
c                        allocation for work space requirements
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              SEPX4 AUTOMATICALLY DISCRETIZES THE SEPARABLE
C                        ELLIPTIC EQUATION WHICH IS THEN SOLVED BY A
C                        GENERALIZED CYCLIC REDUCTION ALGORITHM IN THE
C                        SUBROUTINE POIS.  THE FOURTH ORDER SOLUTION
C                        IS OBTAINED USING THE TECHNIQUE OF DEFFERRED
C                        CORRECTIONS REFERENCED BELOW.
C
C TIMING                 WHEN POSSIBLE, SEPX4 SHOULD BE USED INSTEAD
C                        OF PACKAGE SEPELI.  THE INCREASE IN SPEED
C                        IS AT LEAST A FACTOR OF THREE.
C
C REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
C                        BOUNDARY-VALUE PROBLEMS, BLAISDEL (1968),
C                        WALTHAM, MASS.
C
C                        SWARZTRAUBER, P., AND R. SWEET (1975):
C                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
C                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
C                        EQUATIONS.  NCAR TECHNICAL NOTE
C                          NCAR-TN/IA-109, PP. 135-137.
C***********************************************************************
      SUBROUTINE SEPX4(IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, 
     1   C, D, N, NBDCND, BDC, BDD, COFX, GRHS, USOL, IDMN, PERTRB, 
     2   IERROR)

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IORDER
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDMN
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: ALPHA
      REAL  :: BETA
      REAL  :: C
      REAL  :: D
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: GRHS(IDMN,*)
      REAL  :: USOL(IDMN,*)
      EXTERNAL :: COFX
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      TYPE (fishworkspace) :: W
      INTEGER :: L, K, LOG2N, LENGTH, IRWK, ICWK, I1, I2, I3, I4, I5, I6
     1   , I7, I8, I9, I10, I11, I12, I13
C-----------------------------------------------
C
C     CHECK INPUT PARAMETERS
C
      CALL C4KPRM(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR)
      IF (IERROR /= 0) RETURN 
C
C     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
C
      L = N + 1
      IF (NBDCND == 0) L = N
      K = M + 1
      L = N + 1
C     ESTIMATE LOG BASE 2 OF N
      LOG2N = INT(ALOG(FLOAT(N + 1))/ALOG(2.0) + 0.5)
!     set required work space estimate
      LENGTH = 4*(N + 1) + (10 + LOG2N)*(M + 1)
      IRWK = LENGTH + 6*(K + L) + 1
      ICWK = 0
!     allocate work space
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed
      IF (IERROR == 20) RETURN 
      IERROR = 0
C
C     SET WORK SPACE INDICES
C
      I1 = LENGTH + 1
      I2 = I1 + L
      I3 = I2 + L
      I4 = I3 + L
      I5 = I4 + L
      I6 = I5 + L
      I7 = I6 + L
      I8 = I7 + K
      I9 = I8 + K
      I10 = I9 + K
      I11 = I10 + K
      I12 = I11 + K
      I13 = 1
      CALL S4ELIP(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     +NBDCND,BDC,BDD,COFX,w%rew(I1),w%rew(I2),w%rew(I3),
     +w%rew(i4),w%rew(i5),w%rew(i6),w%rew(i7),w%rew(i8),
     +w%rew(i9),w%rew(i10),w%rew(i11),w%rew(i12),
     +GRHS,USOL,IDMN,w%rew(i13),PERTRB,IERROR)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE SEPX4


      SUBROUTINE S4ELIP(IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, 
     1   C, D, N, NBDCND, BDC, BDD, COFX, AN, BN, CN, DN, UN, ZN, AM, BM
     2   , CM, DM, UM, ZM, GRHS, USOL, IDMN, W, PERTRB, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IORDER
      INTEGER , INTENT(IN) :: M
      INTEGER  :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDMN
      INTEGER , INTENT(INOUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL  :: ALPHA
      REAL  :: BETA
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL  :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL  :: DN(*)
      REAL  :: UN(*)
      REAL  :: ZN(*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL  :: DM(*)
      REAL  :: UM(*)
      REAL  :: ZM(*)
      REAL  :: GRHS(IDMN,*)
      REAL  :: USOL(IDMN,*)
      REAL  :: W(*)
      EXTERNAL:: COFX
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J, I1, MP, NP, IEROR
      REAL :: XI, AI, BI, CI, AXI, BXI, CXI, DYJ, EYJ, FYJ, AX1, CXM, 
     1   DY1, FYN, GAMA, XNU, PRTRB
      LOGICAL :: SINGLR
C-----------------------------------------------
C
C     S4ELIP SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
C     AND COMPUTES A SECOND ORDER SOLUTION IN USOL.  A RETURN JUMP TO
C     SEPELI OCCURRS IF IORDER=2.  IF IORDER=4 A FOURTH ORDER
C     SOLUTION IS GENERATED IN USOL.
C
C
C     SET PARAMETERS INTERNALLY
C
      KSWX = MBDCND + 1
      KSWY = NBDCND + 1
      K = M + 1
      L = N + 1
      AIT = A
      BIT = B
      CIT = C
      DIT = D
      DLY = (DIT - CIT)/FLOAT(N)
C
C     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
C     AND NON-SPECIFIED BOUNDARIES.
C
      USOL(2:M,2:N) = DLY**2*GRHS(2:M,2:N)
      IF (KSWX/=2 .AND. KSWX/=3) THEN
         USOL(1,2:N) = DLY**2*GRHS(1,2:N)
      ENDIF
      IF (KSWX/=2 .AND. KSWX/=5) THEN
         USOL(K,2:N) = DLY**2*GRHS(K,2:N)
      ENDIF
      IF (KSWY/=2 .AND. KSWY/=3) THEN
         USOL(2:M,1) = DLY**2*GRHS(2:M,1)
      ENDIF
      IF (KSWY/=2 .AND. KSWY/=5) THEN
         USOL(2:M,L) = DLY**2*GRHS(2:M,L)
      ENDIF
      IF (KSWX/=2 .AND. KSWX/=3 .AND. KSWY/=2 .AND. KSWY/=3) USOL(1,1)
     1    = DLY**2*GRHS(1,1)
      IF (KSWX/=2 .AND. KSWX/=5 .AND. KSWY/=2 .AND. KSWY/=3) USOL(K,1)
     1    = DLY**2*GRHS(K,1)
      IF (KSWX/=2 .AND. KSWX/=3 .AND. KSWY/=2 .AND. KSWY/=5) USOL(1,L)
     1    = DLY**2*GRHS(1,L)
      IF (KSWX/=2 .AND. KSWX/=5 .AND. KSWY/=2 .AND. KSWY/=5) USOL(K,L)
     1    = DLY**2*GRHS(K,L)
      I1 = 1
C
C     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
C
      MP = 1
      IF (KSWX == 1) MP = 0
      NP = NBDCND
C
C     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
C     IN NINT,MINT
C
      DLX = (BIT - AIT)/FLOAT(M)
      MIT = K - 1
      IF (KSWX == 2) MIT = K - 2
      IF (KSWX == 4) MIT = K
      DLY = (DIT - CIT)/FLOAT(N)
      NIT = L - 1
      IF (KSWY == 2) NIT = L - 2
      IF (KSWY == 4) NIT = L
      TDLX3 = 2.0*DLX**3
      DLX4 = DLX**4
      TDLY3 = 2.0*DLY**3
      DLY4 = DLY**4
C
C     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
C
      IS = 1
      JS = 1
      IF (KSWX==2 .OR. KSWX==3) IS = 2
      IF (KSWY==2 .OR. KSWY==3) JS = 2
      NS = NIT + JS - 1
      MS = MIT + IS - 1
C
C     SET X - DIRECTION
C
      DO I = 1, MIT
         XI = AIT + FLOAT(IS + I - 2)*DLX
         CALL COFX (XI, AI, BI, CI)
         AXI = (AI/DLX - 0.5*BI)/DLX
         BXI = (-2.*AI/DLX**2) + CI
         CXI = (AI/DLX + 0.5*BI)/DLX
         AM(I) = DLY**2*AXI
         BM(I) = DLY**2*BXI
         CM(I) = DLY**2*CXI
      END DO
C
C     SET Y DIRECTION
C
      DYJ = 1.0
      EYJ = -2.0
      FYJ = 1.0
      AN(:NIT) = DYJ
      BN(:NIT) = EYJ
      CN(:NIT) = FYJ
C
C     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
C
      AX1 = AM(1)
      CXM = CM(MIT)
      SELECT CASE (KSWX) 
      CASE (2) 
C
C     DIRICHLET-DIRICHLET IN X DIRECTION
C
         AM(1) = 0.0
         CM(MIT) = 0.0
      CASE (3) 
C
C     DIRICHLET-MIXED IN X DIRECTION
C
         AM(1) = 0.0
         AM(MIT) = AM(MIT) + CXM
         BM(MIT) = BM(MIT) - 2.*BETA*DLX*CXM
         CM(MIT) = 0.0
      CASE (4) 
C
C     MIXED - MIXED IN X DIRECTION
C
         AM(1) = 0.0
         BM(1) = BM(1) + 2.*DLX*ALPHA*AX1
         CM(1) = CM(1) + AX1
         AM(MIT) = AM(MIT) + CXM
         BM(MIT) = BM(MIT) - 2.*DLX*BETA*CXM
         CM(MIT) = 0.0
      CASE (5) 
C
C     MIXED-DIRICHLET IN X DIRECTION
C
         AM(1) = 0.0
         BM(1) = BM(1) + 2.*ALPHA*DLX*AX1
         CM(1) = CM(1) + AX1
         CM(MIT) = 0.0
      END SELECT
c
C     ADJUST IN Y DIRECTION UNLESS PERIODIC
C
      DY1 = AN(1)
      FYN = CN(NIT)
      GAMA = 0.0
      XNU = 0.0
      SELECT CASE (KSWY) 
      CASE (2) 
C
C     DIRICHLET-DIRICHLET IN Y DIRECTION
C
         AN(1) = 0.0
         CN(NIT) = 0.0
      CASE (3) 
C
C     DIRICHLET-MIXED IN Y DIRECTION
C
         AN(1) = 0.0
         AN(NIT) = AN(NIT) + FYN
         BN(NIT) = BN(NIT) - 2.*DLY*XNU*FYN
         CN(NIT) = 0.0
      CASE (4) 
C
C     MIXED - MIXED DIRECTION IN Y DIRECTION
C
         AN(1) = 0.0
         BN(1) = BN(1) + 2.*DLY*GAMA*DY1
         CN(1) = CN(1) + DY1
         AN(NIT) = AN(NIT) + FYN
         BN(NIT) = BN(NIT) - 2.0*DLY*XNU*FYN
         CN(NIT) = 0.0
      CASE (5) 
C
C     MIXED-DIRICHLET IN Y DIRECTION
C
         AN(1) = 0.0
         BN(1) = BN(1) + 2.*DLY*GAMA*DY1
         CN(1) = CN(1) + DY1
         CN(NIT) = 0.0
      END SELECT
      IF (KSWX /= 1) THEN
C
C     ADJUST USOL ALONG X EDGE
C
         IF (KSWX==2 .OR. KSWX==3) THEN
            IF (KSWX==2 .OR. KSWX==5) THEN
               USOL(IS,JS:NS) = USOL(IS,JS:NS) - AX1*USOL(1,JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - CXM*USOL(K,JS:NS)
            ELSE
               USOL(IS,JS:NS) = USOL(IS,JS:NS) - AX1*USOL(1,JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - 2.0*DLX*CXM*BDB(JS:NS)
            ENDIF
         ELSE
            IF (KSWX==2 .OR. KSWX==5) THEN
               USOL(IS,JS:NS) = USOL(IS,JS:NS) + 2.0*DLX*AX1*BDA(JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - CXM*USOL(K,JS:NS)
            ELSE
               USOL(IS,JS:NS) = USOL(IS,JS:NS) + 2.0*DLX*AX1*BDA(JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - 2.0*DLX*CXM*BDB(JS:NS)
            ENDIF
         ENDIF
      ENDIF
      IF (KSWY /= 1) THEN
C
C     ADJUST USOL ALONG Y EDGE
C
         IF (KSWY==2 .OR. KSWY==3) THEN
            IF (KSWY==2 .OR. KSWY==5) THEN
               USOL(IS:MS,JS) = USOL(IS:MS,JS) - DY1*USOL(IS:MS,1)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - FYN*USOL(IS:MS,L)
            ELSE
               USOL(IS:MS,JS) = USOL(IS:MS,JS) - DY1*USOL(IS:MS,1)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - 2.0*DLY*FYN*BDD(IS:MS)
            ENDIF
         ELSE
            IF (KSWY==2 .OR. KSWY==5) THEN
               USOL(IS:MS,JS) = USOL(IS:MS,JS) + 2.0*DLY*DY1*BDC(IS:MS)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - FYN*USOL(IS:MS,L)
            ELSE
               USOL(IS:MS,JS) = USOL(IS:MS,JS) + 2.0*DLY*DY1*BDC(IS:MS)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - 2.0*DLY*FYN*BDD(IS:MS)
            ENDIF
         ENDIF
      ENDIF
C
C     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
C
      IF (IORDER == 4) THEN
         GRHS(IS,JS:NS) = USOL(IS,JS:NS)
         GRHS(MS,JS:NS) = USOL(MS,JS:NS)
         GRHS(IS:MS,JS) = USOL(IS:MS,JS)
         GRHS(IS:MS,NS) = USOL(IS:MS,NS)
      ENDIF
      PERTRB = 0.0
C
C     CHECK IF OPERATOR IS SINGULAR
C
      CALL C4KSNG (MBDCND, NBDCND, ALPHA, BETA, COFX, SINGLR)
C
C     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
C     IF SINGULAR
C
      IF (SINGLR) CALL SEPTRI (MIT, AM, BM, CM, DM, UM, ZM)
      IF (SINGLR) CALL SEPTRI (NIT, AN, BN, CN, DN, UN, ZN)
C
C     ADJUST RIGHT HAND SIDE IF NECESSARY
C
      IF (SINGLR) CALL SEPORT (USOL, IDMN, ZN, ZM, PERTRB)
C
C     COMPUTE SOLUTION
C
C     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
      GRHS(IS:MS,JS:NS) = USOL(IS:MS,JS:NS)
      CALL GENBUN(NP,NIT,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS),IEROR)
C
C     CHECK IF ERROR DETECTED IN POIS
C     THIS CAN ONLY CORRESPOND TO IERROR=12
      IF (IEROR /= 0) THEN
C       SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
         IERROR = 12
         RETURN 
      ENDIF
      IF (IERROR /= 0) RETURN 
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX == 1) THEN
         USOL(K,:L) = USOL(1,:L)
      ENDIF
      IF (KSWY == 1) THEN
         USOL(:K,L) = USOL(:K,1)
      ENDIF
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C
      IF (SINGLR) CALL SEPMIN (USOL, IDMN, ZN, ZM, PRTRB)
C
C     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
C     NOT FLAGGED
C
      IF (IORDER == 2) RETURN 
C
C     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
C
      CALL D4FER (COFX, IDMN, USOL, GRHS)
      IF (SINGLR) CALL SEPORT (USOL, IDMN, ZN, ZM, PERTRB)
C
C     COMPUTE SOLUTION
C
C     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
      GRHS(IS:MS,JS:NS) = USOL(IS:MS,JS:NS)
      CALL GENBUN(NP,NIT,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS),IEROR)
C     CHECK IF ERROR DETECTED IN POIS
C     THIS CAN ONLY CORRESPOND TO IERROR=12
      IF (IEROR /= 0) THEN
C       SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
         IERROR = 12
         RETURN 
      ENDIF
      IF (IERROR /= 0) RETURN 
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX == 1) THEN
         USOL(K,:L) = USOL(1,:L)
      ENDIF
      IF (KSWY == 1) THEN
         USOL(:K,L) = USOL(:K,1)
      ENDIF
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C
      IF (SINGLR) CALL SEPMIN (USOL, IDMN, ZN, ZM, PRTRB)
      RETURN 
      END SUBROUTINE S4ELIP


      SUBROUTINE C4KPRM(IORDER, A, B, M, MBDCND, C, D, N, NBDCND, COFX, 
     1   IDMN, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IORDER
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER , INTENT(IN) :: IDMN
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      EXTERNAL:: COFX
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I
      REAL :: DLX, XI, AI, BI, CI
C-----------------------------------------------
C
C     THIS PROGRAM CHECKS THE INPUT PARAMETERS FOR ERRORS
C
C
C
C     CHECK DEFINITION OF SOLUTION REGION
C
      IF (A>=B .OR. C>=D) THEN
         IERROR = 1
         RETURN 
      ENDIF
C
C     CHECK BOUNDARY SWITCHES
C
      IF (MBDCND<0 .OR. MBDCND>4) THEN
         IERROR = 2
         RETURN 
      ENDIF
      IF (NBDCND<0 .OR. NBDCND>4) THEN
         IERROR = 3
         RETURN 
      ENDIF
C
C     CHECK FIRST DIMENSION IN CALLING ROUTINE
C
      IF (IDMN < 7) THEN
         IERROR = 5
         RETURN 
      ENDIF
C
C     CHECK M
C
      IF (M>IDMN - 1 .OR. M<6) THEN
         IERROR = 6
         RETURN 
      ENDIF
C
C     CHECK N
C
      IF (N < 5) THEN
         IERROR = 7
         RETURN 
      ENDIF
C
C     CHECK IORDER
C
      IF (IORDER/=2 .AND. IORDER/=4) THEN
         IERROR = 8
         RETURN 
      ENDIF
C
C     CHECK THAT EQUATION IS ELLIPTIC
C
      DLX = (B - A)/FLOAT(M)
      DO I = 2, M
         XI = A + FLOAT(I - 1)*DLX
         CALL COFX (XI, AI, BI, CI)
         IF (AI > 0.0) CYCLE 
         IERROR = 10
         RETURN 
      END DO
C
C     NO ERROR FOUND
C
      IERROR = 0
      RETURN 
      END SUBROUTINE C4KPRM


      SUBROUTINE C4KSNG(MBDCND, NBDCND, ALPHA, BETA, COFX, SINGLR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: NBDCND
      REAL , INTENT(IN) :: ALPHA
      REAL , INTENT(IN) :: BETA
      LOGICAL , INTENT(OUT) :: SINGLR
      EXTERNAL:: COFX
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I
      REAL :: XI, AI, BI, CI
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS SUBROUTINE CHECKS IF THE PDE   SEPELI
C     MUST SOLVE IS A SINGULAR OPERATOR
C
      SINGLR = .FALSE.
C
C     CHECK IF THE BOUNDARY CONDITIONS ARE
C     ENTIRELY PERIODIC AND/OR MIXED
C
      IF(MBDCND/=0.AND.MBDCND/=3.OR.NBDCND/=0.AND.NBDCND/=3)RETURN
C
C     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
C
      IF (MBDCND == 3) THEN
         IF (ALPHA/=0.0 .OR. BETA/=0.0) RETURN 
      ENDIF
C
C     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
C     ARE ZERO
C
      DO I = IS, MS
         XI = AIT + FLOAT(I - 1)*DLX
         CALL COFX (XI, AI, BI, CI)
         IF (CI == 0.0) CYCLE 
         RETURN 
      END DO
C
C     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
C
      SINGLR = .TRUE.
      RETURN 
      END SUBROUTINE C4KSNG


      SUBROUTINE D4FER(COFX, IDMN, USOL, GRHS)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IDMN
      REAL  :: USOL(IDMN,*)
      REAL , INTENT(INOUT) :: GRHS(IDMN,*)
      EXTERNAL:: COFX
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /FISH_SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J
      REAL :: XI, AI, BI, CI, UXXX, UXXXX, UYYY, UYYYY, TX, TY
C-----------------------------------------------
C
C     THIS SUBROUTINE FIRST APPROXIMATES THE TRUNCATION ERROR GIVEN BY
C     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY WHERE
C     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 ON THE INTERIOR AND
C     AT THE BOUNDARIES IF PERIODIC(HERE UXXX,UXXXX ARE THE THIRD
C     AND FOURTH PARTIAL DERIVATIVES OF U WITH RESPECT TO X).
C     TX IS OF THE FORM AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
C     AT X=A OR X=B IF THE BOUNDARY CONDITION THERE IS MIXED.
C     TX=0.0 ALONG SPECIFIED BOUNDARIES.  TY HAS SYMMETRIC FORM
C     IN Y WITH X,AFUN(X),BFUN(X) REPLACED BY Y,DFUN(Y),EFUN(Y).
C     THE SECOND ORDER SOLUTION IN USOL IS USED TO APPROXIMATE
C     (VIA SECOND ORDER FINITE DIFFERENCING) THE TRUNCATION ERROR
C     AND THE RESULT IS ADDED TO THE RIGHT HAND SIDE IN GRHS
C     AND THEN TRANSFERRED TO USOL TO BE USED AS A NEW RIGHT
C     HAND SIDE WHEN CALLING BLKTRI FOR A FOURTH ORDER SOLUTION.
C
C
C
C     COMPUTE TRUNCATION ERROR APPROXIMATION OVER THE ENTIRE MESH
C
      DO I = IS, MS
         XI = AIT + FLOAT(I - 1)*DLX
         CALL COFX (XI, AI, BI, CI)
         DO J = JS, NS
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
C
            CALL SEPDX (USOL, IDMN, I, J, UXXX, UXXXX)
            CALL SEPDY (USOL, IDMN, I, J, UYYY, UYYYY)
            TX = AI*UXXXX/12.0 + BI*UXXX/6.0
            TY = UYYYY/12.0
C
C     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
C
            IF (KSWX/=1 .AND. (I==1 .OR. I==K)) TX = AI/3.0*(UXXXX/4.0
     1          + UXXX/DLX)
            IF(KSWY/=1.AND.(J==1.OR.J==L))TY=(UYYYY/4.0+UYYY/DLY)/3.0
            GRHS(I,J) = GRHS(I,J) + DLY**2*(DLX**2*TX + DLY**2*TY)
         END DO
      END DO
C
C     RESET THE RIGHT HAND SIDE IN USOL
C
      USOL(IS:MS,JS:NS) = GRHS(IS:MS,JS:NS)
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE D4FER

      end module FISH
