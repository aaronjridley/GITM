!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!  code will use :
!
! initialisation:
!	call exp_tab8()		! prepare tabulated exponentials
!	call prepCorrUbar()	! or it will be called automatically
!	call prep_projE(hnuGr(0:nbGr),nbGR)


!------- radiom.f90
!

! \
MODULE CRASH_M_projE
  ! /
  !  {U=hnu/Kte,Erad,Brad} will be projected on a working grid , with even log. spacing of U
  !
  implicit none
  PRIVATE !Except
  public:: prep_projE, projSP_E
  !\
  ! Minimal and maximal values for E/T_e ratio
  !/
  real,parameter :: Umin=1d-2, Umax=1.50d2
  public::UMin, UMax
  
  !\
  ! Number of grid points for u=E/T_e (between UMin, UMax)
  !/
  real,parameter :: nbUout =100

  integer,parameter :: mxOut= 300
  public::mxOut
  !\
  ! log(UMax/UMin)/nBOut  - set in prep_proj
  !/
  real,save         :: lgdu
  !\
  ! exp(lgdu)
  !/
  real,save         :: rdu 

  real,save         :: Efirst,Elast
  real,save,dimension(mxOut+1) :: Eout,Uout,SPout
  public::EFirst, ELast,EOut,UOut,SPOut
  !\
  ! Number of energy grouprs in the input 
  ! file of radiation group energies
  !/
  integer,save :: nbIn=0

  !\
  ! Number of grid points in the internal energy grid
  ! between 
  ! = log(Elast/Efirst)/lgdu+1
  !/
  integer,save :: nbOut=0
  public:: nbIn, nBOut

  integer,save :: nbContrib
  integer,parameter :: mxContrib=150000
  integer,save :: to_ctrb(mxContrib),fr_ctrb(mxContrib)
  real,save :: coef_ctrb(mxContrib)
  integer,save :: ctrb_to1(mxContrib),ctrb_to2(mxContrib)


  !-------
contains
  subroutine prep_projE(Ein,nbE)	
    !  Initialization of the projection coefficients 
    !\
    !Inputs
    !/
    !Number of energy groups
    integer,intent(In) :: nbE
    !Photon eergies
    real,intent(In),dimension(0:nbE) :: Ein

    real :: u1,u2,du
    integer :: n,nIn,nOut 
    logical :: b2
    integer :: fr,to
    real :: c
    
    if(nbe<=1)then
       write(*,*)' -E- cannot use RADIOM w/ ng.group=',nbe
       call CON_stop('Stopped')
    end if
    nbIn=nbE

    lgdu=log(Umax/Umin)/nbUout
    rdu=exp(lgdu)
    Efirst=Ein(0)
    Elast=Ein(nbIn)
    nbOut=1+int(log(Elast/Efirst)/lgdu)

    if(nbOut.gt.mxOut)then
       write(*,125) nbOut,mxOut
125    format(// '-E- nbOut(=',i5,') should be LESS THAN mxOut(=',i5,')',//)
       call CON_stop('-R- : prep_projE  nbOut .GT. mxOUT')
    end if
 
    Eout(1)=Efirst

    do n=2,nbOut+1
       Eout(n)=Eout(n-1)*rdu
    end do
 
!    first_ctrb(1:nbOut+1)=0		!! handling 		<<<<<<<
!    last_ctrb(0t:nbOu)=0		!! some exceptions 	<<<<<<<

 
    ! 
    nOut=0
    nbContrib=0
    du=Eout(nbOut+1)-Eout(nbOut)
    LOOP10:	do nIn=1,nbIn	! U(nIn-1)-U(nIn)
       do while (Ein(nIn-1).ge.Eout(nOut+1)) 
          nOut = nOut+1
          du=Eout(nOut+1)-Eout(nOut)
       end do
       if(Ein(nIn).le.Eout(nOut)) then
          cycle LOOP10
       end if
2      u1=max(Ein(nIn-1),Eout(nOut))
       b2=Ein(nIn).lt.Eout(nOut+1)
       if(b2) then
          u2=Ein(nIn)
       else
          u2=Eout(nOut+1)
       end if
       nbContrib=nbContrib+1
       if(nbContrib.le.mxContrib) then
          to_ctrb(nbContrib)=nOut
          fr_ctrb(nbContrib)=nIn	! in group  [ Ein(nIn-1) : Ein(nIn) ]
          coef_ctrb(nbContrib)=(u2-u1)/du
       end if
       !
       !   Ein :	  |      |
       !   Eout:	    |..a   b
       !   u1,u2           |..A B
       !
       !   Ein :	    |       |
       !   Eout:	  a    b..|
       !   u1,u2           A  B..|
       !
       if(.not.b2)then	! input segment overlaps 
          if(nOut.ge.nbOut) exit LOOP10	! more than 1 output segment
          nOut=nOut+1
	!  first_ctrb(nOut)=nIn		!!		 <<<<<<<
          du=Eout(nOut+1)-Eout(nOut)
          goto 2
       end if
    end do LOOP10
    ! 
  
    if(nbContrib.ge.mxContrib) then
       write(*,*)'-D- needs to increase "mxContrib" to ',nbContrib
       call CON_stop('')
    else
       do n=1,nbContrib
          ctrb_to1(n)=0
          ctrb_to2(n)=0
       end do
       do n=1,nbContrib
          ctrb_to2(to_ctrb(n))=n
       end do
       do n=nbContrib,1,-1
          ctrb_to1(to_ctrb(n))=n
       end do
    end if
  end subroutine prep_projE
  !========
  subroutine projSP_E(Te,SPin,nOut,gotoLTE)
    real,intent(In) :: Te,SPin(nbIn)
    integer,intent(Out) :: nOut
    real :: Ufirst,Ulast ,uBef,uAft,c,r,u
    integer :: n1,n2,n3,ctr,fr,to,nBef,nAft,n,tt1,tt2
    real, parameter :: zero=0.d0,one=1.d0
    logical,intent(out) :: gotoLTE
    ! 
    if(nbIn.eq.0) then
       write(*,*)'-P- subroutine "prep_projE" has not been called'
       call CON_stop('')
    end if
    Ufirst=Efirst/Te
    Ulast=Elast/Te
    gotoLTE=.false.
    !
    ! special cases:  range of interest and group span dont overlap:
    if(Ulast.le.Umin) then
       r=zero
       goto 110
    elseif(Ufirst.gt.Umax) then
       gotoLTE=.true.
       n1=nbContrib			!! 	<<<<<<<
       n2=nbContrib-1			!! 	<<<<<<<
       r=one
       goto 110
    end if
    !
    !if(Ulast<10)then
    !   write(*,*) '-W- revise group binning,Ulast(='&
    !        ,Ulast,') found<10 for Te=',Te
    !end if
    !
    Uaft=max(Umin,Ulast)
    Uaft=min(Uaft,Umax)
    Naft=log(Umax/Uaft)/lgdu
    Naft=max(Naft,0)
    !
    Ubef=min(Umax,Ufirst)
    Ubef=max(Ubef,Umin)
    Nbef=log(Ubef/Umin)/lgdu
    Nbef=max(Nbef,0)
    !
    n1=log(Umin/Ufirst)/lgdu
    n1=n1+1
    n2=log(Umax/Ufirst)/lgdu
    n2=n2+1
    n1=max(n1,0)
    n2=max(n1,min(n2,nbOut))
    ! 
    if(n1.gt.1) SPout(1:n1-1)=SPin(1)		!! 	<<<<<<<  
    SPout(max(n1,1):n2)=0			!! 	<<<<<<< 
    if(nBef.eq.0 .and. nAft.eq.0) then
       n1=ctrb_to1(n1+1)
       n3=ctrb_to2(n2)				! 111004
       if(n3.eq.0) n3=ctrb_to2(max(1,n2-1))		! 111004
       n2=n3						! 111004
    elseif(nBef.eq.0)then
       n1=ctrb_to1(n1+1)
       n1=max(n1,1)		! 110827
       n2=nbContrib
    elseif(nAft.eq.0) then
       n1=1
       n2=ctrb_to2(n2)
    else
       n1=1
       n2=nbContrib
    end if

    if(nBef.gt.1) then			!! 	<<<<<<
       SPout(1:nBef)=1			!! 	<<<<<<
       Uout(1:nBef)=Eout(1:nBef)/te 		!! 	<<<<<<
    endif					!! 	<<<<<<
 
    ! 
    nOut=0
    if(nBef.ne.0) then
       u=Ubef
       r=SPin(1)
       if(nBef.gt.mxout) goto 1101
       do n=nBef,1,-1
          u=u/rdu
          SPout(n)=r
          Uout(n)=u
       end do
       nOut=nOut+nBef
       if(nout.gt.mxout) goto 1100
    else				!! <<<<<<
       u=Umin				!! <<<<<<
       Uout(1)=u			!! <<<<<<
       r=SPin(1)			!! <<<<<<
       SPout(1)=r			!! <<<<<<
    end if
    ! 
    do n=nOut+1,min(mxOut,nbOut)
       SPout(n)=zero
    end do
    u=Ubef
    Uout(nOut+1)=u
    if(n2.ge.n1) then
       tt1=n1
       tt2=to_ctrb(n1)
       do ctr=tt1,n2
          to=to_ctrb(ctr)
          u=Eout(to)/te
          if(u.gt.Umin .or. to.ne.tt2) goto 201
          n1=ctr+1
       end do
201    tt1=to_ctrb(n2)
       tt2=n2
       do ctr=tt2,n1,-1
          to=to_ctrb(ctr)
          u=Eout(to)/te
          if(u.gt.Umin .or. to.ne.tt1) goto 202
          n2=ctr
       end do
       ! 
202    if(nBef.eq.0) then
          to=to_ctrb(n1)
          Uout(1)=Eout(to)/te
       end if
       tt1=nOut+1-to_ctrb(n1)
       do ctr=n1,n2
          fr=fr_ctrb(ctr)
          to=to_ctrb(ctr)
          u=Eout(to+1)/te
          to=to+tt1
          if(to.ge.size(Uout)) then
             write(*,*)'to,size(Uout)=',to,size(Uout)
             call CON_stop('to > size(Uout)')
          end if
          Uout(to+1)=u
          c=coef_ctrb(ctr)
          SPout(to)=SPout(to)+c*SPin(fr)
       end do
       nOut=to
    end if
    ! 

    !! ?? que doit-on faire quand  nOut=0 ???

    if(nAft.ne.0) then
       r=zero							!!  <<<<<<
       c=one							!!  <<<<<<
       if(nOut.gt.0) r=Spout(nOut)				!!  <<<<<<
       if(nOut.gt.1) c=min(one,Spout(nOut+1)/Spout(nOut))	!!  <<<<<<
       do n=1,nAft
          if(nOut.ge.nbOut) exit		! 110806
          u=u*rdu
          nOut=nOut+1
!n          SPout(nOut)=zero
	  SPout(nOut)=r						!!  <<<<<<
	  r=r*c							!!  <<<<<<
          Uout(nOut+1)=u
       end do
    end if
    ! 
    return
    ! 
110 nOut=min(mxOut,int(nbUout))
    u=Umin
    do n=1,nOut
       SPout(n)=r
       Uout(n)=u
       u=u*rdu
    end do
    Uout(nOut+1)=u
    ! 
    return
1100 write(*,*)'-R- nOut=',nOut,'  > mxOut=',mxOUt
    call CON_stop('- projSP_E, error=1100 -')
1101 write(*,*)'-R- nbef=',nbef,'  > mxOut=',mxOUt
    call CON_stop('- projSP_E, error=1101 -')
  end subroutine projSP_E  
  ! \
end MODULE CRASH_M_projE
! /
!------- radiom.f90

