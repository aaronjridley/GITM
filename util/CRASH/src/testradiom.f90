!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
 !c
 !c  1000 tables de 3000 point
 !c	a l'ancienne  : 36.4 s   et l'exp. cablée
 !c	a la nouvelle : 38.1 s  et l'exp. cablée
 !c  soit 12 µs par points.
 !c                        9.8 s  avec l'exp/ tabulée.
 !c  soit  3 µs par points.
 !c
	program testMaps
	use M_expTab,only: exp_tab8
	use M_RADIOM	! use M_projE
	implicit none
      include 'psize.h90'
      include 'cntrl.h90'
      include 'garcon.h90'
      include 'radiat.h90'
         real*8 sum(mxOut)
	REAL*8 xne,tablim1,gam
	integer i1,nfgp 
	common/uzb/ 	 & ! um,dum,pop,popE, &
		xne,tablim1,gam,sum,i1,nfgp	
	integer nNe,nTe,mx_N,mx_T
	integer,parameter :: mxNe=170,mxTe=100,mxNe1=60,mxTe1=50
	real*8 NE_1,NE_2,NE_1s(0:mxNe), T_1,T_2,T_1s(0:mxTE)
	parameter (NE_1=1e18,NE_2=1e24, T_1=1e0,T_2=1e5)
	real(8),parameter :: NE_10=1e10,NE_20=1e27,T_10=0.01,T_20=1e8
	real*8 NE_2s(0:mxNe,0:mxTE),TE_2s(0:mxNe,0:mxTE) &
		,TZ_2s(0:mxNe,0:mxTE)  ,r,rn,rt
	real*8 UBAR_2s(0:mxNe,0:mxTE),FRAD_2s(0:mxNe,0:mxTE)
	real*8 BETA_2s(0:mxNe,0:mxTE),ALPHA_2s(0:mxNe,0:mxTE)
		real*8 ubar1
		common /cubar/ubar1
	
	real*8 tdeb, tfin

	real*8 tzl,TE8
	integer nn,is,k,kk  ,doCaltz,doCalte
	real*8 e1,e2,teold,tenew,TZin,TEmid,xs , ts,at32s, IoB(mxOut)
	real*8 te,TZ,ne ,betapm,bu3,frad
	real*8 aSaha,bCoR,zero,one
	parameter (aSaha = 6.02E21, bCoR = 1.34E13*0.2 ,zero=0,one=1)
	character fname*80
		integer pass,npass
		parameter (npass=250)
	integer,parameter :: mxUB=45
	real*8,parameter :: r70=-70d0,r1e3=1e3
	real*8 Uin(301),SPin(301)
	real*8 gr45(46),dgr45(0:45),ones(301)
	data gr45/ &
		 1., 10., 48.75, 50.5, 58.5, 67.5 &
		,68.5, 79.25, 91., 98.5, 105., 127.5 &
		,150., 167.33, 184.67, 202., 219.33, 236.67 &
		,254., 267., 287., 300., 308., 309.5 &
		,320., 335., 345., 361., 364., 369. &
		,376., 378., 403.5, 435.5, 436.5, 446.2 &
		,474., 480., 540., 600., 757., 950. &
		,1250., 2500., 12000., 100000./
	data ones /301*1.d0/
	nn=0
 !c
	call exp_tab8()
	noconw=1
	indx(1)=1
	tablim=2.
	
	projMode=.true.
5	write(0,*) &
		'Te for projection check (0 to skip) ' &
		,',  E/B (<0 -> 1/u),Auger_1oB,Auger_B' &
		,',Auger_A,Auger_1oA' &
		,',Gaunt ?'
	read(*,*,err=5) teold,xs ,Auger_1oB,Auger_B &
		,Auger_A,Auger_1oA &
		,Gaunt
	if(Auger_1oB.le.0) Auger_1oB=20 ! 100
	if(Auger_B.le.0) Auger_B=0
	if(Auger_A  .le.0) Auger_A  =0
	if(Auger_1oA.le.0) Auger_1oA=0
	if(Gaunt.le.0.) Gaunt=0.2
	if(teold.le.0) goto 3
	nfg=45
			print*,'--- xs= ',xs,'  -----'
	if(xs.lt.0) then
	 do is=1,nfg+1
	  Uin(is)=gr45(is)/teold
	  SPin(is)=one/Uin(is)
	 enddo
	else
	 do is=1,nfg+1
	  Uin(is)=gr45(is)/teold
	  SPin(is)=xs
	 enddo
	endif
	nbIn=nfg
	 r=(Umax/Umin)**(one/ mxUB )
	 Uout(1)=Umin
	 do k=2,mxUB+1
	  Uout(k)=Uout(k-1)*r
	 enddo
		print*,'1# call projSP,  nbIn,nbOut=',nbIn,nbOut
	call projSP(Uin,SPin,nbIn)	! ,Uout,SPout,nbOut
		print*,'-> nbIn,nbOut=',nbIn,nbOut
		print*,'1# call prep_projE'
	 dgr45(0:45)=gr45(1:46)
	 nbIn=45
	 	print*,'gr45=',gr45(1),' : ',gr45(46)
	 	print*,'dgr45=',dgr45(0),' : ',dgr45(45)
	 	print*,'144.clls prep_projE, nbIn,nbOut=',nbIn,nbOut
	 call cpu_time(tdeb)
 !c	      ----------
	 call prep_projE(dgr45,nbIn)	! ,nbOut
 !c	      ----------
	 call cpu_time(tfin)
	 print*,'time for prep_projE=',tfin-tdeb,'  nbI,nbO=',nbIn,nbOut
	 TE8=teOld
	goto 5
3	 nbOut=0
	doCaltz=3
	write(0,*) &
	  'nb.group(0=stop), e1, e2, k(0/1/2/3=log/lin/45gr/file) '  !&
	read(*,*,err=3) nfg,e1,e2,k	! ,doCaltz
	if(nfg.le.0) stop
	if(k.eq.0) then
	 if(e2.le.e1) goto 3
	 e2=(e2/e1)**(1./nfg)
	 do is=1,nfg+1
	  eg(is-1)=e1
	  e1=e1*e2
	 enddo
		print*,'geom. prog.'
	elseif(k.eq.1) then
	 if(e2.le.e1) goto 3
	 e2=(e2-e1)/nfg
	 do is=1,nfg+1
	  eg(is-1)=e1
	  e1=e1+e2
	 enddo
		print*,'linear prog.'
	elseif(k.eq.2 .or. k.eq.45) then
	 nfg=45
	 do is=1,nfg+1
	  eg(is-1)=gr45(is)
	 enddo
	else
	 write(0,*) 'fileName ?'
	 read(*,'(a)') fname
	 open(1,file=fname,status='OLD',form='FORMATTED')
	 read(1,*) nfg
	 	print*,'- - found ',nfg,' groups'
	 read(1,*) (eg(is-1),is=1,nfg+1)
	 close(1)
	endif
	print*,'-> eg(1),eg(0,',nfg,')=',eg(0),eg(nfg)
 !c----
 !c----
	if(doNew) then
	 nbIn=nfg
	 call cpu_time(tdeb)
	 	print*,'210.clls prep_projE, nbIn,nbOut=',nbIn,nbOut
 !c	     ----------
	 call prep_projE(eg,nbIn)	! ,nbOut
 !c	     ----------
	 call cpu_time(tfin)
		print*,'prep_projE -> nbOut=',nbOut
	 print*,'dt=',tfin-tdeb,'  nbI,nbO=',nbIn,nbOut
	endif
 !c----
 !c----
1	write(0,*) &
		'Ne(0 for table), TE, Tz_old , I/B (or -Trad)' &
		,',mk-jhg-sumMod-proj-proj+jhg (0/1/2/3/4, + 8 : calte)  ?'
	read(*,*,err=1) deng(1),teold,tzlc(1),xs,k 
	doCalte=k/8
	k=mod(k,8)
	jhg=k.eq.1.or.k.ge.4
	sumMode=k.eq.2
	projMode=k.ge.3
	doNew=doCaltz.ge.2			! 080919
	if(doCalte.ne.0) doCaltz=0
		print*,'testradiom calls prepCorrUbar ..'
		write(0,*)'testradiom calls prepCorrUbar ..'
	call prepCorrUbar()
	TZin=TEold
	TEmid=TZin
			print*,'--- XS= ',xs,'  -----'
	if(xs.lt.0) then
         do is = 1,nfg
          r=0.5*(eg(is-1)+eg(is))/TEold
          if(r.gt.50.) then
           r=r*(-TEold/xs-1)
           r=max(r70,r)
           IoB(is)=exp(-r)
           spint(1,is)=IoB(is)			!!!!!!!
          else
           if(r.gt.1e-3) r=exp(r)-1
           spint(1,is)= r
           r=0.5*(eg(is-1)+eg(is))/(-xs)
           if(r.gt.1e-3) r=exp(r)-1
           spint(1,is)= spint(1,is)/r
           IoB(is)=spint(1,is)			!!!!!!!
          endif
         enddo    !   spint do loop
	else
         do is = 1,nfg
          spint(1,is)= xs
         enddo    !   spint do loop
	endif
	if(TZin.le.0 .or. deng(1).le.0) then
	 if(TZin.lt.0) then
	  mx_N=mxNe
	  mx_T=mxTe
	  rn=(NE_20/NE_10)**(1.d0/mx_N)
	  rt=(T_20/T_10)**(1.d0/mx_T)
	  do nNe=0,mx_N
	   NE_1s(nNe)=NE_10*rn**nNe
	  enddo
	  do nTe=0,mx_T
	   T_1s(nTe)=T_10*rt**nTe
	  enddo
	 else
	  mx_N=mxNe1
	  mx_T=mxTe1
	  rn=(NE_2/NE_1)**(1.d0/mx_N)
	  rt=(T_2/T_1)**(1.d0/mx_T)
	  do nNe=0,mx_N
	   NE_1s(nNe)=NE_1*rn**nNe
	  enddo
	  do nTe=0,mx_T
	   T_1s(nTe)=T_1*rt**nTe
	  enddo
	 endif
	 print *,'Ne=',NE_1s(1),' : ',NE_1s(mxNe)
	 print *,'Te=',T_1s(1),' : ',T_1s(mxTe)
	endif

			print*,'doCaltz=',doCaltz
	if(doCaltz.gt.0) goto 30
 !c
	 if(TZin.le.0 .or. deng(1).le.0) then
	  CALL cpu_time(tdeb)
 !c	       --------
	  	do pass=1,npass
	  do nNe=0,mx_N
	   do nTe=0,mx_T
	    NE_2s(nNe,nTe)=NE_1s(nNe)
	    TE_2s(nNe,nTe)=T_1s(nTe)	    
	    TZ_2s(nNe,nTe)=T_1s(nTe)
	    deng(1)=NE_2s(nNe,nTe)
	    ne=deng(1)
	    teold=TE_2s(nNe,nTe)
	    tenew=teold
	    do k=1,30
	     tzlc(1)=TZ_2s(nNe,nTe)
	     if(xs.lt.0) then
              do is = 1,nfg
               r=0.5*(eg(is-1)+eg(is))/TEold
               if(r.gt.50.) then
                r=r*(-TEold/xs-1)
                r=max(r70,r)
                IoB(is)=exp(-r)
                spint(1,is)=IoB(is)			!!!!!!!
               else
                if(r.gt.1e-3) r=exp(r)-1
                spint(1,is)= r
                r=0.5*(eg(is-1)+eg(is))/(-xs)
                if(r.gt.1e-3) r=exp(r)-1
                spint(1,is)= spint(1,is)/r
                IoB(is)=spint(1,is)			!!!!!!!
               endif
              enddo    !   spint do loop
	     endif
	     call calte(teold,tzlc(1),deng(1),tenew,eg,spint,ones,nfg)
	     if(abs(teold-tenew)/(teold+tenew).lt.1e-3) goto 122
	     teold=tenew
	    enddo
	    tenew=TE_2s(nNe,nTe)	! not converged
		te=TE_2s(nNe,nTe)
		ne=NE_2s(nNe,nTe)
122	    TE_2s(nNe,nTe)=tenew
	    UBAR_2s(nNe,nTe)=ubar1
		te=TE_2s(nNe,nTe)
		ne=NE_2s(nNe,nTe)
		ts=sqrt(te)
		at32s= aSaha*te*ts/ne
		betapm= bCoR/Gaunt*ts*te**3/ne
		bu3=betapm*ubar1**3
                frad= at32s*exp(-ubar1)/(1+bu3)
	    FRAD_2s(nNe,nTe)=frad
	    BETA_2s(nNe,nTe)=betapm
	    ALPHA_2s(nNe,nTe)=at32s
	   enddo
	  enddo
	  	enddo ! pass
	  CALL cpu_time(tfin)
 !c	       --------
	  goto 100
	 endif
 !c
 !c point fixe sur TE avec appel a  CALTE
 !c
20	tenew=teold
	do k=1,30
	 kk=k
	if(xs.lt.0) then
         do is = 1,nfg
          r=0.5*(eg(is-1)+eg(is))/TEold
          if(r.gt.50.) then
           r=r*(-TEold/xs-1)
           r=max(r70,r)
           IoB(is)=exp(-r)
           spint(1,is)=IoB(is)			!!!!!!!
          else
           if(r.gt.1e-3) r=exp(r)-1
           spint(1,is)= r
           r=0.5*(eg(is-1)+eg(is))/(-xs)
           if(r.gt.1e-3) r=exp(r)-1
           spint(1,is)= spint(1,is)/r
           IoB(is)=spint(1,is)			!!!!!!!
          endif
         enddo    !   spint do loop
	endif
	 call calte(teold,tzlc(1),deng(1),tenew,eg,spint,ones,nfg)
	 	print*,'iter#',k,' teold,new=',TEold,TEnew
	 if(abs(teold-tenew)/(teold+tenew).lt.1e-3) goto 2
	  TEmid=teold
	  teold=tenew
	enddo
2	continue
	te=tenew
	ne=deng(1)
		ts=sqrt(te)
		at32s= aSaha*te*ts/ne
		betapm= bCoR/Gaunt*ts*te**3/ne
		bu3=betapm*ubar1**3
                frad= at32s*exp(-ubar1)/(1+bu3)
                print*,'->ubar,frad=',ubar1,frad
        do is = 1,nfg
         IoB(is)= xs
        enddo    !   spint do loop
	call CALTZ(Te,Ne, TZ, eg,IoB,ones,nfg )

	print *,'TZin=',TZin,' -> TEnew=',Tenew
	if(kk.eq.30) then
	 print *,'le point fixe ne converge pas'
	else
	 print * &
		 ,' nb.iter=',kk,'  ubar=',ubar1 &
		 ,'  frad=',frad
	endif
	print *,' TZ=',TZ
	goto 3

 !c
 !c  calcul  TZ  avec appel direct a  CALTZ, puis verif avec calte
 !c
30       do is = 1,nfg
          IoB(is)= xs
         enddo    !   spint do loop
	 if(TZin.le.0 .or. deng(1).le.0) then
	 	print *,'computing Tz ...'
	  CALL cpu_time(tdeb)
 !c	       --------
	  	do pass=1,npass
	  do nNe=0,mx_N
	   do nTe=0,mx_T
	    NE_2s(nNe,nTe)=NE_1s(nNe)
	    TE_2s(nNe,nTe)=T_1s(nTe)	    
	    TZ_2s(nNe,nTe)=T_1s(nTe)	
	    deng(1)=NE_2s(nNe,nTe)
	    ne=deng(1)
	    TE8=TE_2s(nNe,nTe)
	    teold=TE8
	    if(xs.lt.0) then
 !c  Bnu(Trad)/Bnu(Te)=(exp(hnu/Te)-1)/(exp(hnu/Trad)-1)
             do is = 1,nfg
	      r=0.5*(eg(is-1)+eg(is))/teold
              if(r.gt.50.) then
               r=r*(-TEold/xs-1)
               r=max(r70,r)
               IoB(is)=exp(-r)
               spint(1,is)=IoB(is)			!!!!!!!
              else
               if(r.gt.1e-3) r=exp(r)-1
               IoB(is)= r
               r=0.5*(eg(is-1)+eg(is))/(-xs)	! =r*teold/(-xs)
               if(r.gt.1e-3) r=exp(r)-1
               IoB(is)= min(r1e3, IoB(is)/r)
               spint(1,is)=IoB(is)			!!!!!!!
              endif
             enddo    !   spint do loop
	    endif
	    call CALTZ(Te8,Ne, TZl, eg,IoB,ones,nfg )
	    	write(92,*)'ne,te,ubar,tz=',ne,te8,ubar1,tzl
102	    TZ_2s(nNe,nTe)=tzl
	    UBAR_2s(nNe,nTe)=ubar1
		te=TE_2s(nNe,nTe)
		ne=NE_2s(nNe,nTe)
		ts=sqrt(te)
		at32s= aSaha*te*ts/ne
		betapm= bCoR/Gaunt*ts*te**3/ne
		bu3=betapm*ubar1**3
                frad= at32s*exp(-ubar1)/(1+bu3)
	    FRAD_2s(nNe,nTe)=frad
	    BETA_2s(nNe,nTe)=betapm
	    ALPHA_2s(nNe,nTe)=at32s
	   enddo
	  enddo
	  	enddo ! pass
	  CALL cpu_time(tfin)
 !c	       --------
	  goto 100
	 endif
 !c
	ne=deng(1)
		print*,'deng,ne=',deng(1),ne
	TE8=TEold
	tzlc(1)=0
        do is = 1,nfg
         IoB(is)= xs
        enddo    !   spint do loop
	if(xs.lt.0) then
         do is = 1,nfg
          r=0.5*(eg(is-1)+eg(is))/TEold
          if(r.gt.50.) then
           r=r*(-TEold/xs-1)
           r=max(r70,r)
           IoB(is)=exp(-r)
           spint(1,is)=IoB(is)			!!!!!!!
          else
           if(r.gt.1e-3) r=exp(r)-1
           IoB(is)= r
           r=0.5*(eg(is-1)+eg(is))/(-xs)
           if(r.gt.1e-3) r=exp(r)-1
           IoB(is)= min(r1e3, IoB(is)/r)
           spint(1,is)=IoB(is)			!!!!!!!
          endif
         enddo    !   spint do loop
	endif
        ne=deng(1)
	call CALTZ(Te8,Ne, tzl, eg,IoB,ones,nfg )
	    	print*,'caltz : ne,te,ubar,tz=',ne,te8,ubar1,tzl
	    	write(92,*)'caltz : ne,te,ubar,tz=',ne,te8,ubar1,tzl
	tzlc(1)=tzl
	     call calte(teold,tzlc(1),deng(1),tenew, eg,IoB,ones,nfg)
	goto 3
 !c
100	if(doCaltz.ne.0) then
	 write(90,110) '#Ne','#Tz','Ne[cm-3]','Te[eV]','Tz[eV]' &
		 ,'ubar','frad(ubar)','Tz/Te','BT72s','AT32s'
	 write(91,110) '#Ne','#Tz','Ne[cm-3]','Te[eV]','Tz[eV]' &
		 ,'ubar','frad(ubar)','Tz/Te','BT72s','AT32s'
	else
	 write(90,110) '#Ne','#Te','Ne[cm-3]','Te[eV]','Tz[eV]' &
		 ,'ubar','frad(ubar)','Tz/Te','BT72s','AT32s'
	 write(91,110) '#Ne','#Te','Ne[cm-3]','Te[eV]','Tz[eV]' &
		 ,'ubar','frad(ubar)','Tz/Te','BT72s','AT32s'
	endif
		print*,'time used for the whole table=',tfin-tdeb
		print*,'create files fort.90 , fort.91'
110	format(10(1x,a))
111	format(2i5,1p,10e13.5)
	do nNe=0,mx_N,10
	 write(90,*)
	 do nTe=0,mx_T
	  write(90,111) nNe,nTe,NE_2s(nNe,nTe),TE_2s(nNe,nTe) &
		  ,TZ_2s(nNe,nTe),UBAR_2s(nNe,nTe),FRAD_2s(nNe,nTe) &
		  ,TZ_2s(nNe,nTe)/TE_2s(nNe,nTe), BETA_2s(nNe,nTe) &
		  ,ALPHA_2s(nNe,nTe)
	 enddo
	enddo
	do nTe=0,mx_T
	 write(91,*)
	 do nNe=0,mx_N
	  write(91,111) nNe,nTe,NE_2s(nNe,nTe),TE_2s(nNe,nTe) &
		  ,TZ_2s(nNe,nTe),UBAR_2s(nNe,nTe),FRAD_2s(nNe,nTe) &
		  ,TZ_2s(nNe,nTe)/TE_2s(nNe,nTe), BETA_2s(nNe,nTe) &
		  ,ALPHA_2s(nNe,nTe)
	 enddo
	enddo
	close(90)
	close(91)
	stop
	end
 !c
