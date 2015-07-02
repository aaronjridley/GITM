;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
pro harmonics

; Purpose:
;  Convert raw fits file magnetogram into spherical 
;  harmonics ascii file to be used as an input to SWMF 
;
; Inputs:
; 1) Input magnetogram fits file.
; 2) Observatory name.
; 3) Order of harmonics.
; 4) Carrington rotation number.
; 5) what grid to use (Latitude or Sin Latitude).
;
; Outputs:
; 1) File.H - file containing the header of the input fits file.
; 2) File_tec.dat - Tecplot file containing the original magnetogram
;                   map.
; 3) File.dat - file containing the spherical harmonics coefficients.
; This Module reads a raw (RADIAL, not LOS!!!) magnetogram data file and
; generates a magnetogram file in the form of spherical
; harmonics to be use by SWMF. 
  
; ************************ Data Links ********************************
; * MDI:   http://soi.stanford.edu/magnetic/index6.html              *
; * WSO:   http://wso.stanford.edu/forms/prgs.html                   *  
; * GONG:  http://gong.nso.edu/data/magmap/QR/mqs/                   *
; * SOLIS: ftp://solis.nso.edu/synoptic/level3/vsm/merged/carr-rot   *
; * MWO:   ftp://howard.astro.ucla.edu/pub/obs/synoptic_charts       *
; ********************************************************************
; * Field in Gauss: MDI,GONG,SOLIS                                   *
; * Field in microTesla(0.01Gs): WSO, MWO                            *
; ********************************************************************

nMax=180
CR=2029
FileIn='                                          '
Obs='                                          '
UseSineLatGrid=1

; Required input parameters
read,FileIn,prompt='magnetogram fits file:'
FileIn=strtrim(FileIn,2)
read,Obs,prompt='enter observatory name:' 
Obs=strtrim(Obs,2)
read,nHarmonics,prompt='enter order of harmonics (nMax):' 
nHarmonics=strtrim(nHarmonics,2)
read,CR,prompt='enter Carrington Rotation number:' 
CR=strtrim(CR,2)
read,UseSineLatGrid,prompt='Use sine Latitude grid? (yes=1, no=0):' 
UseSineLatGrid=strtrim(UseSineLatGrid,2)

; Defining file names
FileHeader='CR'+CR+'_'+Obs+'.H'
FileOut='CR'+CR+'_'+Obs+'.dat'
FileTec='CR'+CR+'_'+Obs+'_tec.dat'
DataName='Br [G]'

;#######################################################
; Reading data 
;#######################################################
print,''
print,'Reading magnetogram file '
print,''

Data = readfits(FileIn, ImHeader, silent=silent)

print,'Writing header file: ',FileHeader
print,''

openw,lun,FileHeader,/get_lun
printf,lun,ImHeader
free_lun, lun

; Get image dimensions
s=size(Data)
Nx=s(1)
Ny=s(2)

print,'Magnetogram size= '
print,strtrim(Nx,2),'x',strtrim(Ny,2)
print,''

; Removing missing data by multiply B by sin(lat)^8
for i=0L,Ny-1 do begin
   theta=!PI*float(i)/float(Ny)
   for j=0L,Nx-1 do begin
      if(abs(Data(i*Nx+j)) ge 5000.0) then $
         Data(i*Nx+j)=Data(i*Nx+j)*sin(theta)^8
   endfor
endfor

; Writing tecpot file of original magnetogram
print,'Writing TecPlot file: ',FileTec
print,''

openw,lun,FileTec,/get_lun
printf,lun,' TITLE="',FileIn,'"'
printf,lun,'VARIABLES = "',DataName,'"'
printf,lun,'ZONE T="',FileTec,'", I= ',Nx,' J= ',Ny,' , K=1, F=POINT'

for i=0L,Ny-1 do $
   for j=0L,Nx-1 do $
      printf,lun, format = '(1e14.6)',Data(j,i)

free_lun, lun

print,'reading data completed!!!'
print,''

;#######################################################
; Defining arrays size and constants after 
; retriving magnetogram size
;#######################################################

nPhi=Nx
nTheta=Ny
dPhi=2.*!Pi/nPhi
dTheta=!Pi/nTheta
dSinTheta=2.0/nTheta
p_nm=fltarr(nHarmonics+1,nHarmonics+1)
g_nm=fltarr(nHarmonics+1,nHarmonics+1)
h_nm=fltarr(nHarmonics+1,nHarmonics+1)
PNMTheta_III=fltarr(nHarmonics+1,nHarmonics+1,nTheta)
CosMPhi_II=fltarr(nPhi,nHarmonics+1)
SinMPhi_II=fltarr(nPhi,nHarmonics+1)
factorRatio=fltarr(nHarmonics+1)
temp=fltarr(nPhi)
factorRatio(*)=0.0
MaxInt=100001
Sqrt_I=fltarr(MaxInt)
factorRatio=fltarr(nHarmonics+1)
Theta=0.0
CosTheta=0.0
SinTheta=0.0
delta_m0=0
SinThetaM  = 1.0
SinThetaM1 = 1.0
stuff1=0.0
stuff2=0.0
stuff3=0.0
NormalizationFactor=0.0

p_nm(*,*)=0.0
g_nm(*,*)=0.0
h_nm(*,*)=0.0

;#######################################################
; Calculating Associated Legandre polynoms 
;#######################################################
print,'Calculating Legandre polynoms'
print,''

; Calculate sqrt(integer) from 1 to 10000::
for m=0L,MaxInt-1 do $
   Sqrt_I(m) = sqrt(float(m))

; Calculate the ratio sqrt(2m!)/(2^m*m!)
factorRatio(0)=1.0
for m=1,nHarmonics do begin $   
   factorRatio(m)=factorRatio(m-1)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
endfor

for iTheta=0,nTheta-1 do begin
   if (UseSineLatGrid eq 1)then begin
      Theta=!pi*0.5-asin((float(iTheta)+0.5)*dSinTheta-1.0)
   endif else begin
      Theta=!Pi*0.5-((iTheta + 0.50)*dTheta - !Pi*0.50)
   endelse
   CosTheta=cos(Theta)
   if (sin(Theta) ge 1.e-9) then begin
      SinTheta=sin(Theta)
   endif else begin
      SinTheta=1.e-9
   endelse

   SinThetaM  = 1.0
   SinThetaM1 = 1.0

   p_nm(*,*)=0.0
   
   for m=0,nHarmonics do begin
      if (m eq 0) then begin
         delta_m0 = 1
      endif else begin
         delta_m0 = 0
      endelse

      p_nm(m,m)=factorRatio(m)*Sqrt_I((2-delta_m0)*(2*m+1))* $
                SinThetaM
      if (m lt nHarmonics) then $
         p_nm(m+1,m) = p_nm(m,m)*Sqrt_I(2*m+3)*CosTheta
      SinThetaM1 = SinThetaM
      SinThetaM  = SinThetaM*SinTheta
   endfor

   for m=0,nHarmonics-2 do begin
      for n=m+2,nHarmonics do begin
         stuff1         = Sqrt_I(2*n+1)/Sqrt_I(n^2-m^2)
         stuff2         = Sqrt_I(2*n-1)
         stuff3         = Sqrt_I((n-1)^2-m^2)/Sqrt_I(2*n-3)
         p_nm(n,m)  = $
            stuff1*(stuff2*CosTheta*p_nm(n-1,m)-stuff3*p_nm(n-2,m))
      endfor
   endfor

  ;Apply Schmidt normalization::
   for m=0,nHarmonics do begin 
      for n=m,nHarmonics do begin 
         stuff1 = 1.0/Sqrt_I(2*n+1)
         p_nm(n,m)  = p_nm(n,m)*stuff1
         PNMTheta_III(n,m,iTheta)=p_nm(n,m)
      endfor
   endfor
endfor

;#######################################################
; Calculating spherical harmonics coefficiets
;#######################################################

print,'Calculating harmonic coefficients'
print,''

for iPhi=0,nPhi-1 do begin
   for m=0,nHarmonics do begin
      CosMPhi_II(iPhi,m) = cos(float(m)*float(iPhi)*dPhi)
      SinMPhi_II(iPhi,m) = sin(float(m)*float(iPhi)*dPhi)
   endfor
endfor

for n=0L,nHarmonics do begin
   for m=0L,n do begin 
      ; Comment on normalization!!!!
      ;
      ; The analytic normalization factor in the solution of Laplace's eq. is (2n+1)/R_n, 
      ; where R_n=n+1+n(1/Rs)^(2n+1).
      ; However, in this code the coefficients are normalized only with 2n+1 to reproduce 
      ; the coefficients provided by Stanford. The division by R_n is done after
      ; the coefficients are been read in ModMagnetogram.f90.
      
      NormalizationFactor=2.*float(n)+1.0
      for iTheta=0,nTheta-1 do begin
         if (UseSineLatGrid eq 1)then begin
            Theta=!pi*0.5-asin((float(iTheta)+0.5)*dSinTheta-1.0)
         endif else begin
            Theta=!Pi*0.5-((iTheta + 0.50)*dTheta - !Pi*0.50)
         endelse
         CosTheta=cos(Theta)
         if (sin(Theta) ge 1.e-9) then begin
            SinTheta=sin(Theta)
         endif else begin
            SinTheta=1.e-9
         endelse 
                  
         if(UseSineLatGrid eq 1)then begin
            da = dSinTheta * dPhi
         endif else begin
            da=SinTheta*dTheta*dPhi
         endelse
                 
         temp(*)=Data(*,iTheta)*CosMPhi_II(*,m)
         g_nm(n,m)=g_nm(n,m)+ $
                   total(temp)*da*PNMTheta_III(n,m,iTheta)
         temp(*)=Data(*,iTheta)*SinMPhi_II(*,m)
         h_nm(n,m)=h_nm(n,m)+ $
                   total(temp)*da*PNMTheta_III(n,m,iTheta)
      endfor 

      g_nm(n,m)=NormalizationFactor*g_nm(n,m)/(4.*!pi)
      h_nm(n,m)=NormalizationFactor*h_nm(n,m)/(4.*!pi)
      
   endfor 
   print,'finished ',strtrim(fix(float(m)*100./float(nHarmonics)),2),'%'
endfor

g_nm(0,0) = 0.0

;#######################################################
; Writing spherical harmonics file
;#######################################################

print,'writing spherical harmonics file: ',FileOut
print,''

openw,lun,FileOut,/get_lun
printf,lun, format = '(a19,I3,a10,I4,a4)','Coefficients order=',nHarmonics,' center=CT',CR,':180'
printf,lun, format = '(a)','Observation time'
printf,lun, format = '(a45,I3)','B0 angle & Nmax:        0          ',nHarmonics
printf,lun,''
printf,lun,''
printf,lun,''
printf,lun,''
printf,lun,''
printf,lun, format = '(a19,I3,a26)','Max Harmonic Order:',nHarmonics,' Units: Gauss'
printf,lun, format = '(a)',' '
printf,lun, format = '(a)','  l   m      g(G)      h(G)'
printf,lun, format = '(a)',' '

for n=0L,nHarmonics do $ 
   for m=0L,n do $
      printf,lun, format = '(2I3, 2f14.6)',n,m,g_nm(n,m),h_nm(n,m)

free_lun, lun

print, 'Spherical harmonics file is ready to use by SWMF!!!'


end

