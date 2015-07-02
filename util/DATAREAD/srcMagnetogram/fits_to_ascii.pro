;  Copyright (C) 2002 Regents of the University of Michigan, 
;  portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
pro fits_to_ascii, FileIn, DataName, silent=silent

; Purpose:
;  Read fits magnetogram file and write out an ASCII file.
;
; Usage:
;   fits_to_ascii [,FileIn] [, FileOut] [,/silent]
;
; FileIn  - name of the fits file. Default is fitsfile.fits
; FileOut - first part of the names of the output files. Default is fitsfile
;           so the files will by fitstfile.H, fitsfile.dat 

; /silent - suppress verbose information.

if n_elements(FileIn)  eq 0 then FileIn  = 'fitsfile.fits'
if n_elements(FileOut) eq 0 then FileOut = 'fitsfile'

nMax=180

FileHeader= FileOut + '.H'
FileDat   = FileOut + '.dat'
FileTec   = FileOut + '_tec.dat'
FileIdl   = FileOut + '_idl.out' 
DataName  = 'Br [G]'  

Data = read_fits(FileIn, ImHeader, silent=silent)

if not keyword_set(silent) then begin
    print,''
    print,'Writing header file ',FileHeader
    print,''
endif

openw,lun,FileHeader,/get_lun
printf,lun,ImHeader
free_lun, lun

; Get image dimensions
s=size(Data)
nLon=s(1)
nLat=s(2)

;;; This makes no sense...
;;; ; Removing missing data by multiply B by sin(lat)^8
;;; for i=0L,nLat-1 do begin
;;;     theta = !PI*float(i)/float(nLat)
;;;     for j=0L,nLon-1 do begin
;;;         if(abs(Data(i*nLon+j)) ge 5000.0) then $
;;;             Data(i*nLon+j) = Data(i*nLon+j)*sin(theta)^8
;;;     endfor
;;; endfor

if not keyword_set(silent) then begin
    print,''
    print,'Writing simple data file ',FileDat
    print,''
endif

openw,lun,FileDat,/get_lun
printf,lun,'#nMax'
printf,lun,nMax
printf,lun,'#ARRAYSIZE'
printf,lun,strtrim(nLon,2)
printf,lun,strtrim(nLat,2)
printf,lun,'#START'

for i=0L,nLat-1 do begin
    for j=0L,nLon-1 do begin
        printf,lun, format = '(1e14.6)',Data(j,i)
    endfor
endfor

free_lun, lun

if not keyword_set(silent) then begin
    print,''
    print,'Writing TecPlot file ',FileTec
    print,''
endif

openw, lun, FileTec, /get_lun
printf,lun,' TITLE="',FileIn,'"'
printf,lun,'VARIABLES = "',DataName,'"'
printf,lun,'ZONE T="',FileTec,'", I= ',nLon,' J= ',nLat,' , K=1, F=POINT'

for i=0L,nLat-1 do for j=0L,nLon-1 do $
  printf,lun, format = '(1e14.6)',Data(j,i)

free_lun, lun

if not keyword_set(silent) then begin
    print,''
    print,'Writing IDL file ',FileIdl
    print,''
endif

openw, lun, FileIdl, /get_lun
printf,lun,' Longitude [Deg], Latitude [Deg],',DataName
printf,lun, 0, 0.0, 2, 1, 1
printf,lun, nLon,' ',nLat
printf,lun, '0.5'
printf,lun,'Longitude Latitude Br LongitudeShift'

dLon = 360.0/nLon
for i=0L,nLat-1 do begin
   Latitude =  asin((2*i-nLat+1.0)/nLat)/!dtor
   for j=0L,nLon-1 do begin
      printf,lun,format ='(3e14.6)', j*dLon, Latitude, Data(j,i)
   endfor
endfor

free_lun,lun

if not keyword_set(silent) then print,'Conversion done'

end

