PRO rebin_fism, cdate
;if n_elements(cdate) eq 0 then cdate = ''
;cdate = ask('date to rebin (yyyymmdd): ',cdate)

cyear = strmid(cdate,0,4)
cmon = strmid(cdate,4,2)
cday = strmid(cdate,6,2)

print, 'Working on '+cyear+' '+cmon+' '+cday+'...'
idate = [fix(cyear),fix(cmon),fix(cday),0,0,0]
c_a_to_r, idate,rt
cdoy = chopr('00'+tostr(julian_day(rt)),3)

fismfile = file_search('/bigdisk1/Data/FISM/'+cyear+'/FISM*60sec*_'+cyear+cdoy+'_*.sav')
fismfile = fismfile(0)

if strpos(fismfile,'.sav') ge 0 then begin
    restore,fismfile
    newfile = '/bigdisk1/Data/FISM/BinnedFiles/'+cyear+'/fismflux'+cyear+cmon+cday+'.dat'


;read in GITM wavelength information
    
    close,/all
    openr,30, 'wavelow'
    openr,31,'wavehigh'
    nbins = 59

;    nbins = 1
;    wavelow = 260
;    wavehigh = 340
    wavelow=fltarr(nbins)
    wavehigh=fltarr(nbins)
    readf,30,wavelow
    readf,31,wavehigh
    
close,30
close,31
   ;convert to nm
wavelow=wavelow/10.
wavehigh=wavehigh/10.

ntimes = n_elements(utc)
newflux = fltarr(ntimes,nbins)
f2 = newflux
newfluxavg = fltarr(ntimes,nbins)
wavewvg = fltarr(ntimes,nbins)

nbinsold = n_elements(fism_pred(0,*))
wavecount = fltarr(nbinsold)
wavecount(0) = 0.5

;for il = 0, ntimes - 1 do begin
;   fism_pred(il,*) = fism_pred(il,*)*1.24
;endfor

openw,5,newfile
printf, 5,'#START'
for itime = 0, ntimes - 1 do begin
    wavediff = fltarr(nbinsold)
    wavediff(0) = 10000
    wdlow = fltarr(nbinsold+1)
    wdhigh = fltarr(nbinsold+1)
    wdlow(0) = 10000
    wdhigh(0) = 10000

    line = where(wavelow eq wavehigh)
    fism_lines = fltarr(n_elements(line))
    iline = 0
   
    for ibin = 0, nbins - 1 do begin
        fluxcount = 0
            
        ;;;; This is when gitm needs a single wavelength
        if(wavelow(ibin) eq wavehigh(ibin)) then begin
            diff = min(abs(fism_wv - wavelow(ibin)),idiff)
            newflux(itime,ibin) = fism_pred(itime,idiff)
            fism_lines(iline) = fism_wv(idiff)
            iline = iline + 1
            fluxcount = 1
        endif else begin
            fluxcount = wavehigh(ibin) - wavelow(ibin)

            dw = mean(fism_wv(1:nbinsold-1) - fism_wv(0:nbinsold-2))
            dw2 = dw/2.
            ilow = max(where(wavelow(ibin) - (fism_wv-dw2) ge 0))
            ihigh = max(where(wavehigh(ibin) - (fism_wv+dw2) ge 0))
            if ihigh lt 0 then begin
                ihigh = ilow
            endif
            wdiffl = wavelow(ibin) - (fism_wv(ilow)-dw2)
            wdiffh = wavehigh(ibin) - (fism_wv(ihigh)+dw2)

            ;;;; this is when GITM bins are larger than fism bins
            if ihigh - ilow gt 0 then begin

                if wdiffl gt 0 then ilow = ilow + 1
                newflux(itime,ibin) = total(fism_pred(itime,ilow:ihigh))
                
                x = ((fism_wv(ilow)-dw2) - wavelow(ibin))/(dw)
                y = (wavehigh(ibin)-(fism_wv(ihigh)+dw2))/(dw)
                
                newflux(itime,ibin) = newflux(itime,ibin) + $
                  (x * fism_pred(itime,ilow-1)) + $
                  (y * fism_pred(itime,ihigh+1))
                
               subline = where(fism_lines ge wavelow(ibin) and fism_lines lt wavehigh(ibin))
               if subline(0) gt -1 then begin
                   for isubline = 0, n_elements(subline) - 1 do begin
                       newflux(itime,ibin) = newflux(itime,ibin) - $
                         fism_pred(fism_lines(isubline))
                   endfor

                endif
            endif else begin
               ;;;; These are when GITM bins are the same size or
               ;;;; smaller than fism bins

               if wavelow(ibin) ge fism_wv(ilow)-dw2 then begin
                   x = 1 - ((wavelow(ibin) - (fism_wv(ilow)-dw2))/(dw))
               endif else begin
                   x = ((fism_wv(ilow)-dw2) - wavelow(ibin))/(dw)
               endelse
               if wavehigh(ibin) ge fism_wv(ihigh) + dw2 then begin
                   y = (wavehigh(ibin) - (fism_wv(ihigh)+dw2))/(dw)
                   newflux(itime,ibin) = x*fism_pred(itime,ilow) + $
                     y * fism_pred(itime,ilow+1)
               endif else begin
                   y = ((fism_wv(ihigh)+dw2) - wavehigh(ibin))/(dw)
                   newflux(itime,ibin) = (x-y) * fism_pred(itime,ilow)
               endelse
               if subline(0) gt -1 then begin
                   for isubline = 0, n_elements(subline) - 1 do begin
                       newflux(itime,ibin) = newflux(itime,ibin) - $
                         fism_pred(fism_lines(isubline))
                   endfor
               endif
           endelse

        endelse  
      ;  f2(itime,ibin) = newflux(itime,ibin)
      ;  newflux(itime,ibin) = newflux(itime,ibin) / fluxcount
    endfor

hour = fix(fix(utc(itime)/3600.))
min = fix(fix(utc(itime)/60) - fix(hour)*60.)
sec = fix(utc(itime) - fix(hour)*3600. - fix(min)*60.)
printf, 5, fix(cyear),fix(cmon),fix(cday),hour,min,sec,transpose(newflux(itime,*))

endfor
endif

close,5

;openw,1,'fism.txt'
;for i=0,58 do begin
;printf,1,(wavelow(i)+wavehigh(i))/2.,newflux(0,i)
;endfor
;close,1



end

