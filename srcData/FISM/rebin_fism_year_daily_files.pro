
for year = 2016, 2018 do begin

print, year

maindir = '/raid3/Data/FISM/'

dir = maindir+'DailyFiles/'+tostr(year,4)
fismfilelist = findfile(dir+'/FISM_daily_*.sav')

nFiles = n_elements(fismfilelist)

if (nFiles gt 1) then begin

   newfile = maindir+'BinnedFiles/'+tostr(year)+'/fismflux_daily_'+tostr(year)+'.dat'

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
    
   print,wavelow,wavehigh

   close,30
   close,31

   ;convert to nm
   wavelow=wavelow/10.
   wavehigh=wavehigh/10.

   newflux    = fltarr(nFiles,nbins)
   newfluxavg = fltarr(nFiles,nbins)
   wavewvg    = fltarr(nFiles,nbins)

   openw,5,newfile
   printf, 5,'#START'

   for iFile = 0, nFiles-1 do begin

      file = FismFileList(iFile)
      restore, file

      if (iFile eq 0) then nbinsold   = n_elements(fism_pred)


;for il = 0, nFiles - 1 do begin
;   fism_pred(il,*) = fism_pred(il,*)*1.24
;endfor


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
            newflux(iFile,ibin) = fism_pred(idiff)
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
               newflux(iFile,ibin) = total(fism_pred(ilow:ihigh))
                
               x = ((fism_wv(ilow)-dw2) - wavelow(ibin))/(dw)
               y = (wavehigh(ibin)-(fism_wv(ihigh)+dw2))/(dw)
                
               newflux(iFile,ibin) = newflux(iFile,ibin) + $
                                     (x * fism_pred(ilow-1)) + $
                                     (y * fism_pred(ihigh+1))
                
               subline = where(fism_lines ge wavelow(ibin) and $
                               fism_lines lt wavehigh(ibin))
               if subline(0) gt -1 then begin
                  for isubline = 0, n_elements(subline) - 1 do begin
                     stop
                     newflux(iFile,ibin) = newflux(iFile,ibin) - $
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
                  newflux(iFile,ibin) = x*fism_pred(ilow) + $
                                        y * fism_pred(ilow+1)
               endif else begin
                  y = ((fism_wv(ihigh)+dw2) - wavehigh(ibin))/(dw)
                  newflux(iFile,ibin) = (x-y) * fism_pred(ilow)
               endelse
               if subline(0) gt -1 then begin
                  for isubline = 0, n_elements(subline) - 1 do begin
                     stop
                     newflux(iFile,ibin) = newflux(iFile,ibin) - $
                                           fism_pred(fism_lines(isubline))
                  endfor
               endif
            endelse
            
         endelse  
      ;  f2(iFile,ibin) = newflux(iFile,ibin)
      ;  newflux(iFile,ibin) = newflux(iFile,ibin) / fluxcount
      endfor

      doy = long(day_ar) mod 1000

      if (iFile eq 0) then begin
         it = [year,1,doy,00,0,0]
         c_a_to_r, it, rt
         c_r_to_a, it, rt
         printf, 5, it,transpose(newflux(iFile,*)), $
                 format = '(i5,5i3,'+tostr(nBins)+'e12.5)'
      endif
      it = [year,1,doy,12,0,0]
      c_a_to_r, it, rt
      c_r_to_a, it, rt
      printf, 5, it,transpose(newflux(iFile,*)), $
              format = '(i5,5i3,'+tostr(nBins)+'e12.5)'

   endfor

   it = [year,1,doy+1,0,0,0]
   c_a_to_r, it, rt
   c_r_to_a, it, rt
   printf, 5, it,transpose(newflux(iFile-1,*)), $
           format = '(i5,5i3,'+tostr(nBins)+'e12.5)'
   
   close,5

endif

endfor

;openw,1,'fism.txt'
;for i=0,58 do begin
;printf,1,(wavelow(i)+wavehigh(i))/2.,newflux(0,i)
;endfor
;close,1



end

