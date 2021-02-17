
filelist = findfile('3DION*.bin')
nFiles = n_elements(filelist)

for iFile = 0, nFiles-1 do begin

   filename = filelist(iFile)
   print, filename

   read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                           vars, data, rb, cb, bl_cnt, iTime, Version

   if (iFile eq 0) then begin
      all_vertdrift  = fltarr(nLons-4, nFiles)
      all_localtimes = fltarr(nLons-4, nFiles)
   endif

   vertdrift = fltarr(nLons,nAlts)
   alts      = fltarr(nLons,nAlts)
   lons      = fltarr(nLons,nAlts)
   glats     = fltarr(nLons,nAlts)
   mlats     = fltarr(nLons,nAlts)

   for i=0,nLons-1 do begin

      ml = data(22,i,*,2)
      l = where(abs(ml) eq min(abs(ml)))
      l = l(0)
      vertdrift(i,*) = data(17,i,l,*)
      alts(i,*) = data(2,i,l,*)/1000.0
      lons(i,*) = data(0,i,l,*)/!dtor
      glats(i,*) = data(1,i,l,*)/!dtor
      mlats(i,*) = data(22,i,l,*)

   endfor

   p = strpos(filename,'bin')
   psfile = strmid(filename,0,p)+'vd.ps' 

   c_a_to_s, itime, stime

   setdevice, psfile, 'p', 5

   pos = [0.1, 0.55, 0.9, 0.95]

   maxi = max(abs(vertdrift))
   maxi = 100.0
   mini = -maxi

   levels = findgen(31)/30 * 2*maxi - maxi

   xtickname = ['0','90','180','270','360']
   xtickvals = float(xtickname)
   xticks = n_elements(xtickname)-1

   makect, 'mid'
   contour, vertdrift, lons, alts, xrange = [0,360], xstyle = 1, $
            xtitle = 'Geographic Longitude (deg)', $
            ytitle = 'Altitude (km)', $
            nlevels = 31, levels = levels, $
            /fill, pos = pos, $
            title = 'Vertical Drifts at Magnetic Equator at '+stime, $
            yrange = [0,alts(0,nAlts-3)], ystyle = 1, $
            xtickname = xtickname, xtickv = xtickvals, $
            xticks = xticks, xminor = 9

   levels = findgen(5)/4 * 2*maxi - maxi

   contour, vertdrift, lons, alts, xrange = [0,360], xstyle = 1, $
            nlevels = 5, levels = levels, /noerase, $
            color = 0, /follow, pos = pos, $
            yrange = [0,alts(0,nAlts-3)], ystyle = 1, $
            xtickname = xtickname, xtickv = xtickvals, $
            xticks = xticks

   ctpos = [pos(2)+0.01, pos(1), pos(2)+0.03, pos(3)]
   ncolors = 255
   maxmin = [mini,maxi]
   title = 'Vertical Drift (m/s)'
   plotct, ncolors, ctpos, maxmin, title, /right

   pos = [0.1, 0.2, 0.9, 0.50]

   l = where(alts(0,*) gt 350.0)
   l = l(0)

   plot, lons(*,l),vertdrift(*,l), yrange = [mini,maxi], pos = pos, $
         ytitle = 'Vertical Drift @ 350km (m/s)', /noerase, $
         ystyle = 1, xstyle = 1, xrange = [0,360], $
         xtickname = xtickname, xtickv = xtickvals, $
         xticks = xticks, xminor = 9

   oplot, lons(*,l), lons(*,l)*0.0, linestyle = 2

   ut = float(itime(3)) + float(itime(4))/60.0
   lo = reform(lons(*,l))
   localtime = (reform(lons(*,l))/15 + ut + 24.0) mod 24.0

   all_localtimes(*,iFile) = localtime(2:nLons-3)
   all_vertdrift(*,iFile)  = vertdrift(2:nLons-3,l)

   for i=0,xticks do begin
      l = where(lo gt xtickvals(i))
      l = l(0)
      lts = tostr(fix(localtime(l)))+'LT'
      xyouts, lo(l), mini-0.25*maxi, lts, align=0.5
   endfor

   closedevice

endfor

;display,vars

setdevice, 'all.ps', 'p', 5

vd = all_vertdrift(*,nFiles-1)
lo = all_localtimes(*,nFiles-1)
sort_a, lo, vd

plot, lo, vd, yrange = [-50,50], xrange = [0,24], xstyle = 1, $
      pos = [0.1, 0.3, 0.95, 0.7], $
      xtitle = 'Localtime (hours)', ytitle = 'Vertical Drift (m/s)', $
      thick = 5, color = 0, ystyle = 1

oplot, [0.0,24.0], [0.0,0.0], linestyle = 2

for iFile = 0,nFiles-2 do begin

   vd = all_vertdrift(*,iFile)
   lo = all_localtimes(*,iFile)
   sort_a, lo, vd

   c = 255.0 - 255.0*float(iFile)/(nFiles-1)

   oplot, lo, vd, thick = float(iFile)/(nFiles-1)*3.0, color = c

endfor

vd = all_vertdrift(*,nFiles-1)
lo = all_localtimes(*,nFiles-1)
sort_a, lo, vd

oplot, lo, vd, thick = 5


closedevice

end

