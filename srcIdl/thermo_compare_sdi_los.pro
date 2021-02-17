
gitmfiles = findfile('3D*.bin')
gitm_read_header, gitmfiles, GitmTimes, nVars, Vars, $
                  nLons, nLats, nAlts, version

gitm_read_bin, gitmfiles(0), GitmData, firsttime, nVars, Vars, version

GitmLongitudes = reform(GitmData(0,*,0,0))/!dtor
GitmLatitudes  = reform(GitmData(1,0,*,0))/!dtor
GitmAltitudes  = reform(GitmData(2,0,0,*))/1000.0

l = where(GitmAltitudes gt 240.0)
iAlt = l(0)

iTemp   = 15
iEast   = 16
iNorth  = 17
iUp     = 18
iiEast  = 36
iiNorth = 37

dLons = GitmLongitudes(1) - GitmLongitudes(0)
dLats = GitmLatitudes(1)  - GitmLatitudes(0)

list = findfile('*.txt')
nfiles = n_elements(list)
display, list
if (nfiles eq 1) then iFile = 0 $
else begin
   if (n_elements(iFile) eq 0) then iFile = 0
   iFile = ask('sdi file to compare to (-1 for all)',tostr(iFile))
endelse

if (iFile gt -1) then begin
   iStart = fix(iFile)
   iEnd = fix(iFile)
endif else begin
   iStart = 0
   iEnd = nFiles-1
endelse

for iFile = iStart, iEnd do begin

   SdiFilename = list(iFile)

   sdi = sdi_ascii_reader(SdiFilename)

   if (strpos(sdi.header.site,'Mawson') gt -1) then site = 'maw'
   if (strpos(sdi.header.site,'Poker')  gt -1) then site = 'pkr'
   if (strpos(sdi.header.site,'Toolik') gt -1) then site = 'tlk'
   if (strpos(sdi.header.site,'HAARP')  gt -1) then site = 'hrp'

   SdiCenterLongitude = sdi.header.longitude
   SdiCenterLatitude  = sdi.header.latitude

   nPts = sdi.fov.viewing_directions

   SdiLongitudes = (sdi.fov.zone_longitudes+360.0) mod 360.0
   SdiLatitudes  = sdi.fov.zone_latitudes

   SdiZenith  = sdi.fov.zenith_angles
   SdiAzimuth = sdi.fov.azimuth_angles

   x = (SdiZenith) * cos(!pi/2 - SdiAzimuth*!dtor)
   y = (SdiZenith) * sin(!pi/2 - SdiAzimuth*!dtor)

   maxran = max(sqrt(x^2 + y^2))

   iLongitudes = fix((SdiLongitudes-min(gitmlongitudes))/dLons)
   iLatitudes  = fix((SdiLatitudes-min(gitmlatitudes))/dLats)

   rLongitudes = (SdiLongitudes-min(gitmlongitudes))/dLons - iLongitudes
   rLatitudes  = (SdiLatitudes-min(gitmlatitudes))/dLats  - iLatitudes

   nTimesSdi = n_elements(sdi.ALLSKY_TMP_INT)

   SdiTimes = dblarr(nTimesSdi)
   GitmIter = intarr(nTimesSdi)

   north    = fltarr(nTimesSdi,nPts)
   east     = fltarr(nTimesSdi,nPts)
   northI   = fltarr(nTimesSdi,nPts)
   eastI    = fltarr(nTimesSdi,nPts)
   GitmLos  = fltarr(nTimesSdi,nPts)
   GitmLosI = fltarr(nTimesSdi,nPts)
   GitmTemp = fltarr(nTimesSdi,nPts)
   GitmUp   = fltarr(nTimesSdi,nPts)

   gitmlosTotal  = fltarr(4,nTimesSdi)
   gitmtempTotal = fltarr(4,nTimesSdi)
   sdilosTotal   = fltarr(4,nTimesSdi)
   sditempTotal  = fltarr(4,nTimesSdi)

   sditempZ  = fltarr(nTimesSdi)

   gitmtempZ = fltarr(nTimesSdi)
   gitmUpZ   = fltarr(nTimesSdi)

   for iT = 0, nTimesSdi-1 do begin

      itime = [sdi.time.year(iT), $
               sdi.time.month(iT)+1, $
               sdi.time.day(iT), 0, 0, 0]
      c_a_to_r, itime, rtime

      ut = (sdi.time.ut_begin(iT) + sdi.time.ut_end(iT))/2.0 * 3600.0
      SdiTimes(iT) = rtime + ut

      c_r_to_a, itime, SdiTimes(iT)
      c_a_to_ymd, itime, ymd

      dmy = strmid(ymd,6,2)+strmid(ymd,4,2)+strmid(ymd,0,4)

      l = where(GitmTimes gt SdiTimes(iT), c)
      if (c gt 0) then GitmIter(iT) = l(0)-1 $
      else GitmIter(iT) = -1

      if (GitmIter(iT) eq -1) then GitmIter(iT) = 0

      if (iT eq 0 and GitmIter(iT) gt -1) then $
         gitm_read_bin_1var, GitmFiles(GitmIter(iT)), GitmData, $
                             firsttime, nVars, Vars, version, $
                             VarsToGet = [iTemp, iEast, iNorth, iUp, iiEast, iiNorth] $
      else if (GitmIter(iT) ne GitmIter(iT-1)) then $
         gitm_read_bin_1var, gitmfiles(GitmIter(iT)), GitmData, $
                             firsttime, nVars, Vars, $
                             version, VarsToGet = [iTemp, iEast, iNorth, iUp, iiEast, iiNorth]

      GitmTemp2d = reform(GitmData(0,*,*,iAlt))
      GitmEast   = reform(GitmData(1,*,*,iAlt))
      GitmNorth  = reform(GitmData(2,*,*,iAlt))
      GitmUp2d   = reform(GitmData(3,*,*,iAlt))
      GitmEastI  = reform(GitmData(4,*,*,iAlt))
      GitmNorthI = reform(GitmData(5,*,*,iAlt))

      for iP = 0, nPts-1 do begin

         gitmup(iT,iP) = $
            (1.0-rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmUp2d(iLongitudes(iP)  ,iLatitudes(iP)  ) + $
            (1.0-rLongitudes(iP)) * (    rLatitudes(iP)) * GitmUp2d(iLongitudes(iP)  ,iLatitudes(iP)+1) + $
            (    rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmUp2d(iLongitudes(iP)+1,iLatitudes(iP)  ) + $
            (    rLongitudes(iP)) * (    rLatitudes(iP)) * GitmUp2d(iLongitudes(iP)+1,iLatitudes(iP)+1)

         gitmtemp(iT,iP) = $
            (1.0-rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmTemp2d(iLongitudes(iP)  ,iLatitudes(iP)  ) + $
            (1.0-rLongitudes(iP)) * (    rLatitudes(iP)) * GitmTemp2d(iLongitudes(iP)  ,iLatitudes(iP)+1) + $
            (    rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmTemp2d(iLongitudes(iP)+1,iLatitudes(iP)  ) + $
            (    rLongitudes(iP)) * (    rLatitudes(iP)) * GitmTemp2d(iLongitudes(iP)+1,iLatitudes(iP)+1)

         north(iT,iP) = $
            (1.0-rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmNorth(iLongitudes(iP)  ,iLatitudes(iP)  ) + $
            (1.0-rLongitudes(iP)) * (    rLatitudes(iP)) * GitmNorth(iLongitudes(iP)  ,iLatitudes(iP)+1) + $
            (    rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmNorth(iLongitudes(iP)+1,iLatitudes(iP)  ) + $
            (    rLongitudes(iP)) * (    rLatitudes(iP)) * GitmNorth(iLongitudes(iP)+1,iLatitudes(iP)+1)

         east(iT,iP) = $
            (1.0-rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmEast(iLongitudes(iP)  ,iLatitudes(iP)  ) + $
            (1.0-rLongitudes(iP)) * (    rLatitudes(iP)) * GitmEast(iLongitudes(iP)  ,iLatitudes(iP)+1) + $
            (    rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmEast(iLongitudes(iP)+1,iLatitudes(iP)  ) + $
            (    rLongitudes(iP)) * (    rLatitudes(iP)) * GitmEast(iLongitudes(iP)+1,iLatitudes(iP)+1)

         northI(iT,iP) = $
            (1.0-rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmNorthI(iLongitudes(iP)  ,iLatitudes(iP)  ) + $
            (1.0-rLongitudes(iP)) * (    rLatitudes(iP)) * GitmNorthI(iLongitudes(iP)  ,iLatitudes(iP)+1) + $
            (    rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmNorthI(iLongitudes(iP)+1,iLatitudes(iP)  ) + $
            (    rLongitudes(iP)) * (    rLatitudes(iP)) * GitmNorthI(iLongitudes(iP)+1,iLatitudes(iP)+1)

         eastI(iT,iP) = $
            (1.0-rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmEastI(iLongitudes(iP)  ,iLatitudes(iP)  ) + $
            (1.0-rLongitudes(iP)) * (    rLatitudes(iP)) * GitmEastI(iLongitudes(iP)  ,iLatitudes(iP)+1) + $
            (    rLongitudes(iP)) * (1.0-rLatitudes(iP)) * GitmEastI(iLongitudes(iP)+1,iLatitudes(iP)  ) + $
            (    rLongitudes(iP)) * (    rLatitudes(iP)) * GitmEastI(iLongitudes(iP)+1,iLatitudes(iP)+1)

         gitmlos(iT,iP) = north(iT,iP) * cos(SdiAzimuth(iP)*!dtor) * sin(SdiZenith(iP)*!dtor) + $
                          east(iT,iP) *  sin(SdiAzimuth(iP)*!dtor) * sin(SdiZenith(iP)*!dtor)
         gitmlosI(iT,iP) = northI(iT,iP) * cos(SdiAzimuth(iP)*!dtor) * sin(SdiZenith(iP)*!dtor) + $
                           eastI(iT,iP) *  sin(SdiAzimuth(iP)*!dtor) * sin(SdiZenith(iP)*!dtor)

      endfor

      width = 22.5

      l = where(SdiAzimuth gt 360.0-width or SdiAzimuth lt width and $
                SdiZenith gt 30)
      gitmlosTotal(0,iT) = mean(gitmlos(iT,l))
      sdilosTotal(0,iT) = mean(sdi.los_wind_skymap(iT,l))
      gitmtempTotal(0,iT) = mean(gitmtemp(iT,l))
      sditempTotal(0,iT) = mean(sdi.temp_skymap(iT,l))

      for i=1,3 do begin
         l = where(SdiAzimuth gt 90.0*i-width and $
                   SdiAzimuth lt 90*i+width and $
                   SdiZenith gt 30)
         gitmlosTotal(i,iT) = mean(gitmlos(iT,l))
         sdilosTotal(i,iT) = mean(sdi.los_wind_skymap(iT,l))
         gitmtempTotal(i,iT) = mean(gitmtemp(iT,l))
         sditempTotal(i,iT) = mean(sdi.temp_skymap(iT,l))
      endfor

      l = where(SdiZenith eq 0)

      gitmtempZ(iT) = gitmtemp(iT,l(0))
      gitmUpZ(iT)   = gitmUp(iT,l(0))

      l = where(SdiZenith lt 30)
      sditempZ(iT) = mean(sdi.temp_skymap(iT,l))

   endfor

; ---------------------------------------------------------------------------

   UsePkz = 0
   filepkz = findfile('ResultsPKZ'+dmy+'.txt')
   filepkz = filepkz(0)
   if (strlen(filepkz) gt 0) then begin
      read_meriwether_fpi, filepkz, PkzTime, PkzVz, PkzVze, PkzTn, PkzTne, in
      pkztime = rtime + PkzTime - 24.0*60.0*60.0
      UsePkz = 1
   endif

   setdevice, 'compare_temp_'+site+'_'+ymd+'.ps', 'p', 5

   plotdumb

   ppp = 6
   space = 0.01
   pos_space, ppp, space, sizes, ny = ppp

   titles = ['North Temp (K)', 'East Temp (K)', 'South Temp (K)', 'West Temp (K)']

   stime = min(SdiTimes)
   etime = max(SdiTimes)

   time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

   maxi = max([gitmtempTotal, sditempTotal])
   mini = min([gitmtempTotal, sditempTotal])
   range = [0.9*mini, 1.1*maxi]

   for i = 0, 3 do begin

      get_position, ppp, space, sizes, i, pos, /rect
      pos(0) = pos(0)+0.2
      pos(2) = pos(2)-0.2

      if (i eq 3 and UsePkz eq 0) then begin
         xt = xtitle
         xtn = xtickname
      endif else begin
         xt = ' '
         xtn = strarr(10)+' '
      endelse

      plot, SdiTimes-stime, gitmtempTotal(i,*), $
            pos = pos, /noerase, $
            ytitle = titles(i),		$
            xtickname = xtn,			$
            xtitle = xt,			$
            xtickv = xtickv,			$
            xminor = xminor,			$
            xticks = xtickn,   $
            yrange = range, ystyle = 1, $
            thick = 4

      oplot, SdiTimes-stime, sditempTotal(i,*), linestyle = 2, thick = 4

      dx = (etr-btr)/20.0

      oplot, [dx, dx*3], [mini*0.95, mini*0.95], thick = 4
      xyouts, dx*3.1, mini*0.95, 'GITM'

      oplot, [dx*10, dx*10 + dx*2], [mini*0.95, mini*0.95], thick = 4, linestyle = 2
      xyouts, dx*12.1, mini*0.95, 'SDI ('+mkupper(site)+')

      meansdi = mean(sditempTotal(i,*))
      rms  = sqrt(mean((gitmtempTotal(i,*)-sditempTotal(i,*))^2))
      diff = mean(gitmtempTotal(i,*)-sditempTotal(i,*))
      xyouts, dx, maxi*1.03, $
              'RMS : '+tostr(rms)+' K; Diff : '+tostr(diff)+' K; mean SDI : '+$
              tostr(meansdi)+' K'

   endfor

   if (UsePkz eq 1) then begin

      get_position, ppp, space, sizes, 4, pos, /rect
      pos(0) = pos(0)+0.2
      pos(2) = pos(2)-0.2

      xt = xtitle
      xtn = xtickname

      plot, SdiTimes-stime, gitmTempZ, $
            pos = pos, /noerase, $
            ytitle = 'Vz (m/s)',		$
            xtickname = xtn,			$
            xtitle = xt,			$
            xtickv = xtickv,			$
            xminor = xminor,			$
            xticks = xtickn,   $
            yrange = range, ystyle = 1, $
            thick = 4

      oplot, PkzTime-stime, PkzTn, linestyle = 1, thick = 4
      oplot, SdiTimes-stime, SdiTempZ, linestyle = 2, thick = 4

      dx = (etr-btr)/20.0

      oplot, [dx, dx*3], [mini*0.95, 0.95*mini], thick = 4
      xyouts, dx*3.1, mini*0.95, 'GITM'
      
      oplot, [dx*7, dx*7 + dx*2], [mini*0.95, mini*0.95], thick = 4, linestyle = 2
      xyouts, dx*9.1, mini*0.95, 'SDI ('+mkupper(site)+')'

      oplot, [dx*13, dx*13 + dx*2], [mini*0.95, mini*0.95], thick = 4, linestyle = 1
      xyouts, dx*15.1, mini*0.95, 'PKZ'

   endif

   closedevice

; ---------------------------------------------------------------------------

setdevice, 'compare_los_'+site+'_'+ymd+'.ps', 'p', 5

plotdumb

ppp = 6
space = 0.01
pos_space, ppp, space, sizes, ny = ppp

titles = ['North LOS (m/s)', 'East LOS (m/s)', 'South LOS (m/s)', 'West LOS (m/s)']

stime = min(SdiTimes)
etime = max(SdiTimes)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

maxi = max(abs([gitmlosTotal, sdilosTotal]))
range = [-1.1*maxi, 1.1*maxi]

for i = 0, 3 do begin

   get_position, ppp, space, sizes, i, pos, /rect
   pos(0) = pos(0)+0.2
   pos(2) = pos(2)-0.2

   if (i eq 3 and UsePkz eq 0) then begin
      xt = xtitle
      xtn = xtickname
   endif else begin
      xt = ' '
      xtn = strarr(10)+' '
   endelse

   plot, SdiTimes-stime, gitmlosTotal(i,*), $
         pos = pos, /noerase, $
         ytitle = titles(i),		$
         xtickname = xtn,			$
         xtitle = xt,			$
         xtickv = xtickv,			$
         xminor = xminor,			$
         xticks = xtickn,   $
         yrange = range, ystyle = 1, $
         thick = 4

   oplot, SdiTimes-stime, sdilosTotal(i,*), linestyle = 2, thick = 4

   dx = (etr-btr)/20.0

   oplot, [dx, dx*3], [-maxi, -maxi], thick = 4
   xyouts, dx*3.1, -maxi, 'GITM'

   oplot, [dx*10, dx*10 + dx*2], [-maxi, -maxi], thick = 4, linestyle = 2
   xyouts, dx*12.1, -maxi, 'SDI ('+mkupper(site)+')'

   rms  = sqrt(mean((gitmlosTotal(i,*)-sdilosTotal(i,*))^2))
   diff = mean(gitmlosTotal(i,*)-sdilosTotal(i,*))
   xyouts, dx, maxi*0.85,'RMS : '+tostr(rms)+' m/s; Diff : '+tostr(diff)+' m/s'

endfor

if (UsePkz eq 1) then begin

   get_position, ppp, space, sizes, 4, pos, /rect
   pos(0) = pos(0)+0.2
   pos(2) = pos(2)-0.2

   xt = xtitle
   xtn = xtickname

   plot, SdiTimes-stime, gitmUpZ, $
         pos = pos, /noerase, $
         ytitle = 'Vz (m/s)',		$
         xtickname = xtn,			$
         xtitle = xt,			$
         xtickv = xtickv,			$
         xminor = xminor,			$
         xticks = xtickn,   $
         yrange = range/2, ystyle = 1, $
         thick = 4

   oplot, PkzTime-stime, PkzVz, linestyle = 2, thick = 4

   dx = (etr-btr)/20.0

   oplot, [dx, dx*3], [-maxi/2, -maxi/2], thick = 4
   xyouts, dx*3.1, -maxi/2, 'GITM'

   oplot, [dx*10, dx*10 + dx*2], [-maxi/2, -maxi/2], thick = 4, linestyle = 2
   xyouts, dx*12.1, -maxi/2, 'PKZ'

endif

closedevice

combined = [sdi.los_wind_skymap,gitmlos]
maxlos = max(abs(combined))
nTotal = n_elements(sdi.los_wind_skymap)
c = nTotal
while (c gt 0.99*float(nTotal)) do begin
   maxlos = 0.99*maxlos
   l = where(sdi.los_wind_skymap le maxlos, c)
endwhile
LosLevels = (findgen(31)/15-1.0)*maxlos*1.05

combined = [sdi.temp_skymap,gitmtemp]
maxtemp = max(combined)
nTotal = n_elements(combined)
c = nTotal
while (c gt 0.99*float(nTotal)) do begin
   maxtemp = 0.99*maxtemp
   l = where(combined le maxtemp, c)
endwhile

mintemp = min(combined)
c = nTotal
while (c gt 0.99*float(nTotal)) do begin
   mintemp = mintemp+1.0
   l = where(combined ge mintemp, c)
endwhile
mintemp = mintemp*0.95
maxtemp = maxtemp*1.05
TempLevels = findgen(31)/30*(maxtemp-mintemp) + mintemp

; ---------------------------------------------------------------------------

c_r_to_a, itime_start, min(SdiTimes)
c_r_to_a, itime_end, max(SdiTimes)

c_a_to_s, itime_start, stime_start
c_a_to_s, itime_end, stime_end

if (iStart eq iEnd) then begin
   start_time = ask('start time of plotting',strmid(stime_start,0,15))
   if (strlen(start_time) lt 9) then $
      start_time = strmid(stime_start,0,9)+' '+start_time
endif else start_time = stime_start

;--------------------------------------------------------------
; I got sick of typing in the ending date, so if the date is
; the same, assume the user just wants to enter the time
;--------------------------------------------------------------

sdate = strmid(start_time,0,9)
if strpos(stime_end,sdate) gt -1 then 					$
  end_time_default = strmid(stime_end,10,5)				$
else end_time_default = strmid(stime_end,0,15)

if (iStart eq iEnd) then begin

   end_time   = ask('end time of plotting',end_time_default)

;--------------------------------------------------------------
; If the user entered a short string, assume it is just a time
; and add the date on the front
;--------------------------------------------------------------

   if (strlen(end_time) lt 9) then end_time = strmid(start_time,0,9)+' '+end_time

endif else end_time = stime_end

;--------------------------------------------------------------
; Now figure out where in the file these things are, with the
; default to give the user everything
;--------------------------------------------------------------

print, start_time
print, end_time

c_s_to_a, itime_start, start_time
c_a_to_r, itime_start, rtime_start

c_s_to_a, itime_end, end_time
c_a_to_r, itime_end, rtime_end

n_start = where(SdiTimes ge rtime_start)
if n_start(0) ge 0 then n_start = n_start(0) else n_start = 0

n_end = where(SdiTimes ge rtime_end)
if n_end(0) ge 0 then n_end = n_end(0) else n_end = n_elements(Sditimes)-1

;--------------------------------------------------------------
; Now, allow the user to skip times
;--------------------------------------------------------------

print, 'You have selected '+tostr(n_end-n_start+1)+' times to plot.'

if (iStart eq iEnd) then $
   step = ask('number of times to skip between each plot','0')+1 $
else step = (iEnd-iStart)/4

plotnum = -1
page = 1

ppp = 20
space = 0.01

for iT = n_start, n_End, step do begin

   c_r_to_a, itime, SdiTimes(iT)
   c_a_to_s, itime, stime

   plotnum = (plotnum + 1) mod ppp

   if (plotnum eq 0) then begin

      setdevice, 'compare_circle_'+site+'_'+ymd+'_page_'+tostr(page,2)+'.ps', 'l', 5
      page++

      pos_space, ppp, space, sizes, nx = 5

      plotdumb

   endif

   get_position, ppp, space, sizes, plotnum, pos
   pos([1,3]) = pos([1,3]) + space*0.5 + 0.75*space*(float(fix(plotnum/5)))

   if (plotnum eq 0) then begin
      xyouts, 1.0-pos(0), pos(3)+0.03, strmid(stime,0,9), $
              /norm, align=1.0, charsize=1.2
      xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.01, mkupper(site)+' LOS', $
              /norm, align=0.5
   endif

   makect, 'mid'

   sdilos = sdi.los_wind_skymap(iT,*)
   l = where(sdilos lt loslevels(1), c)
   if (c gt 0) then sdilos(l) = loslevels(1)
   l = where(sdilos gt loslevels(29), c)
   if (c gt 0) then sdilos(l) = loslevels(29)

   contour, sdilos, x, y, /irr, $
            pos = pos, /noerase, ystyle = 5, xstyle = 5, $
            levels = LosLevels, /fill, $
            xrange = [-maxran, maxran], $
            yrange = [-maxran, maxran]
   plotmlt, maxran, no00 = 1, no06 = 1, no12 = 1, no18 = 1
   xyouts, pos(0)-0.005, (pos(1)+pos(3))/2, strmid(stime,10,5)+' UT', $
           align=0.5, orient=90, /norm

   if (plotnum+5 ge ppp or iT+step gt nTimesSdi-1) then begin
      ctpos = pos
      ctpos(3) = pos(1)-0.01
      ctpos(1) = ctpos(3)-0.015
      ctpos(0) = ctpos(0)+0.02
      ctpos(2) = ctpos(2)-0.02
      ncolors = 254
      maxmin = mm(LosLevels)
      title = mkupper(site)+' LOS (m/s)'
      plotct, ncolors, ctpos, maxmin, title, /bottom
   endif

   plotnum++
   get_position, ppp, space, sizes, plotnum, pos
   pos([1,3]) = pos([1,3]) + space*0.5 + 0.75*space*(float(fix(plotnum/5)))

   gitmlos1d = reform(gitmlos(iT,*))

   l = where(gitmlos1d lt loslevels(1), c)
   if (c gt 0) then gitmlos1d(l) = loslevels(1)
   l = where(gitmlos1d gt loslevels(29), c)
   if (c gt 0) then gitmlos1d(l) = loslevels(29)

   contour, gitmlos1d, x, y, /irr, $
            pos = pos, /noerase, ystyle = 5, xstyle = 5, $
            levels = LosLevels, /fill, $
            xrange = [-maxran, maxran], $
            yrange = [-maxran, maxran]
   for iP = 0, nPts-1 do begin
      oplot, [x(iP)], [y(iP)], psym = 4
      plots, [x(iP), x(iP) + east(iT,iP)/maxran * 10.0], $
             [y(iP), y(iP) + north(iT,iP)/maxran * 10.0]
   endfor
   plotmlt, maxran, no00 = 1, no06 = 1, no12 = 1, no18 = 1

   if (plotnum-5 lt 0) then $
      xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.01, 'GITM LOS', /norm, align=0.5

   if (plotnum+5 ge ppp or iT+step gt nTimesSdi-1) then begin
      ctpos = pos
      ctpos(3) = pos(1)-0.01
      ctpos(1) = ctpos(3)-0.015
      ctpos(0) = ctpos(0)+0.02
      ctpos(2) = ctpos(2)-0.02
      ncolors = 254
      maxmin = mm(LosLevels)
      title = 'GITM LOS (m/s)'
      plotct, ncolors, ctpos, maxmin, title, /bottom
   endif

   plotnum++
   get_position, ppp, space, sizes, plotnum, pos
   pos([1,3]) = pos([1,3]) + space*0.5 + 0.75*space*(float(fix(plotnum/5)))

   gitmlos1d = reform(gitmlosI(iT,*))

   l = where(gitmlos1d lt 2.5*loslevels(1), c)
   if (c gt 0) then gitmlos1d(l) = 2.5*loslevels(1)
   l = where(gitmlos1d gt 2.5*loslevels(29), c)
   if (c gt 0) then gitmlos1d(l) = 2.5*loslevels(29)

   contour, gitmlos1d, x, y, /irr, $
            pos = pos, /noerase, ystyle = 5, xstyle = 5, $
            levels = 2.5*LosLevels, /fill, $
            xrange = [-maxran, maxran], $
            yrange = [-maxran, maxran]
   for iP = 0, nPts-1 do begin
      oplot, [x(iP)], [y(iP)], psym = 4
      plots, [x(iP), x(iP) + eastI(iT,iP)/maxran * 2.5], $
             [y(iP), y(iP) + northI(iT,iP)/maxran * 2.5]
   endfor
   plotmlt, maxran, no00 = 1, no06 = 1, no12 = 1, no18 = 1

   if (plotnum-5 lt 0) then $
      xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.01, 'GITM Ion LOS', /norm, align=0.5

   if (plotnum+5 ge ppp or iT+step gt nTimesSdi-1) then begin
      ctpos = pos
      ctpos(3) = pos(1)-0.01
      ctpos(1) = ctpos(3)-0.015
      ctpos(0) = ctpos(0)+0.02
      ctpos(2) = ctpos(2)-0.02
      ncolors = 254
      maxmin = mm(2.5*LosLevels)
      title = 'GITM Ion LOS (m/s)'
      plotct, ncolors, ctpos, maxmin, title, /bottom
   endif

   makect,'bry'

   plotnum++
   get_position, ppp, space, sizes, plotnum, pos
   pos([1,3]) = pos([1,3]) + space*0.5 + 0.75*space*(float(fix(plotnum/5)))
   l = where(gitmtemp(iT,*) lt templevels(1), c)
   if (c gt 0) then gitmtemp(iT,l) = templevels(1)
   l = where(gitmtemp(iT,*) gt templevels(29), c)
   if (c gt 0) then gitmtemp(iT,l) = templevels(29)
   contour, gitmtemp(iT,*), x, y, /irr, $
            pos = pos, /noerase, ystyle = 5, xstyle = 5, $
            levels = TempLevels, /fill, $
            xrange = [-maxran, maxran], $
            yrange = [-maxran, maxran]
   plotmlt, maxran, no00 = 1, no06 = 1, no12 = 1, no18 = 1

   if (plotnum-5 lt 0) then $
      xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.01, 'GITM Temp', /norm, align=0.5

   if (plotnum+5 ge ppp or iT+step gt nTimesSdi-1) then begin
      ctpos = pos
      ctpos(3) = pos(1)-0.01
      ctpos(1) = ctpos(3)-0.015
      ctpos(0) = ctpos(0)+0.02
      ctpos(2) = ctpos(2)-0.02
      ncolors = 254
      maxmin = mm(TempLevels)
      title = 'GITM Temp (K)'
      plotct, ncolors, ctpos, maxmin, title, /bottom
   endif

   sditemp = sdi.temp_skymap(iT,*)

   plotnum++
   get_position, ppp, space, sizes, plotnum, pos
   pos([1,3]) = pos([1,3]) + space*0.5 + 0.75*space*(float(fix(plotnum/5)))
   l = where(sditemp lt templevels(1), c)
   if (c gt 0) then sditemp(l) = templevels(1)
   l = where(sditemp gt templevels(29), c)
   if (c gt 0) then sditemp(l) = templevels(29)
   contour, sditemp, x, y, /irr, $
            pos = pos, /noerase, ystyle = 5, xstyle = 5, $
            levels = TempLevels, /fill, $
            xrange = [-maxran, maxran], $
            yrange = [-maxran, maxran]
   plotmlt, maxran, no00 = 1, no06 = 1, no12 = 1, no18 = 1

   if (plotnum-5 lt 0) then $
      xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.01, mkupper(site)+' Temp', /norm, align=0.5

   if (plotnum+5 ge ppp or iT+step gt nTimesSdi-1) then begin
      ctpos = pos
      ctpos(3) = pos(1)-0.01
      ctpos(1) = ctpos(3)-0.015
      ctpos(0) = ctpos(0)+0.02
      ctpos(2) = ctpos(2)-0.02
      ncolors = 254
      maxmin = mm(TempLevels)
      title = mkupper(site)+' Temp (K)'
      plotct, ncolors, ctpos, maxmin, title, /bottom
   endif

endfor

closedevice

endfor

end



