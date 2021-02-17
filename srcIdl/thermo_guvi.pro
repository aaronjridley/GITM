
guvidir = '/raid3/Data/GUVI/'

filelist = findfile('3DALL*.bin')

file = filelist(0)
gitm_read_header, file, StartTime, nVars, Vars, $
                  nLons, nLats, nAlts, version


file = filelist(n_elements(filelist)-1)
gitm_read_header, file, EndTime, nVars, Vars, $
                  nLons, nLats, nAlts, version

nDays = (EndTime-StartTime)/86400 + 1
files = strarr(nDays)

nFiles = 0
for i = 0,nDays-1 do begin
   t = StartTime + 86400.0D * double(i)
   c_r_to_a, itime, t
   sYear = tostr(itime(0))
   sDay = tostr(jday(itime(0),itime(1),itime(2)),3)
   filelist = findfile(guvidir+sYear+'/ON2_'+sYear+'_'+sDay+'m.sav')
   if (strlen(filelist(0)) gt 0) then begin
      files(nFiles) = filelist(0)
      nFiles++
   endif
endfor
filelist = files(0:nFiles-1)
print, filelist

nfiles=n_elements(filelist)
print,'nfiles=',nfiles

for iFile=0,nFiles-1 do begin

   print,'restoring file ... ', filelist(iFile)

   restore,filelist(iFile)
   fdoy = saved_data.fdoy
   lat = saved_data.lat
   lon = saved_data.lon
   on2 = saved_data.on2
   loc = where(fdoy gt 0.0,count)

   if (count gt 0) then begin

      itime = intarr(6)
      itime(0) = fix(saved_data[0].year)
      itime(1) = 1
      timesub = dblarr(count)
      latsub = lat(loc)
      lonsub = lon(loc)
      on2sub = on2(loc)

      for il=0,count-1 do begin
         i = loc(il)
         itime(2) = fix(fdoy(i))
         frac_of_day = fdoy(i)-itime(2)
         itime(3) = fix(frac_of_day*24.0)
         frac_of_hour = frac_of_day*24.0 - itime(3)
         itime(4) = fix(frac_of_hour*60.0)
         frac_of_min = frac_of_hour*60.0-itime(4)
         itime(5) = frac_of_min*60.0
         c_a_to_r, itime, rtime
         timesub(il) = rtime
      endfor

      Guvilats  = latsub
      GuviLons  = lonsub
      GuviOn2   = on2sub
      GuviTimes = timesub

      l = where(GuviTimes ge StartTime and GuviTimes le EndTime,c)

      if (c gt 0) then begin
         GuviTimes = GuviTimes(l)
         GuviLats  = GuviLats(l)
         GuviLons  = GuviLons(l)
         GuviOn2   = GuviOn2(l)
      endif else begin
         print, 'Didnt find any guvi times!'
         stop
      endelse

   endif

   iO_ = 4
   iN2_ = 6
   iNO_ = 8
   varstoget = [iO_,iN2_,iNO_]
   GuviAlts = 0.0*GuviLons - 1.0
   get_gitm_points, GuviTimes, GuviLons, GuviLats, GuviAlts, VarsToGet, GitmData

   get_gitm_alts, gitmalts
   nAlts = n_elements(gitmalts)
   dAlts = fltarr(nAlts)
   dAlts(1:nAlts-2) = (gitmalts(2:nAlts-1)-gitmalts(0:nAlts-3))/2.0
   dAlts(0) = gitmalts(1)-gitmalts(0)
   dAlts(nAlts-1) = gitmalts(nAlts-1)-gitmalts(nAlts-2)
   dAlts = dAlts*1000.0         ; km -> m

   nTimes = n_elements(GuviTimes)
   on2 = fltarr(nTimes)
   no = fltarr(nTimes)

   LocAlts = where(GitmAlts ge 100.0 and GitmAlts lt 150.0)

   for iTime = 0L,nTimes-1 do begin

      o_cell  = reform(GitmData(iTime,0,*))*dAlts
      n2_cell = reform(GitmData(iTime,1,*))*dAlts
      no_cell = reform(GitmData(iTime,2,*))*dAlts
      
      o_int  = o_cell
      n2_int = n2_cell
      o_int(nAlts-1)  = o_cell(nAlts-1)
      n2_int(nAlts-1) = n2_cell(nAlts-1)
      iAlt = nAlts-1
      while (n2_int(iAlt) lt 1.0e21) do begin
         iAlt--
         o_int(iAlt)  = o_int(iAlt+1)  + o_cell(iAlt)
         n2_int(iAlt) = n2_int(iAlt+1) + n2_cell(iAlt)
      endwhile
      x = (n2_int(iAlt) - 1.0e21)/(n2_int(iAlt) - n2_int(iAlt+1))
      n2 = (1.0-x)*n2_int(iAlt) + x*n2_int(iAlt+1)
      o  = (1.0-x)*o_int(iAlt)  + x*o_int(iAlt+1)
      
      on2(iTime) = o/n2
      no(iTime) = total(no_cell(LocAlts))
      
   endfor

;shift = -360.0*float(fix((GuviTimes-GuviTimes(0))/(24.0*3600.0)))
   shift = 0.0
   mini = min(no)
   maxi = max(no)
   mini = 0.0
   maxi = 8.0e18
   levels = findgen(31)/30.0*(maxi-mini)+mini
   
   dLon = 15.0
   dLat =  2.0

   ShiftedLons = guvilons+shift
   minLon = float(fix(min(ShiftedLons)/dLon))*dLon
   maxLon = float(fix(max(ShiftedLons)/dLon)+1)*dLon
   minlat = float(fix(min(guvilats)/dLat)-1)*dLat
   maxlat = float(fix(max(guvilats)/dLat)+1)*dLat

   nLats = (maxLat-minLat)/dlat + 1
   nLons = (maxLon-minLon)/dlon + 1

   GriddedNo = fltarr(nLons,nLats)
   GriddedOn2 = fltarr(nLons,nLats)
   GriddedGuviOn2 = fltarr(nLons,nLats)
   GriddedLat = fltarr(nLons,nLats)
   GriddedLon = fltarr(nLons,nLats)

   dummy = no*0.0 + 1.0

   iLon = round((ShiftedLons-minLon)/dLon)
   iLat = round((guvilats-minlat)/dLat)

   for iLo = 0,nLons-1 do begin
      for iLa = 0,nLats-1 do begin

         GriddedLat(iLo,iLa) = minlat + dLat*iLa
         GriddedLon(iLo,iLa) = minlon + dLon*iLo

         la = where(iLat eq iLa and iLon eq iLo, c)
         if (c gt 0) then begin
            GriddedNo(iLo,iLa) = mean(no(la))
            GriddedOn2(iLo,iLa) = mean(on2(la))
            GriddedGuviOn2(iLo,iLa) = mean(guvion2(la))
         endif

      endfor
   endfor

   for iLo = 0,nLons-1 do begin
      for iLa = 0,nLats-1 do begin
         if (GriddedNo(iLo,iLa) eq 0) then begin
            iLoS = max([0,iLo-1])
            iLoE = min([iLo+1,nLons-1])
            iLaS = max([0,iLa-1])
            iLaE = min([iLa+1,nLats-1])
            small = GriddedNo(iLoS:iLoE,iLaS:iLaE)
            smallon2 = GriddedOn2(iLoS:iLoE,iLaS:iLaE)
            smallguvion2 = GriddedGuviOn2(iLoS:iLoE,iLaS:iLaE)
            l = where(small gt 0.0,c)
            if (c gt 0) then GriddedNo(iLo,iLa) = mean(small(l))
            l = where(smallon2 gt 0.0,c)
            if (c gt 0) then GriddedOn2(iLo,iLa) = mean(smallon2(l))
            l = where(smallGuviOn2 gt 0.0,c)
            if (c gt 0) then GriddedGuviOn2(iLo,iLa) = mean(smallguvion2(l))
         endif
      endfor
   endfor

   stime = min(guvitimes)
   etime = max(guvitimes)
   c_r_to_a, itime, stime
   c_a_to_s, itime, sdate
   c_r_to_a, itime, etime
   c_a_to_s, itime, edate

   c_r_to_a, itime, (stime+etime)/2
   c_a_to_ymd, itime, ymd

   setdevice, 'guvi_no_'+ymd+'.ps','p', 5

   makect, 'rainbow'

   print, mm(GriddedLon)

   pos = [0.05, 0.55, 0.95, 0.95]
   !p.position = pos
   map_set, 0.0, 0.0

   GriddedNo(0,*) = (GriddedNo(0,*) + GriddedNo(nLons-1,*))/2
   GriddedNo(nLons-1,*) = GriddedNo(0,*)
   contour, GriddedNo, GriddedLon, GriddedLat, /cell, levels=levels, $
            pos = pos, xstyle = 1, ystyle = 1, $
            xtitle = 'Shifted Longitude (deg)', ytitle = 'Latitude (deg)',$
            /over

   map_continents
   map_grid, lats = findgen(7)*30-90, glinethick=3

   ncolors = 254
   title = 'NO Column Integral (/m2)'
   ctpos = [0.96, 0.55, 0.98, 0.95]
   plotct, ncolors, ctpos, mm(levels), title, /right

;contour, no, ShiftedLons, guvilats, /irr, /fill, levels=levels, $
;         pos = [0.1, 0.05, 0.95, 0.45], /noerase
;ctpos = [0.96, 0.05, 0.98, 0.45]
;plotct, ncolors, ctpos, mm(levels), title, /right

   xyouts, 0.05, 0.955, strmid(edate,0,15), /norm
   xyouts, 0.95, 0.955, strmid(sdate,0,15), /norm, align=1.0

   closedevice


   setdevice, 'guvi_on2_'+ymd+'.ps','p', 5

   makect, 'rainbow'

   mapsize = 0.6
   xMin = 0.5-mapsize/2.0
   yMax = 0.99
   xMax = xMin + mapsize
   yMin = yMax - mapsize/2.0

   pos = [xMin, yMin, xMax, yMax]
   !p.position = pos
   map_set, 0.0, 0.0

   GriddedGuvion2(0,*) = (GriddedGuvion2(0,*) + GriddedGuvion2(nLons-1,*))/2
   GriddedGuvion2(nLons-1,*) = GriddedGuvion2(0,*)

   GriddedOn2(0,*) = (GriddedOn2(0,*) + GriddedOn2(nLons-1,*))/2
   GriddedOn2(nLons-1,*) = GriddedOn2(0,*)

   rms = sqrt(mean((On2(l) - Guvion2(l))^2))
   dif = mean(On2(l) - Guvion2(l))
   cor = c_correlate(on2,guvion2,0.0)

   maxi = max([GriddedOn2,GriddedGuvion2])*1.05
;   maxi = max(GriddedGuvion2)*1.01
   mini = 0.0

;   maxi = 0.8
;   mini = 0.0

   levels = findgen(31)/30.0*(maxi-mini) + mini

   contour, GriddedGuvion2, GriddedLon, GriddedLat, /cell, levels=levels, $
            pos = pos, xstyle = 1, ystyle = 1, /over, $
            xtitle = 'Shifted Longitude (deg)', ytitle = 'Latitude (deg)'
   map_continents
   map_grid, lats = findgen(7)*30-90, glinethick=3

   xyouts, xMax, yMax+0.005, strmid(sdate,0,15), /norm, align=1.0
   xyouts, xMin, yMax+0.005, strmid(edate,0,15), /norm
   xyouts, (xMax+xMin)/2.0, yMax+0.005, 'GUVI', /norm, align = 0.5

   ncolors = 254
   title = 'O/N2 Ratio'
   ctpos = [xMax+0.01, yMin, xMax+0.03, yMax]
   plotct, ncolors, ctpos, mm(levels), title, /right

   yMax = yMin - 0.03
   yMin = yMax - mapsize/2.0

   pos = [xMin, yMin, xMax, yMax]
   !p.position = pos
   map_set, 0.0, 0.0, /noerase

;   maxi = max([GriddedOn2,GriddedGuvion2])*1.05
;   maxi = max(GriddedOn2)*1.01
;   mini = 0.0

;   maxi = 1.2
;   mini = 0.2

   levels = findgen(31)/30.0*(maxi-mini) + mini

   l = where(GriddedOn2 gt levels(29),c)
   if (c gt 0) then GriddedOn2(l) = levels(29)
   contour, GriddedOn2, GriddedLon, GriddedLat, /cell, levels=levels, $
            pos = pos, /noerase, xstyle = 1, ystyle = 1, /over, $
            xtitle = 'Shifted Longitude (deg)', ytitle = 'Latitude (deg)'
   map_continents
   map_grid, lats = findgen(7)*30-90, glinethick=3

   xyouts, xMax, yMax+0.005, strmid(sdate,0,15), /norm, align=1.0
   xyouts, xMin, yMax+0.005, strmid(edate,0,15), /norm
   xyouts, (xMax+xMin)/2.0, yMax+0.005, 'GITM', /norm, align = 0.5

;linelevels = findgen(11)/10.0*1.5
;contour, GriddedGuvion2, GriddedLon, GriddedLat, /follow, levels=linelevels, $
;         pos = [0.05, 0.05, 0.95, 0.45], xstyle = 5, ystyle = 5, /noerase

   ctpos = [xMax+0.01, yMin, xMax+0.03, yMax]
   plotct, ncolors, ctpos, mm(levels), title, /right

   xyouts, xMin, 0.00, 'RMS  = '+string(rms,format='(f6.3)'), /norm
   xyouts, 0.50, 0.00, 'DIFF = '+string(dif,format='(f6.3)'), /norm,align=0.5
   xyouts, xMax, 0.00, 'COR = '+string(cor,format='(f6.3)'), /norm,align=1.0

;xyouts, 0.05, 0.955, strmid(edate,0,9), /norm

   yMax = yMin - 0.03
   
   space = 0.01
   pos_space, 4, space, sizes, ny = 4
   get_position, 4, space, sizes, 3, pos

   dy = pos(3)-pos(1)
   pos(3) = yMax
   pos(1) = pos(3) - dy

   plot, guvion2, on2, xrange = [0,maxi], yrange = [0,maxi], $
         xtitle = 'GUVI O/N2 Ratio', ytitle = 'GITM O/N2 Ratio', $
         /noerase, pos = pos, psym = 3, thick = 4, xstyle = 1, ystyle = 1

   oplot, [0,maxi], [0,maxi], thick = 4

   closedevice

endfor


end
