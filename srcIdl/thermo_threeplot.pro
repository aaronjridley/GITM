
pro thermo_threeplot, psfile, value, time, lons, lats, mini, maxi, $
                      title, colortitle, maxrange, $
                      vn = vn, ve = ve, apex = apex, lines = lines, $
                      noclosepsfile = noclosepsfile, $
                      noopenpsfile = noopenpsfile, $
                      nlevels = nlevels, uppercolor = uppercolor, $
                      factor = factor

  if (n_elements(apex) gt 0) then PlotApex = 1 else PlotApex = 0
  if (n_elements(lines) gt 0) then PlotLinesOnly = 1 else PlotLinesOnly = 0

  if (n_elements(uppercolor) gt 0) then PlotUpperColor=1 else PlotUpperColor=0
  
  if (n_elements(nLevels) eq 0) then begin
     if (PlotLinesOnly) then nLevels = 11 else nLevels = 31
  endif
  if (mini*maxi lt 0.0) then begin
     maxi = max([maxi,abs(mini)])
     mini = -maxi
  endif
  if (nlevels eq 1) then begin
     levels = [(mini+maxi)/2]
  endif else begin
     levels = findgen(nlevels)*(maxi-mini)/(nlevels-1) + mini
  endelse

  nLons = n_elements(value(*,0))
  nLats = n_elements(value(0,*))
  newvalue = reform(value(1:nLons-2,1:nLats-2))

  newlat = lats(1:nLons-2,1:nLats-2)/!dtor
  newlon = lons(1:nLons-2,1:nLats-2)/!dtor
  nLo  = nLons-2
  nLa  = nLats-2

  newvalue(0,*)       = (newvalue(1,*)+newvalue(nLo-2,*))/2.0
  newvalue(nLo-1,*) = (newvalue(1,*)+newvalue(nLo-2,*))/2.0
  newvalue(*,0)       = mean(newvalue(*,1))
  newvalue(*,nLa-1) = mean(newvalue(*,nLa-2))

  newlon(0,*)       = 0.0
  newlon(nLo-1,*) = 360.0
  newlat(*,0) = -90.0
  newlat(*,nLa-1) =  90.0
  
  c_r_to_a, itime, time
  c_a_to_s, itime, stime

  ut = itime(3) + itime(4)/60.0 + itime(5)/3600.0
  utrot = ut * 15.0

  localtime = (newlon/15.0 + ut) mod 24.0
  angle = 23.0 * !dtor * $
          sin((jday(itime(0),itime(1),itime(2)) - $
               jday(itime(0),3,21))*2*!pi/365.0)
  sza =  acos(sin(angle)*sin(newlat*!dtor) + $
              cos(angle)*cos(newlat*!dtor) * $ 
              cos(!pi*(LocalTime-12.0)/12.0))/!dtor

  if (n_elements(noopenpsfile) eq 0) then begin
     setdevice, psfile, 'p', 5 
     plotdumb
  endif

  if (mini lt 0) then makect,'mid' $
  else makect,'all'

  ppp = 4
  space = 0.03
  pos_space, ppp, space, sizes

;--------------------------------------------------------------
; Figure out where on the page we should be
;--------------------------------------------------------------

  get_position, ppp, space, sizes, 0, pos

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

  no00 = 0
  no06 = 1
  no12 = 0
  no18 = 0
  
  lon = reform(newlon(*,0))+utrot

  contour_circle, newvalue, lon, reform(newlat(0,*)), $
                  mini = mini, maxi = maxi, $
                  nLevels = nLevels, $
                  no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                  pos = pos, $
                  maxrange = maxrange, $
                  nocolor=PlotLinesOnly, nolines=1-PlotLinesOnly, $
                  nomaxmin = PlotLinesOnly

  if (n_elements(vn) gt 0 and n_elements(ve) gt 0) then begin

     if (n_elements(factor) eq 0) then factor = max(sqrt(ve^2+vn^2))/10.0

     step = 2

     for iLat=2,nLats-3,step do begin

        la = newlat(0,iLat)*!dtor
        lonstep = round(step / (cos(la)+0.1))

        for iLon = 2, nLons-3,lonstep do begin
           r = 90.0 - lats(iLon,iLat)/!dtor 
           if (r lt maxrange and r gt 0) then begin
              lo = lons(iLon,iLat) + utrot*!dtor - !pi/2 
              x = r * cos(lo)
              y = r * sin(lo)
              ux = - vn(iLon,iLat) * cos(lo) - ve(iLon,iLat) * sin(lo)
              uy = - vn(iLon,iLat) * sin(lo) + ve(iLon,iLat) * cos(lo)
              ux = ux/factor
              uy = uy/factor
              oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

              u = sqrt(ux^2+uy^2)
              if (u gt 0) then begin
                 t = asin(uy/u)
                 if (ux lt 0) then t = !pi-t
                 t1 = t+!pi/12
                 t2 = t-!pi/12
                 ux1 = 0.6 * u * cos(t1)
                 uy1 = 0.6 * u * sin(t1)
                 ux2 = 0.6 * u * cos(t2)
                 uy2 = 0.6 * u * sin(t2)
                 oplot,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
                 oplot,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0
              endif

           endif
        endfor
     endfor

     x = maxrange*0.95
     y = maxrange*0.95
     ux = 0.0
     uy = -factor*10.0
     ux = ux/factor
     uy = uy/factor
     oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

     u = sqrt(ux^2+uy^2)
     if (u gt 0) then begin
        t = asin(uy/u)
        if (ux lt 0) then t = !pi-t
        t1 = t+!pi/12
        t2 = t-!pi/12
        ux1 = 0.6 * u * cos(t1)
        uy1 = 0.6 * u * sin(t1)
        ux2 = 0.6 * u * cos(t2)
        uy2 = 0.6 * u * sin(t2)
        oplot,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
        oplot,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0
     endif

     xyouts, maxrange, y+0.5*uy, tostr(factor*10)+' m/s'

  endif


  xc = (pos(2)+pos(0))/2
  xr = (pos(2)-pos(0))/2 * 1.01
  yc = (pos(3)+pos(1))/2
  yr = (pos(3)-pos(1))/2 * 1.01
  xp = xc - xr*sin(!pi/4)
  yp = yc + yr*sin(!pi/4)
  xyouts, xp, yp, 'North', $
          /norm, charsize = 0.9, align = 0.5, orient = 45

  xyouts, pos(2)+space/2, pos(3), title, align = 0.5, /norm

  xp = pos(0)
  yp = pos(3)+yr/20.0
  xyouts, xp, yp, strmid(stime,0,9), $
          /norm, charsize = 0.9, align = 0.0

;************************************************************************

  get_position, ppp, space, sizes, 1, pos

  pos(0) = pos(0)-space+0.01
  pos(2) = pos(2)-space+0.01

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

  no00 = 0
  no06 = 1
  no12 = 0
  no18 = 1

  contour_circle, newvalue, lon, -reform(newlat(0,*)), $
                  mini = mini, maxi = maxi, $
                  nLevels = nLevels, $
                  no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                  pos = pos, $
                  maxrange = maxrange, $
                  nocolor=PlotLinesOnly, nolines=1-PlotLinesOnly, $
                  nomaxmin = PlotLinesOnly


  if (n_elements(vn) gt 0 and n_elements(ve) gt 0) then begin

     step = 2

     for iLat=2,nLats-3,step do begin

        la = newlat(0,iLat)*!dtor
        lonstep = round(step / (cos(la)+0.1))

        for iLon = 2, nLons-3,lonstep do begin

           r = 90.0 + lats(0,iLat)/!dtor 
           if (r lt maxrange and r gt 0) then begin
              lo = lons(iLon,iLat) + utrot*!dtor - !pi/2 
              x = r * cos(lo)
              y = r * sin(lo)
              ux =   vn(iLon,iLat) * cos(lo) - ve(iLon,iLat) * sin(lo)
              uy =   vn(iLon,iLat) * sin(lo) + ve(iLon,iLat) * cos(lo)
              ux = ux/factor
              uy = uy/factor
              oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

              u = sqrt(ux^2+uy^2)
              if (u gt 0) then begin
                 t = asin(uy/u)
                 if (ux lt 0) then t = !pi-t
                 t1 = t+!pi/12
                 t2 = t-!pi/12
                 ux1 = 0.6 * u * cos(t1)
                 uy1 = 0.6 * u * sin(t1)
                 ux2 = 0.6 * u * cos(t2)
                 uy2 = 0.6 * u * sin(t2)
                 oplot,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
                 oplot,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0
              endif

           endif

        endfor
     endfor

  endif

  xc = (pos(2)+pos(0))/2
  xr = (pos(2)-pos(0))/2 * 1.01
  yc = (pos(3)+pos(1))/2
  yr = (pos(3)-pos(1))/2 * 1.01
  xp = xc + xr*sin(!pi/4)
  yp = yc + yr*sin(!pi/4)
  xyouts, xp, yp, 'South', $
          /norm, charsize = 0.9, align = 0.5, orient = -45

  xp = pos(2)
  yp = pos(3)+yr/20.0
  xyouts, xp, yp, strmid(stime,10,5)+' UT', $
          /norm, charsize = 0.9, align = 1.0

  if (PlotUpperColor) then begin
     ctpos = [pos(2)+0.005, pos(1), pos(2)+0.02, pos(3)]
     range = [mini,maxi]
     ncolors = 255
     plotct, ncolors, ctpos, range, colortitle, /right
  endif

;************************************************************************

;--------------------------------------------------------------
; Figure out where on the page we should be
;--------------------------------------------------------------

  get_position, ppp, space, sizes, 2, pos1
  get_position, ppp, space, sizes, 3, pos2

  pos = pos1
  pos(2) = pos2(2)-0.04

  !p.position = pos

  utime = itime(3)*3600.0 + $
          itime(4)*60.0 + $
          itime(5)
  utime = utime(0)
  p0lon = utime/3600.0 * 360.0 / 24.0

  map_set, 0.0, 180.0-p0lon, /noerase

  newvalue_limited = newvalue

  if (nlevels gt 1) then begin
     l = where(newvalue_limited lt levels(1),c)
     if (c gt 0) then newvalue_limited(l) = levels(1)
     l = where(newvalue_limited gt levels(nLevels-2),c)
     if (c gt 0) then newvalue_limited(l) = levels(nLevels-2)
  endif

  if (PlotLinesOnly) then begin
     contour, newvalue_limited, newlon, newlat, $
              /follow, /over, $
              levels = levels
  endif else begin
     contour, newvalue_limited, newlon, newlat, $
              /follow, /cell_fill, /over, $
              levels = levels
  endelse

  contour, sza, newlon, newlat, levels = [45,90], thick = 4, /over, $
           c_linestyle = 2

  if (PlotApex) then begin
     contour, apex.alats, apex.glons, apex.glats, $
              levels = [0.0], thick = 4, $
              xrange = [0.0,360.0], xstyle = 1, $
              yrange = [-90.0,90.0], ystyle = 1, /over, $
              c_linestyle = 3
  endif

  if (n_elements(vn) gt 0 and n_elements(ve) gt 0) then begin

     step = 4

     for iLat=4,nLats-5,step do begin

        la = newlat(0,iLat)
;        lonstep = round(step / (cos(la*!dtor)+0.1))
        lonstep = step/2
        for iLon = 2, nLons-3,lonstep do begin

           x = newlon(iLon,iLat)
           y = la

           ; Factor of 2 is due to the width being 2 times the height
           ; I don't know if this is deceptive or not.  I don't think so.

           ; Calculate the true length
           uxS = ve(iLon,iLat)/factor
           uyS = vn(iLon,iLat)/factor
           rS = sqrt(ux^2+uy^2)

           ; Magnitude of the rotated in the proper direction:
           uyP = 2*vn(iLon,iLat)/factor
           rP = sqrt(ux^2+uy^2)

           ux = uxS/rP * rS
           uy = uyP/rP * rS

           oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

           u = sqrt(ux^2+uy^2)

           if (u gt 0) then begin
              t = asin(uy/u)
              if (ux lt 0) then t = !pi-t
              t1 = t+!pi/12
              t2 = t-!pi/12
              ux1 = 0.6 * u * cos(t1)
              uy1 = 0.6 * u * sin(t1)
              ux2 = 0.6 * u * cos(t2)
              uy2 = 0.6 * u * sin(t2)
              oplot,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
              oplot,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0
           endif

        endfor
     endfor

  endif

  map_continents
  map_grid, lats = findgen(19)*10-90, glinethick=3

  !p.position = -1

  mini_tmp = min(newvalue)
  maxi_tmp = max(newvalue)

  if (abs(mini_tmp) gt 1000.0) or (abs(maxi_tmp) gt 1000.0) or        $
     (abs(maxi_tmp) lt 0.1) then begin
     maxs = string(maxi_tmp,format="(e9.2)")
     mins = string(mini_tmp,format="(e9.2)")
  endif else begin
     maxs = string(maxi_tmp,format="(f7.2)")
     mins = string(mini_tmp,format="(f7.2)")
  endelse

  xp = pos(2)
  yp = pos(1)-yr/10.0
  xyouts, xp, yp, 'Max: '+maxs, charsize = 0.9, align = 1.0, /norm

  xp = pos(0)
  yp = pos(1)-yr/10.0
  xyouts, xp, yp, 'Min: '+mins, charsize = 0.9, align = 0.0, /norm

  if (not PlotLinesOnly) then begin
     ctpos = [pos(2)+0.005, pos(1), pos(2)+0.02, pos(3)]
     range = [mini,maxi]
     ncolors = 255
     plotct, ncolors, ctpos, range, colortitle, /right
  endif

  if (n_elements(noclosepsfile) eq 0) then closedevice
     
end
