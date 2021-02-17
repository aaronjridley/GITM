
if (n_elements(filelist) eq 0) then begin
   filelist = findfile("-t /*.save")
   if (strlen(filelist(0)) eq 0) then filelist = findfile("-t *.bin")
endif

filelist = ask('filename to plot',filelist(0))
filelist = findfile(filelist)

nfiles = n_elements(filelist)

if (nFiles eq 1 and strlen(filelist(0)) eq 0) then begin
   print, 'can not find file!'
   stop
endif

if (n_elements(psfile) lt 0) then psfile = filelist(0)+'.ps'
psfile = ask('ps file name',psfile)

iszonalaverage = 0

for iFile = 0, nFiles-1 do begin

   filename = filelist(iFile)

   print, 'Reading file ',filename

   read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                           vars, data, rb, cb, bl_cnt, iTime, Version

   alt = reform(data(2,*,*,*)) / 1000.0
   lat = reform(data(1,*,*,*)) / !dtor
   lon = reform(data(0,*,*,*)) / !dtor

   if (iFile eq 0) then begin

      for i=0,nvars-1 do print, tostr(i)+'. '+vars(i)
      if (n_elements(sel) eq 0) then sel = '3' else sel = tostr(sel)
      sel = fix(ask('which var to plot',sel))

      if (n_elements(smini) eq 0) then smini = '0.0'
      if (n_elements(smaxi) eq 0) then smaxi = '0.0'
      smini = ask('minimum (0.0 for automatic)',smini)
      smaxi = ask('maximum (0.0 for automatic)',smaxi)

      if (n_elements(sel2) eq 0) then sel2 = '33' else sel2 = tostr(sel2)
      sel2 = fix(ask('which var2 to plot',sel2))

      if (n_elements(smini2) eq 0) then smini2 = '0.0'
      if (n_elements(smaxi2) eq 0) then smaxi2 = '0.0'
      smini2 = ask('minimum (0.0 for automatic)',smini2)
      smaxi2 = ask('maximum (0.0 for automatic)',smaxi2)

      for i=0,nalts-1 do print, tostr(i)+'. '+string(alt(2,2,i))
      if (n_elements(iAlt) eq 0) then iAlt='0' else iAlt=tostr(iAlt)
      iAlt = fix(ask('which altitude to plot',iAlt))

      if (n_elements(iAlt2) eq 0) then iAlt2='0' else iAlt2=tostr(iAlt2)
      iAlt2 = fix(ask('which altitude to plot',iAlt2))

      if (n_elements(lat0) eq 0) then lat0 = '65.0' else lat0 = string(lat0)
      lat0 = float(ask('which lat to center on',lat0))

      if (n_elements(dlat) eq 0) then dlat = '45.0' else dlat = string(dlat)
      dlat = float(ask('delta-latitude',dlat))

      if (n_elements(IsLT) eq 0) then IsLT = '1' else IsLT = string(IsLT)
      IsLT = fix(ask('whether to center on a given local time (0=Lon,1=LT)',IsLT))

      if (IsLT) then cLon = 'local time' else cLon = 'longitude'

      if (n_elements(lon0) eq 0) then lon0 = '200.0' else lon0 = string(lon0)
      lon0 = float(ask('which '+cLon+' to center on',lon0))

      if (n_elements(dlon) eq 0) then dlon = '45.0' else dlon = string(dlon)
      dlon = float(ask('delta-longitude',dlon))

   endif

   if (IsLT) then begin
      ut   = float(iTime(3)) + float(iTime(4))/60.0 + float(iTime(5))/3600.0
      diff = (lon0-ut+24.0) mod 24.0
      lo0 = diff * 360.0/24.0
   endif else lo0 = lon0

   if (nFiles gt 1) then begin
      p = strpos(psfile,'.ps')
      if (p gt -1) then psfile = strmid(psfile,0,p)
      psfile_final = psfile+'_'+chopr('000'+tostr(iFile),4)+'.ps'
   endif else begin
      psfile_final = psfile
   endelse

   print, psfile_final
   setdevice,psfile_final,'p',4,0.95

   limit = [lat0-dlat/2, lo0-dlon/2, $
            lat0+dlat/2, lo0+dlon/2]

   ppp = 2
   space = 0.01
   pos_space, ppp, space, sizes
   get_position, ppp, space, sizes, 0, p2
   get_position, ppp, space, sizes, 1, p1

;   p1 = [0.1, 0.48, 0.9, 0.95]
;   p2 = [0.1, 0.00, 0.9, 0.47]

   !p.position = p1

   map_set, lat0, lo0, limit = limit

   iEast = 16
   ve = reform(data(iEast  ,*,*,iAlt))
   vn = reform(data(iEast+1,*,*,iAlt))

   v = reform(data(sel,*,*,iAlt))
   a = reform(lat(*,*,iAlt))
   o = reform(lon(*,*,iAlt))

   l = where(a ge lat0-dlat/2 and a le lat0+dlat/2 and $
             o ge lo0 -dlon/2 and o le lo0 +dlon/2,c)

   if (c gt 0) then begin

      vl = v(l)
      al = a(l)
      ol = o(l)

      colorbar = 'all'
      if (float(smini) eq 0) then mini = min(vl) else mini = float(smini)
      if (float(smaxi) eq 0) then maxi = max(vl) else maxi = float(smaxi)
      if (float(smini) eq 0 and mini lt 0.0) then begin
         maxi = max([maxi,abs(mini)])
         mini = -maxi
         colorbar = 'mid'
      endif

      levels = findgen(31)/30.0 * (maxi-mini) + mini

      loc2 = where(vl lt levels(1),c2)
      if (c2 gt 0) then vl(loc2) = levels(1)

      loc2 = where(vl gt levels(29),c2)
      if (c2 gt 0) then vl(loc2) = levels(29)

      makect, colorbar
      contour, vl, ol, al, nlevels=31, /over, /cell, /irregular, levels =levels

      for i=0L,c-1 do begin

         x0 = o(l(i))
         y0 = a(l(i))

         dy = vn(l(i))/100.0
         dx = ve(l(i))/100.0/cos((a(l(i))+dy)*!dtor)

         if (y0+dy gt 90.0) then begin
            y = 180.0 - (y0+dy)
            dy = y-y0
            dx = dx + 180.0
         endif

         plots, [x0], [y0], psym = 4, symsize=0.5
         plots, [x0,x0+dx], [y0,y0+dy], thick = 3

         length = (dx^2+dy^2)^0.5
         t2     = asin(dy/length)
         if (dx lt 0.0) then t2 = !pi - t2

         x1 = 0.5*length*cos(t2 + 15.0*!pi/180.0)
         y1 = 0.5*length*sin(t2 + 15.0*!pi/180.0)

         x2 = 0.5*length*cos(t2 - 15.0*!pi/180.0)
         y2 = 0.5*length*sin(t2 - 15.0*!pi/180.0)

         plots, [x0+dx, x0+dx-x1], [y0+dy,y0+dy-y1]
         plots, [x0+dx, x0+dx-x2], [y0+dy,y0+dy-y2]

      endfor

      latmax = max(a(l))
      latmin = min(a(l))
      lonmax = max(o(l))
      lonmin = min(o(l))
      x0 = lonmax + 0.1 * (lonmax-lonmin)
      y0 = latmax - 0.3 * (latmax-latmin)

      dy = -200.0/100.0
      dx = 0.0

      if (y0+dy gt 90.0) then begin
         y = 180.0 - (y0+dy)
         dy = y-y0
         dx = dx + 180.0
      endif

      plots, [x0], [y0], psym = 4, symsize=0.5
      plots, [x0,x0+dx], [y0,y0+dy], thick = 3

      length = (dx^2+dy^2)^0.5
      t2     = asin(dy/length)
      if (dx lt 0.0) then t2 = !pi - t2

      x1 = 0.5*length*cos(t2 + 15.0*!pi/180.0)
      y1 = 0.5*length*sin(t2 + 15.0*!pi/180.0)

      x2 = 0.5*length*cos(t2 - 15.0*!pi/180.0)
      y2 = 0.5*length*sin(t2 - 15.0*!pi/180.0)

      plots, [x0+dx, x0+dx-x1], [y0+dy,y0+dy-y1]
      plots, [x0+dx, x0+dx-x2], [y0+dy,y0+dy-y2]

      xyouts, x0+0.2, y0, '200 m/s', charsize = 0.8

   endif
             
   map_set, lat0, lo0, limit = limit, /grid, /continents, /usa, /hires, /noerase

   centerlats = [ 18.344167, -6.88, -7.38,-11.96, 37.75, 42.619,-14.97, 40.13, 42.27, 35.20]
   centerlons = [-66.752778,-38.56,-36.52,-76.86,-84.29,-71.491,-74.89,-88.20,-83.75,-82.85]+360.0

   oplot, centerlons, centerlats, psym = 5, thick = 4

   salt  = tostr(alt(0,0,iAlt))+' km'
   c_a_to_s, itime, stime

   !p.position = -1

   xyouts, p1(2), p1(3)-0.025, salt, alignment = 1.0, /norm
   xyouts, p2(0), p2(3)+0.01, stime, alignment = 0.0, /norm

   pos = [p1(2)+0.01, p1(1), p1(2)+0.03, p1(3)]
   maxmin = mm(levels)
   title = vars(sel)
   plotct,254,pos,maxmin,title,/right

   !p.position = p2

   map_set, lat0, lo0, limit = limit, /noerase

;   iEast = 36
   ve = reform(data(iEast  ,*,*,iAlt2))
   vn = reform(data(iEast+1,*,*,iAlt2))

   v = reform(data(sel2,*,*,iAlt2))
   a = reform(lat(*,*,iAlt2))
   o = reform(lon(*,*,iAlt2))

   l = where(a ge lat0-dlat/2 and a le lat0+dlat/2 and $
             o ge lo0 -dlon/2 and o le lo0 +dlon/2,c)

   if (c gt 0) then begin

      vl = v(l)
      al = a(l)
      ol = o(l)

      colorbar = 'all'
      if (float(smini2) eq 0) then mini2 = min(vl) else mini2 = float(smini2)
      if (float(smaxi2) eq 0) then maxi2 = max(vl) else maxi2 = float(smaxi2)
      if (float(smini2) eq 0 and mini2 lt 0.0) then begin
         maxi2 = max([maxi2,abs(mini2)])
         mini2 = -maxi2
         colorbar = 'mid'
      endif
      levels = findgen(31)/30.0 * (maxi2-mini2) + mini2

      loc2 = where(vl lt levels(1),c2)
      if (c2 gt 0) then vl(loc2) = levels(1)

      loc2 = where(vl gt levels(29),c2)
      if (c2 gt 0) then vl(loc2) = levels(29)

      makect, colorbar
      contour, vl, ol, al, nlevels=31, /over, /cell, /irregular, levels = levels

      for i=0L,c-1 do begin

         x0 = o(l(i))
         y0 = a(l(i))

         dy = vn(l(i))/100.0
         dx = ve(l(i))/100.0/cos((a(l(i))+dy)*!dtor)

         if (y0+dy gt 90.0) then begin
            y = 180.0 - (y0+dy)
            dy = y-y0
            dx = dx + 180.0
         endif

         plots, [x0], [y0], psym = 4, symsize=0.5
         plots, [x0,x0+dx], [y0,y0+dy], thick = 3

         length = (dx^2+dy^2)^0.5
         t2     = asin(dy/length)
         if (dx lt 0.0) then t2 = !pi - t2

         x1 = 0.5*length*cos(t2 + 15.0*!pi/180.0)
         y1 = 0.5*length*sin(t2 + 15.0*!pi/180.0)

         x2 = 0.5*length*cos(t2 - 15.0*!pi/180.0)
         y2 = 0.5*length*sin(t2 - 15.0*!pi/180.0)

         plots, [x0+dx, x0+dx-x1], [y0+dy,y0+dy-y1]
         plots, [x0+dx, x0+dx-x2], [y0+dy,y0+dy-y2]

      endfor

      latmax = max(a(l))
      latmin = min(a(l))
      lonmax = max(o(l))
      lonmin = min(o(l))
      x0 = lonmax + 0.1 * (lonmax-lonmin)
      y0 = latmax - 0.3 * (latmax-latmin)

      dy = -200.0/100.0
      dx = 0.0

      if (y0+dy gt 90.0) then begin
         y = 180.0 - (y0+dy)
         dy = y-y0
         dx = dx + 180.0
      endif

      plots, [x0], [y0], psym = 4, symsize=0.5
      plots, [x0,x0+dx], [y0,y0+dy], thick = 3

      length = (dx^2+dy^2)^0.5
      t2     = asin(dy/length)
      if (dx lt 0.0) then t2 = !pi - t2

      x1 = 0.5*length*cos(t2 + 15.0*!pi/180.0)
      y1 = 0.5*length*sin(t2 + 15.0*!pi/180.0)

      x2 = 0.5*length*cos(t2 - 15.0*!pi/180.0)
      y2 = 0.5*length*sin(t2 - 15.0*!pi/180.0)

      plots, [x0+dx, x0+dx-x1], [y0+dy,y0+dy-y1]
      plots, [x0+dx, x0+dx-x2], [y0+dy,y0+dy-y2]

      xyouts, x0+0.2, y0, '200 m/s', charsize = 0.8


      ; ---------------------------------
      ; plot ion velocities
      ; ---------------------------------

      iEast = 36
      ve = reform(data(iEast  ,*,*,iAlt2))
      vn = reform(data(iEast+1,*,*,iAlt2))

      v = reform(data(sel2,*,*,iAlt2))
      a = reform(lat(*,*,iAlt2))
      o = reform(lon(*,*,iAlt2))

      l = where(a ge lat0-dlat/2 and a le lat0+dlat/2 and $
                o ge lo0-dlon/2 and o le lo0 +dlon/2,c)

      if (c gt 0) then begin

         vl = v(l)
         al = a(l)
         ol = o(l)

         for i=0L,c-1 do begin

            x0 = o(l(i))
            y0 = a(l(i))

            dy = vn(l(i))/100.0
            dx = ve(l(i))/100.0/cos((a(l(i))+dy)*!dtor)

            if (y0+dy gt 90.0) then begin
               y = 180.0 - (y0+dy)
               dy = y-y0
               dx = dx + 180.0
            endif

            plots, [x0,x0+dx], [y0,y0+dy], thick = 6, color = 100

            length = (dx^2+dy^2)^0.5
            t2     = asin(dy/length)
            if (dx lt 0.0) then t2 = !pi - t2

            x1 = 0.5*length*cos(t2 + 15.0*!pi/180.0)
            y1 = 0.5*length*sin(t2 + 15.0*!pi/180.0)

            x2 = 0.5*length*cos(t2 - 15.0*!pi/180.0)
            y2 = 0.5*length*sin(t2 - 15.0*!pi/180.0)

            plots, [x0+dx, x0+dx-x1], [y0+dy,y0+dy-y1], thick = 6, color = 100
            plots, [x0+dx, x0+dx-x2], [y0+dy,y0+dy-y2], thick = 6, color = 100

         endfor

      endif


   endif
             
   map_set, lat0, lo0, limit = limit, /grid, /continents, /usa, /hires, /noerase

   centerlats = [ 18.344167, -6.88, -7.38,-11.96, 37.75, 42.619,-14.97, 40.13, 42.27, 35.20]
   centerlons = [-66.752778,-38.56,-36.52,-76.86,-84.29,-71.491,-74.89,-88.20,-83.75,-82.85]+360.0

   oplot, centerlons, centerlats, psym = 5, thick = 4

   salt  = tostr(alt(0,0,iAlt2))+' km'
   c_a_to_s, itime, stime

   !p.position = -1

   pos = [p2(2)+0.01, p2(1), p2(2)+0.03, p2(3)]
   maxmin = mm(levels)
   title = vars(sel2)
   plotct,254,pos,maxmin,title,/right

   if (IsLT) then begin
      xyouts, (p2(0)+p2(2))/2.0, p2(1)-0.02, tostr(lon0)+' LT', /norm, align=0.5
      ll = lon0-dlon/15/2
      lu = lon0+dlon/15/2
      xyouts, p2(0), p2(1)-0.02, tostr(ll)+' LT', /norm, align=0.5
      xyouts, p2(2), p2(1)-0.02, tostr(lu)+' LT', /norm, align=0.5
   endif

   xyouts, p2(2), p2(3)-0.025, salt, alignment = 1.0, /norm

   closedevice

endfor

end
