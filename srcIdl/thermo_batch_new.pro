
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

psfile = filelist(0)+'.ps'
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
        if (n_elements(sel) eq 0) then sel = '9' else sel = tostr(sel)
        sel = fix(ask('which var to plot',sel))

        if (n_elements(plotlog) eq 0) then plotlog='n' $
        else if (plotlog) then plotlog='y' else plotlog='n' 
        plotlog = ask('whether you want log or not (y/n)',plotlog)
        if (strpos(plotlog,'y') eq 0) then plotlog = 1 else plotlog = 0

        if (n_elements(IsLinePlot) eq 0) then IsLinePlot='c' $
        else if (IsLinePlot) then IsLinePlot='l' else IsLinePlot='c'
        IsLinePlot = mklower(ask('whether you would like a line plot (l) or contour plot (c)',IsLinePlot))

        if (strpos(IsLinePlot,'c') eq -1) then IsLinePlot = 1 $
        else IsLinePlot = 0

        if (IsLinePlot) then begin

           print, 'What 1D would you like to plot?'
           print, '1. Variable vs Altitude'
           print, '2. Variable vs Longitude'
           print, '3. Variable vs Latitude'
           if (n_elements(iLine) eq 0) then iLine='1' else iLine=tostr(iLine)
           iLine = fix(ask('type of plot to make',iLine))

           if (iLine eq 1) then begin
              if (n_elements(IsGlobalAverage) eq 0) then IsGlobalAverage='0'
              IsGlobalAverage = fix(ask('whether you want a global average (0-no,1-yes)','0'))
           endif

           if (not IsGlobalAverage) then begin

              if (iLine ne 2) then begin
                 for i=0,nlons-1 do print, tostr(i)+'. '+string(lon(i,2,2))
                 if (n_elements(iLon) eq 0) then iLon='0' else iLon=tostr(iLon)
                 iLon = fix(ask('which longitude to plot',iLon))
              endif
              
              if (iLine ne 3) then begin
                 for i=0,nlats-1 do print, tostr(i)+'. '+string(lat(2,i,2))
                 if (n_elements(iLat) eq 0) then iLat='0' else iLat=tostr(iLat)
                 iLat = fix(ask('which latitude to plot',iLat))
              endif

              if (iLine ne 1) then begin
                 for i=0,nalts-1 do print, tostr(i)+'. '+string(alt(2,2,i))
                 if (n_elements(iAlt) eq 0) then iAlt='0' else iAlt=tostr(iAlt)
                 iAlt = fix(ask('which altitude to plot',iLat))
              endif

              title = 'Location : '

              if (iLine eq 1) then begin
                 yLine = reform(data(2,iLon,iLat,*))/1000.0
                 xLine = reform(data(sel,iLon,iLat,*))
                 yTitle = 'Altitude (km)'
                 xTitle = strcompress(vars(sel),/remove)
                 if (plotlog) then begin
                    xLine = alog10(xLine)
                    xtitle = 'log('+xtitle+')'
                 endif
                 title = title+tostr(data(0,iLon,iLat,0)/!dtor)+'(deg Lon)'+$
                         ' and '+tostr(data(1,iLon,iLat,0)/!dtor)+'(deg Lat)'
              endif
              if (iLine eq 2) then begin
                 xLine = reform(data(0,*,iLat,iAlt))/!dtor
                 yLine = reform(data(sel,*,iLat,iAlt))
                 yTitle = strcompress(vars(sel),/remove)
                 xTitle = 'Longitude (deg)'
                 if (plotlog) then begin
                    yLine = alog10(yLine)
                    yTitle = 'log('+ytitle+')'
                 endif
                 title = title+tostr(data(1,0, iLat,iAlt)/!dtor)+'(deg Lat)'+$
                         ' and '+tostr(data(2,0, iLat,iAlt)/1000.0)+'(km Alt)'
              endif
              if (iLine eq 3) then begin
                 xLine = reform(data(1,iLon,*,iAlt))/!dtor
                 yLine = reform(data(sel,iLon,*,iAlt))
                 xTitle = strcompress(vars(sel),/remove)
                 yTitle = 'Latitude (deg)'
                 if (plotlog) then begin
                    yLine = alog10(yLine)
                    yTitle = 'log('+ytitle+')'
                 endif
                 title = title+tostr(data(0,iLon,0,iAlt)/!dtor)+'(deg Lon)'+$
                         ' and '+tostr(data(2,iLon,0,iAlt)/1000.0)+'(km Alt)'
              endif

           endif else begin

              dLon = data(0,1,0,0) - data(0,0,0,0)
              dLat = data(1,0,1,0) - data(1,0,0,0)
              r = cRe_ + reform(data(2,*,*,*))
              area = dLon * dLat * r^2 * cos(reform(data(1,*,*,*)))

              print, "here"
              xLine = fltarr(nAlts)
              for i=0,nAlts-1 do begin
                 Var = reform(data(sel,2:nLons-3,2:nLats-3,i))
                 are = reform(area(2:nLons-3,2:nLats-3,i))
                 xLine(i) = total(Var*Are)/total(Are)
              endfor
              yLine = reform(data(2,0,0,*))/1000.0
              yTitle = 'Altitude (km)'
              xTitle = strcompress(vars(sel),/remove)
              if (plotlog) then begin
                 xLine = alog10(yLine)
                 xTitle = 'log('+xtitle+')'
              endif

              Title = 'Global Average'

           endelse

        endif else begin

           if (n_elements(iSecondVar) eq 0) then iSecondVar = '-1' $
           else iSecondVar = tostr(iSecondVar)
           iSecondVar = $
              fix(ask('second var to plot as line contour (-1 for none)',$
                     iSecondVar))
        
           print, '1. Constant Altitude Plot'
           print, '2. Constant Longitude Plot (or Zonal Average)'
           print, '3. Constant Latitude Plot'
           slice = fix(ask('type of plot to make','1'))

           cnt1 = 0
           cnt2 = 0
           cnt3 = 0

;cnt1 is a lat/lon plot
           if (slice eq 1) then cnt1 = 1

;cnt1 is a lat/alt plot
           if (slice eq 2) then cnt3 = 1

;cnt1 is a lon/alt plot
           if (slice eq 3) then cnt2 = 1

           polar = 0
           npolar = 1
           MinLat = 50.0

           if (slice eq 1) then begin
              for i=0,nalts-1 do print, tostr(i)+'. '+string(alt(2,2,i))
              selset = fix(ask('which altitude to plot','0'))

              polar = fix(ask('polar (1) or non-polar (0)','0'))

              if (polar) then begin
                 npolar = fix(ask('North (1) or South (0)','1'))
                 MinLat = abs(float(ask('minimum latitude to plot','50.0')))
              endif

           endif

           if (slice eq 2) then begin
              if (n_elements(IsZonalAverage) eq 0) then IsZonalAverage='0' $
              else IsZonalAverage = tostr(IsZonalAverage)
              IsZonalAverage = $
                 fix(ask('whether you want a Zonal Average (0=n,1=y)',$
                         IsZonalAverage))
              if (not IsZonalAverage) then begin
                 for i=0,nlons-1 do print, tostr(i)+'. '+string(lon(i,2,2))
                 selset = fix(ask('which longitude to plot','0'))
              endif else selset = 0
           endif

           if (slice eq 3) then begin
              for i=0,nlats-1 do print, tostr(i)+'. '+string(lat(2,i,2))
              selset = fix(ask('which latitude to plot','0'))
           endif

           smini = ask('minimum (0.0 for automatic)','0.0')
           smaxi = ask('maximum (0.0 for automatic)','0.0')

           if (not IsZonalAverage) then begin
              if (n_elements(plotVector) eq 0) then plotvector='y' $
              else if (plotvector) then plotvector='y' else plotvector='n'
              plotvector=ask('whether you want vectors or not (y/n)',plotvector)
              if strpos(plotvector,'y') eq 0 then plotvector=1 $
              else plotvector = 0
           endif else plotvector = 0

           if (plotvector) then begin

              PlotNeutrals = fix(ask('plot neutral winds (1) or ions (0)','1'))
; vi_cnt is whether to plot vectors of Vi
              vi_cnt = 1-PlotNeutrals

; vn_cnt is whether to plot vectors of Vn
              vn_cnt = PlotNeutrals

              print,'-1  : automatic selection'
              factors = [1.0, 5.0, 10.0, 20.0, 25.0, $
                         50.0, 75.0, 100.0, 150.0, 200.0,300.0]
              nfacs = n_elements(factors)
              for i=0,nfacs-1 do print, tostr(i+1)+'. '+string(factors(i)*10.0)
              vector_factor = fix(ask('velocity factor','-1'))
           endif else vector_factor = 0

; cursor position variables, which don't matter at this point
           cursor_x = 0.0
           cursor_y = 0.0
           strx = '0.0'
           stry = '0.0'

; yes is whether ghostcells are plotted or not:
           yes = 0
           no  = 1

; yeslog is whether variable should be logged or not:
           if (plotlog) then begin 
              yeslog = 1
              nolog  = 0
           endif else begin
              yeslog = 0
              nolog = 1
           endelse

; yeswrite_cnt is whether we have to output to a ps file or not.
           yeswrite_cnt = 1

; showgridyes says whether to plot the grid or not.
           showgridyes = 0

;plotvectoryes says whether to plot vectors or not
           plotvectoryes = plotvector

; number of points to skip when plotting vectors:
           step = 1

           cursor_cnt = 0

           xrange = [0.0,0.0]

           yrange = [0.0,0.0]

        endelse

     endif

    if (nFiles gt 1) then begin
       p = strpos(psfile,'.ps')
       if (p gt -1) then psfile = strmid(psfile,0,p)
       psfile_final = psfile+'_'+chopr('000'+tostr(iFile),4)+'.ps'
    endif else begin
       psfile_final = psfile
    endelse

    if (IsLinePlot) then begin

       if (not IsGlobalAverage) then begin

          title = 'Location : '

          if (iLine eq 1) then begin
             yLine = reform(data(2,iLon,iLat,*))/1000.0
             xLine = reform(data(sel,iLon,iLat,*))
             yTitle = 'Altitude (km)'
             xTitle = strcompress(vars(sel),/remove)
             if (plotlog) then begin
                xLine = alog10(xLine)
                xtitle = 'log('+xtitle+')'
             endif
             title = title+tostr(data(0,iLon,iLat,0)/!dtor)+'(deg Lon)'+$
                     ' and '+tostr(data(1,iLon,iLat,0)/!dtor)+'(deg Lat)'
          endif
          if (iLine eq 2) then begin
             xLine = reform(data(0,*,iLat,iAlt))/!dtor
             yLine = reform(data(sel,*,iLat,iAlt))
             yTitle = strcompress(vars(sel),/remove)
             xTitle = 'Longitude (deg)'
             if (plotlog) then begin
                yLine = alog10(yLine)
                yTitle = 'log('+ytitle+')'
             endif
             title = title+tostr(data(1,0, iLat,iAlt)/!dtor)+'(deg Lat)'+$
                     ' and '+tostr(data(2,0, iLat,iAlt)/1000.0)+'(km Alt)'
          endif
          if (iLine eq 3) then begin
             xLine = reform(data(1,iLon,*,iAlt))/!dtor
             yLine = reform(data(sel,iLon,*,iAlt))
             xTitle = strcompress(vars(sel),/remove)
             yTitle = 'Latitude (deg)'
             if (plotlog) then begin
                yLine = alog10(yLine)
                yTitle = 'log('+ytitle+')'
             endif
             title = title+tostr(data(0,iLon,0,iAlt)/!dtor)+'(deg Lon)'+$
                     ' and '+tostr(data(2,iLon,0,iAlt)/1000.0)+'(km Alt)'
          endif

       endif else begin

          dLon = data(0,1,0,0) - data(0,0,0,0)
          dLat = data(1,0,1,0) - data(1,0,0,0)
          r = cRe_ + reform(data(2,*,*,*))
          area = dLon * dLat * r^2 * cos(reform(data(1,*,*,*)))

          xLine = fltarr(nAlts)
          for i=0,nAlts-1 do begin
             Var = reform(data(sel,2:nLons-3,2:nLats-3,i))
             are = reform(area(2:nLons-3,2:nLats-3,i))
             xLine(i) = total(Var*Are)/total(Are)
          endfor
          yLine = reform(data(2,0,0,*))/1000.0
          yTitle = 'Altitude (km)'
          xTitle = strcompress(vars(sel),/remove)
          if (plotlog) then begin
             xLine = alog10(yLine)
             xTitle = 'log('+xtitle+')'
          endif

          Title = 'Global Average'

       endelse

       setdevice,psfile_final,'p',4,0.95

       r = max(xline)-min(xline)
       if (r eq 0) then r = 1.0
       xrange = [min(xline)-0.1*r,max(xline)+0.1*r]
       r = max(yline)-min(yline)
       if (r eq 0) then r = 1.0
       yrange = [min(yline)-0.1*r,max(yline)+0.1*r]

       c_a_to_s, itime, stime

       plot, xline, yline, $
             thick = 3, xstyle = 1, ystyle = 1, $
             xtitle = xtitle, ytitle = ytitle, $
             pos = [0.1,0.1,0.9,0.7], $
             xrange = xrange, yrange = yrange, title = stime
       xyouts, 0.1, 0.71, /norm, title
       closedevice

    endif else begin

       smini_final = smini
       smaxi_final = smaxi

       thermo_plot_new,cursor_x,cursor_y,strx,stry,step,nvars,sel,nfiles,$
                   cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
                   1-yeslog,nalts,nlats,nlons,yeswrite_cnt,$
                   polar,npolar,MinLat,showgridyes,	  $
                   plotvectoryes,vi_cnt,vn_cnt,vector_factor,	  $
                   cursor_cnt,data,alt,lat,lon,	  $
                   xrange,yrange,selset,smini_final,smaxi_final,	  $
                   filename,vars, psfile_final, 0, 'all', itime, $
                   iSecondVar=iSecondVar, IsZonalAverage=IsZonalAverage

    endelse

endfor


end
