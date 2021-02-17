pro thermo_plot,cursor_x,cursor_y,strx,stry,step,nvars,sel,nfiles, $
                cnt1,cnt2,cnt3,ghostcells,no,yeslog,  	  $
                nolog,nalts,nlats,nlons,yeswrite_cnt,$
                polar,npolar,MinLat,showgridyes,	  $
                plotvectoryes,vi_cnt,vn_cnt,cf,	  $
                cursor_cnt,data,alt,lat,lon,	  $
                xrange,yrange,selset,smini,smaxi,	  $
                filename,vars, psfile, mars, colortable, itime, $
                plotpole = plotpole

if (n_elements(colortable) eq 0) then colortable = 'mid'
if (strlen(colortable) eq 0) then colortable = 'mid'

if (n_elements(logplot) eq 0) then logplot = yeslog

if (min(data(sel,*,*,*)) lt 0.0) then begin
  logplot = 0
  yeslog = 0
  nolog = 1
endif

if (n_elements(iTime) eq 0) then begin

    if (strpos(filename,"save") gt 0) then begin

        fn = findfile(filename)
        if (strlen(fn(0)) eq 0) then begin
            print, "Bad filename : ", filename
            stop
        endif else filename = fn(0)
        
        l1 = strpos(filename,'.save')
        fn2 = strmid(filename,0,l1)
        len = strlen(fn2)
        l2 = l1-1
        while (strpos(strmid(fn2,l2,len),'.') eq -1) do l2 = l2 - 1
        l = l2 - 13
        year = fix(strmid(filename,l, 2))
        mont = fix(strmid(filename,l+2, 2))
        day  = fix(strmid(filename,l+4, 2))
        hour = float(strmid(filename, l+7, 2))
        minu = float(strmid(filename,l+9, 2))
        seco = float(strmid(filename,l+11, 2))
    endif else begin
        year = fix(strmid(filename,07, 2))
        mont = fix(strmid(filename,09, 2))
        day  = fix(strmid(filename,11, 2))
        hour = float(strmid(filename,14, 2))
        minu = float(strmid(filename,16, 2))
        seco = float(strmid(filename,18, 2))
    endelse

    itime = [year,mont,day,fix(hour),fix(minu),fix(seco)]

endif

c_a_to_s, itime, stime

ut = itime(3) + itime(4)/60.0 + itime(5)/3600.0
if (polar) then utrot = ut * 15.0 else utrot = 0.0

if (cnt1 eq 0 and polar) then polar = 0

;According to the button user selected, get the values for plotting
if cnt1 eq 1 then begin
    if (polar) then MinLat = MinLat else MinLat = -1000.0
    if not (polar) then mr = 1090
    if (polar) then begin
        mr = 90.0 - abs(MinLat)
    endif
;    if (polar) and not (npolar) then begin
;        mr = 90.0 + MinLat
;    endif
    
    nLons = n_elements(lon(*,0,0))

    if (polar) then begin
        if (npolar) then begin
            loc = where(lat(0,*,0) ge abs(MinLat) and abs(lat(0,*,0)) lt 90.0)
        endif else begin
            loc = where(lat(0,*,0) le -abs(MinLat) and abs(lat(0,*,0)) lt 90.0)
        endelse
    endif else begin
        loc = where(lat(0,*,0) ge -200 and abs(lat(0,*,0)) lt 200.0)
    endelse

    if (polar) then begin
       if (lon(0,0,0) lt 0.0) then begin
          datatoplot=reform(data(sel,2:nLons-2,loc,selset))
          datatoplot(nLons-4,*) = datatoplot(0,*)
          iLonS = 2
          iLonE = nLons-2
       endif else begin
          iLonS = 0
          iLonE = nLons-1
          datatoplot=reform(data(sel,0:nLons-1,loc,selset))
          ;datatoplot(nLons-1,*) = datatoplot(0,*)
       endelse
    endif else datatoplot=reform(data(sel,*,loc,selset))

    maxi=max(datatoplot)
    mini=min(datatoplot)

    if (polar) then begin

        if (npolar) then begin
            x = reform( (90.0 - lat(iLonS:iLonE,loc,selset)) * $
                        cos((lon(iLonS:iLonE,loc,selset)+utrot)*!pi/180.0 - !pi/2.0))
            y = reform( (90.0 - lat(iLonS:iLonE,loc,selset)) * $
                        sin((lon(iLonS:iLonE,loc,selset)+utrot)*!pi/180.0 - !pi/2.0))
        endif else begin
            x = reform( (90.0 + lat(iLonS:iLonE,loc,selset)) * $
                        cos((lon(iLonS:iLonE,loc,selset)+utrot)*!pi/180.0 - !pi/2.0))
            y = reform( (90.0 + lat(iLonS:iLonE,loc,selset)) * $
                        sin((lon(iLonS:iLonE,loc,selset)+utrot)*!pi/180.0 - !pi/2.0))
        endelse

        xrange = [-mr,mr]
        yrange = [-mr,mr]
        xtitle=' '
        ytitle=' '

    endif else begin

        x=reform(lon(*,loc,selset))
        y=reform(lat(*,loc,selset))
        if ghostcells eq 1 then begin
            xrange=mm(lon)
            yrange=mm(lat)
        endif
        if ghostcells eq 0 then begin
            xrange=[0,360]
            yrange=[-90,90]
        endif
        xtitle='Longitude (deg)'
        ytitle='Latitude (deg)'
    endelse
    ygrids = n_elements(loc)
    xgrids = nlons
    location = string(alt(0,0,selset),format='(f5.1)')+' km Altitude'
endif

if cnt2 eq 1 then begin
    datatoplot=reform(data(sel,*,selset,*))
    maxi=max(datatoplot)
    mini=min(datatoplot)
    x=reform(lon(*,selset,*))
    y=reform(alt(*,selset,*))
    location = string(lat(0,selset,0),format='(f5.1)')+' deg Latitude'
    xtitle='Longitude (deg)'
    ytitle='Altitude (km)'
    if ghostcells eq 1 then begin
        xrange=mm(lon)
        yrange=mm(alt)
    endif
    if ghostcells eq 0 then begin
        backup_xrange=mm(lon)
        backup_yrange=mm(alt)
        default_xrange=[0,360]
        default_yrange=mm(alt)
;If out of range then use 'mm' to set xrange and yrange values.
;Else use default values.
        if (backup_xrange[0] lt default_xrange[0]) $
          and (backup_xrange[1] gt default_xrange[1]) then begin
            xrange=mm(lon)
            yrange=mm(alt)
        endif else begin
            xrange=[0,360]
            yrange=mm(alt)
        endelse
    endif

    if (min(lon) gt 0 and max(lon) lt 0) then xrange = mm(lon) $
    else xrange = [0,360]

    ygrids=nalts
    xgrids=nlons
endif

if cnt3 eq 1 then begin
    datatoplot=reform(data(sel,selset,*,*))
    maxi=max(datatoplot)
    mini=min(datatoplot)
    x=reform(lat(selset,*,*))
    y=reform(alt(selset,*,*))
    location = string(lon(selset,0,0),format='(f5.1)')+' deg Longitude'
    xtitle='Latitude (deg)'
    ytitle='Altitude (km)'
    if ghostcells eq 1 then begin
        xrange=mm(lat)
        yrange=mm(alt)
    endif
    if ghostcells eq 0 then begin
        backup_xrange=mm(lat)
        backup_yrange=mm(alt)
        default_xrange=[-90,90]
        default_yrange=mm(alt)
;If out of range then use 'mm' to set xrange and yrange values.
;Else use default values.
        if (backup_xrange[0] lt default_xrange[0]) $
          and (backup_xrange[1] gt default_xrange[1]) then begin
            xrange=mm(lat)
            yrange=mm(alt)
        endif else begin
            xrange=[-90,90]
            yrange=mm(alt)
        endelse
    endif
    ygrids=nalts
    xgrids=nlats

    if (min(lat) gt -90.0 and max(lat) lt 90.0) then lat = mm(lat) $
    else xrange = [-90,90]

endif

  ;Calculate the xld, yld according to the cursor position user set.
  ;Calculate and get the array will be plotted.  

;xld = x(*,0)
;yld = y(0,*)
;dist_x=abs(xld-cursor_x)
;dist_y=abs(yld-cursor_y)
;locx=where(dist_x eq min(dist_x))
;locy=where(dist_y eq min(dist_y))
;datald=reform(data(sel,*,locx,locy))

IsManual = 0

print, mini, maxi

if n_elements(smini) eq 0 then smini = '0.0'
if n_elements(smaxi) eq 0 then smaxi = '0.0'

if (float(smini) ne 0 or float(smaxi) ne 0) then begin
    mini = float(smini)
    maxi = float(smaxi)
    mini = mini(0)
    maxi = maxi(0)
    IsManual = 1
endif else begin

    if (logplot) then begin
        if (maxi gt 0.0) then maxi = alog10(maxi)
        if (mini gt 0.0) then mini = alog10(mini)
        if (maxi-mini gt 8) then begin
            mini = maxi-8
            print, "Limiting minimum..."
        endif
    endif 

    mini = mini(0)
    maxi = maxi(0)
    r = (maxi-mini)*0.05
    mini = mini - r
    maxi = maxi + r

endelse

if (abs(mini-maxi) lt 0.001*abs(maxi)) or (abs(mini-maxi) eq 0.0) then begin
   if (not IsManual) then maxi=mini*1.01+1.0
endif

if (mini lt -0.1*maxi and maxi gt 0.0 and not IsManual) then begin

  maxi = max([abs(mini),maxi])
  mini = -maxi

endif

levels = findgen(31)/30.0*(maxi-mini) + mini

loc = where(datatoplot lt levels(1), count)
if (count gt 0) then datatoplot(loc) = levels(1)

 ; Check if user wants to write the result to a file
 ;If user wanted then setdevice. 

if yeswrite_cnt eq 1 then begin

    if (strlen(psfile) eq 0) then psfile = filename+'.ps'
    setdevice,psfile,'l',4,0.95

endif

plotdumb

variable = strcompress(vars(sel),/remove)

;  makect,'wyr'

; makect,'mid'
makect, colortable

clevels = findgen(31)/30 * 253.0 + 1.0

if (polar) then begin
    xstyle = 5 
    ystyle = 5 
endif else begin
    xstyle = 1
    ystyle = 1
endelse

if (not polar and cnt1) then ppp = 2 else ppp = 1

space = 0.075
pos_space, ppp, space, sizes, ny = 1
get_position, ppp, space, sizes, 0, pos

if (not polar) then begin
    if (cnt1) then begin
        r = pos(2) - pos(0)
        pos(2) = pos(0) + r*2.0
    endif else begin
        get_position, ppp, space, sizes, 0, pos, /rect
        pos(1) = pos(1) + space
        pos(3) = pos(3) - space*2.0
        pos(0) = pos(0) + space*1.0
        pos(2) = pos(2) - space*1.0
    endelse
endif

if (n_elements(plotpole) gt 0) then begin
    apexfile = getenv('IDL_EXTRAS')+'apex'+tostr(itime(0))+'.dat'
    print, 'reading file : ',apexfile
    readapex, apexfile, apex
    ln = where(apex.alats eq max(apex.alats))
    ls = where(apex.alats eq min(apex.alats))
    northpolelat = apex.glats(ln)
    northpolelon = apex.glons(ln)
    southpolelat = apex.glats(ls)
    southpolelon = apex.glons(ls)
endif

;If user DOES NOT do the Ghost Cells Contour. 
if ghostcells eq 0 then begin

    ;If user DO NOT do plot log. 

    if (logplot) then begin
        loc = where(datatoplot lt max(datatoplot)*1e-8,count)
        if (count gt 0) then datatoplot(loc) = max(datatoplot)*1e-8
        datatoplot = alog10(datatoplot)
    endif

    if (cnt1) then begin 

        if (not polar) then begin
            locx = where(x(*,0) ge   0.0 and x(*,0) le 360.0,nx)
            locy = where(y(0,*) ge -90.0 and y(0,*) le  90.0,ny)
            d2 = fltarr(nx,ny)
            x2 = fltarr(nx,ny)
            y2 = fltarr(nx,ny)
            for i=nx/2, nx-1 do begin
                d2(i-nx/2,0:ny-1)  = datatoplot(locx(i),locy)
                x2(i-nx/2,0:ny-1)  = x(locx(i),locy)
                y2(i-nx/2,0:ny-1)  = y(locx(i),locy)
                d2(i,0:ny-1)  = datatoplot(locx(i-nx/2),locy)
                x2(i,0:ny-1)  = x(locx(i-nx/2),locy)
                y2(i,0:ny-1)  = y(locx(i-nx/2),locy)
            endfor
        endif else begin
            d2 = datatoplot
            x2 = x
            y2 = y
            ny = n_elements(y2(0,*))
            nx = n_elements(x2(0,*))
        endelse

        if (not polar) then begin

            !p.position = pos
            map_set, title=variable+' at '+location+' at '+$
              strmid(stime,0,15)+' UT'

        endif else begin

            ;-------------------------------------
            ; polar plot
            ;-------------------------------------

            plot, [-mr, mr], [-mr, mr], pos=pos, $
              xstyle=5, ystyle=5,/nodata,/noerase, $
              title=variable+' at '+location+' '+stime

        endelse


    endif else begin

        d2 = datatoplot
        x2 = x
        y2 = y

        plot, mm(x2), mm(y2), /nodata, xstyle = 1, ystyle=1,$
          /noerase, pos = pos, $
          title=variable+' at '+location+' '+stime, $
          xrange=xrange,yrange=yrange

    endelse

;    endif else begin 

;    endelse

    if (cnt1) then begin

        loc = where(d2 gt max(levels(n_elements(levels)-2)),count)
        if (count gt 0) then d2(loc) = max(levels(n_elements(levels)-2))

        contour,d2, x2, y2,POS=pos,$
          levels=levels,xstyle=xstyle,ystyle=ystyle,$
          xrange=xrange,yrange=yrange,$
          c_colors=clevels,$
          xtitle=xtitle,ytitle=ytitle,/cell_fill,/over

    endif else begin
        contour,d2, x2, y2,POS=pos,$
          levels=levels,xstyle=xstyle,ystyle=ystyle,$
          xrange=xrange,yrange=yrange,$
          c_colors=clevels,$
          xtitle=xtitle,ytitle=ytitle,/cell_fill,/noerase
    endelse

    if (not polar and cnt1 and not mars) then map_continents, color = 0

    if (cnt1 and not polar and mars) then begin
        file = '/remotehome/ridley/idl/extras/marsflat.jpg'
        read_jpeg, file, image
        contour, image(2,*,*), levels = [150], pos = pos, /noerase, $
          xstyle =5, ystyle=5, color = 0, thick=1.5
        plot, [-180,180],[-90,90],xstyle=1,ystyle=1,/noerase,/nodata,pos=pos, $
          xtitle = 'Longitude', ytitle = 'Latitude'
    endif

    if (polar) then begin
        plotmlt, mr, /no06, /no12

        if (n_elements(plotpole) gt 0) then begin
            if (npolar) then begin
                polerange = 90.0-northpolelat
                polelon = northpolelon
            endif else begin
                polerange = 90.0+southpolelat
                polelon = southpolelon
            endelse
            npx = polerange * $
              cos((polelon+utrot)*!pi/180.0 - !pi/2.0)
            npy = polerange * $
              sin((polelon+utrot)*!pi/180.0 - !pi/2.0)
            oplot, [npx],[npy], psym = 4, symsize = 3, thick = 10
        endif

    endif

    if not (polar) then mr = 1090

    ;Draw grid.
    if (showgridyes eq 1) then begin
        for i=0,ygrids-1 do begin
            oplot,x(*,i),y(*,i)
        endfor
        for j=0,xgrids-1 do begin
            oplot,x(j,*),y(j,*)
        endfor
    endif
    ;If user set cursor position to plot, then do plot of datald.
    ;Else clean up the cursor text fields on the interface. 

    if cursor_cnt eq 1 then begin
        if cnt1 eq 1 then begin
            plot,datald,alt(0,0,*),xtitle='datald',ytitle='Alt. (deg)'
        endif
        if cnt2 eq 1 then begin
            plot,datald,lat(0,*,0),ystyle=1,xtitle='datald',ytitle='Lat. (deg)'
        endif
        if cnt3 eq 1 then begin
            plot,datald,lon(*,0,0),xtitle='datald',ytitle='Long. (deg)'
        endif
    endif else begin	
        txt=''
        ;widget_control,(*ptr).curx_txt,set_value=txt
        ;widget_control,(*ptr).cury_txt,set_value=txt
    endelse
endif


                                ;If user DOES the Ghost Cells Contour.
if ghostcells eq 1 then begin
    x1=min(lon)
    x2=max(lon)
    y1=min(lat)
    y2=min(lat)
                                ;If user does not want to do the Plot Log.
    if nolog eq 1 then begin
       if (not IsManual) then begin
          maxi=max(datatoplot)
          mini=min(datatoplot)
          if mini eq maxi then maxi=mini+1
       endif
       levels = findgen(31)/30.0*(maxi-mini) + mini

        contour,datatoplot(*,*),x(*,*),y(*,*),POS=pos,$
          levels=levels,xstyle=1,ystyle=1,$
          xrange=xrange,yrange=yrange,$
          title='Contour Plot Thermosphere',c_colors=clevels,$
          xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE
        
;        if cnt1 eq 1 then begin	
;            if (plotvectoryes eq 1) then begin
;                thermo_plotvectors,vars,selset,data,lat,lon,nlats,nlons,cf,vi_cnt,vn_cnt,step,$
;                  polar, mr								
;            endif
;        endif
        
                                ;Draw grid.
        if (showgridyes eq 1) then begin
            for i=0,ygrids-1 do begin
                oplot,mm(x),[y(0,i),y(0,i)]
            endfor
            for j=0,xgrids-1 do begin
                oplot,[x(j,0),x(j,0)],mm(y)
            endfor
        endif

                                ;If user set cursor position to plot, then do plot of datald.
                                ;Else clean up the cursor text fields on the interface. 
        if cursor_cnt eq 1 then begin
            if cnt1 eq 1 then begin
                plot,datald,alt,xtitle='datald',ytitle='Alt. (deg)'
            endif
            if cnt2 eq 1 then begin
                plot,datald,lat,xtitle='datald',ytitle='Lat. (deg)'
            endif
            if cnt3 eq 1 then begin
                plot,datald,lon,xtitle='datald',ytitle='Long. (deg)'
            endif
        endif else begin	
            txt=''
                                ;widget_control,(*ptr).curx_txt,set_value=txt
                                ;widget_control,(*ptr).cury_txt,set_value=txt
        endelse
    endif                       ;End of if nolog eq 1

                                ;If user does want to do the Plot Log. 
    if yeslog eq 1 then begin
        nzsubs=where(datatoplot gt 0, cnt)
        datatoplot(nzsubs)=datatoplot(nzsubs)
        datatoplot=ALOG10(datatoplot)
        maxi=max(datatoplot)
        mini=min(datatoplot)
        if mini eq maxi then maxi=mini+1
        levels = findgen(31)/30.0*(maxi-mini) + mini	
        
        contour,datatoplot(*,*),x(*,*),y(*,*),POS=pos,$
          levels=levels,xstyle=1,ystyle=1,$
          xrange=xrange,yrange=yrange,$
          title='Contour Plot Thermosphere',c_colors=clevels,$
          xtitle=xtitle,ytitle=ytitle,/cell_fill,/NOERASE

;        if cnt1 eq 1 then begin 
;            if (plotvectoryes eq 1) then begin
;                utrot = 180.0
;                thermo_plotvectors,vars,selset,data,lat,lon,nlats,nlons,cf,vi_cnt,vn_cnt,step,$
;                  polar, mr
;            endif
;        endif
        
                                ;Draw grid.
        if (showgridyes eq 1) then begin
            for i=0,ygrids-1 do begin
                oplot,mm(x),[y(0,i),y(0,i)]
            endfor
            for j=0,xgrids-1 do begin
                oplot,[x(j,0),x(j,0)],mm(y)
            endfor
        endif
        
                                ;If user set cursor position to plot, then do plot of datald.
                                ;Else clean up the cursor text fields on the interface. 
        if cursor_cnt eq 1 then begin
            if cnt1 eq 1 then begin
                plot,datald,alt,xtitle='datald',ytitle='Alt. (deg)'
            endif
            if cnt2 eq 1 then begin
                plot,datald,lat,xtitle='datald',ytitle='Lat. (deg)'
            endif
            if cnt3 eq 1 then begin
                plot,datald,lon,xtitle='datald',ytitle='Long. (deg)'
            endif
        endif else begin	
            txt=''
                                ;widget_control,(*ptr).curx_txt,set_value=txt
                                ;widget_control,(*ptr).cury_txt,set_value=txt
        endelse
    endif                       ;End of if yeslog eq 1
endif                           ;End of if yes eq 1


if (cnt1) then plane = 1
if (cnt2) then plane = 2
if (cnt3) then plane = 3

if (plotvectoryes eq 1) then begin
    if (not polar) then $
      plot,xrange,yrange, /nodata,xstyle=5,ystyle=5,/noerase, pos = pos
;    utrot = 180.0

    if (plane eq 1 and not ghostcells and not polar) then begin
        lon = lon + 180.0
        loc = where(lon gt 360.0,count)
        if (count gt 0) then lon(loc) = lon(loc) - 360.0
    endif
    thermo_plotvectors,vars,selset,data,lat, $
      lon, utrot, alt, $
      nlats,nlons,nalts, cf,vi_cnt,vn_cnt,step, polar, mr, plane, npolar
endif 

;stfr = 90.0-67.0
;stft = 309.0 + utrot
;
;cirr = 2.0
;cirt = findgen(17)*2*!pi/16
;cirx = cirr*cos(cirt)
;ciry = cirr*sin(cirt)
;stfx = stfr*cos(stft*!pi/180.0-!pi/2) + cirx
;stfy = stfr*sin(stft*!pi/180.0-!pi/2) + ciry
;polyfill, stfx, stfy, color = 0


                                ;Draw color bar.

;	   pos = [0.82,0.05,0.87,0.96]
pos(0) = pos(2)+0.025
pos(2) = pos(0)+0.03
maxmin = mm(levels)



title = variable
plotct,254,pos,maxmin,title,/right,color=color

maxi=max(datatoplot)
mini=min(datatoplot)

r = (maxmin(1) - maxmin(0)) * 0.03

if (mini gt maxmin(0)-r) then begin
    plots, [0,1], [mini, mini], thick = 3
    plots, [0,0.5], [mini, mini-r], thick = 3
    plots, [0,0.5], [mini, mini+r], thick = 3
    if (abs(mini) lt 10000.0 and abs(mini) gt 0.01) then begin
        smin = strcompress(string(mini, format = '(f10.2)'), /remove)
    endif else begin
        smin = strcompress(string(mini, format = '(e12.3)'), /remove)
    endelse
    xyouts, -0.1, mini, smin, align = 0.5, charsize = 0.75, orient = 90
endif

if (maxi lt maxmin(1)+r) then begin

    plots, [0,1], [maxi, maxi], thick = 3
    plots, [0,0.5], [maxi, maxi-r], thick = 3
    plots, [0,0.5], [maxi, maxi+r], thick = 3
    if (abs(maxi) lt 10000.0 and abs(maxi) gt 0.01) then begin
        smax = strcompress(string(maxi, format = '(f10.2)'), /remove)
    endif else begin
        smax = strcompress(string(maxi, format = '(e12.3)'), /remove)
    endelse
    xyouts, -0.1, maxi, smax, align = 0.5, charsize = 0.75, orient = 90
endif

;If user write the result to a file, then closedevice right now. 
;;if yeswrite_cnt eq 1 then begin
;;closedevice
;;mes=widget_message('Done with writing into file!')
;;endif

if (not IsManual) then begin
   smini = mini
   smaxi = maxi
endif

closedevice

end

