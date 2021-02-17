
IsEarth = 1

filelist = findfile("-t *.bin")

filelist = ask('perturbation filenames to plot',filelist(0))
filelist = findfile(filelist)
nfiles = n_elements(filelist)

if (n_elements(filelist_bg) eq 0) then filelist_bg = findfile("-t *.bin")
filelist_bg = ask('background filenames to plot',filelist_bg(0))
filelist_bg = findfile(filelist_bg)
nfiles_bg = n_elements(filelist_bg)

if (nfiles ne nfiles_bg) then begin

   print, 'nfiles and nfiles_bg not equal to each other. Yikes.'
   stop
   
endif

if (n_elements(psfile) eq 0) then psfile = 'plot.ps'
psfile = ask('ps file name',psfile)

PlotNumber = 0
DoPlotVectors = 0

if (n_elements(MaxRange) eq 0) then MaxRange = 40.0

times = dblarr(nFiles)

for iFile = 0, nFiles-1 do begin

    filename = filelist(iFile)
    print, 'Reading file ',filename
    read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
      vars, data, rb, cb, bl_cnt, iTime, Version

    filename_bg = filelist_bg(iFile)
    read_thermosphere_file, filename_bg, nvars_bg, nalts_bg, nlats_bg, nlons_bg, $
      vars_bg, data_bg, rb, cb, bl_cnt, iTime_bg, Version

    c_a_to_r, itime, rtime
    Times(iFile) = rtime

    iAlt_ = -1
    for iVar = 0, nVars-1 do begin
        p = strpos(mklower(vars(iVar)),"alt")
        if (p eq 0) then begin
            iAlt_ = iVar
            alt = reform(data(iAlt_,*,*,*)) / 1000.0
        endif
    endfor

    iLat_ = -1
    for iVar = 0, nVars-1 do begin
        if (strpos(mklower(vars(iVar)),"lat") gt -1 and iLat_ eq -1) then begin
            iLat_ = iVar
            lat = reform(data(iLat_,*,*,*)) / !dtor
        endif
    endfor

    iLon_ = -1
    IsLocalTime = 0
    for iVar = 0, nVars-1 do begin
        if (strpos(mklower(vars(iVar)),"lon") gt -1 and iLon_ eq -1) then begin
            iLon_ = iVar
            lon = reform(data(iLon_,*,*,*)) / !dtor
        endif
        if (strpos(mklower(vars(iVar)),"local") gt -1) then begin
            iLon_ = iVar
            lon = reform(data(iLon_,*,*,*))
            IsLocalTime = 1
        endif
    endfor

    iVEast_ = -1
    for iVar = 0, nVars-1 do begin
        if (strpos(mklower(vars(iVar)),"v!dn!n(east)") gt -1) then begin
            iVEast_ = iVar
        endif
    endfor

    Is2D = 0
    if (iAlt_ eq -1) then Is2D = 1

    if (iFile eq 0) then begin

        for iVar = 0, nVars-1 do print, tostr(iVar),'. ',vars(iVar)

        VarToPlot = 0
        nVarsToPlot = 0
        while VarToPlot gt -1 do begin
            VarToPlot = fix(ask('Variable Number to plot (-1 to exit)','-1'))
            if (VarToPlot ge 0) and (VarToPlot lt nVars) then begin
                if nVarsToPlot eq 0 then $
                  VarPlotList = [VarToPlot]          $
                else VarPlotList = [VarPlotList,VarToPlot]
                nVarsToPlot = nVarsToPlot + 1
            endif
        endwhile

        if (not(Is2D)) then begin

            print, '1. Altitude Slice'
            print, '2. Longitude Slice'
            print, '3. Latitude Slice'
            if (n_elements(iPlotType) eq 0) then iPlotType = 1
            iPlotType = fix(ask('what type of plot do you want',$
                                tostr(iPlotType)))

            if (iPlotType eq 1) then begin
                nSlicesToPlot = nAlts
                if (nAlts gt 1) then $
                  for iAlt = 0, nAlts-1 do $
                  print, tostr(iAlt),'. ',data(iAlt_,0,0,iAlt)
            endif

            if (iPlotType eq 2) then begin
                nSlicesToPlot = nLats
                if (nLats gt 1) then $
                  for iLat = 0, nLats-1 do $
                  print, tostr(iLat),'. ',data(iLat_,0,iLat,0)/!dtor
            endif

            if (iPlotType eq 3) then begin
                nSlicesToPlot = nLons
                if (nLons gt 1) then $
                  for iLon = 0, nLons-1 do $
                  print, tostr(iLon),'. ',data(iLon_,iLon,0,0)/!dtor
            endif

            if (nSlicesToPlot gt 1) then begin
                SliceToPlot = 0
                nSlicesToPlot = 0
                while SliceToPlot gt -1 do begin
                    SliceToPlot = $
                      fix(ask('Slice number to plot (-1 to exit)','-1'))
                    if (SliceToPlot ge 0) then begin
                        if nSlicesToPlot eq 0 then $
                          SlicePlotList = [SliceToPlot]          $
                        else SlicePlotList = [SlicePlotList,SliceToPlot]
                        nSlicesToPlot = nSlicesToPlot + 1
                    endif
                endwhile
            endif else begin
                SlicePlotList = [0]
            endelse

        endif else begin
            iPlotType = 1
            nSlicesToPlot = 1
            SlicePlotList = [0]
        endelse

        if (iPlotType eq 1) then begin

            if (n_elements(IsPolar) eq 0) then IsPolar = 1
            IsPolar = $
              fix(ask('whether you want polar coordinates (1) or not (0)',$
                      tostr(IsPolar)))

            if (IsPolar) then begin
                if (n_elements(IsNorth) eq 0) then IsNorth = 1
                IsNorth = $
                  fix(ask('whether you want northern hemisphere(1) or not(0)',$
                          tostr(IsNorth)))
                MaxRange = float(ask('maximum range to plot',string(MaxRange)))
            endif

            if (iVEast_ gt -1) then begin
                if (n_elements(DoPlotVectors) eq 0) then DoPlotVectors=1
                DoPlotVectors = fix(ask('whether you want vectors (1=yes)',$
                                        tostr(DoPlotVectors)))
                if (DoPlotVectors) then begin
                    if (n_elements(DoPlotIons) eq 0) then DoPlotIons=0
                    DoPlotIons=fix(ask('whether you want ion vectors (1=yes)',$
                                            tostr(DoPlotIons)))
                    AllData = fltarr(nFiles,nVars,nLons,nLats,nAlts)
                endif
            endif

        endif else IsPolar = 1

        DataToPlotAll = $
          fltarr(nFiles,nVarsToPlot,nLons,nLats,nAlts)

    endif

    for iVar = 0, nVarsToPlot-1 do begin
        iV_ = VarPlotList(iVar)
        DataToPlotAll(iFile,iVar,*,*,*) = data(iV_,*,*,*)-data_bg(iV_,*,*,*)
    endfor

    if (DoPlotVectors) then begin
        AllData(iFile,*,*,*,*) = Data-Data_bg
    endif

endfor

nLevels = 31

if (nVarsToPlot gt 1) then begin
    nX = nVarsToPlot
    if (nSlicesToPlot gt 1) then nY = nSlicesToPlot $
    else nY = nFiles
endif else begin
;    if (nSlicesToPlot gt 1) then begin
        nX = nSlicesToPlot
        nY = nFiles
;    endif
endelse


if (nSlicesToPlot eq 1 and nVarsToPlot eq 1) then begin
   nX = min([nFiles, 3])
   nY = nFiles/nX
endif

print, nx, ny 
if (nY gt 4) then nY = 4

ppp = nX * nY
space = 0.01
if (nY eq 1 or nX eq 1) then begin
    if (iPlotType eq 1 and not IsPolar) then begin
        if (nY eq 1) then $
          setdevice, psfile, 'l', 5, 0.95 $
        else setdevice, psfile, 'p', 5, 0.95
    endif else begin
        space = 0.03
        setdevice, psfile, 'p', 5, 0.95
    endelse
    pos_space, ppp, space, sizes
endif else begin
    if (nX gt nY or (iPlotType eq 1 and not IsPolar)) then $
      setdevice, psfile, 'l', 5, 0.95 $
    else setdevice, psfile, 'p', 5, 0.95
    pos_space, ppp, space, sizes, nx = nX
endelse
plotdumb
makect,'mid'

; Figure out whether we need individual color-tables or one color 
; table for each column.

if (nFiles eq 1 and nSlicesToPlot eq 1) then Individualct = 1
if (nFiles gt 1 and nSlicesToPlot eq 1) then Individualct = 0
if (nSlicesToPlot gt 1) then Individualct = 1

if (nFiles gt 1) then begin
    dt = Times(1)-Times(0)
    if (dt ge 60.0) then IncludeSeconds = 0 else IncludeSeconds = 1
    dt = Times(nFiles-1)-Times(0)
    if (dt gt 86400.0) then IncludeDay = 1 else IncludeDay = 0
endif

BottomCt = 0

for iFile = 0, nFiles-1 do begin

    c_r_to_a, itime, Times(iFile)
    c_a_to_s, itime, sTime, /nice, IncludeSeconds=IncludeSeconds, $
      IncludeDay=IncludeDay

    ut = float(itime(3)) + float(itime(4))/60.0 + float(itime(5))/3600.0
    if (IsPolar and not IsLocalTime) then utrot = ut * 15.0 else utrot = 0.0

    if (nFiles eq 1 and nSlicesToPlot eq 1) then IsOnlyVars = 1 $
    else IsOnlyVars = 0

    for iSlice = 0, nSlicesToPlot-1 do begin

        iS_ = SlicePlotList(iSlice)

        for iVar = 0, nVarsToPlot-1 do begin

            iV_ = VarPlotList(iVar)

            plotctnow = 0

            if (IsPolar) then $
               get_position, ppp, space, sizes, PlotNumber, pos $
            else begin
                get_position, ppp, space, sizes, PlotNumber, pos

                if (IndividualCt) then pos(2) = pos(2)-0.06

                ; This basically tries to make the plots to have the
                ; correct height to width ratio (which should be 1:2)

                if (iPlotType eq 1) then begin
;                    xs = (pos(2) - pos(0))/sizes.xf
;                    ys = (pos(3) - pos(1))/sizes.yf
                    xs = (pos(2) - pos(0))
                    ys = (pos(3) - pos(1))

if (xs*2 gt 1.0/float(sizes.nbx)) then begin
   xsnew = (1.0-space*3)/float(sizes.nbx) / 2
   fac = xsnew / xs
   ys = ys * fac
   xs = xs * fac
   offset = 0.0
endif else begin
   offset = (1.0 - 2.0*xs* float(sizes.nbx))/4
endelse

pos(0) = float(PlotNumber mod sizes.nbx) / float(sizes.nbx) + offset

pos(2) = pos(0)+xs*2
pos(3) = pos(1)+ys

;                    if (xs gt ys) then begin
;                        if (ys*2 lt xs) then begin
;                            xMid = (pos(0)+pos(2))/2.0
;                            pos(0) = xMid - ys*sizes.xf
;                            pos(2) = xMid + ys*sizes.xf
;                        endif else begin
;                            yMid = (pos(1)+pos(3))/2.0
;                            pos(1) = yMid - xs*sizes.yf/4.0
;                            pos(3) = yMid + xs*sizes.yf/4.0
;                        endelse
;                    endif else begin
;                        if (xs*2 lt ys) then begin
;                            yMid = (pos(1)+pos(3))/2.0
;                            pos(1) = yMid - xs*sizes.yf
;                            pos(3) = yMid + xs*sizes.yf
;                        endif else begin
;                            xMid = (pos(0)+pos(2))/2.0
;                            pos(0) = xMid - ys*sizes.xf/4.0
;                            pos(2) = xMid + ys*sizes.xf/4.0
;                        endelse
;                    endelse
                endif
            endelse

            no00 = 1
            no06 = 1
            no12 = 1
            no18 = 1
            ix = fix(PlotNumber) mod sizes.nbx
            iy = fix(PlotNumber)/sizes.nbx

            if (ix eq 0) then no18 = 0
            if (ix eq sizes.nbx-1) then no06 = 0
            if (iy eq 0) then no12 = 0

            if (iy eq 0) then begin
                xp = (pos(0) + pos(2))/2.0
                yp = pos(3) + 0.02
                location = ''
                if (not Is2D) then begin
                    location = ' (at '
                    if (iPlotType eq 1) then $
                      location = $
                      location+string(alt(0,0,iS_),format='(f6.1)')+' km)'
                    if (iPlotType eq 2) then $
                      location = $
                      location+string(lon(iS_,0,0),format='(f5.1)')+' deg)'
                    if (iPlotType eq 3) then $
                      location = $
                      location+string(lat(0,iS_,0),format='(f5.1)')+' deg)'
                endif
                if (not IndividualCt) then $
                  xyouts, xp, yp, /norm, vars(iV_)+location, align = 0.5

                if (IsPolar) then begin
                    if (IsNorth) then $
                      xyouts, 1.0, yp+0.02, /norm, "North", align = 1.0 $
                    else $
                      xyouts, 1.0, yp+0.02, /norm, "South", align = 1.0
                endif

            endif

            pos([1,3]) = pos([1,3]) + iy * space/2

            xtickname = strarr(10)
            ytickname = strarr(10)
            ytitle = ''
            if (no00) then xtickname = xtickname + ' '
            if (no18) then ytickname = ytickname + ' ' 

            if (iPlotType eq 1) then begin

                DataToPlot = reform(DataToPlotAll(iFile,iVar,*,*,iS_))

                Lon1D = reform(lon(*,0,iS_))
                Lat1D = reform(lat(0,*,iS_))

                if (Individualct) then begin
                    mini = min(DataToPlotAll(iFile,iVar,*,*,iS_))
                    maxi = max(DataToPlotAll(iFile,iVar,*,*,iS_))
                endif else begin
                    ; Want all of the plots to have the same color scale
                    mini = min(DataToPlotAll(*,iVar,*,*,iS_))
                    maxi = max(DataToPlotAll(*,iVar,*,*,iS_))
                endelse

                if (mini lt 0.0 and maxi gt 0.0) then begin
                   rmaxi = max([abs(mini),maxi])
                   if (maxi gt 0.3*rmaxi and mini lt -0.3*rmaxi) then begin
                      mini = -rmaxi
                      maxi = rmaxi
                   endif
                endif

                if (IsPolar) then begin

                    ; Simple way of plotting the southern hemisphere is
                    ; to reverse the sign on the latitudes.....
                    if (not IsNorth) then Lat1D = -Lat1D

                    lo = [212.51] + utrot
                    la = [65.13]

                    if (Individualct) then begin

                        contour_circle, DataToPlot, Lon1D+utrot, Lat1D, $
                          no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                          pos = pos, $
                          maxrange = MaxRange, $
                          colorbar = vars(iV_), $
                          mini = mini, maxi = maxi
;, latExtra = la, lonExtra = lo

                        if (PlotNumber eq 0) then begin
                            xyouts, 0.0, pos(3)+0.05, /norm, stime
                        endif

                    endif else begin

                        contour_circle, DataToPlot, Lon1D+utrot, Lat1D, $
                          no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                          pos = pos, $
                          maxrange = MaxRange, $
                          mini = mini, maxi = maxi, $
                          title = sTime
;, latExtra = la, lonExtra = lo
                        
                        BottomCt = 1
                        if (iy eq sizes.nby-1) then plotctnow = 1

                    endelse

                    if (DoPlotVectors) then begin
                        SubData = reform(AllData(iFile,*,*,*,*))
                        LatV = reform(AllData(iFile,iLat_,*,*,*))/!dtor
                        LonV = reform(AllData(iFile,iLon_,*,*,*))/!dtor
                        AltV = reform(AllData(iFile,iAlt_,*,*,*))
                        cf = -1
                        step = 2
                        plot, [-maxrange,maxrange], $
                          [-maxrange,maxrange], xstyle = 5, ystyle = 5, $
                          pos = pos, /noerase, /nodata
                        thermo_plotvectors,vars,iS_,SubData, latV,lonV, $
                          utrot, altV, nlats,nlons, nalts, $
                          cf,DoPlotIons,1-DoPlotIons, $
                          step, IsPolar, MaxRange, iPlotType, IsNorth
                    endif

                                ; , title = title, nolines = nolines
;                    nLevels = nLevels, $

                endif else begin

                    if (not no18) then ytitle = 'Latitude (deg)'

                    lon2d = reform(lon(*,*,iS_))
                    lat2d = reform(lat(*,*,iS_))
                    loc = where(abs(lat2d) le 90.0 and $
                                lon2d ge 0.0 and $
                                lon2d lt 360.0, count)

                    contour, datatoplot(loc), lon2d(loc), lat2d(loc), $
                      /noerase, pos = pos, /fill, $
                      nlevels = 31, xtickname = xtickname, $
                      ytickname = ytickname, xstyle = 1, ystyle = 1, /irr, $
                      levels = findgen(nlevels)*(maxi-mini)/(nlevels-1)+mini, $
                      c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
                      ytitle = ytitle, yrange = [-90,90]

                    !p.position = pos
                    map_set, /noerase
                    map_continents, color = 0

                    if (IndividualCt) then BottomCt = 0 else BottomCt = 1
                    if (iy eq sizes.nby-1 or $
                        iFile eq nFiles-1 or $
                        IndividualCt) then plotctnow = 1
                        
                endelse

            endif

            if (iPlotType eq 2) then $
              DataToPlot = reform(DataToPlotAll(iFile,iVar,iS_,*,*))
            if (iPlotType eq 3) then $
              DataToPlot = reform(DataToPlotAll(iFile,iVar,*,iS_,*))

            if (iPlotType gt 1) then begin
                contour, datatoplot, /noerase, pos = pos, /fill, nlevels = 31
            endif

            if (plotctnow) then begin
                minmax = [mini,maxi]
                ctpos = pos
                if (BottomCt) then begin
                    ctpos(1) = ctpos(1)-0.02
                    ctpos(3) = ctpos(1)+0.01
                    title = vars(iV_)
                endif else begin
                    ctpos(0) = ctpos(2)+0.005
                    ctpos(2) = ctpos(0)+0.01
                    title = ' '
                endelse
                plotct, 255, ctpos, minmax, title, bottom=BottomCt, $
                  right = 1-BottomCt
            endif

            PlotNumber = (PlotNumber + 1) mod ppp

            if (PlotNumber eq 0) then plotdumb 

        endfor

    endfor

endfor

closedevice

end
