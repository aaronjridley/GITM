
IsEarth = 1

filelist = findfile("-t *.bin")

filelist = ask('filename to plot',filelist(0))
filelist = findfile(filelist)
nfiles = n_elements(filelist)

if (n_elements(iDiffPlot) eq 0) then iDiffPlot = 0
iDiffPlot = fix(ask('whether you want difference plots (0=no)',tostr(iDiffPlot)))

sDiff = ''
iPercentDiff = 0
if (iDiffPlot) then begin
   filelistBack = ask('filename(s) to subtract',filelist(0))
   filelistBack = findfile(filelistBack)
   nfilesBack = n_elements(filelistBack)
   if (nFilesBack ne nFiles) then begin
      print, 'nFilesBack ne nFiles!', nFilesBack, nFiles
      stop
   endif
   if (n_elements(iPercentDiff) eq 0) then iPercentDiff = 0
   iPercentDiff = fix(ask('whether to plot percent difference (0=no)',$
                          tostr(iPercentDiff)))
   if (iPercentDiff) then sDiff = ' (Percent Difference)' $
   else sDiff = ' (Difference)'

endif

if (n_elements(psfile) eq 0) then psfile = 'plot.ps'
psfile = ask('ps file name',psfile)

PlotNumber = 0
DoPlotVectors = 0

if (n_elements(MaxRange) eq 0) then MaxRange = 40.0

times = dblarr(nFiles)

if (n_elements(iPlotLocalTime) eq 0) then iPlotLocalTime = 0

nVectors = 0

for iFile = 0, nFiles-1 do begin

    filename = filelist(iFile)

    print, 'Reading file ',filename

    if (iFile eq 0) then begin
       read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                               vars, data, rb, cb, bl_cnt, iTime, Version
       c_a_to_r, itime, rtime
       if (iDiffPlot) then begin
          fileBack = filelistBack(iFile)
          read_thermosphere_file, fileBack, nvars, nalts, nlats, nlons, $
                                  vars, dataBack, rb, cb, bl_cnt,iTime,Version
       endif
    endif else begin
       gitm_read_bin_1var, filename, data, rtime, nVars, Vars, version, $
                          VarsToGet = varplotlist
       if (iDiffPlot) then begin
          fileBack = filelistBack(iFile)
          gitm_read_bin_1var, fileBack, dataBack, rtime, nVars, Vars, $
                              version, VarsToGet = varplotlist
       endif
    endelse

    Times(iFile) = rtime

    if (iFile eq 0) then begin
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

       s = size(data)
       Is2D = 0
       if (iAlt_ eq -1 or s(0) eq 3) then Is2D = 1

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

          if (iPlotType eq 3) then begin
             nSlicesToPlot = nLats
             if (nLats gt 1) then $
                for iLat = 0, nLats-1 do $
                   print, tostr(iLat),'. ',data(iLat_,0,iLat,0)/!dtor
          endif

          if (iPlotType eq 2) then begin
             iPlotLocalTime = $
                fix(ask('whether to plot local time (1) or lon (0)', $
                        tostr(iPlotLocalTime)))
             nSlicesToPlot = nLons
             if (iPlotLocalTime) then begin
                print, 'Enter the local time(s) you would like to plot.'
             endif else begin
                nSlicesToPlot = nLons
                if (nLons gt 1) then $
                   for iLon = 0, nLons-1 do $
                      print, tostr(iLon),'. ',data(iLon_,iLon,0,0)/!dtor
             endelse
          endif

          if (nSlicesToPlot gt 1) then begin
             SliceToPlot = 0
             nSlicesToPlot = 0
             while SliceToPlot gt -1 do begin
                SliceToPlot = $
                   ask('Slice number to plot (-1 to exit)','-1')
                if (iPlotLocalTime) then SliceToPlot = float(SliceToPlot) $
                else SliceToPlot = fix(SliceToPlot)
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

          if (iPlotType eq 2) then begin
             if (n_elements(minX) eq 0) then minX = -90.0
             if (n_elements(maxX) eq 0) then maxX =  90.0
             minX = float(ask('minimum latitude to plot',string(minX)))
             maxX = float(ask('maximum latitude to plot',string(maxX)))
          endif

          if (iPlotType eq 3) then begin
             if (n_elements(minX) eq 0) then minX =   0.0
             if (n_elements(maxX) eq 0) then maxX = 360.0
             minX = float(ask('minimum longitude to plot',string(minX)))
             maxX = float(ask('maximum longitude to plot',string(maxX)))
          endif

          if (iPlotType gt 1) then begin
             if (n_elements(minY) eq 0) then minY = alt(0,0,2)
             if (n_elements(maxY) eq 0) then maxY = alt(0,0,nAlts-3)
             minY = float(ask('minimum Altitude to plot',string(minY)))
             maxY = float(ask('maximum Altitude to plot',string(maxY)))
          endif

       endif else begin
          iPlotType = 1
          nSlicesToPlot = 1
          SlicePlotList = [0]
          nAlts = 1
       endelse

       if (iPlotType eq 1) then begin

          if (n_elements(IsPolar) eq 0) then IsPolar = 1
          IsPolar = $
             fix(ask('whether you want polar coordinates (1) or not (0)',$
                     tostr(IsPolar)))
          
          if (n_elements(IsNorth) eq 0) then IsNorth = 1
          if (IsPolar) then begin
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
                if (DoPlotIons) then begin
                   for iVar = 0, nVars-1 do begin
                      if (strpos(mklower(vars(iVar)),"v!di!n(east)") gt -1) then $
                         iVEast_ = iVar
                   endfor
                endif
                AllData = fltarr(nFiles,nVars,nLons,nLats,nAlts)
                nVectors = 2
                VarPlotList = [VarPlotList, iVEast_, iVEast_+1]
             endif
          endif

       endif else begin
          IsPolar = 0
       endelse

       DataToPlotAll = $
          fltarr(nFiles,nVarsToPlot+nVectors,nLons,nLats,nAlts)

    endif

    for iVar = 0, nVarsToPlot+nVectors-1 do begin

       if (iFile eq 0) then iV_ = VarPlotList(iVar) else iV_ = iVar
       if (iDiffPlot) then begin
          if (Is2D) then BackData = DataBack(iV_,*,*) $
          else BackData = DataBack(iV_,*,*,*)
       endif else BackData = 0.0

       if (Is2D) then DataToPlotAll(iFile,iVar,*,*,0) = data(iV_,*,*) - BackData $
       else DataToPlotAll(iFile,iVar,*,*,*) = data(iV_,*,*,*) - BackData

       if (iPercentDiff and iVar lt nVarsToPlot) then begin
          if (Is2D) then $
             DataToPlotAll(iFile,iVar,*,*,0) = $
             DataToPlotAll(iFile,iVar,*,*,0)/data(iV_,*,*)*100.0 $
          else DataToPlotAll(iFile,iVar,*,*,*) = $
             DataToPlotAll(iFile,iVar,*,*,*)/data(iV_,*,*,*)*100.0
       endif

    endfor

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
   if (iPlotType gt 1) then space = 0.03
    if (nX gt nY or (iPlotType eq 1 and not IsPolar)) then $
      setdevice, psfile, 'l', 5, 0.95 $
    else setdevice, psfile, 'p', 5, 0.95
    pos_space, ppp, space, sizes, nx = nX
endelse
plotdumb

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

cf = -1

iSlices_ = intarr(nFiles,nSlicesToPlot)

if (iPlotLocalTime) then begin

   for iFile = 0, nFiles-1 do begin

      c_r_to_a, itime, Times(iFile)
      ut = float(itime(3)) + float(itime(4))/60.0 + float(itime(5))/3600.0
      localtime = ((reform(lon(*,0,0)) + ut*15.0)/15.0 mod 24.0)

      for iSlice = 0, nSlicesToPlot-1 do begin
         dis = abs(localtime - SlicePlotList(iSlice))
         loc = where(dis eq min(dis))
         iSlices_(iFile, iSlice) = loc(0)
      endfor
   endfor

endif else for iS=0,nSlicesToPlot-1 do iSlices_(*,iS) = SlicePlotList(iS)

for iFile = 0, nFiles-1 do begin

    c_r_to_a, itime, Times(iFile)
    c_a_to_s, itime, sTime, /nice, IncludeSeconds=IncludeSeconds, $
      IncludeDay=IncludeDay

;    if (not includeDay) then begin
       c_a_to_s, itime, sDate
       sDate = strmid(sDate,0,9)
;    endif

    ut = float(itime(3)) + float(itime(4))/60.0 + float(itime(5))/3600.0
    if (IsPolar and not IsLocalTime) then utrot = ut * 15.0 else utrot = 0.0

    if (nFiles eq 1 and nSlicesToPlot eq 1) then IsOnlyVars = 1 $
    else IsOnlyVars = 0

    for iSlice = 0, nSlicesToPlot-1 do begin

       iS_ = iSlices_(iFile,iSlice)
       
       for iVar = 0, nVarsToPlot-1 do begin

          mini = 1.0e32
          maxi = -1.0e32

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
          if (iy eq sizes.nby-1) then begin
             if (not IsPolar) then no00 = 0
          endif

          if (ix gt 0 and iPlotType gt 1) then begin
             pos(0) = pos(0) - space/2.0 * ix
             pos(2) = pos(2) - space/2.0 * ix
          endif

          if (iy eq 0) then begin
             xp = (pos(0) + pos(2))/2.0
             yp = pos(3) + 0.02
             location = ''
             if (not Is2D) then begin
                location = sDiff+' (at '
                if (iPlotType eq 1) then $
                   location = $
                   location+string(alt(0,0,iS_),format='(f6.1)')+' km)'
                if (iPlotType eq 2) then begin
                   if (iPlotLocalTime) then begin
                      location = $
                         location+string(SlicePlotList(iSlice),$
                                         format='(f4.1)')+$
                         ' hours local time)'
                   endif else $
                      location = $
                      location+string(lon(iS_,0,0),format='(f5.1)')+' deg)'
                endif
                if (iPlotType eq 3) then $
                   location = $
                   location+string(lat(0,iS_,0),format='(f5.1)')+' deg)'
             endif
             if (not IndividualCt and iX eq 1) then $
                xyouts, 0.5, yp, /norm, vars(iV_)+location, align = 0.5

             if (ix eq 0) then begin
                ydatepos = pos(3)+0.05*(pos(3)-pos(1))
                xyouts,pos(0), ydatepos, sDate
             endif

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

             cLevel = tostr(alt(0,0,iS_))+' km'

             Lon1D = reform(lon(*,0,iS_))
             Lat1D = reform(lat(0,*,iS_))

             if (Individualct) then begin

                lat2d = reform(lat(*,*,iS_))
                if (not IsNorth and IsPolar) then Lat2D = -Lat2D
                loc = where(lat2d gt  MaxRange, count)

                d = reform(DataToPlotAll(iFile,iVar,*,*,iS_))
                mini = min(d(loc))
                maxi = max(d(loc))

             endif else begin
                  ; Want all of the plots to have the same color scale

                   lat2d = reform(lat(*,*,iS_))
                   if (IsPolar) then begin
                      if (not IsNorth) then Lat2D = -Lat2D
                      loc = where(lat2d gt MaxRange, count)
                   endif else loc = where(abs(lat2d) gt 0.0, count)

                   for iT = 0, nFiles-1 do begin
                      d = reform(DataToPlotAll(iT,iVar,*,*,iS_))
                      mini = min([mini,min(d(loc))])
                      maxi = max([maxi,max(d(loc))])
                   endfor

                endelse

                if (mini lt 0.0 and maxi gt 0.0) then begin
                   rmaxi = max([abs(mini),maxi])
                   if (maxi gt 0.05*rmaxi and mini lt -0.05*rmaxi) then begin
                      mini = -rmaxi
                      maxi = rmaxi
                      makect,'mid'
                   endif else makect,'wyrb'
                endif else makect,'wyrb'

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
                        
                        BottomCt = 1
                        if (iy eq sizes.nby-1) then plotctnow = 1

                    endelse

                    if (DoPlotVectors) then begin

                       if (iFile eq 0) then begin
                          maxV = -1.0e32
                          for iFi=0,nFiles-1 do begin
                             ve = reform(DataToPlotAll(iFi,nVarsToPlot,*,*,iS_))
                             vn = reform(DataToPlotAll(iFi,nVarsToPlot+1,*,*,iS_))
                             vm = sqrt(ve^2+vn^2)
                             if (max(vm) gt maxV) then maxV = max(vm)
                          endfor
                          factors = 10.0*[1.0, 5.0, 10.0, 20.0, 25.0, $
                                          50.0, 75.0, 100.0, 150.0, 200.0]
                          l = where(factors gt maxV)
                          iVectorFactor = l(0)
                       endif

                       SubData = reform(DataToPlotAll(iFile,*,*,*,*))
                       LatV = lat
                       LonV = lon
                       AltV = alt*1000.0
                       step = 2
                       plot, [-maxrange,maxrange], $
                             [-maxrange,maxrange], xstyle = 5, ystyle = 5, $
                             pos = pos, /noerase, /nodata
                       thermo_plotvectors,vars(varplotlist), iS_, SubData, $
                                          latV,lonV, $
                                          utrot, altV, nlats,nlons, nalts, $
                                          iVectorFactor,DoPlotIons,1-DoPlotIons, $
                                          step, IsPolar, MaxRange, $
                                          iPlotType, IsNorth
                    endif

                endif else begin

                    if (not no18) then ytitle = 'Latitude (deg)'

                    lon2d = reform(lon(*,*,iS_))
                    lat2d = reform(lat(*,*,iS_))
                    loc = where(abs(lat2d) le 90.0 and $
                                lon2d ge 0.0 and $
                                lon2d lt 360.0, count)

                    if (mini eq 0 and maxi eq 0) then maxi = 1.0

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

            if (iPlotType eq 2) then begin
               cLevel = tostr(lon(iS_,0,0))+' Deg'
               DataToPlot = reform(DataToPlotAll(iFile,iVar,iS_,*,*))
               xPlot = reform(lat(iS_,*,*))
               yPlot = reform(alt(iS_,*,*))
               xtitle = 'Latitude (deg)'
               if (iVar eq 0) then ytitle = 'Altitude (km)'

               l = where(xPlot ge minX and xPlot le maxX and $
                         yPlot ge minY and yPlot le maxY)
               for iT = 0, nFiles-1 do begin
                  iS = iSlices_(iT,iSlice)
                  d = reform(DataToPlotAll(iT,iVar,iS,*,*))
                  mini = min([mini,min(d(l))])
                  maxi = max([maxi,max(d(l))])
               endfor
               r = (maxi-mini)*0.01
               mini = mini-r
               maxi = maxi+r

            endif
            if (iPlotType eq 3) then begin
               cLevel = tostr(lat(0,iS_,0))+' Deg'
               DataToPlot = reform(DataToPlotAll(iFile,iVar,*,iS_,*))
               xPlot = reform(lon(*,iS_,*))
               yPlot = reform(alt(*,iS_,*))
               xtitle = 'Longitude (deg)'
               if (iVar eq 0) then ytitle = 'Altitude (km)'
               l = where(xPlot ge minX and xPlot le maxX and $
                         yPlot ge minY and yPlot le maxY)
               for iT = 0, nFiles-1 do begin
                  iS = iSlices_(iT,iSlice)
                  d = reform(DataToPlotAll(iT,iVar,*,iS,*))
                  mini = min(d(l))
                  if (mini gt 0) then mini = mini*0.95 else mini = mini*1.05
                  maxi = max(d(l))*1.05
               endfor
            endif

            if (iPlotType gt 1) then begin
               
               ix = fix(PlotNumber) mod sizes.nbx
               iy = fix(PlotNumber)/sizes.nbx

               xticknames = strarr(10)
               yticknames = strarr(10)

               if (ix gt 0) then begin
                  ytitle = ''
                  yticknames = strarr(10)+' '
               endif
               if (iy lt sizes.nby-1) then begin
                  xtitle = ''
                  xticknames = strarr(10)+' '
               endif

               nLevels = 31
               nLevelsLines = 11
               if (mini lt 0.0 and maxi gt 0.0) then begin
                  rmaxi = max([abs(mini),maxi])
                  if (maxi gt 0.05*rmaxi and mini lt -0.05*rmaxi) then begin
                     mini = -rmaxi
                     maxi = rmaxi
                     makect,'mid'
                  endif else makect,'wyrb'
               endif else makect,'wyrb'
               r = (maxi-mini)
               levels = findgen(nLevels)/(nLevels-1)*r + mini

               linelevels = findgen(nLevelsLines)/(nLevelsLines-1)*r + mini

               xrange = [minX,MaxX]
               yrange = [minY,MaxY]

               contour, datatoplot, xPlot, yPlot, /noerase, $
                        pos = pos, /fill, levels = levels, $
                        c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
                        xrange = xrange, yrange = yrange, $
                        xstyle = 1, ystyle = 1, $
                        xtitle = xtitle, ytitle = ytitle, $
                        xtickname = xtickname, $
                        ytickname = ytickname

               contour, datatoplot, xPlot, yPlot, /noerase, $
                        pos = pos, levels = linelevels, $
                        xrange = xrange, yrange = yrange, $
                        xstyle = 5, ystyle = 5, $
                        xtitle = xtitle, ytitle = ytitle, $
                        xtickname = xtickname, $
                        ytickname = ytickname, /follow, $
                        c_linestyle = [3]

               yMean = mean(yrange)
               l = where(yplot(0,*) gt yMean)
               iAlt_ = l(0)

               plot, xplot(*,iAlt_), datatoplot(*,iAlt_), $
                     xrange = xrange, yrange = mm(levels), $
                     xstyle = 5, ystyle = 5, pos = pos, /noerase, $
                     thick = 6

               if (iFile eq 0) then begin
                  nX = n_elements(datatoplot(*,iAlt_))
                  oneD = fltarr(nFiles,nX)
                  xSave = reform(xplot(*,iAlt_))
                  xRangeSave = xRange
                  yRangeSave = mm(levels)
               endif

               oneD(iFile, *) = datatoplot(*,iAlt_)

               if (mini lt 0.0) then begin
                  oplot, xplot(*,iAlt_), datatoplot(*,iAlt_)*0.0, $
                         linestyle = 1, thick = 6
               endif else begin
                  x1 = xrange(0)
                  x2 = x1 + 0.1*(max(xrange(1)-xrange(0)))
                  y1 = (max(levels)+min(levels))/2.0
                  oplot, [x1,x2], [y1,y1], $
                         linestyle = 1, thick = 6
               endelse

               xpos = pos(0) + 0.05*(pos(2)-pos(0))
               ypos = pos(3) - 0.1*(pos(2)-pos(0))
               xyouts, xpos, ypos, stime, /norm

               if (nVarsToPlot eq 1 and iX eq sizes.nbx-1 ) then plotctnow = 1

            endif

            if (plotctnow) then begin
                minmax = [mini,maxi]
                ctpos = pos
                if (BottomCt) then begin
                    ctpos(1) = ctpos(1)-0.03
                    ctpos(3) = ctpos(1)+0.015
                    title = vars(iV_)
                endif else begin
                    ctpos(0) = ctpos(2)+0.005
                    ctpos(2) = ctpos(0)+0.01
                    title = ' '
                    title = vars(iV_)
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

if (nFiles gt 1 and n_elements(oneD) gt 0) then begin

   l = strpos(psfile,'.ps')
   psfileLine = strmid(psfile,0,l)+'_timehistory.ps'

   setdevice, psfileLine, 'p', 5

   makect,'bry'

   l = where(xSave ge xRangeSave(0) and xSave le xRangeSave(1))
   maxi = max(oneD(*,l))
   mini = min(oneD(*,l))

   plot, xSave, oneD(0,*), $
         xRange = xRangeSave, yrange = [mini,maxi], $
         pos = [0.1, 0.3, 0.90, 0.8], $
         xtitle = xtitle, ytitle = vars(iV_)+location, /nodata

   for iFile = 0, nFiles-1 do begin
      c = 245.0*iFile/(nFiles-1) + 5.0
      oplot, xSave, oneD(iFile,*), color = c, thick = 5
   endfor

   ctpos = [0.91, 0.3, 0.95, 0.8]
   maxi = (max(times)-min(times))/3600.0
   minmax = [0.0,maxi]
   title = 'Hours after Start Time'
   plotct, 255, ctpos, minmax, title, /right

   closedevice

endif

end
