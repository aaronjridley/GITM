
wildcard = ask('files to plot tec (e.g., 3DALL*.bin)','3DALL*.bin')

filelist = findfile(wildcard)
file = filelist(0)

gitm_read_header, file, time, nVars, Vars, nLons, nLats, nAlts, version

IsDone = 0
iNO_ = 8
iVar = [0,1,2,iNO_]

nFiles = n_elements(filelist)

TotalNO = fltarr(nFiles)
TotalTime = dblarr(nFiles)

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   print, file
   gitm_read_bin_1var,file, data, time, nVars, Vars, version, VarsToGet = iVar

   TotalTime(iFile) = time

   if (iFile eq 0) then begin

      lats = reform(data(1,2:nLons-3,2:nLats-3,0))
      lons = reform(data(0,2:nLons-3,2:nLats-3,0))
      alts = reform(data(2,0,0,*))
      dLat = lats(0,1) - lats(0,0)
      dlon = lons(1,0) - lons(0,0)
      area = 6372000.0*6372000.0 * cos(lats) * dlon * dlat

      iNO = 3
      iVar = [iNO_]

   endif else iNO = 0

   no = reform(data(iNO,2:nLons-3,2:nLats-3,*))

   for iAlt = 2, nAlts-5 do begin

      da = (alts(iAlt+1) - alts(iAlt-1))
      TotalNO(iFile) = TotalNO(iFile) + total(area*no(*,*,iAlt))*dA

   endfor

endfor

end

