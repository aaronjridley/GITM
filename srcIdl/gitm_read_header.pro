pro gitm_read_header, file, time, nVars, Vars, nLons, nLats, nAlts, version

  if (n_elements(skip) eq 0) then skip = 0
  
  if (n_elements(file) eq 1) then filelist = findfile(file) $
  else filelist = file

  nFiles = n_elements(filelist)

;  if (nFiles gt 1) then Time = dblarr(nFiles) else Time = 0.0D

  DoAppendFile = 0

  OldTime = -1.0d32

  nTimes = 0L

  for iFile = 0L, nFiles-1 do begin

     if (iFile mod 100 eq 0 and nFiles gt 1000) then $
        print, "Progress : ",iFile, nFiles-1

     filein = filelist(iFile)

     close, 1
     openr, 1, filein, /f77

     version = 0.0D

     nLons = 0L
     nLats = 0L
     nAlts = 0L
     nVars = 0L
     
     IsFirstTime = 1

     DoAppendFile = 0

     readu, 1, version
     readu, 1, nLons, nLats, nAlts
     readu, 1, nVars

     Vars = strarr(nVars)
     line = bytarr(40)
     for iVars = 0, nVars-1 do begin
        readu, 1, line
        Vars(iVars) = strcompress(string(line),/remove)
     endfor

     lTime = lonarr(7)
     readu, 1, lTime

     iTime = fix(lTime(0:5))
     c_a_to_r, itime, rtime
     Time = rTime + lTime(6)/1000.0

     if (iFile eq 0) then begin
        if (nFiles gt 1) then begin
           TimeAllNew = dblarr(nFiles)
        endif else begin
           TimeAllNew = dblarr(1)
        endelse
     endif

     TimeAllNew(iFile) = time

     close, 1

  endfor

  time = TimeAllNew
  if (n_elements(time) eq 1) then time = time(0)

end

