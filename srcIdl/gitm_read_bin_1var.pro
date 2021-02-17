

pro gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                        VarsToGet = VarsToGet

  if (n_elements(file) eq 1) then filelist = findfile(file) $
  else filelist = file

  nFiles = n_elements(filelist)

  DoAppendFile = 0

  OldTime = -1.0d32

  nTimes = 0L

  for iFile = 0L, nFiles-1 do begin

     if (iFile mod 100 eq 0 and nFiles gt 100) then $
        print, "Progress : ",iFile, nFiles-1

     filein = filelist(iFile)

     valid = 0
     on_ioerror, skipfile

     close, 1
     openr, 1, filein, /f77

     version = 0.0D

     nLons = 0L
     nLats = 0L
     nAlts = 0L
     nVars = 0L
     
     IsFirstTime = 1

     DoAppendFile = 0

     while (not eof(1)) do begin

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

        iHeaderLength = 8L + 4+4 + $ ; Version + start/stop byte
                        3*4 + 4+4 + $ ; nLons, nLats, nAlts + start/stop byte
                        4 + 4+4 + $ ; nVars + start/stop byte
                        nVars*40 + nVars*(4+4) + $ ;Variable Names + start/stop byte
                        7*4 + 4+4 ; time + start/stop byte

        iDataLength = nLons*nLats*nAlts*8L + 4+4

        tmp = dblarr(nLons, nLats, nAlts)
        
        OnlyGetSome = 0
        if (n_elements(VarsToGet) gt 0) then begin
           nVarsTotal = n_elements(VarsToGet)
           OnlyGetSome = 1
        endif else begin
           nVarsTotal = nVars
        endelse

        if (IsFirstTime and nTimes eq 0) then begin
           Data = dblarr(nVarsTotal, nLons, nLats, nAlts)
           if (nFiles gt 1) then begin
              DataAllNew = dblarr(nFiles, nVarsTotal, nLons, nLats, nAlts)
              TimeAllNew = dblarr(nFiles)
           endif else begin
              DataAllNew = dblarr(1, nVarsTotal, nLons, nLats, nAlts)
              TimeAllNew = dblarr(1)
           endelse
        endif

        iVar = 0

        if (OnlyGetSome) then begin
           for i=0,n_elements(VarsToGet)-1 do begin
              pos = iHeaderLength + VarsToGet(i)*iDataLength
              point_lun, 1, pos
              readu,1,tmp
              data(i,*,*,*) = tmp
           endfor
        endif else begin
           for i=0L,nVars-1 do begin
              readu,1,tmp
              data(i,*,*,*) = tmp
           endfor
        endelse

        if (OnlyGetSome) then begin
           pos = iHeaderLength + nVars*iDataLength
           point_lun, 1, pos
        endif

        ; Figure out whether we have more than one data set in file
        if (not eof(1)) then DoAppendFile = 1

        ; nFiles = 1 and only 1 set of data in the file
        if (nFiles eq 1 and not DoAppendFile) then begin
           ; do nothing.  We have the data exactly as we want it
        endif

        ; nFiles > 1 and only 1 set of data in the file
        if (nFiles gt 1 and not DoAppendFile) then begin
           DataAllNew(nTimes,*,*,*,*) = data
           TimeAllNew(nTimes) = time
           nTimes = nTimes+1
        endif

        ; multiple data sets

        if (DoAppendFile) then begin
           if (IsFirstTime) then begin
              IsFirstTime = 0
              spawn, 'ls -s '+filein,filesize
              filesize = long(filesize(0))*1024.0
              datasize = (nVars*nLons*nLats*nAlts*8L + $
                          8 + 3*4 + 4 + nVars*40 + 7*4) + $
                         2*4*(3+nVars+1+nVars)
              nTimesMax = round(filesize/datasize)
              print, 'Found file with multiple times: ',filein
              print, 'Assuming that file contains nTimes : ',nTimesMax

              if (nFiles gt 1) then begin
                 ; How many data sets could POSSIBLY be there?
                 ; Add together:
                 ; 1. how many we have already read (nTimes)
                 ; 2. How many we think are in this file (nTimesMax)
                 ; 3. How many files we have left (nFiles-iFile)
                 nTimesMax = nTimesMax + nTimes + (nFiles-iFile)
                 print, 'Reallocating to ',nTimesMax

              endif
              
              ; If we already have an array full of stuff, move the
              ; stuff!
              if (nTimes gt 0) then begin
                 DataAllOld = DataAllNew
                 TImeAllOld = TimeAllNew
              endif

              ; Create a new array for all of the data
              DataAllNew = dblarr(nTimesMax, nVarsTotal, nLons, nLats, nAlts)
              TimeAllNew = dblarr(nTimesMax)

              ; If we have old stuff, move it into the new array
              if (nTimes gt 0) then begin
                 for it = 0L, nTimes-1 do begin
                    DataAllNew(it,*,*,*,*) = DataAllOld(it,*,*,*,*)
                    TimeAllNew(iT) = TimeAllOld(iT)
                 endfor
              endif

           endif

           DataAllNew(nTimes,*,*,*,*) = data
           TimeAllNew(nTimes) = time
           nTimes = nTimes+1

        endif

     endwhile

     valid = 1

     skipfile: 
     if (not valid) then begin
        print, "Bad file found : ",filein
;        if (nFiles gt 1) then time(ifile) = 0.0
     endif

     close, 1

  endfor

  if (nTimes gt 1) then begin
     data = DataAllNew
     time = TimeAllNew
  endif

  ; If we had some bad files, we could be running into trouble with
  ; missing data.  Clean it up here.

  if (nTimes gt 1 or nFiles gt 1) then begin
     nTimes = n_elements(time)
     l = where(time gt 0.0,c)
     if (c lt nTimes) then begin
        print, "Found some empty times.  Correcting."
        DataAllNew = dblarr(c, nVarsTotal, nLons, nLats, nAlts)
        TimeAllNew = dblarr(c)
        for i=0L,c-1 do begin
           DataAllNew(i,*,*,*,*) = data(l(i),*,*,*,*)
           TimeAllNew(i) = time(l(i))
        endfor
        data = DataAllNew
        time = TimeAllNew
     endif
  endif

end

