;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf

pro read_thermosphere_file, filelist, nvars, nalts, nlats, nlons, $
                            vars, data, nBLKlat, nBLKlon, nBLK, $
                            iTime, Version

  if (n_elements(nBLKlat) eq 0) then nBLKlat = 0
  if (n_elements(nBLKlon) eq 0) then nBLKlon = 0
  if (n_elements(nBLK)    eq 0) then nBLK    = 0

  filelist = findfile(filelist)

  Version = -1.0

  if (strpos(filelist(0), "save") gt 0) then begin
      restore, filelist(0)
      nBLK = 1
      if (n_elements(iTime) eq 0) then begin
          p = strpos(filelist(0),".save")-1
          while (strpos(strmid(filelist(0),p,1),'.') eq -1) do p = p-1
          iYear   = fix(strmid(filelist(0),p-13,2))
          iMonth  = fix(strmid(filelist(0),p-11,2))
          iDay    = fix(strmid(filelist(0),p-9,2))
          iHour   = fix(strmid(filelist(0),p-6,2))
          iMinute = fix(strmid(filelist(0),p-4,2))
          iSecond = fix(strmid(filelist(0),p-2,2))
          iTime = [iYear, iMonth, iDay, iHour, iMinute, iSecond]
      endif
      return
  endif else begin
      if (strpos(filelist(0), "bin") gt 0) then begin
          gitm_read_bin, filelist, data, time, nVars, Vars, version
          s = size(data)
          if (s(0) eq 4) then begin
              nLons = s(2)
              nLats = s(3)
              nAlts = s(4)
              if (nLons eq 24 and nLats eq 24) then begin
                  if max(data(2,22:23,22:23,0)) eq 0.0 then begin
                      print, "Error in file.  Correcting"
                      data = data(*,0:21,0:21,*)
                      nLons = 22
                      nLats = 22
                  endif
              endif
          endif
          if (s(0) eq 3) then begin
              nLons = s(2)
              nLats = s(3)
              nAlts = 1
          endif
          nBLK = 1
          c_r_to_a, itime, time(0)
          return
      endif
  endelse

  f = filelist(0)

  if (strpos(filelist(0),"b0") eq 0) then begin
      all = findfile('b0*'+strmid(filelist(0),6,18)+'*')
  endif else begin
      all = filelist(0)
  endelse

  nfiles = n_elements(all)

  if (nBLKlat eq 0 and nBLKlon eq 0) then begin

    if (nfiles eq 1) then begin
      nBLKlon = 1
      nBLKlat = 1
    endif

 endif

 if (nBLKlat*nBLKlon eq 0) then begin

     file = all(0)
     print, "Determining Block Information from : ",file
     openr,1,file

      done = 0
      line = ''
      while not done do begin
          readf,1,line
          if strpos(line,'BLOCKS') gt -1 then begin
              nBLKlat = 0L
              nBLKlon = 0L
              nBLKalt = 0L
              readf,1, nBLKalt 
              readf,1, nBLKlat
              readf,1, nBLKlon
              done = 1
          endif

          if (eof(1)) then done = 1

      endwhile

      close,1

      if (nBLKlat*nBLKlon eq 0) then begin
          nBLK = 0
          print, "Could not determine block structure!!!"
          stop
      endif
 endif
  if (nBLKlat*nBLKlon gt nfiles) then begin
    nBLK = 0
    mess = 'There are not enough files to fill the blocks! Blocks:'+ $
           string(nBLKlat*nBLKlon)+'  files:'+string(nfiles)
    print, mess
    stop
  endif else begin
    nBLK = 1
    for n=0,nBLKlat*nBLKlon-1 do begin
      file = all(n)
      print, "reading file : ",file
      openr,1,file
      done = 0
      line = ''
      while not done do begin
          readf,1,line
          if strpos(line,'NUMERICAL') gt -1 then begin
              nvars = 0L
              nlines = 0L
              readf,1, nvars
              readf,1, nalts
              readf,1, nlats
              readf,1, nlons
          endif
          if strpos(line,'BLOCKS') gt -1 then begin
              nBLKlat = 0L
              nBLKlon = 0L
              nBLKalt = 0L
              readf,1, nBLKalt 
              readf,1, nBLKlat
              readf,1, nBLKlon
          endif
          if strpos(line,'TIME') gt -1 then begin
              iYear   = 0
              iMonth  = 0
              iDay    = 0
              iHour   = 0
              iMinute = 0
              iSecond = 0
              readf,1, iYear
              readf,1, iMonth
              readf,1, iDay
              readf,1, iHour
              readf,1, iMinute
              readf,1, iSecond
              iTime = [iYear, iMonth, iDay, iHour, iMinute, iSecond]
          endif
          if strpos(line,'VERSION') gt -1 then begin
              readf,1,Version
          endif
          if strpos(line,'VARIABLE') gt -1 then begin
              if n_elements(nvars) eq 0 then begin
                  print, 'File is in the wrong order, NUMERICAL VALUES must come'
                  print, 'before VARIABLE LIST'
                  stop
              endif else begin
                  vars = strarr(nvars)
                  for i=0,nvars-1 do begin
                      readf,1,format="(I7,a)",j,line
                      vars(i) = line
                  endfor
              endelse
              ;Store the value of vars into the pointer. 
          endif
          if strpos(line,'BEGIN') gt -1 then done = 1
      endwhile
      if (n eq 0) then begin
          nlo = (nlons-4) * nBLKlon + 4
          nla = (nlats-4) * nBLKlat + 4
          data = fltarr(nvars, nlo, nla, nalts)
          tmp = fltarr(nvars)
          format = '('+tostr(nvars)+'E11.3)'
      endif
      
      if (n_elements(iTime) eq 0) then begin
          p = strlen(filelist(0))-1
          while (strpos(strmid(filelist(0),p,1),'.') eq -1) do p = p-1
          iYear   = fix(strmid(filelist(0),p-13,2))
          iMonth  = fix(strmid(filelist(0),p-11,2))
          iDay    = fix(strmid(filelist(0),p-9,2))
          iHour   = fix(strmid(filelist(0),p-6,2))
          iMinute = fix(strmid(filelist(0),p-4,2))
          iSecond = fix(strmid(filelist(0),p-2,2))
          iTime = [iYear, iMonth, iDay, iHour, iMinute, iSecond]
      endif

      line = ''

      for k = 0, nalts-1 do begin
          for j = 0, nlats-1 do begin
              for i = 0, nlons-1 do begin
                  readf,1,tmp, format=format
                  ii = (n mod nBLKlon)*(nlons-4) + i
		  jj = (n / nBLKlon)*(nlats-4) + j
                  if (i ge 2 and i le nlons-3 and $
                      j ge 2 and j le nlats-3) then $
                    data(*,ii,jj,k) = tmp
                  if (jj lt 2) then data(*,ii,jj,k) = tmp
                  if (ii lt 2) then data(*,ii,jj,k) = tmp
                  if (jj gt nla-3) then data(*,ii,jj,k) = tmp
                  if (ii gt nlo-3) then data(*,ii,jj,k) = tmp
              endfor
          endfor	
      endfor
      close,1
    endfor
    nlons = nlo
    nlats = nla

endelse
end

