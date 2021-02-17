
pro thermo_readlog, filelist, data, time

  nFiles = n_elements(filelist)

  nLines = 0L

  for i=0,nFiles-1 do begin

     spawn, 'wc '+filelist(i), nL
     nLines = nLines+long(nL)-2

  endfor

  line = ''

  iLine = 0L

  for i=0,nFiles-1 do begin

     openr, 1, filelist(i)

     while (strpos(line,"START") lt 0) do readf,1,line
     readf,1,line

     if (i eq 0) then begin

        l = strsplit(line)
        nVars = n_elements(l)
        l = [l,strlen(line)]
        Vars = strarr(nVars)
      
        for iVar = 0, nVars-1 do $
           Vars(iVar) = strmid(line,l(iVar),l(iVar+1)-l(iVar))

        data = fltarr(nVars,nLines)
        tmp = fltarr(nVars)

     endif

     while not eof(1) do begin

        readf,1,tmp
        data(*,iLine) = tmp
        iLine = iLine + 1L

     endwhile
     
     close,1

  endfor

  nLines = iLine
  data = data(*,0:nLines-1)
  time = dblarr(nLines)

  for i=0L,nLines-1 do begin

     itime = fix(data(1:6,i))
     c_a_to_r, itime, rtime
     time(i) = rtime

  endfor

end
