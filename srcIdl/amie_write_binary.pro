
pro amie_write_binary, filename, Vars, lats, mlts, time, data, $
                       imf = imf, ae = ae, dst = dst, hpi = hpi, cpcp = cpcp, $
                       Version = Version

  ; IMF should be (nTimes,4) - V, Bx, By, Bz
  ; AE  (nTimes, 4) - AL, AU, AE, AEI?
  ; Dst (nTimes, 2) - Dst, Dsti?
  ; Hpi (nTimes, 2) - HP, Joule Heating
  ; CPCP (nTimes) - CPCP

  nTimes = n_elements(time)

  if (n_elements(imf) lt nTimes*4 and n_elements(imf) gt 0) then begin
     print, "If IMF is provided, it should be (nTimes,4)"
     stop
  endif

  if (n_elements(imf) eq 0) then imf = fltarr(nTimes,4)

  if (n_elements(ae) lt nTimes*4 and n_elements(ae) gt 0) then begin
     print, "If AE is provided, it should be (nTimes,4)"
     stop
  endif

  if (n_elements(ae) eq 0) then ae = fltarr(nTimes,4)

  if (n_elements(dst) lt nTimes*2 and n_elements(dst) gt 0) then begin
     print, "If Dst is provided, it should be (nTimes,2)"
     stop
  endif

  if (n_elements(dst) eq 0) then dst = fltarr(nTimes,2)

  if (n_elements(HPI) lt nTimes*2 and n_elements(HPI) gt 0) then begin
     print, "If HPI is provided, it should be (nTimes,2)"
     stop
  endif

  if (n_elements(hpi) eq 0) then hpi = fltarr(nTimes,2)

  if (n_elements(cpcp) lt nTimes and n_elements(cpcp) gt 0) then begin
     print, "If CPCP is provided, it should be (nTimes)"
     stop
  endif

  if (n_elements(cpcp) eq 0) then cpcp = fltarr(nTimes)

  if (n_elements(version) eq 0) then Version = 0.0

  nVars = n_elements(vars)

  sLats = size(lats)
  if (sLats(0) eq 1) then begin
     nLats = sLats(1)
     colats = 90.0 - lats
  endif
  if (sLats(0) ge 2) then begin
     nLats = sLats(2)
     colats = 90.0 - reform(lats(0,*))
  endif
  sMlts = size(mlts)
  nMlts = sMlts(1)

  print,'writing file : ',filename
  close,1
  openw,1,filename, /f77

  writeu,1, long(nLats), long(nMlts), long(nTimes)

  writeu,1, colats
  writeu,1, mlts(*,0)

  writeu,1, nVars

  tmp = '                                  '
  for i=0,nVars-1 do writeu,1,chopl(Vars(i)+tmp,30)

  print, nTimes

  for i=0,nTimes-1 do begin

     c_r_to_a, itime, time(i)

     writeu, 1, long(i), long(itime(0:4))

     writeu, 1, imf(i,*), reverse(ae(i,*)), $
             reverse(dst(i,*)), hpi(i,*), cpcp(i)
        
     for iVar = 0, nVars-1 do begin
        d = reform(data(i,iVar,*,*))
        writeu, 1, d
     endfor

  endfor

  writeu,1,version

  close,1

end

