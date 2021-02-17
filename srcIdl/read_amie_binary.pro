pro read_amie_binary, amiefile, data, lats, mlts, ut, fields, imf,	$
                      ae, dst, hp, pot, version, date = date,	$
                      ltpos = ltpos, lnpos = lnpos, 			$
                      plotapot = plotapot, speed = speed, by = by, 	$
                      bz = bz, field = field

  if n_elements(plotapot) eq 1 then plotapot = 1 else plotapot = 0
  if n_elements(field) eq 0 then field = -1

  if (strpos(amiefile,"save") gt 0) then begin
      restore, amiefile
      ut = time
  endif else begin

      openr,1,amiefile, /f77

      nlats = 0L
      nmlts = 0L
      ntimes = 0L

      n = 0L
      iyr_tmp = 0L
      imo_tmp = 0L
      ida_tmp = 0L
      ihr_tmp = 0L
      imi_tmp = 0L
      nfields = 0L

      swv_tmp = 0.0
      bx_tmp  = 0.0
      by_tmp  = 0.0
      bz_tmp  = 0.0
      aei_tmp = 0.0
      ae_tmp  = 0.0
      au_tmp  = 0.0
      al_tmp  = 0.0
      dsti_tmp= 0.0
      dst_tmp = 0.0
      hpi_tmp = 0.0
      sjh_tmp = 0.0
      pot_tmp = 0.0

      readu,1,nlats,nmlts,ntimes

      clats = fltarr(nlats)
      mlts = fltarr(nmlts)

      readu,1,clats
      readu,1,mlts

      readu,1,nfields

      ut  = dblarr(ntimes)
      imf = fltarr(ntimes,4)
      ae  = fltarr(ntimes,4)
      dst = fltarr(ntimes,2)
      hp  = fltarr(ntimes,2)
      pot = fltarr(ntimes)
      data  = fltarr(ntimes,nfields,nmlts,nlats)
      dummy = fltarr(nlats,nmlts)
      dummy = fltarr(nmlts,nlats)
      fields = strarr(nfields)

      tmp = bytarr(30)

      for i=0,nfields-1 do begin
          readu,1,tmp
          fields(i) = string(tmp)
      endfor

      for i=0,ntimes-1 do begin

          readu,1,n,iyr_tmp,imo_tmp,ida_tmp,ihr_tmp,imi_tmp

          itime = [fix(iyr_tmp),imo_tmp,ida_tmp,ihr_tmp,imi_tmp,0]

          c_a_to_r, itime, rtime
          ut(i) = rtime

          readu,1,swv_tmp,bx_tmp,by_tmp,bz_tmp, $
            aei_tmp,ae_tmp,au_tmp,al_tmp,    $
            dsti_tmp,dst_tmp,hpi_tmp,sjh_tmp,pot_tmp

          imf(i,0) = bx_tmp
          imf(i,1) = by_tmp
          imf(i,2) = bz_tmp
          imf(i,3) = swv_tmp

          ae(i,0)  = ae_tmp 
          ae(i,1)  = au_tmp 
          ae(i,2)  = al_tmp 
          ae(i,3)  = aei_tmp 

          dst(i,0) = dst_tmp
          dst(i,1) = dsti_tmp

          hp(i,0)  = hpi_tmp
          hp(i,1)  = sjh_tmp

          pot(i)   = pot_tmp

          for j=0,nfields-1 do begin
              readu,1,dummy
              data(i,j,*,*) = dummy
              if (strpos(mklower(fields(j)),'potential') gt -1) then  $
                data(i,j,*,*) = dummy/1000.0
          endfor

      endfor

      version = 2.0
      if (not eof(1)) then readu,1,version

      close,1

      lats = 90.0 - clats

  endelse

  if plotapot then begin

    if (field lt 0) then begin
      for i=0,nfields-1 do print, tostr(i+1)+'. '+fields(i)
      type = fix(ask('field to plot','1'))-1
      if (type lt 0) or (type gt nfields-1) then type = 0
    endif else type = field

    data = reform(data(*,type,*,*))
    fields = fields(type)

    ltpos = fltarr(nmlts,nlats)
    lnpos = fltarr(nmlts,nlats)
    for i=0,nmlts-1 do ltpos(i,*) = lats
    for j=0,nlats-1 do lnpos(*,j) = mlts*360.0/24.0

    date = strarr(ntimes)
    time = strarr(ntimes)

    for i=0,ntimes-1 do begin
      c_r_to_a, itime, ut(i)
      c_a_to_s, itime, stime
      time(i) = strmid(stime,10,5)
      syear = fix(strmid(stime,7,2))
      if syear lt 65 then syear = syear + 2000 else syear = syear + 1900
      date(i) = strmid(stime,3,1)+mklower(strmid(stime,4,2))+' '+	$
                strmid(stime,0,2)+', '+tostr(syear)
    endfor

    ut = time
    speed = reform(imf(*,3))
    by    = reform(imf(*,1))
    bz    = reform(imf(*,2))

    lats = nlats
    mlts = nmlts

  endif else begin

    if field gt -1 then data = reform(data(*,field,*,*))

    nFields = n_elements(fields)
    nLats = n_elements(lats)
    nMlts = n_elements(mlts)
    nTimes = n_elements(ut)

    iEField_ = -1
    iPed_    = -1
    iSJH_    = -1

    for iField = 0, nFields-1 do begin
        fields(iField) = mklower(fields(iField))
        if (strpos(fields(iField),'electric field') gt -1 and $
            iEField_ eq -1) then iEField_ = iField
        if (strpos(fields(iField),'ped') gt -1 and $
            strpos(fields(iField),'total') gt -1 and $
            iPed_ eq -1) then iPed_ = iField
        if (strpos(fields(iField),'joule heating') gt -1 and $
            iSJH_ eq -1) then iSJH_ = iField
    endfor

    if (iEField_ gt -1 and iPed_ gt -1) then begin

        nFields = nFields + 1
        fields = [fields, "field-aligned current (!Mm!XA/m!E2!N)"]

        data_add = fltarr(nTimes, nFields, nMlts, nLats)

        theta = (90.0 - lats)*!dtor
        phi   = mlts * !pi/12.0
        re    = 6372000.0

        fac = fltarr(nMlts, nLats)
        for iTime = 0, nTimes-1 do begin
            data_add(iTime, 0:nFields-2, *, *) = data(iTime, 0:nFields-2, *, *)
            ee = reform(data(iTime, iEField_+0, *, *))
            en = reform(data(iTime, iEField_+1, *, *))
            sp = reform(data(iTime, iPed_,      *, *))
            je = ee * sp
            jn = -en * sp 
            for iMlts=0,nMlts-1 do jn(iMlts,*) = jn(iMlts,*) * sin(theta)
            fac(*,*) = 0.0

                                ; North Component
            for iMlts=0,nMlts-2 do begin
                fac(iMlts,1:nLats-2) = $
                  (jn(iMlts,2:nLats-1) - jn(iMlts,0:nLats-3)) / $
                  (theta(2:nLats-1)-theta(0:nLats-3))/re
                fac(iMlts,nLats-1)   = $
                  (jn(iMlts,nLats-1) - jn(iMlts,nLats-2)) / $
                  (theta(nLats-1)-theta(nLats-2))/re
            endfor

        ; East Component
            for iLat = 1, nLats-1 do begin
                f = max([sin(theta),0.1])
                fac(1:nMlts-2, iLat) = fac(1:nMlts-1, iLat) + $
                  (je (2:nMlts-1,iLat) - je (0:nMlts-3,iLat)) / $
                  (phi(2:nMlts-1)      - phi(0:nMlts-3))/(re*f)
            ; across the 0-24 border
                fac(0, iLat) = fac(0, iLat) + $
                  (je (1,iLat) - je (nMlts-2,iLat)) / $
                  (phi(2)      - phi(0))/(re*f)
            endfor

            fac(nMlts-1, *) = fac(0, *)
            data_add(iTime, nFields-1, *, *) = -fac*1.0e6

            if (version lt 2.2) then $
              data_add(iTime, iSJH_, *, *) = sp * (ee^2 + en^2) * 1000.0

        endfor

        data = data_add

    endif

  endelse

  return

end


