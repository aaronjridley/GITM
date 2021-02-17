
pro thermo_plotvectors,vars,k,data, lat,lon, utrot, alt, nlats,nlons, nalts, $
                       cf,vi_cnt,vn_cnt,step, polar, maxran, plane, npolar

  ;Calculate the factor
  count = 0
  vieast_index = -1
  vneast_index = -1
  if (vi_cnt eq 1) then begin
     count = n_elements(vars)
     test = "Vi (east)"
     for i=0,count-1 do begin
        var=strcompress(vars[i],/remove_all)
        tes = strcompress(test,/remove_all)
        result= STRCMP( var, tes, 8)
        if (result eq 1) then vieast_index = i
     endfor
     if (vieast_index lt 0) then begin ; 3DALL files
        test = "V!Di!N(east)"
        for i=0,count-1 do begin
           var=strcompress(vars[i],/remove_all)
           tes = strcompress(test,/remove_all)
           result= STRCMP( var, tes, 8)
           if (result eq 1) then vieast_index = i
        endfor
     endif
     if (vieast_index lt 0) then begin ; 3DLST files
        test = "Vn(east)(m/s)"
        for i=0,count-1 do begin
           var=strcompress(vars[i],/remove_all)
           tes = strcompress(test,/remove_all)
           result= STRCMP( var, tes, 8)
           if (result eq 1) then vieast_index = i
        endfor
     endif
  endif
  if (vn_cnt eq 1) then begin
     count = n_elements(vars)
     test = "Vn (east)"
     for i=0,count-1 do begin
        var=strcompress(vars[i],/remove_all)
        tes = strcompress(test,/remove_all)
        result= STRCMP( var, tes,8 )
        if (result eq 1) then vneast_index = i
     endfor

     if (vneast_index lt 0) then begin
        test = "V!Dn!N(east)"
        for i=0,count-1 do begin
           var=strcompress(vars[i],/remove_all)
           tes = strcompress(test,/remove_all)
           result= STRCMP( var, tes,8 )
           if (result eq 1) then vneast_index = i
        endfor
     endif
  endif

  factors = [1.0, 5.0, 10.0, 20.0, 25.0, $
             50.0, 75.0, 100.0, 150.0, 200.0]

  if (vi_cnt eq 1) then vindex = vieast_index  $
  else vindex = vneast_index

  if (cf lt 0) then begin
     dist = max(abs(data(vindex,*,*,k))) - factors*10.0
     loc = where(dist lt 0.0, cf)
;     cf = 1
;     mindist = dist(cf-1)
;     for i=1,n_elements(factors)-1 do begin
;        if (dist(i) gt 0 and dist(i) lt mindist) then begin
;           cf = i+1
;           mindist = dist(i)
;        endif
;     endfor
  endif

  factor = factors(cf-1)

  if (plane eq 1) then begin

     for i =0,nlats-1,step do begin

        if (npolar) then la = lat(0,i,k) else la = -lat(0,i,k)

        if (90-la lt maxran) then begin

           for j =0,nlons-1,step do begin
              lo = lon(j,i,k) + utrot

              ux = data(vindex,j,i,k)/factor
              if (polar) then $
                 uy = data((vindex+1),j,i,k)/factor*lat(0,i,k)/abs(lat(0,i,k)) $
              else $
                 uy = data((vindex+1),j,i,k)/factor

              if (polar) then begin
                 x = (90.0 - la) * cos(lo*!pi/180.0 - !pi/2.0)
                 y = (90.0 - la) * sin(lo*!pi/180.0 - !pi/2.0)
                 
                 ulo = ux
                 ula = uy
                      
                 ux = - ula * cos(lo*!pi/180.0 - !pi/2.0)  $
                      - ulo * sin(lo*!pi/180.0 - !pi/2.0)
                 uy = - ula * sin(lo*!pi/180.0 - !pi/2.0) $
                      + ulo * cos(lo*!pi/180.0 - !pi/2.0)

              endif else begin
                 x = lo
                 y = la
              endelse

                  ;ux is the eastward welocity (neutral or ion)
                  ;uy is the northward velocity (neutral or ion)
;                  oplot,[x],[y],psym = 4, color = 0
              oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

              u = sqrt(ux^2+uy^2)
              if (u gt 0) then begin
                 t = asin(uy/u)
                 if (ux lt 0) then t = !pi-t
                 t1 = t+!pi/12
                 t2 = t-!pi/12
                 ux1 = 0.6 * u * cos(t1)
                 uy1 = 0.6 * u * sin(t1)
                 ux2 = 0.6 * u * cos(t2)
                 uy2 = 0.6 * u * sin(t2)
                 oplot,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
                 oplot,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0
              endif
              
           endfor

        endif

     endfor

     if (polar) then begin

        x  =  maxran
        y  =  maxran
        uy = -10.0
        ux =   0.0
        
        plots,[x,x+ux],[y,y+uy]
          
        u = sqrt(ux^2+uy^2)
        t = asin(uy/u)
        if (ux lt 0) then t = !pi-t
        t1 = t+!pi/12
        t2 = t-!pi/12
        ux1 = 0.6 * u * cos(t1)
        uy1 = 0.6 * u * sin(t1)
        ux2 = 0.6 * u * cos(t2)
        uy2 = 0.6 * u * sin(t2)
        plots,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
        plots,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0

        str = string(factor*10.0, format = '(f6.1)') + ' m/s'
        xyouts, x-1.0, y+uy/2.0, str, alignment = 1.0

     endif else begin

        x  =  360.0
        y  =   95.0
        uy =    0.0
        ux =  -10.0

        plots,[x,x+ux],[y,y+uy]

        u = sqrt(ux^2+uy^2)
        t = asin(uy/u)
        if (ux lt 0) then t = !pi-t
        t1 = t+!pi/12
        t2 = t-!pi/12
        ux1 = 0.6 * u * cos(t1)
        uy1 = 0.6 * u * sin(t1)
        ux2 = 0.6 * u * cos(t2)
        uy2 = 0.6 * u * sin(t2)

        plots,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
        plots,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0

        str = string(factor*10.0, format = '(f6.1)') + ' m/s'
        xyouts, x, y+5.0, str, alignment = 1.0

     endelse

  endif

  if (plane eq 2) then begin

     for i =0,nlons-1,step do begin

        lo = lon(i,0,k)
        
        for j =0,nalts-1,step*2 do begin
           al = alt(i,k,j)

           ux = data(vindex,i,k,j)/factor
           uy = data((vindex+2),i,k,j)/factor

           x = lo
           y = al

           ;ux is the eastward welocity (neutral or ion)
           ;uy is the northward velocity (neutral or ion)
           oplot,[x],[y],psym = 4, color = 0
           oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

        endfor

     endfor

     x  =  360.0
     y  =  max(alt)*1.01
     uy =    0.0
     ux =  -10.0

     plots,[x],[y],psym = 4
     plots,[x,x+ux],[y,y+uy]

     str = string(factor*10.0, format = '(f6.1)') + ' m/s'
     xyouts, x, y+5.0, str, alignment = 1.0

  endif

  if (plane eq 3) then begin

     for i =0,nlats-1,step do begin

        la = lat(k,i,0)

        for j =0,nalts-1,step*2 do begin

           al = alt(k,i,j)

           ux = data(vindex+1,k,i,j)/factor
           uy = data((vindex+2),k,i,j)/factor

           x = la
           y = al

           ;ux is the eastward welocity (neutral or ion)
           ;uy is the northward velocity (neutral or ion)
           oplot,[x],[y],psym = 4, color = 0
           oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

        endfor

     endfor

     x  =   90.0
     y  =  max(alt)*1.01
     uy =    0.0
     ux =  -10.0

     plots,[x],[y],psym = 4
     plots,[x,x+ux],[y,y+uy]

     str = string(factor*10.0, format = '(f6.1)') + ' m/s'
     xyouts, x, y+5.0, str, alignment = 1.0

  endif

end

