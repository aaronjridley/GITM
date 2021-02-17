;This program can do vectors. User select factor from the slider bar.
;Select vector step from the slider bar. 
;Polar Plot and select MinLat to do polar plot.

pro thermo_gui_event,event

widget_control,event.top,get_uvalue=ptr,/no_copy
widget_control,event.id,get_uvalue=whichevent

case whichevent of 
    'ROWBLOCK':begin
        widget_control,(*ptr).rb_txt,get_value=txt
        (*ptr).rb = txt[0]
        print,(*ptr).rb
    end
    'COLBLOCK':begin
        widget_control,(*ptr).cb_txt,get_value=txt
        (*ptr).cb = txt[0]
        print,(*ptr).cb
    end
    'FILENAME':begin
        widget_control,(*ptr).file_txt,get_value=txt
        (*ptr).filename=txt[0]
        print,(*ptr).filename
    end
    'SELVARS':begin
        sel=event.index
        (*ptr).sel=sel
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt
    end
    'SEARCH':begin
        (*ptr).cursor_cnt=0
        widget_control,/hourglass
                                ;Get the file
        widget_control,(*ptr).file_txt,get_value=txt
        (*ptr).filename=txt[0]

        widget_control,(*ptr).cb_txt,get_value=txt
        (*ptr).cb = txt[0]
        widget_control,(*ptr).rb_txt,get_value=txt
        (*ptr).rb = txt[0]
        rb = fix((*ptr).rb)
        cb = fix((*ptr).cb)
        
        file=(*ptr).filename
        if file eq '' then begin
            mes=widget_message('Please enter the file to search!')
        endif else begin
            widget_control,/hourglass
;            print,'file is:',file
            filelist = findfile(file)
            nfiles = n_elements(filelist)
            (*ptr).nfiles=nfiles
;            print,'filelist is:',filelist
;            print,'nfiles:',nfiles
            if filelist[0] eq '' or nfiles eq 0 then begin
                me=widget_message('There is no such a file!')
            endif else begin
                read_thermosphere_file, filelist, nvars, nalts, nlats, nlons,vars,data,rb,cb,bl_cnt, iTime
                (*ptr).bl_cnt = bl_cnt
                if (*ptr).bl_cnt eq 1 then begin
                                ;Store the array 'data' to the '*(*ptr).data'
                    *(*ptr).vars=vars
                    *(*ptr).data=data
                    *(*ptr).itime=itime

                                ;Store the nvars,nalts,nlats,nlons to the pointer.
                    (*ptr).nvars=nvars
                    (*ptr).nalts=nalts
                    (*ptr).nlats=nlats
                    (*ptr).nlons=nlons
;                    print,'nvars,nalts,nlats,nlons:'
;                    print,nvars,nalts,nlats,nlons
                endif
            endelse
        endelse	
        if (*ptr).bl_cnt eq 1 then begin
                                ;Set the VARS to the vars_list
            widget_control,(*ptr).vars_list,set_value=*(*ptr).vars
            
                                ;Reset the maximum number of the sliders
            cnt1=(*ptr).lat_lon_cnt
            cnt2=(*ptr).alt_lon_cnt
            cnt3=(*ptr).alt_lat_cnt
            if (n_elements(nalts) eq 0) then nalts = 4
            if (n_elements(nlats) eq 0) then nlats = 4
            if (n_elements(nlons) eq 0) then nlons = 4
            if cnt1 eq 1 then begin
                widget_control,(*ptr).sli,set_slider_max=nalts-1
            endif
            if cnt2 eq 1 then begin
                widget_control,(*ptr).sli,set_slider_max=nlats-1
            endif
            if cnt3 eq 1 then begin
                widget_control,(*ptr).sli,set_slider_max=nlons-1
            endif
        endif
    end
    'MULTIFILES':begin
        (*ptr).cursor_cnt=0
        widget_control,/hourglass
                                ;Get the file
        widget_control,(*ptr).file_txt,get_value=txt
        (*ptr).filename=txt[0]

        widget_control,(*ptr).cb_txt,get_value=txt
        (*ptr).cb = txt[0]
        widget_control,(*ptr).rb_txt,get_value=txt
        (*ptr).rb = txt[0]
        rb = (*ptr).rb
        cb = (*ptr).cb
        
;        print,'rb & cb:',(*ptr).rb,(*ptr).cb
        file=(*ptr).filename
        if file eq '' then begin
            mes=widget_message('Please enter the file to search!')
        endif else begin
            widget_control,/hourglass
            filelist = findfile(file)
            nfiles = n_elements(filelist)
            (*ptr).nfiles=nfiles
            if filelist(0) eq '' or nfiles eq 0 then begin
                me=widget_message('There is no such a file!')
            endif else begin
              ;Find the multiple files according to the filename user entered. 
                if (strpos(filelist(0), "save") gt 0) then begin
                    pos 	 = strpos(file,"_")
                    fil 	 = strtrim(strmid(file,0,pos))
                    filsear  = "*"+".save"
;                    print,'filsear is:',filsear
                    mulfiles = findfile(filsear)
                    n 	 = n_elements(mulfiles)	
;                    print,'mulfiles:',mulfiles
;                    print,'n is:',n
                    
                    bl_cnt = 0
                                ;Read multiple files.
                    read_thermosphere_file, mulfiles(0), nvars, nalts, nlats, $
                      nlons,vars,data,rb,cb,bl_cnt, iTime
                    
                    sel = (*ptr).sel
                    
                    if sel eq -1 then begin

                    endif else begin
                                ;Get cursor_x & cursor_y from widget_text.
                        widget_control,(*ptr).curx_txt,get_value=txt
                        strx	    =	txt[0]
                        cursor_x =   double(strx)
                        widget_control,(*ptr).cury_txt,get_value=txt
                        stry	    =   txt[0]
                        cursor_y =   double(stry)
                        
                                ;Set the value of step. 
                        if ((*ptr).vyes_cnt eq 1) then begin
                            step = (*ptr).vs
                        endif
                        
                                ;Clean up the display window
                        plotdumb

                                ;Initialize values
                        vars	  = *(*ptr).vars ;Get the variables' list
                        sel	  = (*ptr).sel ;Get the selected variable.
                        nfiles = (*ptr).nfiles ;Get the files users searched for 
                        cnt1   = (*ptr).lat_lon_cnt ;Get the counter of selecting lat_lon button
                        cnt2   = (*ptr).alt_lon_cnt ;Get the counter of selecting alt_lon button
                        cnt3   = (*ptr).alt_lat_cnt ;Get the counter of selecting alt_lat button 
                        yes	  = (*ptr).yes_cnt ;Get the counter of selecting yes (Ghost Cells) button 
                        no	  = (*ptr).no_cnt ;Get the counter of selecting no (Ghost Cells) button
                        yeslog = (*ptr).yes_logcnt ;Get the counter of selecting yes Plot Log button
                        nolog  = (*ptr).no_logcnt ;Get the counter of selecting no Plot Log button
                        nalts  = (*ptr).nalts ;Get the number of nalts
                        nlats  = (*ptr).nlats ;Get the number of nlats
                        nlons  = (*ptr).nlons ;Get the number of nlons
                        
                        yeswrite_cnt = (*ptr).yeswrite_cnt
                        
                        polar  = (*ptr).polar
                        npolar = (*ptr).npolar
                        if (polar eq 1) and (npolar eq 1) then begin
                            MinLat = (*ptr).nplat
                        endif
                        if (polar eq 1) and (npolar eq 0) then begin
                            MinLat = (*ptr).splat
                        endif
                        
                        showgridyes   = (*ptr).gridyes_cnt
                        plotvectoryes = (*ptr).vyes_cnt
                        vi_cnt 	 = (*ptr).vi_cnt
                        vn_cnt 	 = (*ptr).vn_cnt
                        cf 		 = (*ptr).cf
                        vs 		 = (*ptr).vs
                        cursor_cnt    = (*ptr).cursor_cnt

                                ;Counter 'cnt' is to count if any of lat_lon,alt_lon,alt_lat buttons selected.
                        cnt=1
                        if (cnt1 eq 0) and (cnt2 eq 0) and (cnt3 eq 0) then begin
                            cnt=0
                        endif
                        
                   ;Check if there are any variables in the text_list. 
                   ;If there is not, then initialize the counter of selecting 
              ;plot button (plot_cnt) to zero and pop out a warning message. 
              ;If there are variables, then start to do the plot.  
                        if (vars(0) eq '') or (sel eq -1) or (cnt eq 0) then begin
                            (*ptr).plot_cnt  =  0
                            me  =  widget_message('Please do search, select a variable and a setting first!')
                        endif else begin
                            (*ptr).plot_cnt  =  1 ;Set the plot counter to 1.
                                ;Initialize the data,alt,lat,lon,k,xrange,yrange
                            data   = *(*ptr).data ;Get the array data which be plotted. 
                            alt    = reform(data(2,*,*,*))/1000.0 ;Calculate alt,lat,lon and k.
                            lat    = reform(data(1,*,*,*))*180.0/!pi
                            lon    = reform(data(0,*,*,*))*180.0/!pi
                            k	   = n_elements(alt(0,0,*))-1
                            xrange = [0,0] ;Initialize xrange and yrange.
                            yrange = [0,0]

                            selset = 0
                            if cnt1 eq 1 then begin
                                selset  = (*ptr).lat_lon_max
                            endif
                            if cnt2 eq 1 then begin
                                selset  = (*ptr).alt_lon_max
                            endif
                            if cnt3 eq 1 then begin
                                selset  = (*ptr).alt_lat_max
                            endif

                            widget_control,(*ptr).min_txt,get_value=smini
                            widget_control,(*ptr).max_txt,get_value=smaxi
                            
                            smini = smini
                            smaxi = smaxi
                            
                            filename = mulfiles(0)

                            widget_control,(*ptr).writefile_txt,get_value=txt
                            psfile=strtrim(string(txt[0]),2)
                            mars = (*ptr).mars
                            colortable = (*ptr).colortable
                            itime = *(*ptr).itime

                                ;Call the subprogram to do the plots.
                            thermo_plot,cursor_x,cursor_y,strx,stry,step,vars,sel,nfiles,	  $
                              cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
                              nolog,nalts,nlats,nlons,yeswrite_cnt,$
                              polar,npolar,MinLat,showgridyes,	  $
                              plotvectoryes,vi_cnt,vn_cnt,cf,	  $
                              cursor_cnt,data,alt,lat,lon,	  $
                              xrange,yrange,selset,smini,smaxi,	  $
                              filename,vars, psfile,mars, colortable, iTime


                            
                        endelse
;                        print,'smini,smaxi:',smini,smaxi
                        mini = smini
                        maxi = smaxi
                                ;Show min and max number to the interface.
                        min=strtrim(string(mini),2)
                        max=strtrim(string(maxi),2)
                        widget_control,(*ptr).min_txt,set_value=min
                        widget_control,(*ptr).max_txt,set_value=max
                                ;Reset set cursor counter to zero.
                        (*ptr).cursor_cnt=0
                    endelse

                    (*ptr).bl_cnt = bl_cnt
                    if (bl_cnt eq 1) then begin
                                ;Store the array 'data' to the '*(*ptr).data'
                        *(*ptr).vars=vars
                        *(*ptr).data=data

                                ;Store the nvars,nalts,nlats,nlons to the pointer.
                        (*ptr).nvars=nvars
                        (*ptr).nalts=nalts
                        (*ptr).nlats=nlats
                        (*ptr).nlons=nlons
;                        print,'nvars,nalts,nlats,nlons:'
;                        print,nvars,nalts,nlats,nlons

                                ;Set the VARS to the vars_list
                        widget_control,(*ptr).vars_list,set_value=*(*ptr).vars
                        
                                ;Reset the maximum number of the sliders
                        cnt1=(*ptr).lat_lon_cnt
                        cnt2=(*ptr).alt_lon_cnt
                        cnt3=(*ptr).alt_lat_cnt

                        if cnt1 eq 1 then begin
                            widget_control,(*ptr).sli,set_slider_max=nalts-1
                        endif
                        if cnt2 eq 1 then begin
                            widget_control,(*ptr).sli,set_slider_max=nlats-1
                        endif
                        if cnt3 eq 1 then begin
                            widget_control,(*ptr).sli,set_slider_max=nlons-1
                        endif
                        if n gt 1 then begin
                            for i = 0,n-1 do begin
                                print,'In the',i,' Loop!'
                                print,'mulfiles(i):',mulfiles(i)
                                read_thermosphere_file, mulfiles(i), nvars, nalts,$
                                  nlats, nlons,vars,data,rb,cb,bl_cnt, iTime
                                ;Store the data to pointer
                                *(*ptr).data = data
                                *(*ptr).itime = itime

                                if sel eq -1 then begin

                                endif else begin						
                                ;Get cursor_x & cursor_y from widget_text.
                                    widget_control,(*ptr).curx_txt,get_value=txt
                                    strx	    =	txt[0]
                                    cursor_x =   double(strx)
                                    widget_control,(*ptr).cury_txt,get_value=txt
                                    stry	    =   txt[0]
                                    cursor_y =   double(stry)
                                    
                                ;Set the value of step. 
                                    if ((*ptr).vyes_cnt eq 1) then begin
                                        step = (*ptr).vs
                                    endif
                                    
                                ;Clean up the display window
                                    plotdumb

                                ;Initialize values
                                    vars	  = *(*ptr).vars
                                    sel	  = (*ptr).sel
                                    nfiles = (*ptr).nfiles
                                    cnt1   = (*ptr).lat_lon_cnt
                                    cnt2   = (*ptr).alt_lon_cnt
                                    cnt3   = (*ptr).alt_lat_cnt
                                    yes	  = (*ptr).yes_cnt
                                    no	  = (*ptr).no_cnt
                                    yeslog = (*ptr).yes_logcnt
                                    nolog  = (*ptr).no_logcnt
                                    nalts  = (*ptr).nalts
                                    nlats  = (*ptr).nlats
                                    nlons  = (*ptr).nlons
                                    
                                    yeswrite_cnt = (*ptr).yeswrite_cnt
                                    
                                    polar  = (*ptr).polar
                                    npolar = (*ptr).npolar
                                    if (polar eq 1) and $
                                      (npolar eq 1) then begin
                                        MinLat = (*ptr).nplat
                                    endif
                                    if (polar eq 1) and $
                                      (npolar eq 0) then begin
                                        MinLat = (*ptr).splat
                                    endif
                                    
                                    showgridyes   = (*ptr).gridyes_cnt
                                    plotvectoryes = (*ptr).vyes_cnt
                                    vi_cnt 	 = (*ptr).vi_cnt
                                    vn_cnt 	 = (*ptr).vn_cnt
                                    cf 		 = (*ptr).cf
                                    vs 		 = (*ptr).vs
                                    cursor_cnt    = (*ptr).cursor_cnt

                                ;Counter 'cnt' is to count if any of lat_lon,alt_lon,alt_lat buttons selected.
                                    cnt=1
                                    if (cnt1 eq 0) and (cnt2 eq 0) and (cnt3 eq 0) then begin
                                        cnt=0
                                    endif
                                    
                                ;Check if there are any variables in the text_list. 
                                ;If there is not, then initialize the counter of selecting 
                                ;plot button (plot_cnt) to zero and pop out a warning message. 
                                ;If there are variables, then start to do the plot.  
                                    if (vars(0) eq '') or (sel eq -1) or (cnt eq 0) then begin
                                        (*ptr).plot_cnt  =  0
                                        me  =  widget_message('Please do search, select a variable and a setting first!')
                                    endif else begin
                                        (*ptr).plot_cnt  =  1 ;Set the plot counter to 1.
                                ;Initialize the data,alt,lat,lon,k,xrange,yrange
                                        data   = *(*ptr).data ;Get the array data which be plotted. 
                                        alt    = reform(data(2,*,*,*))/1000.0 ;Calculate alt,lat,lon and k.
                                        lat    = reform(data(1,*,*,*))*180.0/!pi
                                        lon    = reform(data(0,*,*,*))*180.0/!pi
                                        k	   = n_elements(alt(0,0,*))-1
                                        xrange = [0,0] ;Initialize xrange and yrange.
                                        yrange = [0,0]

                                        selset = 0
                                        if cnt1 eq 1 then begin
                                            selset  = (*ptr).lat_lon_max
                                        endif
                                        if cnt2 eq 1 then begin
                                            selset  = (*ptr).alt_lon_max
                                        endif
                                        if cnt3 eq 1 then begin
                                            selset  = (*ptr).alt_lat_max
                                        endif

                                        widget_control,(*ptr).min_txt,get_value=smini
                                        widget_control,(*ptr).max_txt,get_value=smaxi
                                        
                                        smini = smini
                                        smaxi = smaxi
                                        
                                        filename = mulfiles[i]

                                        widget_control,(*ptr).writefile_txt,get_value=txt
                                        psfile=strtrim(string(txt[0]),2)
                                        mars = (*ptr).mars
                                        colortable = (*ptr).colortable
                                        itime = *(*ptr).itime


                                ;Call the subprogram to do the plots.
                                        thermo_plot,cursor_x,cursor_y,strx,stry,step,vars,sel,nfiles,	  $
                                          cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
                                          nolog,nalts,nlats,nlons,yeswrite_cnt,$
                                          polar,npolar,MinLat,showgridyes,	  $
                                          plotvectoryes,vi_cnt,vn_cnt,cf,	  $
                                          cursor_cnt,data,alt,lat,lon,	  $
                                          xrange,yrange,selset,smini,smaxi,	  $
                                          filename,vars, psfile,mars, $
                                          colortable, iTime


                                        
                                    endelse
;                                    print,'smini,smaxi:',smini,smaxi
                                    mini = smini
                                    maxi = smaxi
                                ;Show min and max number to the interface.
                                    min=strtrim(string(mini),2)
                                    max=strtrim(string(maxi),2)
                                    widget_control,(*ptr).min_txt,set_value=min
                                    widget_control,(*ptr).max_txt,set_value=max
                                ;Reset set cursor counter to zero.
                                    (*ptr).cursor_cnt=0

                                endelse

                            endfor
                        endif	
                    endif			
                endif else begin
                    mes = widget_message('Please enter the file with .save extension!')
                endelse
            endelse
        endelse	
    end
    'LAT_LON':begin
        (*ptr).cursor_cnt=0
        widget_control,(*ptr).sli,sensitive=1,set_slider_max=(*ptr).nalts-1
        (*ptr).lat_lon_cnt=1
        (*ptr).alt_lon_cnt=0
        (*ptr).alt_lat_cnt=0
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt

        widget_control,(*ptr).sli,get_value=max
        data=*(*ptr).data       ;Get the array data which be plotted. 
        mv = data(2,0,0,max)/1000.0
        widget_control,(*ptr).slider_txt,set_value=string(mv,format='(f6.1)')

    end
    'ALT_LON':begin
        (*ptr).cursor_cnt=0
        widget_control,(*ptr).sli,sensitive=1,set_slider_max=(*ptr).nlats-1
        (*ptr).lat_lon_cnt=0
        (*ptr).alt_lon_cnt=1
        (*ptr).alt_lat_cnt=0
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt

        widget_control,(*ptr).sli,get_value=max
        data=*(*ptr).data       ;Get the array data which be plotted. 
        mv = data(1,0,max,0)*180.0/!pi
        widget_control,(*ptr).slider_txt,set_value=string(mv,format='(f6.1)')

    end
    'ALT_LAT':begin
        (*ptr).cursor_cnt=0
        widget_control,(*ptr).sli,sensitive=1,set_slider_max=(*ptr).nlons-1
        (*ptr).lat_lon_cnt=0
        (*ptr).alt_lon_cnt=0
        (*ptr).alt_lat_cnt=1
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt

        widget_control,(*ptr).sli,get_value=max
        data=*(*ptr).data       ;Get the array data which be plotted. 
        mv = data(0,max,0,0)*180.0/!pi
        widget_control,(*ptr).slider_txt,set_value=string(mv,format='(f6.1)')

    end
    'SLI':begin
        (*ptr).cursor_cnt=0
        widget_control,(*ptr).sli,get_value=max
                                ;Get value from user regarding the button they chose.
        cnt1=(*ptr).lat_lon_cnt
        cnt2=(*ptr).alt_lon_cnt
        cnt3=(*ptr).alt_lat_cnt
        data=*(*ptr).data       ;Get the array data which be plotted. 
        if cnt1 eq 1 then begin
            (*ptr).lat_lon_max=max
            mv = data(2,0,0,max)/1000.0
        endif
        if cnt2 eq 1 then begin
            (*ptr).alt_lon_max=max
            mv = data(1,0,max,0)*180.0/!pi
        endif
        if cnt3 eq 1 then begin
            (*ptr).alt_lat_max=max
            mv = data(0,max,0,0)*180.0/!pi
        endif
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt
        widget_control,(*ptr).slider_txt,set_value=string(mv,format='(f6.1)')

                                ;Get cursor_x & cursor_y from widget_text.
        widget_control,(*ptr).curx_txt,get_value=txt
        strx	    =	txt[0]
        cursor_x =   double(strx)
        widget_control,(*ptr).cury_txt,get_value=txt
        stry	    =   txt[0]
        cursor_y =   double(stry)
        
                                ;Set the value of step. 
        if ((*ptr).vyes_cnt eq 1) then begin
            step = (*ptr).vs
        endif
        
                                ;Clean up the display window
        plotdumb

                                ;Initialize values
        vars	  = *(*ptr).vars ;Get the variables' list
        sel	  = (*ptr).sel  ;Get the selected variable.
        nfiles = (*ptr).nfiles  ;Get the files users searched for 
        cnt1   = (*ptr).lat_lon_cnt ;Get the counter of selecting lat_lon button
        cnt2   = (*ptr).alt_lon_cnt ;Get the counter of selecting alt_lon button
        cnt3   = (*ptr).alt_lat_cnt ;Get the counter of selecting alt_lat button 
        yes	  = (*ptr).yes_cnt ;Get the counter of selecting yes (Ghost Cells) button 
        no	  = (*ptr).no_cnt ;Get the counter of selecting no (Ghost Cells) button
        yeslog = (*ptr).yes_logcnt ;Get the counter of selecting yes Plot Log button
        nolog  = (*ptr).no_logcnt ;Get the counter of selecting no Plot Log button
        nalts  = (*ptr).nalts   ;Get the number of nalts
        nlats  = (*ptr).nlats   ;Get the number of nlats
        nlons  = (*ptr).nlons   ;Get the number of nlons
        
        yeswrite_cnt = (*ptr).yeswrite_cnt
        
        polar  = (*ptr).polar
        npolar = (*ptr).npolar
        if (polar eq 1) and (npolar eq 1) then begin
            MinLat = (*ptr).nplat
        endif
        if (polar eq 1) and (npolar eq 0) then begin
            MinLat = (*ptr).splat
        endif
        
        showgridyes   = (*ptr).gridyes_cnt
        plotvectoryes = (*ptr).vyes_cnt
        vi_cnt 	 = (*ptr).vi_cnt
        vn_cnt 	 = (*ptr).vn_cnt
        cf 		 = (*ptr).cf
        vs 		 = (*ptr).vs
        cursor_cnt    = (*ptr).cursor_cnt

                                ;Counter 'cnt' is to count if any of lat_lon,alt_lon,alt_lat buttons selected.
        cnt=1
        if (cnt1 eq 0) and (cnt2 eq 0) and (cnt3 eq 0) then begin
            cnt=0
        endif
        
                                ;Check if there are any variables in the text_list. 
                                ;If there is not, then initialize the counter of selecting 
                                ;plot button (plot_cnt) to zero and pop out a warning message. 
                                ;If there are variables, then start to do the plot.  
        if (vars(0) eq '') or (sel eq -1) or (cnt eq 0) then begin
            (*ptr).plot_cnt  =  0
            me  =  widget_message('Please do search, select a variable and a setting first!')
        endif else begin
            (*ptr).plot_cnt  =  1 ;Set the plot counter to 1.
                                ;Initialize the data,alt,lat,lon,k,xrange,yrange
            data   = *(*ptr).data ;Get the array data which be plotted. 
            alt    = reform(data(2,*,*,*))/1000.0 ;Calculate alt,lat,lon and k.
            lat    = reform(data(1,*,*,*))*180.0/!pi
            lon    = reform(data(0,*,*,*))*180.0/!pi
            k	   = n_elements(alt(0,0,*))-1
            xrange = [0,0]      ;Initialize xrange and yrange.
            yrange = [0,0]

            selset = 0
            if cnt1 eq 1 then begin
                selset  = (*ptr).lat_lon_max
            endif
            if cnt2 eq 1 then begin
                selset  = (*ptr).alt_lon_max
            endif
            if cnt3 eq 1 then begin
                selset  = (*ptr).alt_lat_max
            endif

            widget_control,(*ptr).min_txt,get_value=smini
            widget_control,(*ptr).max_txt,get_value=smaxi
            
            smini = smini
            smaxi = smaxi
            
            filename = (*ptr).filename

            widget_control,(*ptr).writefile_txt,get_value=txt
            psfile=strtrim(string(txt[0]),2)
            mars = (*ptr).mars
            colortable = (*ptr).colortable
            itime = *(*ptr).itime
        

                                ;Call the subprogram to do the plots.
            thermo_plot,cursor_x,cursor_y,strx,stry,step,vars,sel,nfiles,	  $
              cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
              nolog,nalts,nlats,nlons,yeswrite_cnt,$
              polar,npolar,MinLat,showgridyes,	  $
              plotvectoryes,vi_cnt,vn_cnt,cf,	  $
              cursor_cnt,data,alt,lat,lon,	  $
              xrange,yrange,selset,smini,smaxi,	  $
              filename,vars, psfile, mars, colortable, iTime


            
        endelse
;        print,'smini,smaxi:',smini,smaxi
        if (n_elements(smini) eq 0) then smini = '0.0'
        if (n_elements(smaxi) eq 0) then smaxi = '0.0'
        mini = smini
        maxi = smaxi
                                ;Show min and max number to the interface.
        min=strtrim(string(mini),2)
        max=strtrim(string(maxi),2)
        widget_control,(*ptr).min_txt,set_value=min
        widget_control,(*ptr).max_txt,set_value=max
                                ;Reset set cursor counter to zero.
        (*ptr).cursor_cnt=0
        
    end
    'PLOT':begin
                                ;Get cursor_x & cursor_y from widget_text.
        widget_control,(*ptr).curx_txt,get_value=txt
        strx	    =	txt[0]
        cursor_x =   double(strx)
        widget_control,(*ptr).cury_txt,get_value=txt
        stry	    =   txt[0]
        cursor_y =   double(stry)
        
                                ;Set the value of step. 
        if ((*ptr).vyes_cnt eq 1) then begin
            step = (*ptr).vs
        endif
        
                                ;Clean up the display window
        plotdumb

                                ;Initialize values
        vars	  = *(*ptr).vars ;Get the variables' list
        sel	  = (*ptr).sel  ;Get the selected variable.
        nfiles = (*ptr).nfiles  ;Get the files users searched for 
        cnt1   = (*ptr).lat_lon_cnt ;Get the counter of selecting lat_lon button
        cnt2   = (*ptr).alt_lon_cnt ;Get the counter of selecting alt_lon button
        cnt3   = (*ptr).alt_lat_cnt ;Get the counter of selecting alt_lat button 
        yes	  = (*ptr).yes_cnt ;Get the counter of selecting yes (Ghost Cells) button 
        no	  = (*ptr).no_cnt ;Get the counter of selecting no (Ghost Cells) button
        yeslog = (*ptr).yes_logcnt ;Get the counter of selecting yes Plot Log button
        nolog  = (*ptr).no_logcnt ;Get the counter of selecting no Plot Log button
        nalts  = (*ptr).nalts   ;Get the number of nalts
        nlats  = (*ptr).nlats   ;Get the number of nlats
        nlons  = (*ptr).nlons   ;Get the number of nlons
        
        yeswrite_cnt = (*ptr).yeswrite_cnt
        
        polar  = (*ptr).polar
        npolar = (*ptr).npolar
        if (polar eq 1) and (npolar eq 1) then begin
            MinLat = (*ptr).nplat
        endif
        if (polar eq 1) and (npolar eq 0) then begin
            MinLat = (*ptr).splat
        endif

        mars = (*ptr).mars
        colortable = (*ptr).colortable
        
        showgridyes   = (*ptr).gridyes_cnt
        plotvectoryes = (*ptr).vyes_cnt
        vi_cnt 	 = (*ptr).vi_cnt
        vn_cnt 	 = (*ptr).vn_cnt
        cf 		 = (*ptr).cf
        vs 		 = (*ptr).vs
        cursor_cnt    = (*ptr).cursor_cnt

                                ;Counter 'cnt' is to count if any of lat_lon,alt_lon,alt_lat buttons selected.
        cnt=1
        if (cnt1 eq 0) and (cnt2 eq 0) and (cnt3 eq 0) then begin
            cnt=0
        endif
        
                                ;Check if there are any variables in the text_list. 
                                ;If there is not, then initialize the counter of selecting 
                                ;plot button (plot_cnt) to zero and pop out a warning message. 
                                ;If there are variables, then start to do the plot.  
        if (vars(0) eq '') or (sel eq -1) or (cnt eq 0) then begin
            (*ptr).plot_cnt  =  0
            me  =  widget_message('Please do search, select a variable and a setting first!')
        endif else begin
            (*ptr).plot_cnt  =  1 ;Set the plot counter to 1.
                                ;Initialize the data,alt,lat,lon,k,xrange,yrange
            data   = *(*ptr).data ;Get the array data which be plotted. 
            alt    = reform(data(2,*,*,*))/1000.0 ;Calculate alt,lat,lon and k.
            lat    = reform(data(1,*,*,*))*180.0/!pi
            lon    = reform(data(0,*,*,*))*180.0/!pi
            k	   = n_elements(alt(0,0,*))-1
            xrange = [0,0]      ;Initialize xrange and yrange.
            yrange = [0,0]

            selset = 0
            if cnt1 eq 1 then begin
                selset  = (*ptr).lat_lon_max
            endif
            if cnt2 eq 1 then begin
                selset  = (*ptr).alt_lon_max
            endif
            if cnt3 eq 1 then begin
                selset  = (*ptr).alt_lat_max
            endif

            widget_control,(*ptr).min_txt,get_value=smini
            widget_control,(*ptr).max_txt,get_value=smaxi
            
            smini = smini
            smaxi = smaxi
            
            filename = (*ptr).filename

                                ; Check if user wants to write the result to a file
                                ;If user wanted then setdevice. 
            if yeswrite_cnt eq 1 then begin
                widget_control,(*ptr).writefile_txt,get_value=txt
                file_name=strtrim(string(txt[0]),2)
                                ;Check if user has already entered a file name. 
                if (file_name eq '') or (n_elements(file_name) eq 0) then begin
                    mes=widget_message('Please enter the PS file name!')
                endif else begin
                    setdevice,file_name,'l',4,0.95
                endelse
            endif

            widget_control,(*ptr).writefile_txt,get_value=txt
            psfile=strtrim(string(txt[0]),2)
            mars = (*ptr).mars
            colortable = (*ptr).colortable
            itime = *(*ptr).itime
        
                                ;Call the subprogram to do the plots.
            thermo_plot,cursor_x,cursor_y,strx,stry,step,vars,sel,nfiles,	  $
              cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
              nolog,nalts,nlats,nlons,yeswrite_cnt,$
              polar,npolar,MinLat,showgridyes,	  $
              plotvectoryes,vi_cnt,vn_cnt,cf,	  $
              cursor_cnt,data,alt,lat,lon,	  $
              xrange,yrange,selset,smini,smaxi,	  $
              filename,vars, psfile,mars, colortable, iTime

                                ;If user write the result to a file, then closedevice right now. 
            if yeswrite_cnt eq 1 then begin
                closedevice
                mes=widget_message('Done with writing into file!')
            endif
            
        endelse
        if n_elements(smini) eq 0 then smini = 0
        if n_elements(smaxi) eq 0 then smaxi = 0
        mini = smini
        maxi = smaxi
                                ;Show min and max number to the interface.
        min=strtrim(string(mini),2)
        max=strtrim(string(maxi),2)
        widget_control,(*ptr).min_txt,set_value=min
        widget_control,(*ptr).max_txt,set_value=max
                                ;Reset set cursor counter to zero.
        (*ptr).cursor_cnt=0
    end
    'MINI':
    'MAXI':
    'GRIDYES':begin
        (*ptr).gridyes_cnt = 1
    end
    'GRIDNO':begin
        (*ptr).gridyes_cnt = 0
    end
    'MARSYES':begin
        (*ptr).mars = 1
    end
    'MARSNO':begin
        (*ptr).mars = 0
    end
    'wyr':begin
        (*ptr).colortable = 'wyr'
    end
    'mid':begin
        (*ptr).colortable = 'mid'
    end
    'all':begin
        (*ptr).colortable = 'all'
    end
    'rwb':begin
        (*ptr).colortable = 'rwb'
    end
    'bwr':begin
        (*ptr).colortable = 'bwr'
    end
    'red':begin
        (*ptr).colortable = 'red'
    end
    'wb':begin
        (*ptr).colortable = 'wb'
    end
    'bw':begin
        (*ptr).colortable = 'bw'
    end
    'VYES':begin
        widget_control,(*ptr).vector_buts,sensitive = 1
;        widget_control,(*ptr).vector_buts,/set_button
        widget_control,(*ptr).vi_but,sensitive = 1
;        widget_control,(*ptr).cf_sli,/set_button
        widget_control,(*ptr).cf_sli,sensitive=1
;        widget_control,(*ptr).vs_sli,/set_button
        widget_control,(*ptr).vs_sli,sensitive=1
        (*ptr).vyes_cnt = 1
        (*ptr).vs_cnt=1
        (*ptr).cf_cnt=1
        (*ptr).vn_cnt=1
        (*ptr).vi_cnt=0
    end
    'VNO':begin
        widget_control,(*ptr).vector_buts,sensitive = 0
        widget_control,(*ptr).vector_buts,set_button=0
;        widget_control,(*ptr).cf_sli,/set_button
        widget_control,(*ptr).cf_sli,sensitive=0
;        widget_control,(*ptr).vs_sli,/set_button
        widget_control,(*ptr).vs_sli,sensitive=0	
        widget_control,(*ptr).vs_sli,get_value=vs
        widget_control,(*ptr).cf_sli,get_value=cf
        (*ptr).vyes_cnt = 0
        (*ptr).vs_cnt=0
        (*ptr).cf_cnt=0
        (*ptr).vs_cnt=0
        (*ptr).cf_cnt=0
        (*ptr).vs=vs
        (*ptr).cf=cf
    end
    'VI':begin
        (*ptr).vi_cnt=1
        (*ptr).vn_cnt=0
    end
    'VN':begin
        (*ptr).vi_cnt=0
        (*ptr).vn_cnt=1
    end
    'CF':begin
        widget_control,(*ptr).cf_sli,get_value=cf
        (*ptr).cf=cf

                                ;Get cursor_x & cursor_y from widget_text.
        widget_control,(*ptr).curx_txt,get_value=txt
        strx	    =	txt[0]
        cursor_x =   double(strx)
        widget_control,(*ptr).cury_txt,get_value=txt
        stry	    =   txt[0]
        cursor_y =   double(stry)
        
                                ;Set the value of step. 
        if ((*ptr).vyes_cnt eq 1) then begin
            step = (*ptr).vs
        endif
        
                                ;Clean up the display window
        plotdumb

                                ;Initialize values
        vars	  = *(*ptr).vars ;Get the variables' list
        sel	  = (*ptr).sel  ;Get the selected variable.
        nfiles = (*ptr).nfiles  ;Get the files users searched for 
        cnt1   = (*ptr).lat_lon_cnt ;Get the counter of selecting lat_lon button
        cnt2   = (*ptr).alt_lon_cnt ;Get the counter of selecting alt_lon button
        cnt3   = (*ptr).alt_lat_cnt ;Get the counter of selecting alt_lat button 
        yes	  = (*ptr).yes_cnt ;Get the counter of selecting yes (Ghost Cells) button 
        no	  = (*ptr).no_cnt ;Get the counter of selecting no (Ghost Cells) button
        yeslog = (*ptr).yes_logcnt ;Get the counter of selecting yes Plot Log button
        nolog  = (*ptr).no_logcnt ;Get the counter of selecting no Plot Log button
        nalts  = (*ptr).nalts   ;Get the number of nalts
        nlats  = (*ptr).nlats   ;Get the number of nlats
        nlons  = (*ptr).nlons   ;Get the number of nlons
        
        yeswrite_cnt = (*ptr).yeswrite_cnt
        
        polar  = (*ptr).polar
        npolar = (*ptr).npolar
        if (polar eq 1) and (npolar eq 1) then begin
            MinLat = (*ptr).nplat
        endif
        if (polar eq 1) and (npolar eq 0) then begin
            MinLat = (*ptr).splat
        endif
        
        showgridyes   = (*ptr).gridyes_cnt
        plotvectoryes = (*ptr).vyes_cnt
        vi_cnt 	 = (*ptr).vi_cnt
        vn_cnt 	 = (*ptr).vn_cnt
        cf 		 = (*ptr).cf
        vs 		 = (*ptr).vs
        cursor_cnt    = (*ptr).cursor_cnt

                                ;Counter 'cnt' is to count if any of lat_lon,alt_lon,alt_lat buttons selected.
        cnt=1
        if (cnt1 eq 0) and (cnt2 eq 0) and (cnt3 eq 0) then begin
            cnt=0
        endif
        
                                ;Check if there are any variables in the text_list. 
                                ;If there is not, then initialize the counter of selecting 
                                ;plot button (plot_cnt) to zero and pop out a warning message. 
                                ;If there are variables, then start to do the plot.  
        if (vars(0) eq '') or (sel eq -1) or (cnt eq 0) then begin
            (*ptr).plot_cnt  =  0
            me  =  widget_message('Please do search, select a variable and a setting first!')
        endif else begin
            (*ptr).plot_cnt  =  1 ;Set the plot counter to 1.
                                ;Initialize the data,alt,lat,lon,k,xrange,yrange
            data   = *(*ptr).data ;Get the array data which be plotted. 
            alt    = reform(data(2,*,*,*))/1000.0 ;Calculate alt,lat,lon and k.
            lat    = reform(data(1,*,*,*))*180.0/!pi
            lon    = reform(data(0,*,*,*))*180.0/!pi
            k	   = n_elements(alt(0,0,*))-1
            xrange = [0,0]      ;Initialize xrange and yrange.
            yrange = [0,0]

            selset = 0
            if cnt1 eq 1 then begin
                selset  = (*ptr).lat_lon_max
            endif
            if cnt2 eq 1 then begin
                selset  = (*ptr).alt_lon_max
            endif
            if cnt3 eq 1 then begin
                selset  = (*ptr).alt_lat_max
            endif

            widget_control,(*ptr).min_txt,get_value=smini
            widget_control,(*ptr).max_txt,get_value=smaxi
            
            smini = smini
            smaxi = smaxi
            
            filename = (*ptr).filename

            widget_control,(*ptr).writefile_txt,get_value=txt
            psfile=strtrim(string(txt[0]),2)
            mars = (*ptr).mars
            colortable = (*ptr).colortable
            itime = *(*ptr).itime


                                ;Call the subprogram to do the plots.
            thermo_plot,cursor_x,cursor_y,strx,stry,step,vars,sel,nfiles,	  $
              cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
              nolog,nalts,nlats,nlons,yeswrite_cnt,$
              polar,npolar,MinLat,showgridyes,	  $
              plotvectoryes,vi_cnt,vn_cnt,cf,	  $
              cursor_cnt,data,alt,lat,lon,	  $
              xrange,yrange,selset,smini,smaxi,	  $
              filename,vars, psfile, mars, colortable, iTime


            
        endelse
;        print,'smini,smaxi:',smini,smaxi
        mini = smini
        maxi = smaxi
                                ;Show min and max number to the interface.
        min=strtrim(string(mini),2)
        max=strtrim(string(maxi),2)
        widget_control,(*ptr).min_txt,set_value=min
        widget_control,(*ptr).max_txt,set_value=max
                                ;Reset set cursor counter to zero.
        (*ptr).cursor_cnt=0
        
    end
    'VS':begin
        widget_control,(*ptr).vs_sli,get_value=vs
        (*ptr).vs=vs

                                ;Get cursor_x & cursor_y from widget_text.
        widget_control,(*ptr).curx_txt,get_value=txt
        strx	    =	txt[0]
        cursor_x =   double(strx)
        widget_control,(*ptr).cury_txt,get_value=txt
        stry	    =   txt[0]
        cursor_y =   double(stry)
        
                                ;Set the value of step. 
        if ((*ptr).vyes_cnt eq 1) then begin
            step = (*ptr).vs
        endif
        
                                ;Clean up the display window
        plotdumb

                                ;Initialize values
        vars	  = *(*ptr).vars ;Get the variables' list
        sel	  = (*ptr).sel  ;Get the selected variable.
        nfiles = (*ptr).nfiles  ;Get the files users searched for 
        cnt1   = (*ptr).lat_lon_cnt ;Get the counter of selecting lat_lon button
        cnt2   = (*ptr).alt_lon_cnt ;Get the counter of selecting alt_lon button
        cnt3   = (*ptr).alt_lat_cnt ;Get the counter of selecting alt_lat button 
        yes	  = (*ptr).yes_cnt ;Get the counter of selecting yes (Ghost Cells) button 
        no	  = (*ptr).no_cnt ;Get the counter of selecting no (Ghost Cells) button
        yeslog = (*ptr).yes_logcnt ;Get the counter of selecting yes Plot Log button
        nolog  = (*ptr).no_logcnt ;Get the counter of selecting no Plot Log button
        nalts  = (*ptr).nalts   ;Get the number of nalts
        nlats  = (*ptr).nlats   ;Get the number of nlats
        nlons  = (*ptr).nlons   ;Get the number of nlons
        
        yeswrite_cnt = (*ptr).yeswrite_cnt
        
        polar  = (*ptr).polar
        npolar = (*ptr).npolar
        if (polar eq 1) and (npolar eq 1) then begin
            MinLat = (*ptr).nplat
        endif
        if (polar eq 1) and (npolar eq 0) then begin
            MinLat = (*ptr).splat
        endif
        
        showgridyes   = (*ptr).gridyes_cnt
        plotvectoryes = (*ptr).vyes_cnt
        vi_cnt 	 = (*ptr).vi_cnt
        vn_cnt 	 = (*ptr).vn_cnt
        cf 		 = (*ptr).cf
        vs 		 = (*ptr).vs
        cursor_cnt    = (*ptr).cursor_cnt

                                ;Counter 'cnt' is to count if any of lat_lon,alt_lon,alt_lat buttons selected.
        cnt=1
        if (cnt1 eq 0) and (cnt2 eq 0) and (cnt3 eq 0) then begin
            cnt=0
        endif
        
                                ;Check if there are any variables in the text_list. 
                                ;If there is not, then initialize the counter of selecting 
                                ;plot button (plot_cnt) to zero and pop out a warning message. 
                                ;If there are variables, then start to do the plot.  
        if (vars(0) eq '') or (sel eq -1) or (cnt eq 0) then begin
            (*ptr).plot_cnt  =  0
            me  =  widget_message('Please do search, select a variable and a setting first!')
        endif else begin
            (*ptr).plot_cnt  =  1 ;Set the plot counter to 1.
                                ;Initialize the data,alt,lat,lon,k,xrange,yrange
            data   = *(*ptr).data ;Get the array data which be plotted. 
            alt    = reform(data(2,*,*,*))/1000.0 ;Calculate alt,lat,lon and k.
            lat    = reform(data(1,*,*,*))*180.0/!pi
            lon    = reform(data(0,*,*,*))*180.0/!pi
            k	   = n_elements(alt(0,0,*))-1
            xrange = [0,0]      ;Initialize xrange and yrange.
            yrange = [0,0]

            selset = 0
            if cnt1 eq 1 then begin
                selset  = (*ptr).lat_lon_max
            endif
            if cnt2 eq 1 then begin
                selset  = (*ptr).alt_lon_max
            endif
            if cnt3 eq 1 then begin
                selset  = (*ptr).alt_lat_max
            endif

            widget_control,(*ptr).min_txt,get_value=smini
            widget_control,(*ptr).max_txt,get_value=smaxi
            
            smini = smini
            smaxi = smaxi
            
            filename = (*ptr).filename

            widget_control,(*ptr).writefile_txt,get_value=txt
            psfile=strtrim(string(txt[0]),2)
            mars = (*ptr).mars
            colortable = (*ptr).colortable
            itime = *(*ptr).itime
        
                                ;Call the subprogram to do the plots.
            thermo_plot,cursor_x,cursor_y,strx,stry,step,vars,sel,nfiles,	  $
              cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
              nolog,nalts,nlats,nlons,yeswrite_cnt,$
              polar,npolar,MinLat,showgridyes,	  $
              plotvectoryes,vi_cnt,vn_cnt,cf,	  $
              cursor_cnt,data,alt,lat,lon,	  $
              xrange,yrange,selset,smini,smaxi,	  $
              filename,vars, psfile,mars, colortable, iTime


            
        endelse
;        print,'smini,smaxi:',smini,smaxi
        mini = smini
        maxi = smaxi
                                ;Show min and max number to the interface.
        min=strtrim(string(mini),2)
        max=strtrim(string(maxi),2)
        widget_control,(*ptr).min_txt,set_value=min
        widget_control,(*ptr).max_txt,set_value=max
                                ;Reset set cursor counter to zero.
        (*ptr).cursor_cnt=0
        
    end
    'YES':begin
        (*ptr).cursor_cnt=0
        (*ptr).yes_cnt=1
        (*ptr).no_cnt=0
        widget_control,(*ptr).no_but,set_button=0
        widget_control,(*ptr).yes_but,sensitive=1
        widget_control,(*ptr).no_but,sensitive=1
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt
    end
    'NO':begin
        (*ptr).cursor_cnt=0
        (*ptr).no_cnt=1
        (*ptr).yes_cnt=0
        widget_control,(*ptr).yes_but,set_button=0
        widget_control,(*ptr).yes_but,sensitive=1
        widget_control,(*ptr).no_but,sensitive=1
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt
    end

    'POLAR' : begin
;        print, "Polar"
        (*ptr).polar = 1
        widget_control,(*ptr).subpolar_buts,sensitive = 1
;        widget_control,(*ptr).subpolar_buts,/set_button
;        widget_control,(*ptr).npolar_but,/set_button
        widget_control,(*ptr).npolar_but,sensitive=1
        widget_control,(*ptr).polar_sli,sensitive=1
    end
    'POLARSLI': begin
        widget_control,(*ptr).polar_sli,get_value = value
        polar = (*ptr).polar
        npolar= (*ptr).npolar
        if (polar eq 1) and (npolar eq 1) then begin
            (*ptr).nplat = value	
;            print,'nplat is:',(*ptr).nplat
        endif
        if (polar eq 1) and (npolar eq 0) then begin
            (*ptr).splat = value
;            print,'splat is:',(*ptr).splat
        endif

                                ;Get cursor_x & cursor_y from widget_text.
        widget_control,(*ptr).curx_txt,get_value=txt
        strx	    =	txt[0]
        cursor_x =   double(strx)
        widget_control,(*ptr).cury_txt,get_value=txt
        stry	    =   txt[0]
        cursor_y =   double(stry)
        
                                ;Set the value of step. 
        if ((*ptr).vyes_cnt eq 1) then begin
            step = (*ptr).vs
        endif
        
                                ;Clean up the display window
        plotdumb

                                ;Initialize values
        vars	  = *(*ptr).vars ;Get the variables' list
        sel	  = (*ptr).sel  ;Get the selected variable.
        nfiles = (*ptr).nfiles  ;Get the files users searched for 
        cnt1   = (*ptr).lat_lon_cnt ;Get the counter of selecting lat_lon button
        cnt2   = (*ptr).alt_lon_cnt ;Get the counter of selecting alt_lon button
        cnt3   = (*ptr).alt_lat_cnt ;Get the counter of selecting alt_lat button 
        yes	  = (*ptr).yes_cnt ;Get the counter of selecting yes (Ghost Cells) button 
        no	  = (*ptr).no_cnt ;Get the counter of selecting no (Ghost Cells) button
        yeslog = (*ptr).yes_logcnt ;Get the counter of selecting yes Plot Log button
        nolog  = (*ptr).no_logcnt ;Get the counter of selecting no Plot Log button
        nalts  = (*ptr).nalts   ;Get the number of nalts
        nlats  = (*ptr).nlats   ;Get the number of nlats
        nlons  = (*ptr).nlons   ;Get the number of nlons
        
        yeswrite_cnt = (*ptr).yeswrite_cnt
        
        polar  = (*ptr).polar
        npolar = (*ptr).npolar
        if (polar eq 1) and (npolar eq 1) then begin
            MinLat = (*ptr).nplat
        endif
        if (polar eq 1) and (npolar eq 0) then begin
            MinLat = (*ptr).splat
        endif
        
        showgridyes   = (*ptr).gridyes_cnt
        plotvectoryes = (*ptr).vyes_cnt
        vi_cnt 	 = (*ptr).vi_cnt
        vn_cnt 	 = (*ptr).vn_cnt
        cf 		 = (*ptr).cf
        vs 		 = (*ptr).vs
        cursor_cnt    = (*ptr).cursor_cnt

                                ;Counter 'cnt' is to count if any of lat_lon,alt_lon,alt_lat buttons selected.
        cnt=1
        if (cnt1 eq 0) and (cnt2 eq 0) and (cnt3 eq 0) then begin
            cnt=0
        endif
        
                                ;Check if there are any variables in the text_list. 
                                ;If there is not, then initialize the counter of selecting 
                                ;plot button (plot_cnt) to zero and pop out a warning message. 
                                ;If there are variables, then start to do the plot.  
        if (vars(0) eq '') or (sel eq -1) or (cnt eq 0) then begin
            (*ptr).plot_cnt  =  0
            me  =  widget_message('Please do search, select a variable and a setting first!')
        endif else begin
            (*ptr).plot_cnt  =  1 ;Set the plot counter to 1.
                                ;Initialize the data,alt,lat,lon,k,xrange,yrange
            data   = *(*ptr).data ;Get the array data which be plotted. 
            alt    = reform(data(2,*,*,*))/1000.0 ;Calculate alt,lat,lon and k.
            lat    = reform(data(1,*,*,*))*180.0/!pi
            lon    = reform(data(0,*,*,*))*180.0/!pi
            k	   = n_elements(alt(0,0,*))-1
            xrange = [0,0]      ;Initialize xrange and yrange.
            yrange = [0,0]

            selset = 0
            if cnt1 eq 1 then begin
                selset  = (*ptr).lat_lon_max
            endif
            if cnt2 eq 1 then begin
                selset  = (*ptr).alt_lon_max
            endif
            if cnt3 eq 1 then begin
                selset  = (*ptr).alt_lat_max
            endif

            widget_control,(*ptr).min_txt,get_value=smini
            widget_control,(*ptr).max_txt,get_value=smaxi
            
            smini = smini
            smaxi = smaxi
            
            filename = (*ptr).filename

            widget_control,(*ptr).writefile_txt,get_value=txt
            psfile=strtrim(string(txt[0]),2)
            mars = (*ptr).mars
            colortable = (*ptr).colortable
            itime = *(*ptr).itime
        

                                ;Call the subprogram to do the plots.
            thermo_plot,cursor_x,cursor_y,strx,stry,step,vars,sel,nfiles,	  $
              cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
              nolog,nalts,nlats,nlons,yeswrite_cnt,$
              polar,npolar,MinLat,showgridyes,	  $
              plotvectoryes,vi_cnt,vn_cnt,cf,	  $
              cursor_cnt,data,alt,lat,lon,	  $
              xrange,yrange,selset,smini,smaxi,	  $
              filename,vars, psfile,mars, colortable, iTime


            
        endelse
;        print,'smini,smaxi:',smini,smaxi
        mini = smini
        maxi = smaxi
                                ;Show min and max number to the interface.
        min=strtrim(string(mini),2)
        max=strtrim(string(maxi),2)
        widget_control,(*ptr).min_txt,set_value=min
        widget_control,(*ptr).max_txt,set_value=max
                                ;Reset set cursor counter to zero.
        (*ptr).cursor_cnt=0
        
    end
    'CART'  : begin
;        print, "Cartesean"
        (*ptr).polar = 0
        widget_control,(*ptr).subpolar_buts,sensitive = 0
        widget_control,(*ptr).polar_sli,sensitive=0
    end
    'NP':begin
        (*ptr).npolar = 1
        widget_control,(*ptr).polar_sli,set_slider_min = 0
        widget_control,(*ptr).polar_sli,set_slider_max = 70
        widget_control,(*ptr).polar_sli,get_value=nplat
        (*ptr).nplat=nplat
    end
    'SP':begin
        (*ptr).npolar = 0
        widget_control,(*ptr).polar_sli,set_slider_min = -70
        widget_control,(*ptr).polar_sli,set_slider_max = 0
        widget_control,(*ptr).polar_sli,get_value=splat
        (*ptr).splat=splat
    end

    'LOGYES':begin
;        print,'LOG YES.'
        (*ptr).cursor_cnt=0
        (*ptr).yes_logcnt=1
        (*ptr).no_logcnt=0
        widget_control,(*ptr).no_logbut,set_button=0
        widget_control,(*ptr).no_logbut,sensitive=1
        widget_control,(*ptr).yes_logbut,sensitive=1
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt
    end
    'LOGNO':begin
        (*ptr).cursor_cnt=0
        (*ptr).yes_logcnt=0
        (*ptr).no_logcnt=1
        widget_control,(*ptr).yes_logbut,set_button=0
        widget_control,(*ptr).yes_logbut,sensitive=1
        widget_control,(*ptr).no_logbut,sensitive=1
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        widget_control,(*ptr).min_txt,set_value=txt
        widget_control,(*ptr).max_txt,set_value=txt
    end
    'SETCURSOR':begin
        (*ptr).cursor_cnt=1
        plot_cnt=(*ptr).plot_cnt
                                ;Reset text widgets
        txt=''
        widget_control,(*ptr).curx_txt,set_value=txt
        widget_control,(*ptr).cury_txt,set_value=txt
        if plot_cnt eq 0 then begin
            mes=widget_message("Please do plot first!")
        endif else begin
            cursor,x,y,/wait,/DATA
            (*ptr).cursor_x=x
            (*ptr).cursor_y=y
;            print,'x,y:',x,y
            x=strtrim(string(x),2)
            y=strtrim(string(y),2)
            widget_control,(*ptr).curx_txt,set_value=x
            widget_control,(*ptr).cury_txt,set_value=y
            cnt1=(*ptr).lat_lon_cnt
            cnt2=(*ptr).alt_lon_cnt
            cnt3=(*ptr).alt_lat_cnt
            datatoplot=*(*ptr).data
        endelse
    end
    'CURX':
    'CURY':
    'YESWRITE':begin
        (*ptr).cursor_cnt=0
        (*ptr).yeswrite_cnt=1
        (*ptr).nowrite_cnt=0
        widget_control,(*ptr).nowrite_but,set_button=0
        widget_control,(*ptr).yeswrite_but,sensitive=1
        widget_control,(*ptr).nowrite_but,sensitive=1
        widget_control,(*ptr).writefile_txt,sensitive=1
        
    end
    'NOWRITE':begin
        (*ptr).cursor_cnt=0
        (*ptr).yeswrite_cnt=0
        (*ptr).nowrite_cnt=1
        widget_control,(*ptr).yeswrite_but,set_button=0
        widget_control,(*ptr).yeswrite_but,sensitive=1
        widget_control,(*ptr).nowrite_but,sensitive=1
        widget_control,(*ptr).writefile_txt,sensitive=0
    end
    'WRITEFILE':
    'CAN':begin
        widget_control,event.top,/destroy
    end
endcase
if (whichevent ne 'CAN') and (whichevent ne 'OK') then begin
    widget_control,event.top,set_uvalue=ptr,/no_copy
endif
return
NOTDONE:
widget_control,(*ptr).list,scr_xsize=120,set_value=files
return
end

pro thermo_gui

title='Plot Thermoshpere'
mainbase=widget_base(title=title,/column)
base = widget_base(mainbase,/column)
txt=''

block_base=widget_base(base,/column,/frame)
blockn_base=widget_base(block_base,/row)
lab = widget_label(blockn_base,value='Enter block number:',/align_left)
rb_txt=widget_text(blockn_base,value=txt,uvalue='ROWBLOCK',xsize=6,/editable,/align_left)
lab = widget_label(blockn_base,value='X',/align_left)
cb_txt=widget_text(blockn_base,value=txt,uvalue='COLBLOCK',xsize=6,/editable,/align_left)

dirlab=widget_base(base,/column,/frame)
lab=widget_label(dirlab,value='Enter file name (can use *) [p0000b0001i0200.dat]:',/align_left)
filename_base=widget_base(dirlab,/row)
filelist = findfile("-t *.save")
;if (strlen(filelist(0)) eq 0) then filelist = findfile("-t *.dat")
if (strlen(filelist(0)) eq 0) then filelist = findfile("-t *.bin")
txt=filelist(0)

file_txt=widget_text(filename_base,value=txt,uvalue='FILENAME',xsize=46,/editable,/align_left)
sear_but=widget_button(filename_base,value='Search',uvalue='SEARCH',xsize=60)
multifiles_but = widget_button(filename_base,value='Multiple',uvalue='MULTIFILES',xsize=60)

lists_base=widget_base(base,/row,/frame)
vars_base=widget_base(lists_base,/column)
sel=-1
lab=widget_label(vars_base,value='Variables',/align_center)
vars=''

ys = 100
if (float(!version.release) gt 5.55) then ys = 20
vars_list=widget_list(vars_base,value=vars,uvalue='SELVARS',xsize=20,$
	ysize=ys,/align_right)
nalts=1
nlats=1
nlons=1
butt_base=widget_base(lists_base,/column)
lab=widget_label(butt_base,value='Select Settings',/align_center)
but_base=widget_base(butt_base,/column,/align_center,/frame)
buts_base=widget_base(but_base,/row,/exclusive,/align_center)
slis_base=widget_base(but_base,/row,/align_center)
lat_lon_but=widget_button(buts_base,value='Lat & Lon',uvalue='LAT_LON',/no_release)
alt_lon_but=widget_button(buts_base,value='Alt & Lon',uvalue='ALT_LON',/no_release)
alt_lat_but=widget_button(buts_base,value='ALT & LAT',uvalue='ALT_LAT',/no_release)
sli=widget_slider(slis_base,/drag,maximum=nlons,minimum=0,uvalue='SLI',xsize=260)

txt_base=widget_base(but_base,/row,/align_center)
lab=widget_label(txt_base,value='Location')
txt=''
slider_txt=widget_text(txt_base,value=txt,uvalue='LOC',xsize=12)

grid_base     = widget_base(but_base,/row,/align_left)
lab           = widget_label(grid_base,value='Show Grid ?       ')

grid_buts     = widget_base(grid_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
grid_yes_but  = widget_button(grid_buts,value='Yes',uvalue='GRIDYES',/no_release)
grid_no_but   = widget_button(grid_buts,value='No',uvalue='GRIDNO',/no_release)

mars_base     = widget_base(but_base,/row,/align_left)
lab           = widget_label(mars_base,value='Use Mars Map ?    ')
mars_buts     = widget_base(mars_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
mars_yes_but  = widget_button(mars_buts,value='Yes',uvalue='MARSYES',/no_release)
mars_no_but   = widget_button(mars_buts,value='No',uvalue='MARSNO',/no_release)
mars = 0

ghost_base = widget_base(but_base,/row,/align_left)
lab        = widget_label(ghost_base,value='Show Ghost Cells? ')
gh_buts    = widget_base(ghost_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
yes_but    = widget_button(gh_buts,value='Yes',uvalue='YES',/no_release)
no_but     = widget_button(gh_buts,value='No',uvalue='NO',/no_release)

polar_base     = widget_base(but_base,/row,/align_left)
lab            = widget_label(polar_base,value='Polar Plot ?')
polar_buts     = widget_base(polar_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
yes_polarbut   = widget_button(polar_buts,value='Yes', uvalue='POLAR', /no_release)
no_polarbut    = widget_button(polar_buts,value='No',  uvalue='CART',  /no_release)
subpolar_buts  = widget_base(polar_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
npolar_but     = widget_button(subpolar_buts,value='N.P',uvalue='NP',/no_release)
spolar_but     = widget_button(subpolar_buts,value='S.P',uvalue='SP',/no_release)
polarsli_base  = widget_base(but_base,/row,/align_center)
lab	       = widget_label(polarsli_base,value = 'MinLat:')
polar_sli      = widget_slider(polarsli_base,/drag,maximum=70,minimum=0,uvalue='POLARSLI',xsize=260)

;widget_control,polar_buts,/set_button
polar = 0
widget_control,no_polarbut, sensitive = 1
;widget_control,subpolar_buts,/set_button
widget_control,subpolar_buts,sensitive = 0
widget_control,polar_sli,sensitive = 0

log_base   = widget_base(but_base,/row,/align_left)
lab        = widget_label(log_base,value='Plot Log (Log10) ?')
log_buts   = widget_base(log_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
yes_logbut = widget_button(log_buts,value='Yes',uvalue='LOGYES',/no_release)
no_logbut  = widget_button(log_buts,value='No',uvalue='LOGNO',/no_release)

vector_base=widget_base(but_base,/row,/align_left)
lab=widget_label(vector_base,value='Plot Vectors ?    ')
vector_yesno_buts=widget_base(vector_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
vyes_but=widget_button(vector_yesno_buts,value='Yes',uvalue='VYES',/no_release)
vno_but=widget_button(vector_yesno_buts,value='No',uvalue='VNO',/no_release)
vector_buts=widget_base(vector_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
vi_but=widget_button(vector_buts,value='Vi',uvalue='VI',/no_release)
vn_but=widget_button(vector_buts,value='Vn',uvalue='VN',/no_release)

vector_sli_base=widget_base(but_base,/column,/align_left,GRID_LAYOUT=1)
vector_cf_base = widget_base(vector_sli_base,/row)
cf_lab=widget_label(vector_cf_base,value ='Calculate Factor: ')
cf_sli=widget_slider(vector_cf_base,/drag,maximum=10,minimum=1,uvalue='CF',xsize=200)
vector_vs_base = widget_base(vector_sli_base,/row)
vs_lab = widget_label(vector_vs_base,value='Vectors Step:     ')
vs_sli=widget_slider(vector_vs_base,/drag,maximum=5,minimum=1,uvalue='VS',xsize=200)

vyes_cnt=0
vno_cnt=1
;widget_control,vector_yesno_buts,/set_button
widget_control,vno_but,sensitive=1
widget_control,vector_buts,sensitive=0

cursor_base=widget_base(but_base,/row,/align_left)
cur_but=widget_button(cursor_base,value='Set Cursor',uvalue='SETCURSOR',/no_release)
lab=widget_label(cursor_base,value='X:')
txt=''
curx_txt=widget_text(cursor_base,value=txt,uvalue='CURX',xsize=12)
lab=widget_label(cursor_base,value='Y:')
txt=''
cury_txt=widget_text(cursor_base,value=txt,uvalue='CURY',xsize=12)

third_base=widget_base(lists_base,/column)

txt_base=widget_base(third_base,/row,/align_left)
lab=widget_label(txt_base,value='Mini')
txt=''
min_txt=widget_text(txt_base,value=txt,uvalue='MINI',xsize=12,/editable)
lab=widget_label(txt_base,value='Maxi')
max_txt=widget_text(txt_base,value=txt,uvalue='MAXI',xsize=12,/editable)
widget_control,sli,sensitive=0

base=widget_base(third_base,/column)
buts=widget_base(base,/row)
pc_base=widget_base(base,/row,/align_right)
lab=widget_label(buts,value='Write to a file?')
yn_buts=widget_base(buts,/row,/exclusive)
nowrite_but=widget_button(yn_buts,value='No',uvalue='NOWRITE',/no_release)
nowrite_cnt=1
;widget_control,buts,/set_button
widget_control,nowrite_but,sensitive=0
yeswrite_but=widget_button(yn_buts,value='Yes',uvalue='YESWRITE',/no_release)
yeswrite_cnt=0
txt=''
writefile_txt=widget_text(buts,value=txt,uvalue='WRITEFILE',xsize=36,/editable)
okID=widget_button(pc_base,value='PLOT',uvalue='PLOT',xsize=100)
cancelID=widget_button(pc_base,value='CANCEL',uvalue='CAN',xsize=100)

;main_color_base=widget_base(third_base,/row,/align_left)

color_base     = widget_base(third_base,/row,/align_left)
lab            = widget_label(color_base, value='Color Table : ')
color_buts     = widget_base(color_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
color_mid_but  = widget_button(color_buts,value='mid',uvalue='mid',/no_release)
color_all_but  = widget_button(color_buts,value='all',uvalue='all',/no_release)
color_wyr_but  = widget_button(color_buts,value='wyr',uvalue='wyr',/no_release)
color_rwb_but  = widget_button(color_buts,value='rwb',uvalue='rwb',/no_release)
color_bwr_but  = widget_button(color_buts,value='bwr',uvalue='bwr',/no_release)
;color_red_but  = widget_button(color_buts,value='red',uvalue='red',/no_release)
;color_wb_but   = widget_button(color_buts,value='wb',uvalue='wb',/no_release)
;color_bw_but   = widget_button(color_buts,value='bw',uvalue='bw',/no_release)

color2_base     = widget_base(third_base,/row,/align_left)
lab            = widget_label(color2_base,value='              ')
color2_buts     = widget_base(color2_base,/row,/exclusive,/align_left,GRID_LAYOUT=1)
;color2_mid_but  = widget_button(color2_buts,value='mid',uvalue='mid',/no_release)
;color2_wyr_but  = widget_button(color2_buts,value='wyr',uvalue='wyr',/no_release)
;color2_rwb_but  = widget_button(color2_buts,value='rwb',uvalue='rwb',/no_release)
;color2_bwr_but  = widget_button(color2_buts,value='bwr',uvalue='bwr',/no_release)
color2_red_but  = widget_button(color2_buts,value='red',uvalue='red',/no_release)
color2_wb_but   = widget_button(color2_buts,value='wb',uvalue='wb',/no_release)
color2_bw_but   = widget_button(color2_buts,value='bw',uvalue='bw',/no_release)
colortable = 'mid'


no_cnt=1
;widget_control,gh_buts,/set_button
widget_control,no_but,sensitive=0

no_logcnt=1
;widget_control,log_buts,/set_button
widget_control,no_logbut,sensitive=0
plot_cnt=0
cursor_cnt=0

gridyes_cnt=0
;widget_control,grid_buts,/set_button
widget_control,grid_no_but,sensitive=1



cf_cnt=0
vs_cnt=0
;widget_control,cf_sli,/set_button
widget_control,cf_sli,sensitive=0
;widget_control,vs_sli,/set_button
widget_control,vs_sli,sensitive=0


if nowrite_cnt eq 1 then begin
	widget_control,writefile_txt,sensitive=0
    endif
itime = [1965,0,0,0,0,0]

ptr = ptr_new({status:0, $
               lists_base:lists_base,$
               file_txt:file_txt,$
               rb_txt:rb_txt,$
               cb_txt:cb_txt,$
               rb:'',cb:'',$
               vars_list:vars_list,$
               filename:'',$
               vars:ptr_new(vars),$
               data:ptr_new(data),$
               nvars:0,$
               nalts:0,nlats:0,nlons:0,$
               lat_lon_cnt:0,alt_lon_cnt:0,alt_lat_cnt:0,$
	sli:sli,lat_lon_but:lat_lon_but,alt_lon_but:alt_lon_but,lat_lon_max:0,$
	alt_lon_max:0,alt_lat_max:0,alt_lat_but:alt_lat_but,sel:-1,nfiles:0,$
	slider_txt   : slider_txt,   $
        min_txt      : min_txt,      $
        max_txt      : max_txt,      $
        itime        : ptr_new(itime), $
        gh_buts      : gh_buts,      $
        no_but       : no_but,       $
	yes_but      : yes_but,      $
        no_cnt       : 1,            $
        yes_cnt      : 0,            $
        yes_logbut   : yes_logbut,   $
        bl_cnt       : 0,            $
	no_logbut    : no_logbut,    $
	yes_polarbut : yes_polarbut, $
	no_polarbut  : no_polarbut,  $
        polar        : 0,            $
	npolar	     : 1, 	     $
        mars         : 0,            $
        colortable   : 'mid',        $
	subpolar_buts: subpolar_buts,$
	npolar_but   : npolar_but,   $
	polarsli_base: polarsli_base,$
	polar_sli    : polar_sli,    $
	splat	     : 0,	     $
	nplat        : 0,	     $
        no_logcnt:1,yes_logcnt:0,log_buts:log_buts,$
	plot_cnt:0,cursor_x:0,cursor_y:0,curx_txt:curx_txt,cury_txt:cury_txt,$
	cursor_cnt:0,nowrite_cnt:1,yeswrite_cnt:0,writefile_txt:writefile_txt,$
	yeswrite_but:yeswrite_but,nowrite_but:nowrite_but,yn_buts:yn_buts,$
	grid_yes_but:grid_yes_but,grid_no_but:grid_no_but,vyes_but:vyes_but,vector_buts:vector_buts,$
	vno_but:vno_but,vi_but:vi_but,vn_but:vn_but,gridyes_cnt:0,vyes_cnt:0,$
	vi_cnt:0,vn_cnt:1,cf_sli:cf_sli,vs_sli:vs_sli,cf_cnt:0,vs_cnt:0,vs:1,cf:1})

widget_control,mainbase,set_uvalue=ptr
widget_control,mainbase,/realize
Xmanager,'thermo_gui',mainbase
ptr_free,ptr
end







