year = 2009
for imonth = 1, 12 do begin
    ndays = d_in_m(year,imonth)

    for iday = 1, ndays do begin
        date = tostr(year)+chopr('0'+tostr(imonth),2)+chopr('0'+tostr(iday),2)

        rebin_fism,date

    endfor

endfor

end
