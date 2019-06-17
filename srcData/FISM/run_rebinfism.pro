year = 2013
d_in_m = [31,28,31,30,31,30,31,31,30,31,30,31]
if (year mod 4 eq 0) then d_in_m(1)=29
for imonth = 1, 12 do begin
    ndays = d_in_m(imonth-1)

    for iday = 1, ndays do begin
        date = tostr(year)+chopr('0'+tostr(imonth),2)+chopr('0'+tostr(iday),2)

        rebin_fism,date

    endfor

endfor

end
