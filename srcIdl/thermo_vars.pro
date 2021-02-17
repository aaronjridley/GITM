
filelist = findfile('3D???*.bin')
file = filelist(0)

gitm_read_header, file, GitmTimes, nVars, Vars, $
                  nLons, nLats, nAlts, version

display, vars
c_r_to_a, itime, GitmTimes
c_a_to_s, itime, stime
print, 'Start Time : ',stime

file = filelist(n_elements(filelist)-1)

gitm_read_header, file, GitmTimes, nVars, Vars, $
                  nLons, nLats, nAlts, version

c_r_to_a, itime, GitmTimes
c_a_to_s, itime, stime
print, 'End Time : ',stime

end

