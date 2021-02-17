
start = ask('starting characters of the satellite',start)

filelist_new = findfile(start+"*.bin")

gitm_read_bin, filelist_new, data, time, nVars, Vars, version

save, file = start+'.save', filelist_new, data, time, nVars, Vars, version

end

