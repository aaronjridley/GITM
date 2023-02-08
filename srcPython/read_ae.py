#!/usr/bin/env python

import numpy as np
import re
from datetime import datetime

#-----------------------------------------------------------------------------
# read in an SME file from JHU/APL 
#-----------------------------------------------------------------------------

def read_ae(filein):

    fpin = open(filein, 'r')

    # Read header:
    for line in fpin:

        m = re.match(r'<year>',line)
        if m:
            break

    times = []
    ae = []
    al = []
    au = []
        
    # Read main data
    for line in fpin:
        cols = line.split()

        times.append(datetime(int(cols[0]),
                              int(cols[1]),
                              int(cols[2]),
                              int(cols[3]),
                              int(cols[4]),
                              int(cols[5]),
                              0))
        ae.append(float(cols[6]))
        al.append(float(cols[7]))
        au.append(float(cols[8]))

    data = {'time' : times,
            'ae' : ae,
            'al' : al,
            'au' : au}
        
    return data
        
