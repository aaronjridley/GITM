#!/usr/bin/env python
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

'''
Open a GITM 3D file adn create a plot similar to the example given by Aaron.
Note that as pybats.gitm is more developed, a plot like this should be made
using syntax like,
>>>a=gitm.GitmBin('filename')
>>>a.add_alt_slice(0, 'Rho', add_cbar=True)
That's how most pybats stuff works right now.
'''

# Import shit.  I needed a lot of shit this time.  
import numpy as np
from spacepy.pybats import gitm
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

# Open file.
a=gitm.GitmBin('./3DALL_t061213_000000.bin')

# Make contour of rho at lowest altitude (index 0).
# Convert lat lon from rad to degrees.
p=180.0/np.pi
f=plt.figure() #make a fig.
ax=f.add_subplot(111) #make an ax.

# Create the contour for an altitude slice and call it 'cnt' (no jokes, please.)
# The '61' is the number of contours; you could use a vector of values to set
# levels manually if you wish.  get_cmap accepts any of the color map names
# from the colormap demo pic from the Matplotlib gallery; adding '_r' 
# reverses the colormap.
cnt=ax.contourf(a['Longitude'][:,:,0]*p,
                 p*a['Latitude'][:,:,0],
                 a['Rho'][:,:,0], 61, cmap=get_cmap('Spectral_r'))

# Configure axis.
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(r'$\rho$ at Altitude=%5.2f$km$' % (a['Altitude'][0,0,0]/1000.0))
f.suptitle('File=%s'%(a.attrs['file']))

# Add a colorbar and set the tick format to exponential notation.
cb=plt.colorbar(cnt)
cb.formatter=FormatStrFormatter('%7.2E')
cb.update_ticks()

# Add the quivers.
ax.quiver(a['Longitude'][:,:,0]*p, p*a['Latitude'][:,:,0],
          a['V!Dn!N (east)'][:,:,0],a['V!Dn!N (north)'][:,:,0])

# Draw to screen.
if plt.isinteractive():
    plt.draw() #In interactive mode, you just "draw".
else:
    # W/o interactive mode, "show" stops the user from typing more 
    # at the terminal until plots are drawn.
    plt.show()

