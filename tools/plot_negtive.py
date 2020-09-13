#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import sys
import os
from astropy.io import fits

hdu = fits.open( sys.argv[1] )[0]
dd = hdu.data

fig, axs = plt.subplots( 1, 3,  figsize=(3*5,5) )

d = dd
vmin = d.min()
vmax = d.max()
norm = mplc.SymLogNorm( 1e-9, vmin=vmin, vmax=vmax  )
ax = axs[0]
img = ax.imshow( d, norm=norm, cmap=cm.jet  )
ax.set_title( "all value" )
plt.colorbar( img, ax=ax, shrink=0.8 )

d = dd.copy()
ax = axs[1]
ax.set_title( "postive value" )
img = ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.jet  )
plt.colorbar( img, ax=ax, shrink=0.8 )

d = -dd.copy()
ax = axs[2]
ax.set_title( "negtive value" )
img = ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.jet  )
plt.colorbar( img, ax=ax, shrink=0.8 )

fig.savefig( sys.argv[1][:-5] + '_value.png', dpi=300 )

