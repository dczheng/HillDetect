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
d = hdu.data

hdu = fits.open( sys.argv[2] )[0]
d2 = hdu.data

diff = d - d2

fig, axs = plt.subplots( 1, 3,  figsize=(3*5,5) )

vmin = d.min()
vmax = d.max()
print( "vmin: %e, vmax: %e"%( vmin, vmax ) )
norm = mplc.SymLogNorm( 1e-9, vmin=vmin, vmax=vmax  )

ax = axs[0]
img = ax.imshow( d, norm=norm  )
plt.colorbar( img, ax=ax, shrink=0.8 )
ax.set_title( os.path.basename( sys.argv[1] ) )

ax = axs[1]
img = ax.imshow( d2, norm=norm  )
plt.colorbar( img, ax=ax, shrink=0.8 )
ax.set_title( 'bkg' )

ax = axs[2]
img = ax.imshow( diff, norm=norm  )
plt.colorbar( img, ax=ax, shrink=0.8 )
ax.set_title( 'res' )

fig.savefig( sys.argv[3], dpi=300 )

