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

fig, ax = plt.subplots( 1, 1,  figsize=(5,5) )

vmin = d.min()
vmax = d.max()
norm = mplc.SymLogNorm( 1e-9, vmin=vmin, vmax=vmax  )
img = ax.imshow( d, norm=norm  )
plt.colorbar( img, ax=ax, shrink=0.8 )

fig.savefig( sys.argv[1][:-4] + 'png' )

