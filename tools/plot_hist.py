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
from mpl_toolkits.mplot3d import Axes3D

hdu = fits.open( sys.argv[1] )[0]
d = hdu.data
m, n = d.shape
d = d.reshape( m*n, 1 )
d1 = d[d>0]
d2 = -d[d<0]
d = np.hstack( [d1, d2] )
print( d.shape )

v = 2e-4
print( d[d>v].shape )
d = d[d<v]

plt.hist( d, bins=1000 )
plt.savefig( sys.argv[1][:-5] + '_hist.png', dpi=300 )

