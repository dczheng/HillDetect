#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import h5py
import sys
import os
from astropy.io import fits

fn = sys.argv[1]
f = h5py.File( fn, 'r' )
data = f['data'][()]
bkg_s = f['bkg_s'][()]
bkg = f['bkg'][()]
vmin = data[data>0].min()
vmax = data.max()

fig, axs = plt.subplots( 2, 2, figsize=(2*5, 2*5) )

norm = mplc.LogNorm( vmin=vmin, vmax=vmax )
ax = axs[0,0]
img = ax.imshow( data, norm=norm )
ax.set_title( "data" )
plt.colorbar( img, ax=ax  )

ax = axs[0,1]
img = ax.imshow( bkg_s, norm=norm )
ax.set_title( "bkg_s" )
plt.colorbar( img, ax=ax  )

ax = axs[1,0]
img = ax.imshow( bkg, norm=norm )
ax.set_title( "bkg" )
plt.colorbar( img, ax=ax  )

ax = axs[1,1]
img = ax.imshow( data-bkg, norm=norm )
ax.set_title( "data - bkg" )
plt.colorbar( img, ax=ax  )

fig.savefig( fn[:-4] + 'png', dpi=300  )
