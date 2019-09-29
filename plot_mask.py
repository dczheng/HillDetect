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
from matplotlib.patches import Rectangle

fits_file = sys.argv[1]
outdir = sys.argv[2]

bname = os.path.basename( fits_file )
print( bname )

fig, axs = plt.subplots( 2, 1, figsize=(5, 5*2) )
#fits_data = fits.open( fits_file )[0].data[0,0,:,:]
hdu = fits.open( fits_file )[0]
fits_data = hdu.data 
fits_header = hdu.header


axs[0].imshow( fits_data, norm=mplc.LogNorm(), cmap=cm.jet )

f = h5py.File( "%s/%s_lset0_edges_regs.hdf5"%( outdir, bname ), 'r' )
g = f[ 'Regs' ]
NRegs = g.attrs[ 'NRegs' ]
print( "NRegs: %i"%NRegs )

for i in range( NRegs ):
    r = g[ 'reg-%i'%i ]
    y = r[0]
    x = r[1]
    fits_data[y, x] = 0

axs[1].imshow( hdu.data, norm=mplc.LogNorm(), cmap=cm.jet )
for i in range(2):
    axs[i].invert_yaxis()
hdu.writeto( fits_file[:-5] + '_mask.fits', overwrite=True )
fig.savefig( 'mask.png' )
