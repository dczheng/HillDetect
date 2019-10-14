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

fig, axs = plt.subplots( 3, 1, figsize=(5, 5*3) )
#fits_data = fits.open( fits_file )[0].data[0,0,:,:]
hdu = fits.open( fits_file )[0]
print( hdu.data.shape )
fits_data = hdu.data 
fits_header = hdu.header



f = h5py.File( "%s/%s_lset0_edges_regs.hdf5"%( outdir, bname ), 'r' )
g = f[ 'Regs' ]
NRegs = g.attrs[ 'NRegs' ]
CRPIX = f.attrs[ 'CRPIX' ]
NAXIS = f.attrs[ 'NAXIS' ]

r = Rectangle( [CRPIX[1], CRPIX[0]], width=NAXIS[1], height=NAXIS[0],\
                fill=False )
axs[0].imshow( fits_data[0,0,:,:], norm=mplc.LogNorm(), cmap=cm.jet )
axs[0].add_artist( r )

print( "CRPIX:", CRPIX )
print( "NAXIS:", NAXIS )
print( "NRegs: %i"%NRegs )

for i in range( NRegs ):
    r = g[ 'reg-%i'%i ]
    y = r[0] + CRPIX[0]
    x = r[1] + CRPIX[1]
    #print( y, x )

    if CRPIX[0] in y or ( CRPIX[0]+NAXIS[0]-1 ) in y \
       or CRPIX[1] in x or ( CRPIX[1]+NAXIS[1]-1 ) in x:
       continue
    #print( y.min(), y.max(), x.min(), x.max(), (y.max()-y.min()) * (x.max()-x.min()) )
    #if (y.max()-y.min()) * (x.max()-x.min()) > 1000:
    #    print( y, x )
    fits_data[0, 0, y.min():y.max()+1, x.min():x.max()+1] = 0
    #if  i==100:
    #    break

#axs[1].imshow( hdu.data[0,0,CRPIX[0]:CRPIX[0]+NAXIS[0],CRPIX[1]:CRPIX[1]+NAXIS[1]],\
#        norm=mplc.LogNorm(), cmap=cm.jet )
axs[1].imshow( hdu.data[0,0,:,:],\
        norm=mplc.LogNorm(), cmap=cm.jet )
r = Rectangle( [CRPIX[1], CRPIX[0]], width=NAXIS[1], height=NAXIS[0],\
                fill=False )
axs[1].add_artist( r )

axs[2].imshow( hdu.data[0,0,CRPIX[0]:CRPIX[0]+NAXIS[0],CRPIX[1]:CRPIX[1]+NAXIS[1]],\
        norm=mplc.LogNorm(), cmap=cm.jet )
for i in range(3):
    axs[i].invert_yaxis()
hdu.writeto( fits_file[:-5] + '_mask.fits', overwrite=True )
fig.savefig( 'mask.png', dpi=300)


