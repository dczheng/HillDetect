#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import matplotlib.colors as mplc

fits_file = sys.argv[1]
h5_file1 = sys.argv[2]
h5_file2 = sys.argv[3]
out_fits_file = fits_file[:-5] + "_add.fits"

hdu = fits.open( fits_file )[0]
fits_data = hdu.data[0,0,:,:]
fits_header = hdu.header

f_h5_1 = h5py.File( h5_file1, 'r' )
f_h5_2 = h5py.File( h5_file2, 'r' )

fig, axs = plt.subplots( 2, 2, figsize=(10,10) )

d1 = np.zeros( fits_data.shape )
d2 = np.zeros( fits_data.shape )
d3 = np.zeros( fits_data.shape )

ds = [ d1, d2 ]
f_h5s = [ f_h5_1, f_h5_2 ]

for fi in range(2):
    f = f_h5s[fi]
    Ngroup = f.attrs[ 'Ngroup' ]
    for i in range( Ngroup ):
        g = f[ 'Group%i'%i ] 
        Nreg = g.attrs[ 'NReg' ]
        for j in range( Nreg ):
            reg = g[ 'Reg%i'%j ]
            xy = reg[ 'region' ]
            x = xy[1,:]
            y = xy[0,:]
            ds[fi][ y, x ] = 1
        d3[ y, x] += fits_data[ y, x ]
axs[0,0].imshow( fits_data, norm=mplc.LogNorm() )
axs[0,1].imshow( d3, norm=mplc.LogNorm() )
axs[1,0].imshow( d1 )
axs[1,1].imshow( d1 )

idx = [ 3, 4 ]
for i in idx:
    del( fits_header['NAXIS%i'%i] )
    del( fits_header['CTYPE%i'%i] )
    del( fits_header['CUNIT%i'%i] )
    del( fits_header['CRPIX%i'%i] )
    del( fits_header['CDELT%i'%i] )
    del( fits_header['CRVAL%i'%i] )


fig.savefig( 'add_hills.png', dpi=300 )
fits.writeto( out_fits_file, data=d3, header=fits_header,\
        overwrite=True )

f_h5_1.close()
f_h5_2.close()
