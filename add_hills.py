#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import matplotlib.colors as mplc
import os

param_file = sys.argv[1]

ls = open( param_file ).readlines()
for l in ls:

    ll = l.split()
    if ll == []:
        continue

    if "FileName" == ll[0]:
        fits_file = ll[1]

    if "OutputDir" == ll[0]:
        OutputDir = ll[1]

    if "CuttingXStart" == ll[0]:
        CutXStart = float(ll[1])

    if "CuttingXEnd" == ll[0]:
        CutXEnd = float(ll[1])

    if "CuttingYStart" == ll[0]:
        CutYStart = float(ll[1])

    if "CuttingYEnd" == ll[0]:
        CutYEnd = float(ll[1])

bname = os.path.basename( fits_file )
bname_mask = '%s_mask.fits'%(bname[:-5])
fmt = '%-25s : %s'
print( fmt%( "fits file", fits_file ) )
print( fmt%( "Output dir", OutputDir ) )
print( fmt%( "base name", bname ) )
print( fmt%( "base name [mask]", bname_mask ) )
print( fmt%( "CutXStart", str(CutXStart) ) )
print( fmt%( "CutXEnd", str(CutXEnd) ) )
print( fmt%( "CutYStart", str(CutYStart) ) )
print( fmt%( "CutYEnd", str(CutYEnd) ) )

h5_file1 = "%s/%s/fof_regs.hdf5"%( OutputDir, bname )
h5_file2 = "%s/%s/fof_regs.hdf5"%( OutputDir, bname_mask )
out_fits_file = fits_file[:-5] + "_add.fits"

hdu = fits.open( fits_file )[0]
fits_data = hdu.data
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
m, n = fits_data.shape
m0 = int( m * CutYStart )
m1 = int( m * CutYEnd )
n0 = int( n * CutXStart )
n1 = int( n * CutXEnd )
axs[0,0].imshow( fits_data, norm=mplc.LogNorm() )
axs[0,1].imshow( d3[m0:m1, n0:n1], norm=mplc.LogNorm() )
axs[1,0].imshow( d1[m0:m1, n0:n1] )
axs[1,1].imshow( d1[m0:m1, n0:n1] )

fig.savefig( 'add_hills.png', dpi=300 )
fits.writeto( out_fits_file, data=d3, header=fits_header,\
        overwrite=True )

f_h5_1.close()
f_h5_2.close()
