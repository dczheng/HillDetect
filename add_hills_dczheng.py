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
from matplotlib import cm

param_file = sys.argv[1]
N = int( sys.argv[2] )

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

f_h5s = []
t = '' 
for i in range(N):
    fn  ="%s/%s%s.fits/fof_regs.hdf5"%(OutputDir, bname[:-5], t)
    print( 'load %s'%fn )
    f_h5s.append( h5py.File( fn, 'r' ) )
    t = t + '_mask'

hdu = fits.open( fits_file )[0]
fits_data = hdu.data
fits_header = hdu.header

f_h5_1 = h5py.File( h5_file1, 'r' )
f_h5_2 = h5py.File( h5_file2, 'r' )

t = 0.1
fs = 5
fx = 1 / ( 2 + 2*t )
fy = 1 / ( 1 + 2*t )
fig = plt.figure( figsize=(fs/fx, fs/fy) )
axs = [fig.add_axes([i*fx+fx*t, fy*t, fx, fy]) for i in range(2)]
for i in range(2):
    axs[i].set_xticks([])
    axs[i].set_yticks([])

ds = []
for i in range(N):
    ds.append( np.zeros( fits_data.shape ) )
d_add = np.zeros( fits_data.shape )

for fi in range(N):
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
        d_add[y, x] += fits_data[ y, x ]

m, n = fits_data.shape
m0 = int( m * CutYStart )
m1 = int( m * CutYEnd )
n0 = int( n * CutXStart )
n1 = int( n * CutXEnd )

#t = 'source'
#for i in range(N):
#    axs[0,i+1].imshow( ds[i], norm=mplc.LogNorm(), cmap=cm.jet )
#    axs[1,i+1].imshow( ds[i][m0:m1, n0:n1], norm=mplc.LogNorm(), cmap=cm.jet )
#    axs[0,i+1].set_title( t )
#    t = t + '_mask'

d = fits_data[m0:m1, n0:n1]
norm = mplc.LogNorm( vmin = d[d>0].min(), vmax = d.max() )
axs[0].imshow( fits_data[m0:m1, n0:n1], norm=norm, cmap=cm.jet )
axs[0].set_title( 'map', fontsize=20 )

axs[1].imshow( d_add[m0:m1, n0:n1], norm=norm, cmap=cm.jet )
axs[1].set_title( 'source', fontsize=20 )

fig.savefig( 'add_hills_dczheng.pdf', dpi=100 )

f_h5_1.close()
f_h5_2.close()
