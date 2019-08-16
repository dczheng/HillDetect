#!/home/hkli/anaconda3/bin/python3
#!/usr/bin/env python3

import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from astropy.io import fits
from matplotlib import cm
import sys

fits_fn = sys.argv[1]
edge_fn = sys.argv[2]
region_fn = sys.argv[3]
catalog_fn = sys.argv[4]
out_fn = sys.argv[5]

fits_data = fits.open( fits_fn )[0].data[0,0,:,:]
#fits_data = fits.open( fits_fn )[0].data
m, n = fits_data.shape
m = int( 0.4 * m )
n = int( 0.4 * n )
fits_data = fits_data[ :m, :n ]

cmap = cm.Spectral
cmap = cm.seismic
cmap = cm.viridis
cmap = cm.magma
cmap = cm.ocean
img = plt.imshow( fits_data, norm=mplc.LogNorm(), cmap=cmap )
plt.colorbar( img )
ax = plt.gca()

lines = open( edge_fn ).readlines()
for i in range( len(lines) ):
    lines[i] = [float(k) for k in lines[i][:-1].split()]
    lines[i] = lines[i][1:]
m = len( lines ) // 2
for i in range( m ):
    plt.plot( lines[i*2], lines[i*2+1], 'r.', ms=0.5 )

'''
lines = open( region_fn ).readlines()
for i in range( len(lines) ):
    lines[i] = [float(k) for k in lines[i][:-1].split()]
    lines[i] = lines[i][1:]
m = len( lines ) // 2
for i in range( m ):
    plt.plot( lines[i*2], lines[i*2+1], 'r.', ms=0.5 )
'''

cata = np.loadtxt( catalog_fn )
plt.plot( cata[:,1], cata[:,2], 'y*', markersize=3 )

ax = plt.gca()
ax.invert_yaxis()

plt.savefig( out_fn )

