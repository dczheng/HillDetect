#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import h5py
import sys

outdir = sys.argv[1]
filepre = sys.argv[2] 
N = int( sys.argv[3] )

n = 5
m = N // n
if ( N % n != 0 ):
    m += 1

fig1, axs1 = plt.subplots( m, n, figsize=(n*5, m*5) )
fig2, axs2 = plt.subplots( m, n, figsize=(n*5, m*5) )

fmap = h5py.File( "%s/%s_lset1_map.hdf5"%(outdir,  filepre), 'r' )
fmap_after = h5py.File( "%s/%s_lset1_map_after.hdf5"%(outdir,  filepre), 'r' )
flines = h5py.File( "%s/%s_lset1_lines.hdf5"%(outdir,  filepre), 'r' )
for i in range(N):
    ax1 = axs1[ i//n, i%n ]
    img1 = fmap[ '/Group%i/map'%i ][()]
    ax1.imshow( img1, norm=mplc.LogNorm(), cmap=cm.jet )

    ax2 = axs2[ i//n, i%n ]
    img2 = fmap_after[ '/Group%i/map'%i ][()]
    ax2.imshow( img2, norm=mplc.LogNorm(), cmap=cm.jet )

    iters = flines['/Group%i'%i].attrs[ 'iters' ]
    l = flines['/Group%i/lines-%i'%(i, iters//3)][()]
    ax1.plot( l[0,:], l[1,:], 'o')
    ax2.plot( l[0,:], l[1,:], 'o')


fmap.close()
fmap_after.close()
flines.close()
fig1.savefig( 'lset1-before.png' )
fig2.savefig( 'lset1-after.png' )
