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

fig, axs = plt.subplots( m, n, figsize=(n*10, m*10) )
fig1, axs1 = plt.subplots( m, n, figsize=(n*10, m*10) )

fmap = h5py.File( "%s/%s_map1.hdf5"%(outdir,  filepre), 'r' )
fmap1 = h5py.File( "%s/%s_map1_after.hdf5"%(outdir,  filepre), 'r' )
fregs = h5py.File( "%s/%s_fof_regs.hdf5"%(outdir,  filepre), 'r' )
for i in range(N):
    ax = axs[ i//n, i%n ]
    ax1 = axs1[ i//n, i%n ]

    img = fmap[ '/Group%i/map'%i ][()]
    ax.imshow( img, norm=mplc.LogNorm(), cmap=cm.jet )

    img1 = fmap1[ '/Group%i/map'%i ][()]
    ax1.imshow( img1, norm=mplc.LogNorm(), cmap=cm.jet )

    NRegs = fregs['/Group%i'%i].attrs[ 'NRegs' ]
    #print( NRegs )
    f = 0
    for j in range(NRegs):
        r = fregs[ '/Group%i/Reg-%i'%(i,j) ][()]
        ax.plot( r[1,:], r[0,:], 'k*' )
        ax1.plot( r[1,:], r[0,:], 'k*' )
        f = f + img[ r[0,:], r[1,:] ].sum()
    print( "flux_fof: %.3e,  flux_tot: %.3e, percent: %.2f%%"%(f, img.sum(),\
    f / img.sum() * 100) )
    mm, nn = img.shape
    ax.text( 0.1*nn, 0.1*mm, "%.2f%%"%(f/img.sum()*100), fontsize=60 )
    ax1.text( 0.1*nn, 0.1*mm, "%.2f%%"%(f/img.sum()*100), fontsize=60 )


fmap.close()
fmap1.close()
fregs.close()
fig.savefig( 'fof.png')
fig1.savefig( 'fof_after.png')
