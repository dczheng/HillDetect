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

fig, axs = plt.subplots( m, n, figsize=(n*5, m*5) )

fmap = h5py.File( "%s/%s_map1.hdf5"%(outdir,  filepre), 'r' )
fregs = h5py.File( "%s/%s_fof_regs.hdf5"%(outdir,  filepre), 'r' )
for i in range(N):
    ax = axs[ i//n, i%n ]
    img = fmap[ '/Group%i/map'%i ][()]
    ax.imshow( img, norm=mplc.LogNorm(), cmap=cm.jet )

    NRegs = fregs['/Group%i'%i].attrs[ 'NRegs' ]
    #print( NRegs )
    for j in range(NRegs):
        r = fregs[ '/Group%i/Reg-%i'%(i,j) ][()]
        ax.plot( r[1,:], r[0,:], '*' )


fmap.close()
fregs.close()
fig.savefig( 'fof.png' )
