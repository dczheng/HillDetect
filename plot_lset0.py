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

fig, axs = plt.subplots( 2,1, figsize=(5, 5*2) )

f = h5py.File( "%s/%s_lset0_map.hdf5"%(outdir,  filepre), 'r' )
map0 = f[ '/map'][()]
f.close()

f = h5py.File( "%s/%s_lset0_lines.hdf5"%(outdir,  filepre), 'r' )
iters = f[ '/Lines' ].attrs[ 'iters' ]
lines = f[ '/Lines/lines-%i'%iters ][()]
print( iters )
f.close()

axs[0].imshow( map0, norm=mplc.LogNorm(), cmap=cm.jet )
axs[0].plot( lines[0,:], lines[1,:], 'b.', ms=0.5  )

img = np.zeros( map0.shape )
f = h5py.File( "%s/%s_map1.hdf5"%(outdir,  filepre), 'r' )
for i in range(N):
    print( "group: %i"%i )
    g = f[ '/Group%i'%i ]
    m0, n0 = g.attrs[ 'REPIXS' ]
    print( 'REPIXS: %i %i'%(m0, n0) )
    map1 = g[ 'map' ]
    m, n = map1.shape
    img[ m0:m0+m, n0:n0+n ] = map1
f.close()


axs[1].imshow( img, norm=mplc.LogNorm(), cmap=cm.jet )

for i in range(2):
    axs[i].invert_yaxis()

fig.savefig( 'lset0.png' )
