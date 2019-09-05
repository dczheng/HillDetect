#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import sys

filepre = sys.argv[1] 
data_raw = np.loadtxt( "./test_outputs/%s_map.dat"%filepre )
img = np.zeros( data_raw.shape )
N = 25
ds = [np.loadtxt("test_outputs/%s_%05i_map.dat"%(filepre,i)) for i in range(N)]
ll = [open( './test_outputs/%s_%05i_lset_lines.dat'%(filepre, i) ).readlines() for i in range(N) ]

for i in range(N):
    d = ds[i]
    xmin = int(d[0,0])
    ymin = int(d[0,1])
    xmax = int( d[0,2] )
    ymax = int( d[0,3] )
    d = d[1:,:]
    m, n = d.shape
    print( d.shape )
    print( ymax-ymin, xmax-xmin )
    img[ ymin:ymax, xmin:xmax ] += d

plt.imshow( img, norm=mplc.LogNorm(), cmap=cm.jet )
plt.savefig( 'lset1_lset-pos.png' )

plt.close()

T = 5
fig, axs = plt.subplots( T, T, figsize=(10,10) )
for i in range(T*T):
    d = ds[i]
    d = d[1:,:]
    ax = axs[i//T][i%T]
    ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.jet )

for i in range(T*T):
    ax = axs[i//T][i%T]
    d = ds[i]
    xmin = int(d[0,0])
    ymin = int(d[0,1])

    l = ll[i]
    m = len( l )
    x = np.array([float(i) for i in l[m-2].split()[4:]])
    y = np.array([float(i) for i in l[m-1].split()[4:]])
    x = x - xmin
    y = y - ymin
    print( len(x) )

    ax.plot( x, y, 'b.', ms=1  )

plt.savefig( 'lset1_lset.png' )
