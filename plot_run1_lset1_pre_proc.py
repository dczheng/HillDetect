#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import sys

filepre = sys.argv[1] 
data_raw = np.loadtxt( "run1_raw_map.dat" )
img = np.zeros( data_raw.shape )
N = 25
ds = [np.loadtxt("run1_%03i_before.dat"%i) for i in range(N)]

T = 5
fig, axs = plt.subplots( T, T, figsize=(10,10) )
for i in range(T*T):
    d = ds[i]
    d = d[1:,:]
    ax = axs[i//T][i%T]
    ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.jet )

plt.savefig( 'lset1_lset-before.png' )

plt.close()

ds = [np.loadtxt("run1_%03i_after.dat"%i) for i in range(N)]

T = 5
fig, axs = plt.subplots( T, T, figsize=(10,10) )
for i in range(T*T):
    d = ds[i]
    d = d[1:,:]
    ax = axs[i//T][i%T]
    ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.jet )

plt.savefig( 'lset1_lset-after.png' )
