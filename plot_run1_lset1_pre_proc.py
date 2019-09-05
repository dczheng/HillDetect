#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import sys

filepre = sys.argv[1] 
data_raw = np.loadtxt( "test_outputs/%s_map.dat"%filepre )
img = np.zeros( data_raw.shape )
N = 25
ds = [np.loadtxt("test_outputs/%s_%05i_before.dat"%(filepre,i)) for i in range(N)]

T = 5
fig, axs = plt.subplots( T, T, figsize=(10,10) )
for i in range(T*T):
    d = ds[i]
    d = d[1:,:]
    ax = axs[i//T][i%T]
    ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.jet )

plt.savefig( 'lset1_lset-before.png' )

plt.close()

ds = [np.loadtxt("test_outputs/%s_%05i_after.dat"%(filepre,i)) for i in range(N)]

T = 5
fig, axs = plt.subplots( T, T, figsize=(10,10) )
for i in range(T*T):
    d = ds[i]
    d = d[1:,:]
    ax = axs[i//T][i%T]
    ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.jet )

plt.savefig( 'lset1_lset-after.png' )
