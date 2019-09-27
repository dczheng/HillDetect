#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import sys

filepre = sys.argv[1] 
outdir = "new_outputs"
data_raw = np.loadtxt( "%s/%s_lset0_map.dat"%(outdir, filepre) )
print( data_raw.shape )

fig, axs = plt.subplots( 2,1, figsize=(5, 5*2) )

lines    = open( './%s/%s_lset0_lines.dat'%(outdir, filepre) ).readlines()
m = len( lines )
print( m )
if m % 2 != 0:
    print( "error" )
    exit()

axs[0].imshow( data_raw, norm=mplc.LogNorm(), cmap=cm.jet )
axs[0].set_title( 'lset lines' )

x = np.array([float(i) for i in lines[m-2].split()[4:]])
y = np.array([float(i) for i in lines[m-1].split()[4:]])
#print( len(x), x.min(), x.max(), y.min(), y.max() )
axs[0].plot( x, y, 'b.', ms=0.5 )

lines    = open( './%s/%s_lset0_edges.dat'%(outdir, filepre) ).readlines()
m = len( lines )
print( m )
if m % 2 != 0:
    print( "error" )
    exit()

axs[1].imshow( data_raw, norm=mplc.LogNorm(), cmap=cm.jet )
axs[1].set_title( 'fof edges' )

i = 0
index = 0
while( index<25 ):
#while( i<m ):
    x = np.array([float(i) for i in lines[i].split()[1:]])
    y = np.array([float(i) for i in lines[i+1].split()[1:]])
    #print( len(x), x.min(), x.max(), y.min(), y.max() )
    axs[1].plot( x, y, 'b.', ms=0.5 )
    if 0 not in x and 0 not in y:
        index += 1
    i += 2

fig.savefig( 'lset0.png' )

