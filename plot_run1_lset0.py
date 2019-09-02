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
print( data_raw.shape )

lines    = open( './test_outputs/%s_lset_lines.dat'%filepre ).readlines()
m = len( lines )
print( m )
if m % 2 != 0:
    print( "error" )
    exit()

plt.imshow( data_raw, norm=mplc.LogNorm(), cmap=cm.jet )

x = np.array([float(i) for i in lines[m-2].split()[4:]])
y = np.array([float(i) for i in lines[m-1].split()[4:]])
#print( len(x), x.min(), x.max(), y.min(), y.max() )
plt.plot( x, y, 'b.', ms=0.5 )

plt.savefig( 'lset0_lset.png' )

plt.close()

lines    = open( './test_outputs/%s_edge.dat'%filepre ).readlines()
m = len( lines )
print( m )
if m % 2 != 0:
    print( "error" )
    exit()

plt.imshow( data_raw, norm=mplc.LogNorm(), cmap=cm.jet )

i = 0
index = 0
while( index<25 ):
#while( i<m ):
    x = np.array([float(i) for i in lines[i].split()[1:]])
    y = np.array([float(i) for i in lines[i+1].split()[1:]])
    #print( len(x), x.min(), x.max(), y.min(), y.max() )
    plt.plot( x, y, 'b.', ms=0.5 )
    if 0 not in x and 0 not in y:
        index += 1
    i += 2

plt.savefig( 'lset0_fof.png' )

