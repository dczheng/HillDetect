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
data_raw = np.loadtxt( "./%s/%s_lset0_map.dat"%(outdir, filepre) )
img = np.zeros( data_raw.shape )
N = 25 

for i in range(N):
    d = open( "./%s/%s_lset1_map_%03i.dat"%(outdir, filepre, i) ).readlines()
    print( d )

exit()
vmax = ds[0].max()
for i in range(N):
    vmax = np.max( [vmax, ds[i].max()] )
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
    img[ymin, xmin:xmax] = vmax 
    img[ymax, xmin:xmax] = vmax 
    img[ymin:ymax, xmin] = vmax 
    img[ymin:ymax, xmax] = vmax 

plt.imshow( img, norm=mplc.LogNorm(), cmap=cm.jet )
plt.savefig( 'pos.png', dpi=200 )

