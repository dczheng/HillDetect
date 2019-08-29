#!/usr/bin/env python3
from my_work_env import *

ds = [ np.loadtxt( f ) for f in sys.argv[1:] ]
m, n = ds[0].shape
ds[0] = ds[0][ m//4:m//4*3, n//4:n//4*3 ]
fig, axs = plt.subplots( 2, 1, figsize=(5,2*5) )

for i in range(2):
    axs[i].imshow( ds[i], norm=mplc.LogNorm(), cmap=cm.jet )

fig.savefig( "comp.png" )

