#!/home/hkli/anaconda3/bin/python3

from matplotlib import use
use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import h5py
import sys
import os
from astropy.io import fits
from matplotlib.patches import Rectangle

fits_file = sys.argv[1]
outdir = sys.argv[2]

bname = os.path.basename( fits_file )
print( bname )

fig, axs = plt.subplots( 2, 1, figsize=(5, 5*2) )
fits_data = fits.open( fits_file )[0].data[0,0,:,:]


f = h5py.File( "%s/%s_fof_regs.hdf5"%( outdir, bname ), 'r' )
NGroups = f.attrs[ 'NGroups' ]
Gcrpix = f.attrs[ 'GCRPIX' ] 
Gnaxis = f.attrs[ 'GNAXIS' ] 
print( "GCRPIX: %i %i"%(Gcrpix[0], Gcrpix[1]) )
print( "Gnaxis: %i %i"%(Gnaxis[0], Gnaxis[1]) )
print( NGroups )


fits_data[Gcrpix[0]+Gnaxis[0], Gcrpix[1]:Gcrpix[1]+Gnaxis[1]] = 0 
fits_data[Gcrpix[0]:Gcrpix[0]+Gnaxis[0], Gcrpix[1]] = 0 
fits_data[Gcrpix[0]:Gcrpix[0]+Gnaxis[0], Gcrpix[1]+Gnaxis[1]] = 0 

axs[0].imshow( fits_data, norm=mplc.LogNorm(), cmap=cm.jet )

r = Rectangle( [Gcrpix[1], Gcrpix[0]], width=Gnaxis[1], height=Gnaxis[0],\
                fill=False )
axs[0].add_artist( r )

axs[1].imshow( fits_data[Gcrpix[0]:Gcrpix[0]+Gnaxis[0],\
                         Gcrpix[1]:Gcrpix[1]+Gnaxis[1]],\
                         norm=mplc.LogNorm(), \
                         cmap = cm.jet )

index = 0
for gidx in range( NGroups ):
    g = f[ '/Group%i'%gidx ]
    crpixy, crpixx = g.attrs[ 'CRPIX' ]
    Nregs = g.attrs[ 'NRegs' ]
    print( "CRPIX: %i %i"%( crpixy, crpixx ) )
    print( "Nregs: %i"%Nregs )

    for j in range( Nregs ):
        c = g.attrs[ 'Center-%i'%j ]
        flux = g.attrs[ 'Flux-%i'%j ]
        y = c[0] + crpixy
        x = c[1] + crpixx
        axs[1].plot( [x], [y], 'b*', ms=2  )
        #axs[1].text( x, y, "%i [%.2f]"%(index, flux), fontsize=8 )
        axs[1].text( x, y, "%i"%(index), fontsize=8 )
        index += 1
for i in range(2):
    axs[i].invert_yaxis()
print( index )

fig.savefig( 'final.png' )

