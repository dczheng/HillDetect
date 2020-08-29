#!/usr/bin/env python3

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

fits_file = sys.argv[1]
output_dir = sys.argv[2]

bname = os.path.basename( fits_file )
fmt = "%-25s : %s"
print( fmt%( "fits file", fits_file  ) )
print( fmt%("base name", bname) )
print( fmt%("HillDetect output dir", output_dir) )
output_dir = "%s/%s"%(output_dir, bname)
print( fmt%("data output dir", output_dir) )
plot_output_dir = "%s_plot"%bname
print( fmt%("plot output dir", plot_output_dir) )

if not os.path.exists( plot_output_dir ):
    os.mkdir( plot_output_dir )

fs =  [\
"%s/lset_regs.hdf5"%output_dir,\
"%s/lset_lines.hdf5"%output_dir,\
]
h5_fs = []
for f in fs:
    print( "load '%s'"%f )
    h5_fs.append( h5py.File( f, 'r' ) )
lset_regs, lset_lines = h5_fs

hdu = fits.open( fits_file )[0]
fits_h = hdu.header
fits_d = hdu.data[0,0,:,:]

Nreg = lset_regs.attrs[ "NReg" ]
print( fmt%("Nreg", str(Nreg)) )

fig, axs = plt.subplots( 2, 1, figsize=(5,10) )
norm = mplc.LogNorm()
axs[0].imshow( fits_d, norm=norm )
axs[1].imshow( fits_d, norm=norm )

cutstart = lset_regs.attrs[ 'CutStart'  ]
cutend = lset_regs.attrs[ 'CutEnd'  ]
print( 'cut start: ', cutstart  )
print( 'cut end: ', cutend  )

iters = lset_lines.attrs['iters']
line = lset_lines[ 'line-%i'%iters ][()]
y = line[0,:]
x = line[1,:]
axs[0].plot( x, y, 'b.', ms=0.3 )

for i in range(Nreg):
    r = lset_regs[ '/Reg-%i'%i ][()]
    y = r[0,:]
    x = r[1,:]
    axs[1].plot( x, y, 'k.', ms=0.3 )

    c = lset_regs.attrs[ 'Center-%i'%i ]
    xyerr = lset_regs.attrs[ 'XYerr-%i'%i ]
    sig = lset_regs.attrs[ 'Sigma-%i'%i ]
    print( "center:", c )
    print( "err:", xyerr )
    print( "sigma:", sig )
    axs[0].plot( [c[1]], [c[0]], '*' )
    axs[0].errorbar( [c[1]], [c[0]], yerr=[xyerr[0]], xerr=[xyerr[1]] )
    axs[1].plot( [c[1]], [c[0]], '*' )
    axs[1].text( c[1], c[0], '%.2e'%sig )


fig.savefig( "%s/lset_single.png"%plot_output_dir )

for f in h5_fs:
    f.close()
