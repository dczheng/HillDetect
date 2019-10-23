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
Ngroup = int( sys.argv[3] )

bname = os.path.basename( fits_file )
fmt = "%-25s : %s"
print( fmt%( "fits file", fits_file  ) )
print( fmt%("base name", bname) )
print( fmt%("Ngroup", str(Ngroup)) )
print( fmt%("HillDetect output dir", output_dir) )
output_dir = "%s/%s"%(output_dir, bname)
print( fmt%("data output dir", output_dir) )
plot_output_dir = "%s_plot"%bname
print( fmt%("plot output dir", plot_output_dir) )

if not os.path.exists( plot_output_dir ):
    os.mkdir( plot_output_dir )

fs =  [\
"%s/map.hdf5"%output_dir,\
"%s/map_after.hdf5"%output_dir,\
"%s/fof_regs.hdf5"%output_dir,\
"%s/lset_edges_regs.hdf5"%output_dir,\
"%s/lset_lines.hdf5"%output_dir,\
]
h5_fs = []
for f in fs:
    print( "load '%s'"%f )
    h5_fs.append( h5py.File( f, 'r' ) )
maps, maps_after, fof_regs, lset_er, lset_lines = h5_fs

hdu = fits.open( fits_file )[0]
fits_h = hdu.header
fits_d = hdu.data[0,0,:,:]

def plot_maps():

    print( "plot maps ..." )
    n = 5
    m = Ngroup // n
    if ( Ngroup % n != 0 ):
        m += 1

    fig, axs = plt.subplots( m, n, figsize=(n*10, m*10) )
    fig1, axs1 = plt.subplots( m, n, figsize=(n*10, m*10) )

    for i in range(Ngroup):
        ax = axs[ i//n, i%n ]
        ax1 = axs1[ i//n, i%n ]

        crpix = maps[ '/Group%i'%i ].attrs[ 'CRPIX' ]
        #print( crpix )
        img = maps[ '/Group%i/map'%i ][()]
        ax.imshow( img, norm=mplc.LogNorm(), cmap=cm.jet )

        img1 = maps_after[ '/Group%i/map'%i ][()]
        ax1.imshow( img1, norm=mplc.LogNorm(), cmap=cm.jet )

        NRegs = fof_regs['/Group%i'%i].attrs[ 'NReg' ]
        #print( NRegs )
        f = 0
        for j in range(NRegs):
            r = fof_regs[ '/Group%i/Reg-%i'%(i,j) ][()]
            x = r[1,:]
            y = r[0,:]
            y = y - crpix[0]
            x = x - crpix[1]
            ax.plot( x, y, 'k*' )
            ax1.plot( x, y, 'k*' )
            f = f + img[ y, x ].sum()
        print( "flux_fof: %.3e,  flux_tot: %.3e, percent: %.2f%%"%(f, img.sum(),\
        f / img.sum() * 100) )
        mm, nn = img.shape
        ax.text( 0.1*nn, 0.1*mm, "%.2f%%"%(f/img.sum()*100), fontsize=60 )
        ax1.text( 0.1*nn, 0.1*mm, "%.2f%%"%(f/img.sum()*100), fontsize=60 )


    fig.savefig( '%s/map.png'%plot_output_dir)
    fig1.savefig( '%s/map_after.png'%plot_output_dir )

def plot_lset():
    
    print( "plot lset ..." )

    fig, axs = plt.subplots( 2,2, figsize=(2*5, 2*5) )
    cutstart = maps.attrs[ 'CutStart' ]
    cutend = maps.attrs[ 'CutEnd' ]
    print( 'cut start: ', cutstart )
    print( 'cut end: ', cutend )
    fits_dd = fits_d[ cutstart[0]:cutend[0], cutstart[1]:cutend[1] ]
    vmin = fits_d[fits_d>0].min()
    vmax = fits_d.max()
    norm = mplc.LogNorm( vmin=vmin, vmax=vmax )

    axs[0,0].imshow( fits_d, norm=norm, cmap=cm.jet )
    axs[0,1].imshow( fits_dd, norm=norm, cmap=cm.jet )

    iters = lset_lines.attrs[ 'iters' ]
    print( "iters: %i"%iters )
    line = lset_lines[ 'lines-%i'%iters ][()]
    y = line[0,:] - cutstart[0]
    x = line[1,:] - cutstart[1]
    axs[0,1].plot( x, y, 'b.', ms=0.3  )

    img = np.zeros( fits_dd.shape )
    for i in range(Ngroup):
        g = maps[ '/Group%i'%i ]
        m0, n0 = g.attrs[ 'CRPIX' ]
        mm = g[ 'map' ]
        m, n = mm.shape
        m0 = m0 - cutstart[0]
        n0 = n0 - cutstart[1]
        img[ m0:m0+m, n0:n0+n ] = mm
    axs[1,0].imshow( img, norm=norm, cmap=cm.jet )

    fits_mask = fits_d.copy()
    NN = maps.attrs['NGroup']
    for i in range(NN):
        g = maps[ '/Group%i'%i ]
        m0, n0 = g.attrs[ 'CRPIX' ]
        mm = g[ 'map' ]
        m, n = mm.shape
        fits_mask[m0:m0+m, n0:n0+n] = -1 

    axs[1,1].imshow( fits_mask, norm=norm, cmap=cm.jet )

    fig.savefig( '%s/lset.png'%plot_output_dir,dpi=300 )
    fits.writeto( '%s/%s_mask.fits'%(plot_output_dir, bname[:-5]),\
                data=fits_mask, header=fits_h, overwrite=True )

plot_maps()
plot_lset()

for f in h5_fs:
    f.close()
