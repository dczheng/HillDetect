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

sep = '-' * 50
params = {} 
data = {}

def load_parameters( fn ):
    lines = open( fn ).readlines()
    global parmas
    for l in lines:
        ll = l.split()
        if ll == []:
            continue
        params[ll[0]] = ll[1]

def plot_lset():
   
    print( "plot lset ..." )

    img = data['fits_img']

    vmin = img[img>0].min()
    vmax = img.max()
    norm = mplc.LogNorm( vmin=vmin, vmax=vmax )

    fig, axs = plt.subplots( 2,2, figsize=(2*5, 2*5) )

    axs[0,0].imshow( img, norm=norm, cmap=cm.jet )
    axs[0,1].imshow( img[data['mcut0']: data['mcut1'], data['ncut0']: data['ncut1']],\
                    norm=norm, cmap=cm.jet )

    lset_lines = h5py.File( "%s/lset_lines.hdf5"%params['data_dir'], 'r' )
    iters = lset_lines.attrs[ 'iters' ]
    print( "iters: %i"%iters )
    line = lset_lines[ 'iter%i/line'%iters ][()]
    y = line[0,:] - data['mcut0'] 
    x = line[1,:] - data['ncut0'] 
    axs[0,1].plot( x, y, 'b.', ms=0.1  )

    fig.savefig( '%s/lset.png'%params['plot_dir'], dpi=300 )

    return

    iters = lset_lines.attrs[ 'iters' ]
    print( "iters: %i"%iters )
    line = lset_lines[ 'iter%i/line'%iters ][()]
    y = line[0,:] - mcut0 
    x = line[1,:] - ncut0
    axs[0,1].plot( x, y, 'b.', ms=0.1  )

    img = np.zeros( fits_dd.shape )
    global Ngroup
    for i in range(Ngroup):
        g = maps[ '/Group%i'%i ]
        m0, n0 = g.attrs[ 'CRPIX' ]
        mm = g[ 'map' ]
        m, n = mm.shape
        m0 = m0 - mcut0
        n0 = n0 - ncut0
        img[ m0:m0+m, n0:n0+n ] = mm
    axs[1,0].imshow( img, norm=norm, cmap=cm.jet )

    fits_mask = fits_d.copy()
    NN = fof_regs.attrs['Ngroup']
    Ntot = 0
    for i in range( NN ):
        g = fof_regs[ 'Group%i'%i ]
        Nreg = g.attrs[ 'NReg' ]
        Ntot += Nreg
        for j in range( Nreg ):
            reg = g[ 'Reg%i'%j ]
            xy = reg[ 'region' ]
            fits_mask[ xy[0,:], xy[1,:] ] = fits_d.min() 

    print( "Ngroup: %i Ntot: %i"%(NN, Ntot) )

    #print( fits_mask.max(), fits_mask.min(), len(fits_mask[fits_mask>0]) )
    img = axs[1,1].imshow( fits_mask, norm=norm, cmap=cm.jet )
    #plt.colorbar( img, ax=axs[1,1] )

    fig.savefig( '%s/lset.png'%plot_output_dir,dpi=300 )
    fits.writeto( '%s_mask.fits'%(fits_file[:-5]),\
                data=fits_mask, header=fits_h, overwrite=True )

def show_params():

    global params

    n = 0
    for k in params:
        if len(k) > n: n = len(k)
    fmt = '%%-%is   %%s'%n
    print( sep )
    print( "HILLDETECT PARAMETERS:" )
    for k in params:
        print( fmt%(k, params[k]) )
    print( sep )

def main():
   
    global params, data

    load_parameters( sys.argv[1] )
    bname = os.path.basename( params['FileName'] )
    params['plot_dir'] = os.path.join( params['OutputDir'], 'plot_' + bname )
    params['data_dir'] = os.path.join( params['OutputDir'], bname )
    show_params()

    if not os.path.exists( params['plot_dir'] ):
        os.makedirs( params['plot_dir'] )

    data['fits_img'] = fits.open( params['FileName'] )[0].data
    m, n = data['fits_img'].shape
    data['mcut0'] = int( m * float(params['CuttingYStart']) )
    data['mcut1'] = int( m * float(params['CuttingYEnd']) )
    data['ncut0'] = int( n * float(params['CuttingXStart']) )
    data['ncut1'] = int( n * float(params['CuttingXEnd']) )

    plot_lset()



if __name__ == '__main__':
    main()

exit()
    
ls = open( param_file ).readlines()
for l in ls:
    ll = l.split()
    if ll == []:
        continue

    if "FileName" == ll[0]:
        fits_file = ll[1]
    if "OutputDir" == ll[0]:
        output_dir = ll[1]
    if "OnlyFoF" == ll[0]:
        onlyfof = int(ll[1])

    if "CuttingXStart" == ll[0]:
        CutXStart = float(ll[1])
    if "CuttingXEnd" == ll[0]:
        CutXEnd = float(ll[1])
    if "CuttingYStart" == ll[0]:
        CutYStart = float(ll[1])
    if "CuttingYEnd" == ll[0]:
        CutYEnd = float(ll[1])
bname = os.path.basename( fits_file )
fmt = "%-25s : %s"

if not onlyfof:
    Ngroup = int( sys.argv[2] )
    print( fmt%("Ngroup", str(Ngroup)) )
else:
    fits_file_ori = sys.argv[2]
    print( fmt%("fits file orig", fits_file_ori) )

 
print( fmt%( "fits file", fits_file  ) )
print( fmt%("base name", bname) )
print( fmt%("HillDetect output dir", output_dir) )
output_dir = "%s/%s"%(output_dir, bname)
print( fmt%("data output dir", output_dir) )
plot_output_dir = "%s"%bname
print( fmt%("plot output dir", plot_output_dir) )
print( fmt%("onlyfof", str(onlyfof)) )
print( fmt%( "CutXStart", str(CutXStart) ) )
print( fmt%( "CutXEnd", str(CutXEnd) ) )
print( fmt%( "CutYStart", str(CutYStart) ) )
print( fmt%( "CutYEnd", str(CutYEnd) ) )



if not onlyfof:
    fs =  [\
    "%s/map.hdf5"%output_dir,\
    "%s/map_after.hdf5"%output_dir,\
    "%s/fof_regs.hdf5"%output_dir,\
    "%s/lset_regs.hdf5"%output_dir,\
    "%s/lset_lines.hdf5"%output_dir,\
    ]
    h5_fs = []
    for f in fs:
        print( "load '%s'"%f )
        h5_fs.append( h5py.File( f, 'r' ) )
    maps, maps_after, fof_regs, lset_regs, lset_lines = h5_fs

    if Ngroup == 0:
        Ngroup = maps.attrs[ "Ngroup" ]
        print( fmt%("set Ngroup to", str(Ngroup)) )
else:
    fn = "%s/only_fof_regs.hdf5"%output_dir
    print( "load '%s'"%fn )
    h5_f = h5py.File( fn, 'r' )

def plot_only_fof():
    
    hdu = fits.open( sys.argv[2] )[0]
    fits_d_ori = hdu.data

    vmin = fits_d_ori[fits_d_ori>0].min()
    vmax = fits_d_ori.max()
    norm = mplc.LogNorm( vmin=vmin, vmax=vmax )
    fig, axs = plt.subplots( 2, 2, figsize=(20, 20) )
    axs[0,0].imshow( fits_d_ori, norm=norm, cmap=cm.jet )
    axs[0,1].imshow( fits_d_ori[mcut0:mcut1, ncut0:ncut1], norm=norm, cmap=cm.jet )

    d = np.zeros( fits_d.shape )
    g = h5_f[ 'Group0' ]
    Nreg = g.attrs[ 'NReg' ]
    print( "Nreg: %i"%Nreg )
    for i in range( Nreg ):
        reg = g[ 'Reg%i'%i ]
        r = reg[ 'region'][()]
        x = r[1,:]
        y = r[0,:]
        d[ y, x ] = 1
    axs[1,0].imshow( fits_d[mcut0:mcut1, ncut0:ncut1], norm=norm, cmap=cm.jet )


    '''
    t = [ -1, 1 ]
    for i in range(Nreg):
        reg = g[ 'Reg%i'%i ]
        r = reg[ 'region'][()]
        x = r[1,:]
        y = r[0,:]
        for k in range(len(x)):

            if d[ y[k], x[k] ] == 0:
                continue

            for xx in t:
                for yy in t:
                    if y[k]+yy < mcut0 or y[k]+yy > mcut1 or \
                        x[k]+xx < ncut0 or x[k]+xx > ncut1:
                        continue

                    if d[ y[k]+yy, x[k]+xx ] == 0:
                        #print( 'ok')
                        axs[1,0].plot( [x[k]-ncut0], [y[k]-mcut0], 'k.', ms=0.1 )
                        break
    '''
    axs[1,1].imshow( d[mcut0:mcut1, ncut0:ncut1], norm=mplc.LogNorm() )



    fig.savefig( '%s/only_fof.png'%plot_output_dir,dpi=300 )

def plot_maps():

    print( "plot maps ..." )

    n = 5

    if n > Ngroup:
        n = Ngroup // 2

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
        img1[img1==img1.min()] =  -1
        ax1.imshow( img1, norm=mplc.LogNorm(), cmap=cm.jet )

        group = fof_regs['/Group%i'%i]
        NRegs = group.attrs[ 'NReg' ]
        #print( NRegs )
        f = 0
        for j in range(NRegs):
            reg = group[ 'Reg%i'%j ]
            r = reg[ 'region'][()]
            c = reg.attrs[ 'center' ]
            xyerr = reg.attrs[ 'xyerr' ]
            x = r[1,:]
            y = r[0,:]
            y = y - crpix[0]
            x = x - crpix[1]
            c = c - crpix
            ax.plot( x, y, 'k*' )
            ax1.plot( x, y, 'k*' )

            ax.plot( [c[1]], [c[0]], '*', ms=40 )
            ax1.plot( [c[1]], [c[0]], '*', ms=40 )


            ax.errorbar( [c[1]], [c[0]], yerr=[xyerr[0]], xerr=[xyerr[1]], linewidth=10 )
            ax1.errorbar( [c[1]], [c[0]], yerr=[xyerr[0]], xerr=[xyerr[1]], linewidth=10 )

    fig.savefig( '%s/map.png'%plot_output_dir)
    fig1.savefig( '%s/map_after.png'%plot_output_dir )

if not onlyfof:
    plot_maps()
    plot_lset()
else:
    plot_only_fof()

if not onlyfof:
    for f in h5_fs:
        f.close()
else:
    h5_f.close()
