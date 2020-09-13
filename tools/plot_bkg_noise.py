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
from load_parameters import load_parameters 

params = load_parameters(sys.argv[1])
#print( params )
output_dir = os.path.join(params['OutputDir'], os.path.basename( params['FileName'] ))
print( 'output_dir:', output_dir )
interp_name = { '0': 'bilinear', '1': 'bicubic' }
fn_bkg = os.path.join( output_dir, "bkg_%s.hdf5"%interp_name[params['BkgEstInterpMethod']] )
print( 'fn_bkg:', fn_bkg )
fn_noi = os.path.join( output_dir, "noise_%s.hdf5"%interp_name[params['NoiseEstInterpMethod']] )
print( 'fn_noise:', fn_noi )

bkg_data_name = [ 'src_noise_bkg', 'bkg_s', 'bkg', 'src_noise' ]
noi_data_name = [ 'src_noise', 'noise_s', 'noise', 'src' ]

data = {}

print( "load " + fn_bkg )
f = h5py.File( fn_bkg, 'r' )
for n in bkg_data_name:
    print( "load %s"%n )
    data[n] = f[n][()]
f.close()

print( "load noise " + fn_noi )
f = h5py.File( fn_noi, 'r' )
for n in noi_data_name:
    print( "load %s"%n )
    data[n] = f[n][()]
f.close()

fig, axs = plt.subplots( 2, 4, figsize=(4*5, 2*5) )
data_name = bkg_data_name + noi_data_name

d = data['src_noise_bkg']
vmin = d.min()
vmax = d.max()

for i in range(len(data_name)):
    print( 'plot %s ...'%data_name[i] )
    ax = axs[ i//4, i%4 ]
    vmin = data[data_name[i]].min()
    vmax = data[data_name[i]].max()
    norm = mplc.SymLogNorm( 1e-9, vmin=vmin, vmax=vmax )
    img = ax.imshow( data[data_name[i]], norm=norm )
    plt.colorbar( img, ax=ax, shrink=0.8  )

    if i < 4:
        ax.set_title( "%s(%s)"%(data_name[i], interp_name[params['BkgEstInterpMethod']]) )
    else:
        ax.set_title( "%s(%s)"%(data_name[i], interp_name[params['NoiseEstInterpMethod']]) )

fig.savefig( os.path.join( output_dir, "bkg_noise.png" ), dpi=300  )
