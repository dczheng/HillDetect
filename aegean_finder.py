#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )

import sys
import os
import logging
from astropy.io import fits
from astropy.wcs import WCS
from AegeanTools.fits_image import Beam
from AegeanTools.catalogs import save_catalog
from AegeanTools.source_finder import SourceFinder
from AegeanTools.AeRes import make_model
from AegeanTools.wcs_helpers import WCSHelper

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mplc

import tools
import pandas as pd

import time

param_file = sys.argv[1]

ls = open( param_file ).readlines()

for l in ls:
    ll = l.split()
    if ll == []:
        continue

    if "FileName" == ll[0]:
        fits_file = ll[1]

    if "OutputDir" == ll[0]:
        OutputDir = ll[1]

    #if "Beam" == ll[0]:
    #    beam = float(ll[1])

    if "InnerClip" == ll[0]:
        innerclip = float(ll[1])

    if "OuterClip" == ll[0]:
        outerclip = float(ll[1])

    if "Cores" == ll[0]:
        cores = int(ll[1])

    if "CuttingXStart" == ll[0]:
        cutx0 = float( ll[1] )

    if "CuttingXEnd" == ll[0]:
        cutx1 = float( ll[1] )

    if "CuttingYStart" == ll[0]:
        cuty0 = float( ll[1] )

    if "CuttingYEnd" == ll[0]:
        cuty1 = float( ll[1] )

    if "FluxMin" == ll[0]:
        fluxmin = float( ll[1] )

bname = os.path.basename( fits_file )
aegean_outputdir = "%s/%s_aegean"%( OutputDir, bname )
catalog_file = "%s/catalog.csv"%aegean_outputdir
catalog_file_point = "%s/catalog_point.csv"%aegean_outputdir
plot_file = "%s/aegean.png"%aegean_outputdir

if not os.path.exists( aegean_outputdir ):
    os.makedirs( aegean_outputdir )

fmt = '%-25s : %s'
print( fmt%( "fits file", fits_file ) )
print( fmt%( "base name", bname ) )
print( fmt%( "Output dir", OutputDir ) )
#print( fmt%( "Beam", str(beam) ) )
print( fmt%( "Aegean OutputDir", aegean_outputdir ) )
print( fmt%( "InnerClip", str(innerclip) ) )
print( fmt%( "OuterClip", str(outerclip) ) )
print( fmt%( "Cores", str(cores) ) )
print( fmt%( "catalog file", catalog_file ) )
print( fmt%( "catalog file [point]", catalog_file_point ) )
print( fmt%( "plot file", plot_file ) )
print( fmt%( "cutx0", cutx0 ) )
print( fmt%( "cutx1", cutx1 ) )
print( fmt%( "cuty0", cuty0 ) )
print( fmt%( "cuty1", cuty1 ) )
print( fmt%( "fluxmin", fluxmin ) )

hdu = fits.open( fits_file )[0]
fits_data = hdu.data
fits_header = hdu.header

mywcs_helper = WCSHelper.from_file( fits_file )

log = logging.getLogger( bname )
sf = SourceFinder( log=log )
found = sf.find_sources_in_image( filename = fits_file,\
                                  innerclip = innerclip,\
                                  outerclip = outerclip,\
                                  cores = 5,\
                                  nonegative=True \
                                  )
save_catalog( catalog_file, found )

shape = fits_data.shape
m, n = shape
m0 = int(m * cuty0)
m1 = int(m * cuty1)
n0 = int(n * cutx0)
n1 = int(n * cutx1)

img =  make_model( sf.sources, shape, mywcs_helper )

fig, axs = plt.subplots( 2, 2, figsize=(10,10) )
axs[0,0].imshow( fits_data, norm=mplc.LogNorm(), cmap=cm.jet )
#axs[0,1].imshow( fits_data, norm=mplc.LogNorm(), cmap=cm.jet )
axs[0,1].imshow( img, norm=mplc.LogNorm(), cmap=cm.jet )
axs[1,0].imshow( fits_data[m0:m1, n0:n1],\
                norm=mplc.LogNorm(), cmap=cm.jet )
axs[1,1].imshow( fits_data[m0:m1, n0:n1],\
                norm=mplc.LogNorm(), cmap=cm.jet )
axs[1,1].imshow( img[m0:m1, n0:n1],\
                norm=mplc.LogNorm(), cmap=cm.jet )
fig.savefig( plot_file )

data = { \
       'x': [],\
       'y': [],\
       'flux': [],\
       'ra':[],\
       'dec':[],\
       'ra[deg]':[],\
       'dec[deg]':[]\
       }

for i in range(m):
    for j in range(n):
        if img[i, j] > fluxmin:
            data[ 'y' ].append(i)
            data[ 'x' ].append(j)
            data['flux'].append(img[i,j])

ras_deg, decs_deg, ras, decs = tools.pix2world(\
                data['x'], data['y'], fits_file )

data['ra'] = ras
data['dec'] = decs
data['ra[deg]'] = ras_deg
data['dec[deg]'] = decs_deg

df = pd.DataFrame(data)
df.to_csv( catalog_file_point, index=False )
