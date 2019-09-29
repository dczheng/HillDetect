#!/home/hkli/anaconda3/bin/python3

from matplotlib import use
use( 'agg' )

import numpy as np
import h5py
import sys
import os
from astropy.io import fits
from astropy.wcs import WCS

fits_file = sys.argv[1]
lset_file = sys.argv[2]
fn_osm = sys.argv[3]

mywcs = WCS( fits_file )
#hdu = fits.open( fits_file )[0]
#fits_header = hdu.header
#crpix1 = fits_header['CRPIX1']
#crpix2 = fits_header['CRPIX2']
#crval1 = fits_header['CRVAL1']
#crval2 = fits_header['CRVAL2']
#cdelt1 = fits_header['CDELT1']
#cdelt2 = fits_header['CDELT2']

f = h5py.File( lset_file, 'r' )

NGroups = f.attrs[ 'NGroups' ]
Gcrpix = f.attrs[ 'GCRPIX' ] 
Gnaxis = f.attrs[ 'GNAXIS' ] 
print( "GCRPIX: %i %i"%(Gcrpix[0], Gcrpix[1]) )
print( "Gnaxis: %i %i"%(Gnaxis[0], Gnaxis[1]) )
print( NGroups )

def mypix2world( rai, deci ):
    t = [ [ rai[i], deci[i], 0, 0 ] for i in range( len(rai) ) ]
    r = mywcs.wcs_pix2world( t, 0 )
    return (r[:,0], r[:,1])

ddd = []
for gidx in range( NGroups ):
    g = f[ '/Group%i'%gidx ]
    crpixy, crpixx = g.attrs[ 'CRPIX' ]
    Nregs = g.attrs[ 'NRegs' ]

    for j in range( Nregs ):
        pix = g[ 'Reg-%i'%j ][()]
        flux = g[ 'Flux-%i'%j ][()]
        #c = g.attrs[ 'Center-%i'%j ]
        #print( pix )
        #print(flux)
        #print(c)
        deci = pix[0,:] + crpixy + Gcrpix[0]
        rai = pix[0,:] + crpixx + Gcrpix[1]
        #print(len(deci))
        ra, dec = mypix2world(rai, deci)
        #print(ra, dec)
        for k in range( len(ra) ):
            ddd.append( ( ra[k], dec[k], flux[k] ) )

#print( ddd )

f_osm = open( fn_osm, "w" )
f_osm.write( "# Frequency = 158.000 [MHz]\n" )
f_osm.write( "# Pixel size = 2.00 [arcsec]\n" )
f_osm.write( "# K2JyPixel = 7.211e-8\n" )
f_osm.write( "# RA0 = 0.0000 [deg]\n" )
f_osm.write( "# Dec0 = -27.0000 [deg]\n" )
f_osm.write( "# Minimum value = 1.0000e-04 [K]\n" )
f_osm.write( "# Maximum value = inf [K]\n" )
f_osm.write( "# Source count = %i (0.0%%)\n"%(len(ddd)) )
f_osm.write( "#\n" )
fmt = "%-15s%-15s%-15s\n"
f_osm.write(fmt%( "# R.A[deg]", "Dec.[deg]", "flux[Jy]" ))
fmt = "%e,%e,%e\n"
for l in ddd:
    f_osm.write( fmt%l )

f_osm.close()
