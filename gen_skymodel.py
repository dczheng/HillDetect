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
outdir = sys.argv[2]

bname = os.path.basename( fits_file )
print( bname )

mywcs = WCS( fits_file )
hdu = fits.open( fits_file )[0]
fits_header = hdu.header
print( hdu.data.shape )

crpix1 = fits_header['CRPIX1']
crpix2 = fits_header['CRPIX2']
crval1 = fits_header['CRVAL1']
crval2 = fits_header['CRVAL2']
cdelt1 = fits_header['CDELT1']
cdelt2 = fits_header['CDELT2']
freq   = fits_header['CRVAL3'] / 1e6
print( "crpix: %i %i"%(crpix1, crpix2) )
print( "crval: %g %g"%(crval1, crval2) )
print( "cdelt: %g %g"%(cdelt1, cdelt2) )

f = h5py.File( "%s/%s_fof_regs.hdf5"%( outdir, bname ), 'r' )
NGroups = f.attrs[ 'NGroups' ]
Gcrpix = f.attrs[ 'GCRPIX' ] 
Gnaxis = f.attrs[ 'GNAXIS' ] 
print( "GCRPIX: %i %i"%(Gcrpix[0], Gcrpix[1]) )
print( "Gnaxis: %i %i"%(Gnaxis[0], Gnaxis[1]) )
print( NGroups )


index = 0
ddd = []
def skymodel_write_head( fd, name ):
    fd.write( "source {\n" )
    fd.write( "  name \"%s\"\n"%name )
def skymodel_write_tail( fd ):
    fd.write( "}\n" )

def skymodel_write_comp( fd, ra, dec, freq, flux ):
    fd.write( "    component {\n" )
    fd.write( "      type point\n" )
    fd.write( "      position %s %s\n"%(ra, dec) )
    fd.write( "      measurement {\n" )
    fd.write( "        frequency %g MHz\n"%freq )
    fd.write( "        fluxdensity Jy %g 0 0 0 \n"%flux )
    fd.write( "    }\n" )
    fd.write( "  }\n" )

f_skymdel = open( "skymodel_regs.txt", "w" )
f_skymdel.write( "skymodel fileformat 1.1\n" )

idx = 0
for gidx in range( NGroups ):
    g = f[ '/Group%i'%gidx ]
    crpixy, crpixx = g.attrs[ 'CRPIX' ]
    Nregs = g.attrs[ 'NRegs' ]
    #print( "CRPIX: %i %i"%( crpixy, crpixx ) )
    #print( "Nregs: %i"%Nregs )


    for j in range( Nregs ):

        c = g['Reg-%i'%j ][()]
        flux = g[ 'Flux-%i'%j ][()]
        skymodel_write_head( f_skymdel, "Group-%i-Reg-%i-%i"%(gidx, j, idx) )

        for k in range(len(flux)):
            deci = c[0, k] + crpixy + Gcrpix[0]
            rai =   c[1, k] + crpixx + Gcrpix[1]
            r = mywcs.wcs_pix2world( [[rai, deci, 0, 0], [0, 0, 0, 0]], 0 )[0]
            ra1 = r[0]
            dec1 = r[1]
            if ra1 < 0:
                ra1 += 360
            ra = ra1
            h = int(ra/(15.0))
            ra = ra - h * 15 
            m = int( ra/(15)*60 )
            s = (ra - m*15/60) / 360 * 24 * 60 * 60
            ra2 = "%ih%im%is"%( h, m, s )
    
            dec = dec1
            sign = dec / np.abs(dec)
            dec = np.abs(dec)
            d = int(dec) * sign
            dec = (dec - int(dec)) * 60
            m = int(dec)
            s = (dec - int(dec)) * 60
            dec2 = "%id%im%is"%( d, m, s )
            skymodel_write_comp( f_skymdel, ra2, dec2, freq, flux[k] )

        skymodel_write_tail( f_skymdel )
        idx += 1

f_skymdel.close()

