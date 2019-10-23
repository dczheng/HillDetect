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


fout = open( 'catalog.txt', 'w' )
fout_csv = open( 'catalog.csv', 'w' )
fmt = "%6s %26s %26s %16s %6s %6s %16s %16s\n"
fmt_csv = "%s,%s,%s,%s,%s,%s,%s,%s\n"
fout.write( fmt%( 'index', 'ra', 'dec', 'flux', 'rai', 'deci', 'ra[deg]', 'dec[deg]' ) )
fout_csv.write( fmt_csv%( 'index', 'ra', 'dec', 'flux', 'rai', 'deci', 'ra[deg]', 'dec[deg]' ) )
fmt = "%6i %26s %26s %.10e %6i %6i %.10e %.10e\n"
fmt_csv = "%i,%s,%s,%.10e,%i,%i,%.10e,%.10e\n"

index = 0
ddd = []
for gidx in range( NGroups ):
    g = f[ '/Group%i'%gidx ]
    crpixy, crpixx = g.attrs[ 'CRPIX' ]
    Nregs = g.attrs[ 'NRegs' ]
    #print( "CRPIX: %i %i"%( crpixy, crpixx ) )
    #print( "Nregs: %i"%Nregs )

    for j in range( Nregs ):
        c = g.attrs[ 'Center-%i'%j ]
        flux = g.attrs[ 'FluxTot-%i'%j ]
        deci = c[0] + crpixy + Gcrpix[0]
        rai = c[1] + crpixx + Gcrpix[1]

        r = mywcs.wcs_pix2world( [[rai, deci, 0, 0], [0, 0, 0, 0]], 0 )[0]
        
        ra1 = r[0]
        dec1 = r[1]
        #print( ra1, dec1 )
        #exit()
        #ra1 =  ( rai-crpix1 ) * cdelt1 + crval1
        #dec1 = ( deci-crpix2 ) * cdelt2 + crval2
        if ra1 < 0:
            ra1 += 360
        ra = ra1
        h = int(ra/(15.0))
        ra = ra - h * 15 
        m = int( ra/(15)*60 )
        s = (ra - m*15/60) / 360 * 24 * 60 * 60
        ra2 = "%ih%im%.10fs"%( h, m, s )

        dec = dec1
        sign = dec / np.abs(dec)
        dec = np.abs(dec)
        d = int(dec) * sign
        dec = (dec - int(dec)) * 60
        m = int(dec)
        s = (dec - int(dec)) * 60
        dec2 = "%id%im%.10fs"%( d, m, s )

        ddd.append( (index, ra2, dec2, flux, rai, deci, ra1, dec1) )
        index += 1
for i in range(index-1):
    for j in range(i, index):
        if ddd[i][3] < ddd[j][3]:
            t = ddd[i]
            ddd[i] = ddd[j]
            ddd[j] = t
def skymodel_write( fd, name, ra, dec, freq, flux ):
    fd.write( "source {\n" )
    fd.write( "  name \"%s\"\n"%name )
    fd.write( "    component {\n" )
    fd.write( "      type point\n" )
    fd.write( "      position %s %s\n"%(ra, dec) )
    fd.write( "      measurement {\n" )
    fd.write( "        frequency %g MHz\n"%freq )
    fd.write( "        fluxdensity Jy %g 0 0 0 \n"%flux )
    fd.write( "    }\n" )
    fd.write( "  }\n" )
    fd.write( "}\n" )

f_skymdel = open( "skymodel.txt", "w" )
f_skymdel.write( "skymodel fileformat 1.1\n" )
for l in ddd:
    skymodel_write( f_skymdel, "group%i"%(l[0]), l[1], l[2], freq, l[3] )
for l in ddd:
    fout.write( fmt%l )
    fout_csv.write( fmt_csv%l )

f_skymdel.close()
fout.close()
fout_csv.close()
