#!/usr/bin/env python3

import h5py
import sys
import os
import pandas as pd
from astropy.wcs import WCS
import numpy as np
from astropy.io import fits

param_file  = sys.argv[1]
ls = open( param_file ).readlines()

for l in ls:
    if "FileName" in l and l[0] != '%':
        fits_file = l.split()[1]
    if "OutputDir" in l and l[0] != '%':
        output_dir = l.split()[1]
    if "DisableSecondFinder" in l and l[0] != '%':
        disable_second_find = int(l.split()[1])

bname = os.path.basename( fits_file )

fmt = "%-25s : %s"
print( fmt%("fits file", fits_file))
print( fmt%("bname", bname))
print( fmt%("SecondFinder", str(disable_second_find)))

output_dir = "%s/%s"%( output_dir, bname )
print( fmt%("output dir", output_dir) )

hdu = fits.open( fits_file )[0]
cdelt1 = hdu.header[ 'CDELT1' ]
cdelt2 = hdu.header[ 'CDELT2' ]


def pix2world( x, y, fits_file ):
    
    mywcs = WCS( fits_file )

    ras = []
    decs = []
    ras_deg = []
    decs_deg = []
    for i in range( len(x) ):
        rai = x[i]
        deci = y[i]
        r = mywcs.wcs_pix2world( [ [rai, deci, 0, 0], [0, 0, 0, 0] ], 0 )[0]
    
        ra1 = r[0]
        dec1 = r[1]
        if ra1 < 0:
            ra1 += 360
        ra = ra1
        h = int(ra/(15.0))
        ra = ra - h * 15 
        m = int( ra/(15)*60 )
        s = (ra - m*15/60) / 360 * 24 * 60 * 60
        ra2 = "%ih%im%.10es"%( h, m, s )
    
        dec = dec1
        sign = dec / np.abs(dec)
        dec = np.abs(dec)
        d = int(dec) * sign
        dec = (dec - int(dec)) * 60
        m = int(dec)
        s = (dec - int(dec)) * 60
        dec2 = "%id%im%.10es"%( d, m, s )
        ras_deg.append( ra1 )
        decs_deg.append( dec1 )
        ras.append( ra2 )
        decs.append( dec2 )
    return ras_deg, decs_deg, ras, decs

def gen_catalog( fn_prefix ):

    fn = "%s/%s_regs.hdf5"%(output_dir, fn_prefix)
    print( "generate catalog from `%s` ..."%fn )

    f = h5py.File( fn, 'r' )

    Ngroup = f.attrs[ 'Ngroup' ]
    print( "Ngroup: %i"%Ngroup )

    data = { 
        'gidx':[],
        'ridx':[],
        'mean_inner': [],
        'mean_outer': [],
        'sigma_inner': [],
        'sigma_outer': [],
        'rms_inner': [],
        'rms_outer': [],
        'center-x':[],
        'center-y':[],
        'flux_tot':[],
        'peak_flux':[],
        'xerr':[],
        'yerr': [],
        'ra_err':[],
        'dec_err': [],
        'ra': [],
        'dec': [],
        'ra[deg]': [],
        'dec[deg]': []
    }

    data2 = { 
        'gidx':[],
        'ridx':[],
        'x':[],
        'y': [],
        'flux': [],
        'ra':[],
        'dec':[],
        'ra[deg]':[],
        'dec[deg]':[]
    }

    Ntot = 0
    for i in range(Ngroup):
        g = f[ 'Group%i'%i ]
        Nreg = g.attrs[ 'NReg' ]
        Ntot = Ntot + Nreg
        if Nreg < 0:
            print( "Error" )
            exit()
        for j in range( Nreg ):
            reg = g[ 'Reg%i'%j ]
            data['gidx'].append( i )
            data['ridx'].append( j )
            data['mean_inner'].append( g.attrs[ "mean_inner" ] )
            data['mean_outer'].append( g.attrs[ "mean_outer" ] )
            data['rms_inner'].append( g.attrs[ "rms_inner" ] )
            data['rms_outer'].append( g.attrs[ "rms_outer" ] )
            data['sigma_inner'].append( g.attrs[ "sigma_inner" ] )
            data['sigma_outer'].append( g.attrs[ "sigma_outer" ] )
            c = reg.attrs[ 'center' ]
            data['center-y'].append( c[0] )
            data['center-x'].append( c[1] )
            xyerr = reg.attrs[ 'xyerr' ]
            data['yerr'].append( xyerr[0] )
            data['xerr'].append( xyerr[1] )
            data['peak_flux'].append( reg.attrs[ 'peak_flux' ] )
            data['flux_tot'].append( reg.attrs[ 'flux_tot' ] )

            flux = reg[ 'flux' ]
            xy = reg[ 'region' ]

            for k in range( len(flux) ):
                data2[ 'y' ].append( xy[0, k] )
                data2[ 'x' ].append( xy[1, k] )
                data2[ 'gidx' ].append( i )
                data2[ 'ridx' ].append( j )
                data2[ 'flux' ].append( flux[k] )


    print( "Ntot: %i"%Ntot )
    ra_deg, dec_deg, ra, dec = pix2world( data['center-x'], data['center-y'], fits_file )
    data['ra'] = ra
    data['dec'] = dec
    data['ra[deg]'] = ra_deg
    data['dec[deg]'] = dec_deg
    data['ra_err'] = list(np.array( data['xerr'] ) * np.abs(cdelt1))
    data['dec_err'] = list(np.array( data['yerr'] ) * np.abs(cdelt2))

    #for k in data.keys():
    #    print( k, len(data[k]) )

    ra_deg, dec_deg, ra, dec = pix2world( data2['x'], data2['y'], fits_file )
    data2['ra'] = ra
    data2['dec'] = dec
    data2['ra[deg]'] = ra_deg
    data2['dec[deg]'] = dec_deg

    df = pd.DataFrame(data)
    df2 = pd.DataFrame(data2)

    if not os.path.exists( bname ):
        os.mkdir( bname )

    output_file  = '%s/%s_catalog_%s.csv'%(bname,bname, fn_prefix)
    output_file2 = '%s/%s_catalog_%s_flux.csv'%(bname,bname, fn_prefix)
    df.to_csv( output_file ,index=False)
    df2.to_csv( output_file2 ,index=False)
    f.close()

gen_catalog( 'lset' )
if not disable_second_find:
    gen_catalog( 'fof' )
