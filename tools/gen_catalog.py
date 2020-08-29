#!/usr/bin/env python3

import h5py
import sys
import os
import pandas as pd
from astropy.wcs import WCS
import numpy as np
from astropy.io import fits
import tools

param_file  = sys.argv[1]
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
    if "DisableSecondFinder" == ll[0]:
        disable_second_find = int(ll[1])

bname = os.path.basename( fits_file )

fmt = "%-25s : %s"
print( fmt%("fits file", fits_file))
print( fmt%("bname", bname))
print( fmt%("SecondFinder", str(disable_second_find)))
print( fmt%("onlyfof", str(onlyfof)))

output_dir = "%s/%s"%( output_dir, bname )
print( fmt%("output dir", output_dir) )

hdu = fits.open( fits_file )[0]
cdelt1 = hdu.header[ 'CDELT1' ]
cdelt2 = hdu.header[ 'CDELT2' ]

def gen_catalog( fn_prefix ):

    fn = "%s/%s_regs.hdf5"%(output_dir, fn_prefix)
    print( "generate catalog from `%s` ..."%fn )

    f = h5py.File( fn, 'r' )

    Ngroup = f.attrs[ 'Ngroup' ]
    print( "Ngroup: %i"%Ngroup )

    if onlyfof:
        data = { 
            'gidx':[],
        'ridx':[],
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
    else:
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
            c = reg.attrs[ 'center' ]
            data['center-y'].append( c[0] )
            data['center-x'].append( c[1] )
            xyerr = reg.attrs[ 'xyerr' ]
            data['yerr'].append( xyerr[0] )
            data['xerr'].append( xyerr[1] )
            data['peak_flux'].append( reg.attrs[ 'peak_flux' ] )
            data['flux_tot'].append( reg.attrs[ 'flux_tot' ] )

            if not onlyfof:
                data['mean_inner'].append( g.attrs[ "mean_inner" ] )
                data['mean_outer'].append( g.attrs[ "mean_outer" ] )
                data['rms_inner'].append( g.attrs[ "rms_inner" ] )
                data['rms_outer'].append( g.attrs[ "rms_outer" ] )
                data['sigma_inner'].append( g.attrs[ "sigma_inner" ] )
                data['sigma_outer'].append( g.attrs[ "sigma_outer" ] )

            flux = reg[ 'flux' ]
            xy = reg[ 'region' ]

            for k in range( len(flux) ):
                data2[ 'y' ].append( xy[0, k] )
                data2[ 'x' ].append( xy[1, k] )
                data2[ 'gidx' ].append( i )
                data2[ 'ridx' ].append( j )
                data2[ 'flux' ].append( flux[k] )


    print( "Ntot: %i"%Ntot )
    ra_deg, dec_deg, ra, dec = tools.pix2world( data['center-x'], data['center-y'], fits_file )
    data['ra'] = ra
    data['dec'] = dec
    data['ra[deg]'] = ra_deg
    data['dec[deg]'] = dec_deg
    data['ra_err'] = list(np.array( data['xerr'] ) * np.abs(cdelt1))
    data['dec_err'] = list(np.array( data['yerr'] ) * np.abs(cdelt2))

    #for k in data.keys():
    #    print( k, len(data[k]) )

    ra_deg, dec_deg, ra, dec = tools.pix2world( data2['x'], data2['y'], fits_file )
    data2['ra'] = ra
    data2['dec'] = dec
    data2['ra[deg]'] = ra_deg
    data2['dec[deg]'] = dec_deg

    df = pd.DataFrame(data)
    df2 = pd.DataFrame(data2)

    if not os.path.exists( bname ):
        os.mkdir( bname )

    output_file  = '%s/%s_catalog_%s.csv'%(bname,bname, fn_prefix)
    output_file2 = '%s/%s_catalog_%s_point.csv'%(bname,bname, fn_prefix)
    df.to_csv( output_file ,index=False)
    df2.to_csv( output_file2 ,index=False)
    f.close()

if onlyfof:
    gen_catalog( 'only_fof' )
else:
    gen_catalog( 'lset' )
    if not disable_second_find:
        gen_catalog( 'fof' )

