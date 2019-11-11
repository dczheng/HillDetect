#!/usr/bin/env python3

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from astropy.coordinates import Angle
import astropy.units as au

def pix2world( x, y, fits_file ):
    

    mywcs = WCS( fits_file )
    r = mywcs.wcs_pix2world(\
            np.array([ [x[i], y[i]] for i in range(len(x)) ]),\
            0 )

    ra_deg = r[:,0]
    dec_deg = r[:,1]

    t = Angle( ra_deg * au.deg ).hms
    ra = []
    for i in range(len(t[0])):
        ra.append( "%ih%im%g"%(t[0][i], t[1][i], t[2][i]) )

    dec = Angle( dec_deg * au.deg )

    return ra_deg, dec_deg, ra, dec
