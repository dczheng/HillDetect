#!/home/hkli/anaconda3/bin/python3
#!/usr/bin/env python3

import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from astropy.io import fits
import sys

fits_name = sys.argv[1]
fits_dir  = sys.argv[2]
lines_dir = sys.argv[3]

lines = open( lines_dir + '/' + fits_name + '.out').readlines()
fits_data = fits.open( fits_dir + '/' + fits_name )[0].data
png_dir = './pngs/'

N = len(lines)
dn = 10
i = 0
while( i<N-1 ):
    x = [ int(v) for v in lines[i].split() ]
    y = [ int(v) for v in lines[i+1].split() ]
    ii = x[0]
    if x[1] == 0:
        continue
    print( 'iter: %i, n: %i'%(x[0], x[1]) )

    x = x[2:]
    y = y[2:]
    plt.imshow( fits_data, norm=mplc.LogNorm() )
    ax = plt.gca()
    ax.invert_yaxis()
    plt.plot( x, y, '.' )

    plt.savefig( png_dir + '/%s_%05i.png'%(fits_name, ii) )
    plt.close()

    i += 2*dn
