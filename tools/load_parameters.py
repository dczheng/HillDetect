#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )


def load_parameters( fn ):
    params = {}
    lines = open( fn ).readlines()
    for l in lines:
        ll = l.split()
        if ll == []:
            continue
        params[ll[0]] = ll[1]
    return params

