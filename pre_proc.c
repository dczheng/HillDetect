/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

void sigma_clipping() {

    long i;
    double sigma, mu, vmin;

    vmin = 1e100;
    for ( i=0; i<Npixs; i++ ) {
        if ( isnan( Data[i] ) )
            continue;
        vmin = ( vmin < Data[i] ) ? vmin : Data[i];
    }

    for( i=0, mu=0; i<Npixs; i++ ) {
        if ( isnan( Data[i] ) )
            Data[i] = vmin;

        mu += Data[i];
    }

    mu /= Npixs;

    for( i=0, sigma=0; i<Npixs; i++ ) {
        sigma += SQR( Data[i] - mu );
    }
    sigma = sqrt( sigma / Npixs );
    vmin = mu + All.RSigma * sigma;  

    for( i=0; i<Npixs; i++ ) {
        if ( Data[i] < vmin )
            Data[i] = vmin; 
    }

}

void normalize() {

    double vmin, vmax, dv;
    long i;

    vmax = -1e100;
    vmin = 1e100;

    for( i=0; i<Npixs; i++ ) {

        if ( isnan( Data[i] ) )
            continue;

        if ( All.LogNorm ) {
            if ( Data[i] > 0 )
                vmin = ( vmin < Data[i] ) ? vmin : Data[i];
        }
        else {
            vmin = ( vmin < Data[i] ) ? vmin : Data[i];
        } 

        vmax = ( vmax > Data[i] ) ? vmax : Data[i];
    }


    for( i=0; i<Npixs; i++ ) {

        if ( isnan( Data[i] ) )
            Data[i] = vmin;

        if ( All.LogNorm )
            if ( Data[i] < vmin )
                Data[i] = vmin;
    }

    if ( All.LogNorm )
        dv =log( vmax / vmin );
    else
        dv = vmax - vmin;

    for ( i=0; i<Npixs; i++ ) {
        if ( All.LogNorm )
            Data[i] = log( Data[i]/vmin ) / dv;
        else
            Data[i] = ( Data[i] - vmin ) / dv;
    }

}
