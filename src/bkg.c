/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

double *bkg_fitting( double *data, int M, int N, void *params ) {

    int m, n, pn, w, h, i, j, ii, jj, w0, h0, iis, iie, jjs, jje, sig_idx, sig, pad, idx;
    double *z, *x, *y, *bkg_local, *bkg, t, R, *p;
    p = params;
    pn = p[0];
    m = p[1];
    n = p[2];
    pad = p[3];
    R = p[4];
    w0 = N / n;
    h0 = M / m;

    if ( M % m != 0 || N % n != 0  )
        endrun( "error parameters\n" );

    z = malloc( sizeof( double ) * (w0+2*pad) * (h0+2*pad) );
    x = malloc( sizeof( double ) * (w0+2*pad) * (h0+2*pad) );
    y = malloc( sizeof( double ) * (w0+2*pad) * (h0+2*pad) );
    bkg = malloc( sizeof(double)  * M * N);
    for( i=0; i<M*N; i++ )
        bkg[i] = 0;

    sig = m*n / 10;
    sig_idx = 0;

    for( i=0; i<m; i++ ) {

        iis = i == 0  ? 0 : -pad;
        iie = i == m-1 ? h0 : h0+pad; 
        h = iie - iis;

        for( j=0; j<n; j++ ) {

            if ( sig_idx % sig == 0 )
                printf( "[bkg fitting] %3.0f%%\n", ((double)sig_idx)/(m*n) * 100 );
            sig_idx ++;

            jjs = j == 0 ? 0 : -pad;
            jje = j == n-1 ? w0 : w0+pad;
            w = jje - jjs;

            for( ii=0; ii<h; ii++ ) {
                for( jj=0; jj<w; jj++ ) {
                    t = data[ (i*h0+ii+iis) * N + (j*w0+jj+jjs) ];
                    idx = ii*w+jj;
                    z[idx] = t;
                    y[idx] = ii;
                    x[idx] = jj;
                }
            }

            t = sigma_clipping( z, h*w, R );
            idx = 0;
            for( ii=0; ii<h*w; ii++  ) {
                if ( z[ii] == t )
                    continue;
                x[idx] = x[ii];
                y[idx] = y[ii];
                z[idx] = z[ii];
                idx++;
            }

            if ( ((double)idx / (h*w)) < 0.3 ) {
                printf( "WARNING [SigmaClipping %.2f] Nfitting[%i] / Ntot[%i] < 0.3\n", R, idx, h*w );
            }

            bkg_local = poly_2d_fitting( x, y, z, idx, h, w, pn );
            for( ii=0; ii<h0; ii++ ) {
                for( jj=0; jj<w0; jj++ ) {
                    bkg[ (i*h0+ii) * N + (j*w0 + jj) ] = bkg_local[(ii-iis)*w+jj-jjs];
                }
            }
            free( bkg_local );
        }
    }

    free(z);
    free(x);
    free(y);

    return bkg;

}

double * background( double *data, int M, int N, int mode, void *params ) {

    if ( mode == 0 ) {
        return bkg_fitting( data, M, N, params );
    }
    return NULL;
}
