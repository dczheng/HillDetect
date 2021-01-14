/*
    dczheng
    created 2019-11-06
*/

#include "allvars.h"

int *fof_map;

void get_mean_sigma_rms( double *data, int *flag, int N,
                     double *mean, double *sigma, double *rms ) {
    int i, flags, *flag_n;

    if ( NULL == flag ) {
        *mean = *sigma = *rms = 0;
        for( i=0; i<N; i++ ) {
            *mean += data[i];
            *rms  += SQR( data[i] );
        }

        *mean /= N;
        *rms = sqrt( *rms / N ); 

        for( i=0; i<N; i++ )
            *sigma += SQR( data[i] - (*mean) ) ;
        *sigma = sqrt( *sigma / N );
        return;
    }

    //put_start;
    flag_n = malloc( sizeof(int) * N );
    memset( flag_n, 0, sizeof(int) * N );

    flags = 0;
    for( i=0; i<N; i++ ) {
        flag_n[flag[i]]++;
        if ( flag[i] > flags )
            flags = flag[i];
    }
    flags ++;

    for(i=0; i<flags; i++)
        mean[i] = sigma[i] = rms[i] = 0;

    for( i=0; i<N; i++ ) {
        mean[flag[i]] += data[i];
        rms[flag[i]] += SQR(data[i]);
    }

    for( i=0; i<flags; i++ )
        if ( flag_n[i] ) {
            mean[i] /= flag_n[i];
            rms[i] = sqrt(rms[i]/flag_n[i]);
        }

    for( i=0; i<N; i++ )
        sigma[flag[i]] += SQR( data[i] - mean[flag[i]] ) ;

    for( i=0; i<flags; i++ ) 
        if ( flag_n[i] )
            sigma[i] = sqrt( sigma[i]/flag_n[i] );

    free( flag_n );

}

void find_vmin_vmax( double *d, int N, double *vmin, double *vmax ) {
    int i;
    *vmin = DBL_MAX;
    *vmax = -DBL_MAX;
    for( i=0; i<N; i++ ) {
        *vmin = ( d[i] < *vmin ) ? d[i] : *vmin;
        *vmax = ( d[i] > *vmax ) ? d[i] : *vmax;
    }
}

double calc_mean( double *data, int N ) {

    int i;
    double m;
    m = 0;
    for( i=0; i<N; i++ )
        m += data[i];
    return m / N;
}

double calc_sigma( double * data, int N, double *mean ) {
    int i;
    double m, s;
    if ( NULL == mean )
        m = calc_mean( data, N );
    else 
        m =  *mean;
    s = 0;
    for( i=0; i<N; i++ ) {
        s += SQR( data[i] - m ) ;
    }
    s = sqrt( s / N );
    return s;
}

double sigma_clipping( double *data, int N, double R ) { 
    int i;
    double m, s, vmin, vmax, v1, v2, vinv;
    find_vmin_vmax( data, N, &vmin, &vmax );
    m = calc_mean( data, N );
    s = calc_sigma( data, N, &m );
    v1 = m - R*s;
    v2 = m + R*s;
    vinv = vmin / 10;

    //printf( "mean: %f, sigma: %f\n", m, s );
    for( i=0; i<N; i++ )
        if ( data[i] < v1 || data[i] > v2 )
            data[i] = vinv;
    return vinv;
}

double *poly_2d_fitting( double *x, double *y, double *z, int nz, int M, int N, int pn) {

    int an, i, j, ii, jj, k, m, n, sign;
    gsl_matrix *XTX, *XTX_inv, *tmp;
    double *xs, *ys, s, t, a, b, *xt, *xty, *as, *fitting, *data, vmin, vmax, z1, z2, b_a, z2_z1, vmax_vmin;
    gsl_permutation *pm;

    if ( nz == 0 ) {
       endrun("error parameters for fitting\n" );
    }

    for( i=0; i<nz; i++ )
        if ( x[i] > N-1 || y[i] > M-1 || x[i] < 0 || y[i] < 0 ) {
            printf( "%f %f %i %i\n", x[i], y[i], M, N );
            endrun("error parameters for fitting\n" );
        }

    an = ( pn + 2 ) * ( pn + 1 ) / 2;

    fitting = malloc( sizeof(double) * M * N );
    data = malloc( sizeof(double) * nz );
    xt = malloc( sizeof(double) * nz * an );
    xty = malloc( sizeof(double) * an );
    as = malloc( sizeof(double) * an );

    XTX = gsl_matrix_alloc( an, an );
    XTX_inv = gsl_matrix_alloc( an, an );
    tmp = gsl_matrix_alloc( an, an );
    pm = gsl_permutation_alloc(an);

    xs = (double*) malloc( sizeof( double ) * (pn+1) );
    ys = (double*) malloc( sizeof( double ) * (pn+1) );

    a = 1;
    b = 2;
    b_a = b-a;
    z1 = 1;
    z2 = 2;
    z2_z1 = z2-z1;

    find_vmin_vmax( z, M*N, &vmin, &vmax );
    vmax_vmin = vmax - vmin;

    for( i=0; i<nz; i++ ) {
        data[i] = (z[i] - vmin) / vmax_vmin * z2_z1 + z1;
    }

    xs[0] = ys[0] = 1;
    for( i=0; i<nz; i++ ) {

        t = y[i] / M * b_a + a;
        for( k=1; k<pn+1; k++ )
            ys[k] = ys[k-1] * t;

        t = x[i] / N * b_a + a;
        for( k=1; k<pn+1; k++ )
            xs[k] = xs[k-1] * t;

        n = 0;
        for( ii=0; ii<pn+1; ii++ )
            for( jj=0; jj<=ii; jj++ ) {
                xt[ n * nz + i ] = xs[jj] * ys[ii-jj];
                n++;
            }

    }

    for( i=0; i<an; i++ )
        for( j=0; j<an; j++ ){
            s = 0;
            for( k=0; k<nz; k++ )
                s += xt[ i * nz + k ] * xt[ j * nz + k ];
            gsl_matrix_set( XTX, i, j, s );
        }

    gsl_matrix_memcpy(tmp, XTX);
    sign = 0;
    gsl_linalg_LU_decomp(tmp, pm, &sign);
    gsl_linalg_LU_invert(tmp, pm, XTX_inv);

    for( i=0; i<an; i++ ) {
        s = 0;
        for( j=0; j<nz; j++ )
            s += xt[ i * nz + j ] * data[j];
        xty[i] = s;
    }

    for( i=0; i<an; i++ ) {
        s = 0;
        for( j=0; j<an; j++ ) 
            s += gsl_matrix_get( XTX_inv, i, j ) * xty[j];
        as[i] = s;
    }

    xs[0] = ys[0] = 1;
    for( i=0; i<M; i++ ) {

        t = ((double)i) / M * b_a + a;
        for( k=1; k<pn+1; k++ )
            ys[k] = ys[k-1] * t;

        for( j=0; j<N; j++ ) {

            t = ((double)j) / N * b_a + a;
            for( k=1; k<pn+1; k++ )
                xs[k] = xs[k-1] * t;

            m = i * N + j;
            s = 0;
            n = 0;
            for( ii=0; ii<pn+1; ii++ )
                for( jj=0; jj<=ii; jj++ ) {
                    s += xs[jj] * ys[ii-jj] * as[n];
                    n++;
                }
            fitting[m] = s;
        }
    }

    free( xt );
    free( xty );
    free( as );
    gsl_matrix_free( XTX );
    gsl_matrix_free( XTX_inv );
    gsl_matrix_free(tmp);
    gsl_permutation_free(pm);
    free(data);

    for( i=0; i<M*N; i++ ) {
        fitting[i] = (fitting[i] - z1) / z2_z1 * vmax_vmin + vmin;
    }
    return fitting;

}
