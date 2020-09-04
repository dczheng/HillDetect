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

void find_vmin_vmax( double *d, int N, double *up, double *down, double *vmin, double *vmax ) {
    int i;
    *vmin = DBL_MAX;
    *vmax = -DBL_MAX;
    for( i=0; i<N; i++ ) {
        if ( NULL != up && d[i] > *up ) continue;
        if ( NULL != down && d[i] < *down ) continue;
        *vmin = ( d[i] < *vmin ) ? d[i] : *vmin;
        *vmax = ( d[i] > *vmax ) ? d[i] : *vmax;
    }
}

void get_mean_sigma( double *data, int N, double *skip, double *mean, double *sigma ) {
    int i, n;
    n = *mean = *sigma = 0;
    for( i=0; i<N; i++ ) {
        if ( skip != NULL && data[i] == *skip ) continue;
        *mean += data[i];
        n++;
    }
    *mean /= n;

    for( i=0; i<N; i++ ) {
        if ( skip != NULL && data[i] == *skip ) continue;
        *sigma += SQR( data[i] - *mean ) ;
    }
    *sigma = sqrt( *sigma / n );
}
