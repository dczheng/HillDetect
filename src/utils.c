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

int get_mean_sigma_comp( const void *a, const void *b ) {
    if ( *((double*)a) > *((double*)b) )
        return -1;
    else
        return 1;
}

void get_mean_sigma( double *data, int N, double *skip, double *mean, double *sigma, double *median, double *buf ) {
    int i, n;
    n = *mean = *sigma = 0;
    *median = *skip;
    for( i=0; i<N; i++ ) {
        if ( skip != NULL && data[i] == *skip ) continue;
        *mean += data[i];
        if ( NULL != buf )
            buf[n] = data[i];
        n++;
    }

    if ( n == 0 ) {
        *mean = *skip;
        *sigma = 0;
        return;
    }

    *mean /= n;

    //printf( "[1]" );
    //for( i=0; i<n; i++ )
    //    printf( "%g ", buf[i] );
    //printf( "\n[2]" );
    if ( NULL != buf ) {
        qsort(buf, n, sizeof(double), &get_mean_sigma_comp);
        *median = buf[n/2];
    }

    //for( i=0; i<n; i++ ) {
    //    if ( i > 0 ) {
    //        if ( buf[i] >= buf[i-1] )
    //            printf( "1" );
    //        else
    //            printf( "0" );
    //    }
    //    //printf( "%g ", buf[i] );
    //}
    //printf( "\n%g\n", *median );
    //
    /*
    endrun("1");
     */

    for( i=0; i<N; i++ ) {
        if ( skip != NULL && data[i] == *skip ) continue;
        *sigma += SQR( data[i] - *mean ) ;
    }
    *sigma = sqrt( *sigma / n );
}
