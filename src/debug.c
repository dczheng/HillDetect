#include "allvars.h"
// those functions used for debug
void print_data( double *d, int x0, int x1, int y0, int y1, int flag ) {

    put_sep(0);
    printf( "x0: %i, x1: %i, y0: %i, y1: %i\n", x0, x1, y0, y1 );
   int ti, tj;
    int tindex;
    double tmin, tmax;
    tmin = 1e20;
    tmax = -tmin;
    for( ti=y0; ti<y1; ti++ ) {
        for( tj=x0; tj<x1; tj++ ) {
            tindex = ti * Width + tj;
            printf( "%g ", d[ tindex ] );
            if ( d[tindex] > tmax )
                tmax = d[tindex];
            if ( d[tindex] < tmin )
                 tmin = d[tindex];
        }
        if ( d[tj] )
        printf( "\n" );
    }
    printf( "tmin: %g, tmax: %g\n", tmin, tmax );
    if ( flag )
        endrun("");
    put_sep(0);

}

void output_data( char *fn ) {

    int i, j;
    FILE *fd;

    fd = fopen( fn, "w" );
    for( i=0; i<Height; i++ ) {
        for ( j=0; j<Width; j++ )
            fprintf( fd, "%g ", Data[i*Width+j] );
        fprintf( fd, "\n" );
    }
    fclose( fd );

}
