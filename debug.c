#include "allvars.h"
void print_data( double *d, int x0, int x1, int y0, int y1, int flag ) {
    put_sep;
    printf( "x0: %i, x1: %i, y0: %i, y1: %i\n", x0, x1, y0, y1 );
    int ti, tj;
    long tindex;
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
        endrun( 1212121 );
    put_sep;

}
