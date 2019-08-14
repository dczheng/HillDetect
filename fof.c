#include "allvars.h"

int *fof_map;

//#define FOF_SINGLE_DEBUG
void fof_single_finder( long p0 ) {

    int  i, j, x0, y0;
    long p1, l, s, ss, p_next;


    if ( fof_map[p0] == 0 )
        return;

    x0 = p0 % Width;
    y0 = p0 / Width;
    p_next = -1;

#ifdef FOF_SINGLE_DEBUG
    printf( "start\n" );
    fflush( stdout );
#endif
    for( i=-1; i<2; i++ )
        for( j=-1; j<2; j++ ) {
            p1 = (y0+i) * Width + (x0+j);
#ifdef FOF_SINGLE_DEBUG
            printf( "p0: %li, p1:  %li, NPixs: %li (%i*%i)\n",
            p0, p1, Npixs, Width, Height );
            fflush( stdout );
#endif
            if ( p1 == p0     ||
                 p1 < 0       ||
                 p1 > Npixs-1 ||
                 fof_map[p1]  == 0 )
                continue;

#ifdef FOF_SINGLE_DEBUG
            printf( "1\n" );
            fflush( stdout );
#endif
            /*
            fof_map[p1] = 0;
            */

            if ( Head[p1] == Head[p0] )
                continue;

            if ( Len[Head[p1]] > Len[Head[p0]] ) {
                l = p1;
                s = p0;
             }
             else {
                l = p0;
                s = p1;
             }

#ifdef FOF_SINGLE_DEBUG
            printf( "2\n" );
            fflush( stdout );
#endif
            Next[Tail[Head[l]]] = Head[s];
            Tail[Head[l]] = Tail[Head[s]];
            Len[Head[l]] += Len[Head[s]];
            Len[Head[s]] = 1;

            ss = Head[s];
            do{
                /*
                printf( "p0: %li, p1: %li l: %li, s: %li, ss: %li\n"
                    "Head[l]: %li, Head[s]: %li\n",
                p0, p1, l, s, ss, Head[l], Head[s] );
                */
                Head[ss] = Head[l];
            } while( (ss=Next[ss]) >= 0 );

#ifdef FOF_SINGLE_DEBUG
            printf( "3\n" );
            fflush( stdout );
#endif
            p_next = p1;

    }

#ifdef FOF_SINGLE_DEBUG
    printf( "end\n\n" );
    fflush( stdout );
#endif

    if ( p_next >= 0 )
        fof_single_finder( p_next );

}

void fof_init() {

    long i;
    Next = malloc( sizeof(long) * Npixs );
    Head = malloc( sizeof(long) * Npixs );
    Len = malloc( sizeof(long) * Npixs );
    Tail = malloc( sizeof(long) * Npixs );
    for( i=0; i<Npixs; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i] = 1;
    }
}

void fof_reset() {

    long i;
    for( i=0; i<Npixs; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i] = 1;
    }
}

void fof_free() {

    free( Next );
    free( Head );
    free( Tail );
    free( Len );

}

int fof_compare_len( const void *a, const void *b ) {
    return ( (*(long*)a) < *((long*)b) ) ? 1 : -1;
}

void fof_sort () {

    long *tmp, p;
    tmp = malloc( sizeof(long) * Npixs * 3 );
    for ( p=0; p<Npixs; p++ ) {
        tmp[ 3*p ] = Len[p];
        tmp[ 3*p+1 ] = Head[p];
        tmp[ 3*p+2 ] = Tail[p];
    }
    
    qsort( tmp, Npixs, sizeof(long)*3, fof_compare_len );

    for ( p=0; p<Npixs; p++ ) {
        Len[p] = tmp[ 3*p ];
        Head[p] = tmp[ 3*p+1 ];
        Tail[p] = tmp[ 3*p+2 ];
    }
    free( tmp );
}

void fof_edge() {

    long i;
    fof_map = malloc( sizeof(int) * Npixs );

    for( i=0; i<Npixs; i++ ) {
        fof_map[i] = 0;
    }

    for( i=0; i<edgen; i++ ) {
        fof_map[ edgey[i] * Width + edgex[i] ] = 1;
    }

    for( i=0; i<Npixs; i++ ) {
         fof_single_finder( i );
    }

    free( fof_map );
    fof_sort();
    for( i=0; i<Npixs; i++ )
        if ( Len[i] == 1 )
            break;
    NfofEdge = i;

}

void fof_region() {

    long p;
    int np;
    fof_map = malloc( sizeof(int) * Npixs );

    np = 0;
    for( p=0; p<Npixs; p++ ) {
        fof_map[p] = 0;
        if ( Phi[p] > 0 ) {
           fof_map[p] = 1;
           np ++;
        }

    }

/*
    int i, j;
    FILE *fd;
    fd = fopen( "fof_map.dat", "w" );
    for( i=0; i<Height; i++ ) {
        for( j=0; j<Width; j++ )
            fprintf( fd, "%i ", fof_map[i*Width+j] );
        fprintf( fd, "\n" );
    }
    fclose( fd );

    fd = fopen( "phi_map.dat", "w" );
    for( i=0; i<Height; i++ ) {
        for( j=0; j<Width; j++ )
            fprintf( fd, "%g ", Phi[i*Width+j] );
        fprintf( fd, "\n" );
    }
    fclose( fd );
*/


    printf( "NPixs: %li, np; %i, nm: %li\n", Npixs, np, Npixs - np );

    //long sigs = Npixs / 100;
    for( p=0; p<Npixs; p++ ) {
        //if ( sigs != 0 && p % sigs == 0 )
        //    printf( "%.2f%%\n", ((double)p)/Npixs * 100 );
        //if ( ((double)i) / Npixs > 0.8 )
        //    break;
       fof_single_finder( p );
    }

    free( fof_map );
    fof_sort();
    for( p=0; p<Npixs; p++ )
        if ( Len[p] == 1 )
            break;
    NfofRegion = p;
}

void find_region_init() {
    fof_init();
}

void find_region_free() {
    fof_free();
}

void fof_edge_save( int iter ) {

    char buf[100];
    FILE *fd;
    long i, p;
    int xi, yi;

    if ( ThisTask == 0 )
        printf( "save edge ...\n" );
    sprintf( buf, "%s/%s_edge_%04i.dat", All.OutputDir, FileName, iter );
    fd = fopen( buf, "w" );
    for( i=0; i<NfofEdge; i++ ) {
        p = Head[i];
        fprintf( fd, "%4li ", Len[i] );
        while( p>=0 ) {
            xi = p % Width;
            fprintf( fd, "%4i ",  xi );
            p = Next[p];
        }
        fprintf( fd, "\n" );

        p = Head[i];
        fprintf( fd, "%4li ", Len[i] );
        while( p>=0 ) {
            yi = p / Width;
            fprintf( fd, "%4i ",  yi );
            p = Next[p];
        }
        fprintf( fd, "\n" );
     }
    fclose( fd );

}

void fof_region_save( int iter ) {

    char buf[100];
    FILE *fd;
    long i, p;
    int xi, yi;

    if ( ThisTask == 0 )
        printf( "save region ...\n" );
    sprintf( buf, "%s/%s_region_%04i.dat", All.OutputDir, FileName, iter );
    fd = fopen( buf, "w" );
    
    for( i=0; i<NfofRegion; i++ ) {
        p = Head[i];
        fprintf( fd, "%4li ", Len[i] );
        while( p>=0 ) {
            xi = p % Width;
            fprintf( fd, "%4i ",  xi );
            p = Next[p];
        }
        fprintf( fd, "\n" );

        p = Head[i];
        fprintf( fd, "%4li ", Len[i] );
        while( p>=0 ) {
            yi = p / Width;
            fprintf( fd, "%4i ",   yi );
            p = Next[p];
        }
        fprintf( fd, "\n" );
     }

    fclose( fd );

}

void fof_catalog_save( int iter ) {

    char buf[100];
    FILE *fd;
    long i, p;
    int xi, yi;
    double flux, x, y;

    if ( ThisTask == 0 )
        printf( "save catalog ...\n" );

    sprintf( buf, "%s/%s_catalog_%04i.dat", All.OutputDir, FileName, iter );
    fd = fopen( buf, "w" );

    for( i=0; i<NfofRegion; i++ ) {
        p = Head[i];
        x = y = flux = 0;
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            x += xi;
            y += yi;
            /*
            if ( i == 0 )
                printf( "%i %i %g\n", xi, yi, DataRaw[ yi*Width + xi ] );
            */
            flux += DataRaw[ yi * Width + xi ];
            p = Next[p];
        }

        xi = x / Len[i];
        yi = y / Len[i];

        fprintf( fd, "%li %i %i %g\n", Len[i], xi, yi, flux );
    }
        
    fclose(fd);
}

void find_region( int iter ) {

    if ( ThisTask == 0 )
        printf( "find region ...\n" );
    fof_reset();
    fof_edge();
    fof_edge_save( iter );

    fof_reset();
    fof_region();
    fof_region_save( iter );
    fof_catalog_save( iter );

    if ( ThisTask == 0 )
        printf( "[fof], edge: %i, region: %i\n", NfofEdge, NfofRegion );

}
