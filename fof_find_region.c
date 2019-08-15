/*
    dczheng
    created 2019-08-14
*/

#include "allvars.h"

int *fof_map;

void fof_edge() {

    long i;

    for( i=0; i<Npixs; i++ ) {
        fof_map[i] = 0;
    }

    for( i=0; i<edgen; i++ ) {
        fof_map[ edgey[i] * Width + edgex[i] ] = 1;
    }

    fof();

    for( i=0; i<Npixs; i++ )
        if ( Len[i] == 1 )
            break;
    NfofEdge = i;

}

void fof_region() {

    long p;
    int np;

    np = 0;
    for( p=0; p<Npixs; p++ ) {
        fof_map[p] = 0;
        if ( Phi[p] > 0 ) {
           fof_map[p] = 1;
           np ++;
        }

    }

    if ( np > Npixs / 2 ) {
        np = 0;
        for( p=0; p<Npixs; p++ ) {
            fof_map[p] = 0;
            if ( Phi[p] < 0 ) {
            fof_map[p] = 1;
            np ++;
            }
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

    fof();

    for( p=0; p<Npixs; p++ )
        if ( Len[p] == 1 )
            break;
    NfofRegion = p;
}

void find_region_init() {
    fof_map = malloc( sizeof(int) * Npixs );
    fof_init( fof_map, Width, Height );
}

void find_region_free() {
    fof_free();
    free(fof_map);
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

void fof_ds9_region_save( int iter ) {

    char buf[100];
    FILE *fd;
    long i, p, l;
    int xi, yi;

    char *h_ds9 ="# Region file format: DS9 version 4.1 \n" \
        "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" " \
        "select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n" \
        "image\n\n";

    if ( ThisTask == 0 )
        printf( "save ds9 region ...\n" );

    sprintf( buf, "%s/%s_ds9_region_%04i.reg", All.OutputDir, FileName, iter );
    fd = fopen( buf, "w" );
    fprintf( fd, "%s", h_ds9 );

    for( i=0; i<NfofEdge; i++ ) {
        fprintf( fd, "polygon(" );
        p = Head[i];
        for( l=0; l<Len[i]; l++ ) {
            xi = p % Width;
            yi = p / Width;
            if ( l == Len[i]-1 )
                fprintf( fd, "%i,%i",  xi, yi );
            else 
                fprintf( fd, "%i,%i,",  xi, yi );
            p = Next[p];
        }
        fprintf( fd, ")\n" );
    }
        
    fclose(fd);
}

void find_region( int iter ) {

    if ( ThisTask == 0 )
        printf( "find region ...\n" );
    fof_reset();
    fof_edge();
    fof_edge_save( iter );
    fof_ds9_region_save( iter );

    fof_reset();
    fof_region();
    fof_region_save( iter );
    fof_catalog_save( iter );

    if ( ThisTask == 0 )
        printf( "[fof], edge: %i, region: %i\n", NfofEdge, NfofRegion );

}
