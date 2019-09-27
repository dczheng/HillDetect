/*
    dczheng
    created 2019-08-14
*/

#include "allvars.h"

int *fof_map;

void fof_edge() {

    int i;

    for( i=0; i<Npixs; i++ ) {
        fof_map[i] = 0;
    }

    for( i=0; i<edgen; i++ ) {
        fof_map[ edgey[i] * Width + edgex[i] ] = 1;
    }

    fof();

    for( i=0; i<Npixs; i++ ) {
//        printf( "%i \n", Len[i] );
        if ( Len[i] == 1 )
            break;
    }
    NfofEdge = i;

}

void fof_region() {

    int p;
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


    //printf( "NPixs: %i, np; %i, nm: %i\n", Npixs, np, Npixs - np );

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

void fof_edge_save( int iter, int mode ) {

    int i, p;
    int xi, yi;

    if ( mode == 0 )
        if ( ThisTask )
            return;

    //writelog( "save edge ...\n" );

    for( i=0; i<NfofEdge; i++ ) {
        p = Head[i];
        fprintf( EdgesFd, "%4i ", Len[i] );
        while( p>=0 ) {
            xi = p % Width;
            fprintf( EdgesFd, "%4i ",  xi + XShift );
            p = Next[p];
        }
        fprintf( EdgesFd, "\n" );

        p = Head[i];
        fprintf( EdgesFd, "%4i ", Len[i] );
        while( p>=0 ) {
            yi = p / Width;
            fprintf( EdgesFd, "%4i ",  yi + YShift );
            p = Next[p];
        }
        fprintf( EdgesFd, "\n" );
     }

}

void fof_region_save( int iter ) {

    int i, p;
    int xi, yi;
    
    for( i=0; i<NfofRegion; i++ ) {
        p = Head[i];
        fprintf( RegsFd, "%4i ", Len[i] );
        while( p>=0 ) {
            xi = p % Width;
            fprintf( RegsFd, "%4i ",  xi + XShift );
            p = Next[p];
        }
        fprintf( RegsFd, "\n" );

        p = Head[i];
        fprintf( RegsFd, "%4i ", Len[i] );
        while( p>=0 ) {
            yi = p / Width;
            fprintf( RegsFd, "%4i ",   yi + YShift );
            p = Next[p];
        }
        fprintf( RegsFd, "\n" );
     }

}

void fof_catalog_save( int iter ) {

}

void fof_ds9_region_save( int iter ) {

    char buf[100];
    FILE *fd;
    int i, p, l;
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

void find_region( int iter, int mode ) {

    if ( mode == 0 ) {
        writelog( "find region ...\n" );
    }
    else{
        fprintf( LogFileFd, "find region ...\n" );
    }

    fof_reset();
    fof_edge();
    fof_edge_save( iter, mode );

    if ( mode == 0 )
        return;

    //fof_ds9_region_save( iter );

    fof_reset();
    fof_region();
    fof_region_save( iter );
    fof_catalog_save( iter );

    fprintf( LogFileFd, "[fof], edge: %i, region: %i\n", NfofEdge, NfofRegion );
    put_end();

}
