/*
    dczheng
    created 2019-08-14
*/

#include "allvars.h"

int *fof_map;

void lset_find_edge() {

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
    Nfof = i;

}

void lset_find_region() {

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
    fof();
    for( p=0; p<Npixs; p++ )
        if ( Len[p] == 1 )
            break;
    Nfof = p;
}

void group_finder() {

    int p, i, xi, yi, *data, index, h5_ndims, crpix[2], c[2], j;
    hsize_t h5_dims[2];
    hid_t h5_dsp, h5_ds, h5_attr;
    double flux_tot, *flux, f, fmax;
    char buf[100];

    group_finder_init();

    for( p=0; p<Npixs; p++ ) {
        fof_map[p] = 0;
        if ( Data[p] > 0 ) {
           fof_map[p] = 1;
        } 

    }

    fof();

    //printf( "%i %i %i\n", Npixs, Width, Height );
    for( p=0; p<Npixs; p++ ) {
        //printf( "%i %i %i\n", Head[p], Next[Head[p]], Len[p] );
        if ( Len[p] == 1 )
            break;
    }
    Nfof = p;
    //printf( "Nfof: %i\n", Nfof );
    /*
    for( p=0; p<Npixs; p++ ) {
        printf( "%i ", Next[p] );
    }
    printf( "\n" );
    */

    h5_dsp = H5Screate( H5S_SCALAR );
    h5_attr = H5Acreate( h5_RegsGroup, "NRegs", H5T_NATIVE_INT, h5_dsp,
        H5P_DEFAULT);
    H5Awrite( h5_attr, H5T_NATIVE_INT, &Nfof );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp  );

    crpix[0] = YShift;
    crpix[1] = XShift;

    h5_ndims = 1;
    h5_dims[0] = 2;
    h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
    h5_attr = H5Acreate( h5_RegsGroup, "CRPIX", H5T_NATIVE_INT, h5_dsp,
        H5P_DEFAULT);
    H5Awrite( h5_attr, H5T_NATIVE_INT, crpix );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp  );

    data = malloc( sizeof(int) * Npixs * 2 );
    flux = malloc( sizeof(double) * Npixs );
    for( i=0; i<Nfof; i++ ) {
        p = Head[i];
        index = 0;
        if ( All.PeakCenterFlag ) {
            fmax = -1e10;
        }
        else {
            c[0] = c[1] = 0;
        }
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            data[index] = yi;
            data[index+Len[i]] = xi;
            if  ( index >= Npixs  * 2 || index+Len[i] >= Npixs * 2  ) {
                printf( "%i %i %i [%i %i]\n", index, Len[i], Npixs, Width, Height );
                endrun("can't be!");
            }
                
            f = DataRaw[ (yi + YShift) * WidthGlobal + ( xi + XShift ) ];

            if ( All.PeakCenterFlag ) {
                if ( f>fmax ) {
                    fmax = f;
                    c[0] = yi;
                    c[1] = xi;
                }
            }
            else {
                c[0] += yi;
                c[1] += xi;
            }
            /*
            if ( ThisTask == 0 )
                printf( "%i [%i %i]\n", p, xi, yi );
            */
            flux[index] = f;
            index ++;
            p = Next[p];
        }

        for( j=0,flux_tot=0; j<Len[i]; j++ )
            flux_tot += flux[j];

        if ( !All.PeakCenterFlag ) {
            c[0] /= (double)Len[i];
            c[1] /= (double)Len[i];
        }

        sprintf( buf, "Center-%i", i );
        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_attr = H5Acreate( h5_RegsGroup, buf, H5T_NATIVE_INT, h5_dsp,
                H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_INT, c );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        sprintf( buf, "Reg-%i", i );
        h5_ndims = 2;
        h5_dims[0] = 2;
        h5_dims[1] = Len[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_ds = H5Dcreate( h5_RegsGroup, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        H5Dwrite( h5_ds, H5T_NATIVE_INT, h5_dsp, H5S_ALL, H5P_DEFAULT, data );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );

        sprintf( buf, "Flux-%i", i );
        h5_ndims = 1;
        h5_dims[0] = Len[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_ds = H5Dcreate( h5_RegsGroup, buf, H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT  );
        H5Dwrite( h5_ds, H5T_NATIVE_DOUBLE, h5_dsp, H5S_ALL, H5P_DEFAULT, flux );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );

        sprintf( buf, "FluxTot-%i", i );
        h5_dsp = H5Screate( H5S_SCALAR );
        h5_attr = H5Acreate( h5_RegsGroup, buf, H5T_NATIVE_DOUBLE, h5_dsp,
            H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, &flux_tot );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

     }

    free( data );
    free( flux );
    group_finder_free();
}

void group_finder_init() {
    fof_map = malloc( sizeof(int) * Npixs );
    fof_init( fof_map, Width, Height );
}

void group_finder_free() {
    fof_free();
    free(fof_map);
}

void lset_edge_region_save( int iter, int mode, int flag ) {

    char buf[100];
    int i, p, index, xi, yi, *data, h5_ndims;
    hsize_t h5_dims[2];
    hid_t h5_dsp, h5_ds, h5_attr;

    if ( mode == 0 )
        if ( ThisTask )
            return;

    data = malloc( sizeof(data) * Npixs * 2 );


        h5_dsp = H5Screate( H5S_SCALAR );
        if ( flag == 0 ) {
            h5_attr = H5Acreate( h5_EdgesGroup, "NEdges", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        }
        else {
            h5_attr = H5Acreate( h5_RegsGroup, "NRegs", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        }
        H5Awrite( h5_attr, H5T_NATIVE_INT, &Nfof  );
        H5Aclose( h5_attr  );
        H5Sclose( h5_dsp  );

    for( i=0; i<Nfof; i++ ) {
        p = Head[i];
        index = 0;
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            data[index] = yi;
            data[index+Len[i]] = xi;
            index ++;
            p = Next[p];
        }


        h5_ndims = 2;
        h5_dims[0] = 2;
        h5_dims[1] = Len[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL   );

        if ( flag == 0 ) {
            sprintf( buf, "edge-%i", i );
            h5_ds = H5Dcreate( h5_EdgesGroup, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        }
        else {
            sprintf( buf, "reg-%i", i );
            h5_ds = H5Dcreate( h5_RegsGroup, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        }

        H5Dwrite( h5_ds, H5T_NATIVE_INT, h5_dsp, H5S_ALL, H5P_DEFAULT, data );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );

     }

     free( data );

}

void ds9_region_save( int iter ) {

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

    for( i=0; i<Nfof; i++ ) {
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

void lset_group_finder( int iter, int mode ) {

    if ( mode == 0 ) {
        writelog( "find region ...\n" );
    }
    else{
        fprintf( LogFileFd, "find region ...\n" );
    }

    fof_reset();
    lset_find_edge();
    lset_edge_region_save( iter, mode, 0 );

    if ( mode == 0 && All.Lset1 )
        return;

    //fof_ds9_region_save( iter );

    fof_reset();
    lset_find_region();
    lset_edge_region_save( iter, mode, 1 );

    fprintf( LogFileFd, "[fof], edge: %i, region: %i\n", Nfof, Nfof );
    put_end();

}
