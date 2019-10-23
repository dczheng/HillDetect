/*
    dczheng
    created 2019-08-14
*/

#include "allvars.h"

int *fof_map;

void group_finder() {

    int p, i, xi, yi, *data, index, h5_ndims, c[2], j, Nfof;
    hsize_t h5_dims[2];
    hid_t h5_dsp, h5_ds, h5_attr, h5_g;
    double flux_tot, *flux, f, fmax;
    char buf[100];

    put_start;
    group_finder_init();

    for( p=0; p<Npixs; p++ ) {
        fof_map[p] = 0;
        if ( Data[p] > SigmaClippingVmin ) {
           fof_map[p] = 1;
        } 

    }

    fof();

    for( p=0; p<Npixs; p++ ) {
        if ( Len[p] == 1 )
            break;
    }
    Nfof = p;

    sprintf( buf, "Group%i", CurGroup  );
    h5_g = H5Gcreate( h5_fof, buf, 0  );


    h5_dsp = H5Screate( H5S_SCALAR );
    h5_attr = H5Acreate( h5_g, "NReg", H5T_NATIVE_INT, h5_dsp,
        H5P_DEFAULT);
    H5Awrite( h5_attr, H5T_NATIVE_INT, &Nfof );
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
            data[index] = yi + YShift + HStartCut;
            data[index+Len[i]] = xi + XShift + WStartCut;
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
            flux[index] = f;
            index ++;
            p = Next[p];
        }

        for( j=0,flux_tot=0; j<Len[i]; j++ ) {
            flux[j] *= fabs( CDELT1*CDELT2 ) / All.Beam;
            flux_tot += flux[j];
        }

        if ( !All.PeakCenterFlag ) {
            c[0] /= (double)Len[i];
            c[1] /= (double)Len[i];
        }

        sprintf( buf, "Center-%i", i );
        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_attr = H5Acreate( h5_g, buf, H5T_NATIVE_INT, h5_dsp,
                H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_INT, c );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        sprintf( buf, "Reg-%i", i );
        h5_ndims = 2;
        h5_dims[0] = 2;
        h5_dims[1] = Len[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_ds = H5Dcreate( h5_g, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        H5Dwrite( h5_ds, H5T_NATIVE_INT, h5_dsp, H5S_ALL, H5P_DEFAULT, data );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );

        sprintf( buf, "Flux-%i", i );
        h5_ndims = 1;
        h5_dims[0] = Len[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_ds = H5Dcreate( h5_g, buf, H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT  );
        H5Dwrite( h5_ds, H5T_NATIVE_DOUBLE, h5_dsp, H5S_ALL, H5P_DEFAULT, flux );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );

        sprintf( buf, "FluxTot-%i", i );
        h5_dsp = H5Screate( H5S_SCALAR );
        h5_attr = H5Acreate( h5_g, buf, H5T_NATIVE_DOUBLE, h5_dsp,
            H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, &flux_tot );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        sprintf( buf, "PeakFlux-%i", i );
        h5_dsp = H5Screate( H5S_SCALAR );
        h5_attr = H5Acreate( h5_g, buf, H5T_NATIVE_DOUBLE, h5_dsp,
            H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, &fmax );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

     }

    free( data );
    free( flux );
    group_finder_free();
    put_end;
}

void group_finder_init() {
    fof_map = malloc( sizeof(int) * Npixs );
    fof_init( fof_map, Width, Height );
}

void group_finder_free() {
    fof_free();
    free(fof_map);
}


void lset_group_finder_save() {

    put_start;
    char buf[100];
    int i, p, index, xi, yi, *data, h5_ndims, Xs[2];
    hsize_t h5_dims[2];
    hid_t h5_dsp, h5_ds, h5_attr, h5_f, h5_g;

    data = malloc( sizeof(data) * Npixs * 2 );

    sprintf( buf, "%s/lset_edges_regs.hdf5", All.OutputDir  );
    h5_f = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT  );

    h5_ndims = 1;
    h5_dims[0] = 2;
    Xs[0] = HStartCut;
    Xs[1] = WStartCut;
    h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
    h5_attr = H5Acreate( h5_f, "CRPIX", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp );

    Xs[0] = HEndCut-HStartCut;
    Xs[1] = WEndCut-WStartCut;
    h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
    h5_attr = H5Acreate( h5_f, "NAXIS", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp );

    // save edge
    h5_g = H5Gcreate( h5_f, "Edges", 0  );

    h5_dsp = H5Screate( H5S_SCALAR );
    h5_attr = H5Acreate( h5_g, "NEdge", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
    H5Awrite( h5_attr, H5T_NATIVE_INT, &lset_Nedge  );
    H5Aclose( h5_attr  );
    H5Sclose( h5_dsp  );

    for( i=0; i<lset_Nedge; i++ ) {
        p = lset_Head_edge[i];
        index = 0;
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            /*
            printf( "%i %i %i, %i %i, %i\n",
                    Width, Height, Npixs, xi, yi, lset_Len_edge[i] );
            */
            data[index] = yi;
            data[index+lset_Len_edge[i]] = xi;
            index ++;
            p = lset_Next_edge[p];
        }

        h5_ndims = 2;
        h5_dims[0] = 2;
        h5_dims[1] = lset_Len_edge[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL   );

        sprintf( buf, "edge-%i", i );
        h5_ds = H5Dcreate( h5_g, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        H5Dwrite( h5_ds, H5T_NATIVE_INT, h5_dsp, H5S_ALL, H5P_DEFAULT, data );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );
     }
    H5Gclose( h5_g );

    // save regions
    h5_g = H5Gcreate( h5_f, "Regs", 0  );

    h5_dsp = H5Screate( H5S_SCALAR );
    h5_attr = H5Acreate( h5_g, "NReg", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
    H5Awrite( h5_attr, H5T_NATIVE_INT, &lset_Nreg  );
    H5Aclose( h5_attr  );
    H5Sclose( h5_dsp  );

    for( i=0; i<lset_Nreg; i++ ) {
        p = lset_Head_reg[i];
        index = 0;
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            data[index] = yi;
            data[index+lset_Len_reg[i]] = xi;
            index ++;
            p = lset_Next_reg[p];
        }

        h5_ndims = 2;
        h5_dims[0] = 2;
        h5_dims[1] = lset_Len_reg[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL   );

        sprintf( buf, "Reg-%i", i );
        h5_ds = H5Dcreate( h5_g, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        H5Dwrite( h5_ds, H5T_NATIVE_INT, h5_dsp, H5S_ALL, H5P_DEFAULT, data );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );

     }

     H5Gclose( h5_g );

     H5Fclose( h5_f );

     free( data );
     put_end;

}

void lset_group_finder() {

    int p, np;

    put_start;
    mytimer_start;
    group_finder_init();

    // find edge
    printf( "find edges\n" );
    for( p=0; p<Npixs; p++ ) {
        fof_map[p] = 0;
    }
    for( p=0; p<edgen; p++ ) {
        fof_map[ edgey[p] * Width + edgex[p] ] = 1;
    }
    fof();
    for( p=0; p<Npixs; p++ ) {
        if ( Len[p] == 1 )
            break;
    }

    lset_Nedge = p;
    memcpy( lset_Next_edge, Next, sizeof(int)*Npixs );
    memcpy( lset_Head_edge, Head, sizeof(int)*Npixs );
    memcpy( lset_Len_edge, Len, sizeof(int)*Npixs );

    // find region
    printf( "find regions\n" );
    fof_reset();
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
    lset_Nreg = p;
    memcpy( lset_Next_reg, Next, sizeof(int)*Npixs );
    memcpy( lset_Head_reg, Head, sizeof(int)*Npixs );
    memcpy( lset_Len_reg, Len, sizeof(int)*Npixs );

    group_finder_free();

    printf( "lset, Nedge: %i, Nreg: %i\n",\
                lset_Nedge, lset_Nreg );

    lset_group_finder_save();

    mytimer_end;
    put_end;

}
