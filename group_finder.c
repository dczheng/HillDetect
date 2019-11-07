/*
    dczheng
    created 2019-08-14
*/

#include "allvars.h"

int *fof_map;

void group_finder() {

    int p, i, xi, yi, xig, yig, *data, index,
            j, Nfof, *flag;
    double flux_tot, *flux, f, fmax, xyerr[2], c[2], *img,
            mean[2], sigma[2], rms[2], N;
    char buf[100];
    hsize_t h5_dims[2];
    hid_t h5_g, h5_gg;

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

    sprintf( buf, "Group%i", CurGroup );
    h5_g = H5Gcreate( h5_fof, buf, 0 );
    hdf5_write_attr_scalar( h5_g, H5T_NATIVE_INT, "NReg", &Nfof );

    data = malloc( sizeof(int) * Npixs * 2 );
    flux = malloc( sizeof(double) * Npixs );
    img = malloc( sizeof(double) * Npixs );
    flag = malloc( sizeof(int) * Npixs );
    memset( flag, 0, sizeof(int) * Npixs );

    N = 0;
    for( i=0; i<Nfof; i++ ) {
        sprintf( buf, "Reg%i", i );
        h5_gg = H5Gcreate( h5_g, buf, 0 );

        p = Head[i];
        index = 0;
        fmax = -1e10;
        c[0] = c[1] = 0;
        N += Len[i];

        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            xig = xi + XShift + WStartCut;
            yig = yi + YShift + HStartCut;
            flag[p] = 1;
            data[index] = yig;
            data[index+Len[i]] = xig;
            if  ( index >= Npixs  * 2 || index+Len[i] >= Npixs * 2  ) {
                printf( "%i %i %i [%i %i]\n", index, Len[i], Npixs, Width, Height );
                endrun("can't be!");
            }
                
            f = DataRaw[ (yi + YShift) * WidthGlobal + ( xi + XShift ) ] *
                    fabs( CDELT1*CDELT2 ) / All.Beam;

            if ( All.PeakCenterFlag ) {
                if ( f>fmax ) {
                    c[0] = yig;
                    c[1] = xig;
                }
            }
            else {
                c[0] += yig * f;
                c[1] += xig * f;
            }

            if ( f>fmax )
               fmax = f;

            flux[index] = f;
            index ++;
            p = Next[p];
        }

        for( j=0,flux_tot=0; j<Len[i]; j++ ) {
            flux_tot += flux[j];
        }

        if ( !All.PeakCenterFlag ) {
            c[0] /= flux_tot;
            c[1] /= flux_tot;
        }

        p = Head[i];
        xyerr[0] = xyerr[1] = 0;
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            xig = xi + XShift + WStartCut;
            yig = yi + YShift + HStartCut;
            xyerr[0] += SQR( yig-c[0] );
            xyerr[1] += SQR( xig-c[1] );
            p = Next[p];
        }

        xyerr[0] = sqrt( xyerr[0] / Len[i] );
        xyerr[1] = sqrt( xyerr[1] / Len[i] );

        h5_dims[0] = 2;

        hdf5_write_attr_nd( h5_gg, H5T_NATIVE_DOUBLE, h5_dims, 1, "center", c );

        hdf5_write_attr_nd( h5_gg, H5T_NATIVE_DOUBLE, h5_dims, 1, "xyerr", xyerr );

        h5_dims[0] = 2;
        h5_dims[1] = Len[i];
        hdf5_write_data( h5_gg, H5T_NATIVE_INT, h5_dims, 2, "region", data );

        h5_dims[0] = Len[i];
        hdf5_write_data( h5_gg, H5T_NATIVE_DOUBLE, h5_dims, 1, "flux", flux );

        hdf5_write_attr_scalar( h5_gg, H5T_NATIVE_DOUBLE, "flux_tot", &flux_tot );

        hdf5_write_attr_scalar( h5_gg, H5T_NATIVE_DOUBLE, "peak_flux", &fmax );

        H5Gclose( h5_gg );
    }

    for( i=0; i<2; i++ )
        mean[i] = sigma[i] = rms[i] = 0;

    for( p=0; p<Npixs; p++ ){
        xi = p % Width;
        yi = p / Width;
        img[p] = DataRaw[ (yi + YShift) * WidthGlobal + ( xi + XShift ) ] *
               fabs( CDELT1*CDELT2 ) / All.Beam;
    }


    get_mean_sigma_rms( img, flag, Npixs, mean, sigma, rms );

    hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "mean_outer", mean );
    hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "mean_inner", mean+1 );
    hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "rms_outer", rms );
    hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "rms_inner", rms+1 );
    hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "sigma_outer", sigma );
    hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "sigma_inner", sigma+1 );

    H5Gclose( h5_g  );

    free( data );
    free( flux );
    free( flag );
    free( img );

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
    int i, p, index, xi, yi, xig, yig, *data, Xs[2], 
            xmin, xmax, ymin, ymax, N, j, *flag, NN, t;
    double mean[2], sigma[2], rms[2],
                *img, f, fmax, *flux, flux_tot, xyerr[2], c[2];
    hsize_t h5_dims[2];
    hid_t h5_f, h5_g, h5_gg;

    data = malloc( sizeof(int) * Npixs * 2 );
    flux = malloc( sizeof(double) * Npixs );

    sprintf( buf, "%s/lset_regs.hdf5", All.OutputDir  );
    h5_f = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT  );

    h5_dims[0] = 2;
    Xs[0] = HStartCut;
    Xs[1] = WStartCut;

    hdf5_write_attr_nd( h5_f, H5T_NATIVE_INT, h5_dims, 1, "CRPIX", Xs );

    h5_dims[0] = 2;
    Xs[0] = HEndCut-HStartCut;
    Xs[1] = WEndCut-WStartCut;
    hdf5_write_attr_nd( h5_f, H5T_NATIVE_INT, h5_dims, 1, "NAXIS", Xs );

    hdf5_write_attr_scalar( h5_f, H5T_NATIVE_INT, "Ngroup", &lset_Nreg );

    h5_dims[0] = 2;
    Xs[0] = HStartCut;
    Xs[1] = WStartCut;
    hdf5_write_attr_nd( h5_f, H5T_NATIVE_INT, h5_dims, 1, "CutStart", Xs );

    h5_dims[0] = 2;
    Xs[0] = HEndCut;
    Xs[1] = WEndCut;
    hdf5_write_attr_nd( h5_f, H5T_NATIVE_INT, h5_dims, 1, "CutEnd", Xs );

    img = malloc( sizeof(double) * Npixs);
    flag = malloc( sizeof(int) * Npixs);

    for( i=0; i<lset_Nreg; i++ ) {

        sprintf( buf, "Group%i", i );
        h5_g = H5Gcreate( h5_f, buf, 0 );
        t = 1;
        hdf5_write_attr_scalar( h5_g, H5T_NATIVE_INT, "NReg", &t );
        h5_gg = H5Gcreate( h5_g, "Reg0", 0 );

        p = lset_Head[i];
        index = 0;
        xmin = ymin = INT_MAX;
        xmax = ymax = -INT_MAX;

        fmax = -1e10;
        c[0] = c[1] = 0;

        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            xig = xi + WStartCut;
            yig = yi + HStartCut;
            data[index] = yig;
            data[index+lset_Len[i]] = xig;

            f = DataRaw[ yi * Width + xi ] * 
                    fabs( CDELT1*CDELT2 ) / All.Beam;

            flux[index] = f;

            if ( All.PeakCenterFlag ) {
                if ( f>fmax ) {
                    c[0] = yig;
                    c[1] = xig;
                }
            }
            else {
                c[0] += yig * f;
                c[1] += xig * f;
            }

            if ( f>fmax )
                  fmax = f;

            index ++;

            xmin = ( xi<xmin ) ? xi : xmin;
            xmax = ( xi>xmax ) ? xi : xmax;
            ymin = ( yi<ymin ) ? yi : ymin;
            ymax = ( yi>ymax ) ? yi : ymax;
            p = lset_Next[p];
        }


        for( j=0,flux_tot=0; j<lset_Len[i]; j++ ) {
            flux_tot += flux[j];
        }

        if ( !All.PeakCenterFlag ) {
            c[0] /= flux_tot;
            c[1] /= flux_tot;
        }

        p = lset_Head[i];
        xyerr[0] = xyerr[1] = 0;
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            xig = xi + WStartCut;
            yig = yi + HStartCut;
            xyerr[0] += SQR( yig-c[0] );
            xyerr[1] += SQR( xig-c[1] );
            p = lset_Next[p];
        }

        xyerr[0] = sqrt( xyerr[0] / lset_Len[i] );
        xyerr[1] = sqrt( xyerr[1] / lset_Len[i] );

        N = xmax-xmin+1;
        NN = N * (ymax-ymin+1);

        for( yi=ymin; yi<=ymax; yi++ ) {
            for( xi=xmin; xi<=xmax; xi++ )
                img[ yi*N + xi ] = DataRaw[ yi * Width + xi ] * 
                    fabs( CDELT1*CDELT2 ) / All.Beam;
            }

        memset( flag, 0, sizeof(int) * NN);
        p = lset_Head[i];
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            flag[ (yi-ymin)*N + (xi-xmin) ] = 1;
            p = lset_Next[p];
        }


        get_mean_sigma_rms( img, flag, NN, mean, sigma, rms );

        hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "mean_outer", mean );

        hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "mean_inner", mean+1 );

        hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "rms_outer", rms );

        hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "rms_inner", rms+1 );

        hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "sigma_outer", sigma );

        hdf5_write_attr_scalar( h5_g, H5T_NATIVE_DOUBLE,
                        "sigma_inner", sigma+1 );

        h5_dims[0] = 2;
        hdf5_write_attr_nd( h5_gg, H5T_NATIVE_DOUBLE, h5_dims, 1, "center", c );

        h5_dims[0] = 2;
        hdf5_write_attr_nd( h5_gg, H5T_NATIVE_DOUBLE, h5_dims, 1, "xyerr", xyerr );

        h5_dims[0] = lset_Len[i];
        hdf5_write_data( h5_gg, H5T_NATIVE_DOUBLE, h5_dims, 1, "flux", flux );

        hdf5_write_attr_scalar( h5_gg, H5T_NATIVE_DOUBLE, "flux_tot", &flux_tot );

        hdf5_write_attr_scalar( h5_gg, H5T_NATIVE_DOUBLE, "peak_flux", &fmax );

        h5_dims[0] = 2;
        h5_dims[1] = lset_Len[i];
        hdf5_write_data( h5_gg, H5T_NATIVE_DOUBLE, h5_dims, 2, "region", data );

        H5Gclose( h5_g );
        H5Gclose( h5_gg );

     }

     H5Fclose( h5_f );
    
     free(data);
     free(flux);
     free(img);
     free(flag);
     put_end;

}

void lset_group_finder() {

    int p, np, pc, i, x, y, flag;

    put_start;
    mytimer_start;
    group_finder_init();

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

    memcpy( lset_Next, Next, sizeof(int)*Npixs );
    memcpy( lset_Head, Head, sizeof(int)*Npixs );
    memcpy( lset_Len, Len, sizeof(int)*Npixs );

    group_finder_free();

    printf( "lset, Nreg: %i\n", lset_Nreg );

    pc = 0;
    for (i=0; i<lset_Nreg; i++) {
        p = lset_Head[i];
        flag = 0;
        while( p>=0 ) {
            x = p % Width;
            y = p / Width;
            p = lset_Next[p];
            if ( x == 0 || y == 0 ||
                 x == Width-1 || y == Height-1
                 ) {
                flag = 1;
                break;
            }
        }
        if ( flag == 0 ) {
            lset_Head[pc] = lset_Head[i];
            lset_Len[pc] = lset_Len[i];
            pc++;
        }
    }

    lset_Nreg = pc;
    printf( "Nreg [remove edge]: %i\n", lset_Nreg );

    for (i=0; i<lset_Nreg; i++)
        if ( lset_Len[i] < All.MinReg )
            break;
    lset_Nreg = i;
    printf( "Nreg [>%i]: %i\n", All.MinReg, lset_Nreg );

    lset_group_finder_save();

    mytimer_end;
    put_end;

}
