/*
    dczheng
    created 2019-08-14
*/

#include "allvars.h"

int *fof_map;

void group_finder() {

    int p, i, xi, yi, xig, yig, *data, index, h5_ndims, j, Nfof;
    hsize_t h5_dims[2];
    hid_t h5_dsp, h5_ds, h5_attr, h5_g;
    double flux_tot, *flux, f, fmax, xyerr[2], c[2], *img, mean, sigma, N, vmin;
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
    img = malloc( sizeof(double) * Npixs );

    vmin = 1e100;
    for( yi=0; yi<Height; yi++ )
        for( xi=0; xi<Width; xi++ ) {
            f = DataRaw[ (yi + YShift) * WidthGlobal + ( xi + XShift ) ] *
                    fabs( CDELT1*CDELT2 ) / All.Beam;
            img[yi*Width+xi] = f;
            if ( f<vmin )
                vmin = f;
        }

    N = Npixs;
    for( i=0; i<Nfof; i++ ) {
        p = Head[i];
        index = 0;
        fmax = -1e10;
        c[0] = c[1] = 0;
        N -= Len[i];

        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            xig = xi + XShift + WStartCut;
            yig = yi + YShift + HStartCut;
            img[ yi*Width+xi ] = vmin/10;
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

        sprintf( buf, "Center-%i", i );
        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_attr = H5Acreate( h5_g, buf, H5T_NATIVE_DOUBLE, h5_dsp,
                H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, c );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        sprintf( buf, "XYerr-%i", i );
        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_attr = H5Acreate( h5_g, buf, H5T_NATIVE_DOUBLE, h5_dsp,
                H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, xyerr );
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

     if ( N <= 0)  {
         printf( "WARNING: no pixs to estimate sigma!\n" );
         sigma = -1;
     }
     else {
         mean = 0;
         for( p=0; p<Npixs; p++ ) {
             if ( img[p]<vmin )
                 continue;
             mean += img[p];
         }
         mean /= N;

         sigma = 0;
         for( p=0; p<Npixs; p++ ) {
             if ( img[p]<vmin )
                 continue;
             sigma += SQR(img[p]-mean);
         }
         sigma /= N;
         sigma = sqrt( sigma );
     }

     sprintf( buf, "Sigma" );
     h5_dsp = H5Screate( H5S_SCALAR  );
     h5_attr = H5Acreate( h5_g, buf, H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT  );
     H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, &sigma  );
     H5Aclose( h5_attr );
     H5Sclose( h5_dsp  );

    H5Gclose( h5_g  );

    free( data );
    free( flux );
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
    int i, p, index, xi, yi, xig, yig, *data, h5_ndims, Xs[2], 
            xmin, xmax, ymin, ymax, N, j;
    double mean, sigma, *img, f, fmax, *flux, flux_tot, xyerr[2], c[2], vmin;
    hsize_t h5_dims[2];
    hid_t h5_dsp, h5_ds, h5_attr, h5_f;

    data = malloc( sizeof(int) * Npixs * 2 );
    flux = malloc( sizeof(double) * Npixs );

    sprintf( buf, "%s/lset_regs.hdf5", All.OutputDir  );
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

    h5_dsp = H5Screate( H5S_SCALAR );
    h5_attr = H5Acreate( h5_f, "NReg", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
    H5Awrite( h5_attr, H5T_NATIVE_INT, &lset_Nreg  );
    H5Aclose( h5_attr  );
    H5Sclose( h5_dsp  );

    Xs[0] = HStartCut;
    Xs[1] = WStartCut;
    h5_ndims = 1;
    h5_dims[0] = 2;
    h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
    h5_attr = H5Acreate( h5_f, "CutStart", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
    H5Aclose( h5_attr );

    Xs[0] = HEndCut;
    Xs[1] = WEndCut;
    h5_ndims = 1;
    h5_dims[0] = 2;
    h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
    h5_attr = H5Acreate( h5_f, "CutEnd", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
    H5Aclose( h5_attr );

    img = malloc( sizeof(double) * Npixs);
    vmin = 1e100;
    for( yi=0; yi<Height; yi++ )
        for( xi=0; xi<Width; xi++ ) {
         xig = xi + WStartCut;
         yig = yi + HStartCut;
         f = DataRaw[ yig * Width + xig ] * 
                    fabs( CDELT1*CDELT2 ) / All.Beam;
         if ( f<vmin )
             vmin = f;
    }

    for( i=0; i<lset_Nreg; i++ ) {

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

            f = DataRaw[ yig * Width + xig ] * 
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
            //printf( "%i %i %i %i, %g\n", xig, yig, c[0], c[1], f );

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
        //printf( "%i, %i\n", c[0], c[1] );

        p = lset_Head[i];
        xyerr[0] = xyerr[1] = 0;
        while( p>=0 ) {
            xi = p % Width;
            yi = p / Width;
            xig = xi + WStartCut;
            yig = yi + HStartCut;
            xyerr[0] += SQR( yig-c[0] );
            xyerr[1] += SQR( xig-c[1] );
            //printf( "%i %i %g %g\n", yig, xig, c[0], c[1] );
            p = lset_Next[p];
        }

        xyerr[0] = sqrt( xyerr[0] / lset_Len[i] );
        xyerr[1] = sqrt( xyerr[1] / lset_Len[i] );
        //printf( "%g %g\n", xyerr[0], xyerr[1] );

        N = (xmax-xmin+1)*(ymax-ymin+1) - lset_Len[i];
        if ( N <= 0)  {
            printf( "WARNING: no pixs to estimate sigma!\n" );
            sigma = -1;
        }
        else {
            memset( img, 0, sizeof(double)*Npixs );
            for( yi=ymin; yi<=ymax; yi++ )
                for( xi=xmin; xi<=xmax; xi++ ) {
                    img[yi*Width+xi] = DataRaw[yi*Width+xi] * 
                                    fabs( CDELT1*CDELT2 ) / All.Beam;
                }

            p = lset_Head[i];
            while( p>=0 ) {
                xi = p % Width;
                yi = p / Width;
                img[yi*Width+xi] =  vmin/10;
                p = lset_Next[p];
            }

            mean = 0;
            for( yi=ymin; yi<=ymax; yi++ )
                for( xi=xmin; xi<=xmax; xi++ ) {
                    if ( img[yi*Width+xi] < vmin )
                        continue;
                    mean += img[yi*Width+xi];
                }
    
            mean /= N;
            sigma = 0;
            for( yi=ymin; yi<=ymax; yi++ )
                for( xi=xmin; xi<=xmax; xi++ ) {
                    if ( img[yi*Width+xi] < vmin )
                        continue;
                    sigma += SQR(img[yi*Width+xi]-mean);
                }
    
            sigma /=  N;
            //printf( "mean: %g, sigma: %g\n", mean, sigma );

            if ( sigma<0 ) {
                printf( "%i %i %i %i %i\n",
                    xmin+WStartCut,
                    xmax+WStartCut,
                    ymin+HStartCut,
                    ymax+HStartCut,
                    lset_Len[i]
                    );
                printf( "%g ", sigma );
    
                endrun( "error" );
            }
            sigma = sqrt(sigma);
        }

        sprintf( buf, "Center-%i", i );
        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_attr = H5Acreate( h5_f, buf, H5T_NATIVE_DOUBLE, h5_dsp,
                H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, c );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        sprintf( buf, "XYerr-%i", i );
        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_attr = H5Acreate( h5_f, buf, H5T_NATIVE_DOUBLE, h5_dsp,
                H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, xyerr );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        sprintf( buf, "Flux-%i", i );
        h5_ndims = 1;
        h5_dims[0] = lset_Len[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
        h5_ds = H5Dcreate( h5_f, buf, H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT  );
        H5Dwrite( h5_ds, H5T_NATIVE_DOUBLE, h5_dsp, H5S_ALL, H5P_DEFAULT, flux );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );

        sprintf( buf, "FluxTot-%i", i );
        h5_dsp = H5Screate( H5S_SCALAR );
        h5_attr = H5Acreate( h5_f, buf, H5T_NATIVE_DOUBLE, h5_dsp,
            H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, &flux_tot );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        sprintf( buf, "PeakFlux-%i", i );
        h5_dsp = H5Screate( H5S_SCALAR );
        h5_attr = H5Acreate( h5_f, buf, H5T_NATIVE_DOUBLE, h5_dsp,
            H5P_DEFAULT);
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, &fmax );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        h5_ndims = 2;
        h5_dims[0] = 2;
        h5_dims[1] = lset_Len[i];
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL   );
        sprintf( buf, "Reg-%i", i );
        h5_ds = H5Dcreate( h5_f, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        H5Dwrite( h5_ds, H5T_NATIVE_INT, h5_dsp, H5S_ALL, H5P_DEFAULT, data );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp  );

        sprintf( buf, "Sigma-%i", i );
        h5_dsp = H5Screate( H5S_SCALAR  );
        h5_attr = H5Acreate( h5_f, buf, H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT  );
        H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, &sigma  );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

        sprintf( buf, "Sigma-Npixs-%i", i );
        h5_dsp = H5Screate( H5S_SCALAR  );
        h5_attr = H5Acreate( h5_f, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT  );
        H5Awrite( h5_attr, H5T_NATIVE_INT, &N  );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp  );

     }

     H5Fclose( h5_f );
    
     free(data);
     free(flux);
     free(img);
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
            if ( x == 0 || y == 0 || x == Width-1 || y == Height-1 ) {
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
