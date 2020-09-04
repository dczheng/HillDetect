/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

double RSigma;
int LogNorm;

void sigma_clipping( double sigma, double mean, int flag ) {

    int i;
    double vmin, vcut, vmax, rms, vmin_p;

    vmin = 1e30;
    vmin_p = 1e30;
    vmax = -vmin;
    for( i=0; i<Npixs; i++ ) {
        vmin = ( Data[i] < vmin ) ? Data[i] : vmin;
        vmax = ( Data[i] > vmax ) ? Data[i] : vmax;

        if ( Data[i] > 0 )
            vmin_p = ( Data[i] < vmin_p ) ? Data[i] : vmin_p;
    }

    if ( sigma == 0 ) {
        get_mean_sigma_rms( Data, NULL, Npixs, &mean, &sigma, &rms );
    }

    vcut = mean + RSigma * sigma;
    //printf( "%g %g %g %g\n", mean, RSigma, sigma, vcut );
    //vcut = 1e-11;

    for( i=0; i<Npixs; i++ ) {
        if ( Data[i] < vcut )
            Data[i] = vcut; 
    }

    writelog( flag, "vmin: %g, vmin[+]: %g, vmax: %g\n"
            "Rsigma: %g, mean: %g, sigma: %g, vcut: %g\n",
            vmin, vmin_p, vmax, RSigma, mean, sigma, vcut );
    SigmaClippingVmin = vcut;
}

void background_estimation() {

    int i, j, ii, jj, w, h, m, n, N, GridN, GridM, i0, j0, k, flag, l, r;
    double *d, RS, sigma, mean, v_invalid, vmin, vmax;
    hsize_t h5_dims[2];
    hid_t h5_f;
    char buf[200], buf2[20];

    m = All.BkgEstm;
    n = All.BkgEstn;
    N = All.BkgEstN;
    RS = All.BkgEstRSigma;
    GridN =  All.BkgEstGridN;
    GridM =  All.BkgEstGridM;

    if ( GridN == 0 )
        GridN = n;
    if ( GridM == 0 )
        GridM = m;



    printf( "BkgEstm: %i, BkgEstn: %i, BkgEstN: %i, BkgGridM: %i, BkgGridN: %i, BkgEstRSigma: %g\n",
                m, n, N, GridM, GridN, RSigma);
    w = WidthGlobal / GridM + 1;
    h = HeightGlobal / GridN + 1;

    printf( "Bkg Dim: %i x %i\n", w, h );

    d = malloc( sizeof(double) * m * n );
    Bkg_s = malloc( sizeof(double) * w * h );
    Bkg = malloc( sizeof(double) * Npixs );

    flag = 0;
    v_invalid = DataRawMin / 100;

    mytimer_start;
    writelog( 0, "do clipping ...\n" );

    for( i=0; i<h; i++ ) {
        ii = i * GridM; 
        for( j=0; j<w; j++ ) {
            jj = j * GridN;
            for( i0 = 0; i0 < m; i0++ )
                for( j0 = 0; j0 < n; j0++ ) {
                    l = ii + i0 - m/2;
                    r = jj + j0 - n/2;
                    if ( l < 0 || l > HeightGlobal-1 || r < 0 || r > WidthGlobal-1 )
                        d[i0 * n + j0] = v_invalid;
                    else
                        d[i0 * n + j0] = DataRaw[ l * Width + r  ];
                }
            for( k=0; k<N; k++ ) {
                get_mean_sigma( d, m*n, &v_invalid, &mean, &sigma );
                find_vmin_vmax( d, m*n, NULL, NULL, &vmin, &vmax );
                //printf( "[%i] vmin: %.3e, vmax: %.3e, mean: %.3e, sigma: %.3e, [ %.3e, %.3e ]\n", k, vmin, vmax, mean, sigma,
                 //       mean - RS * sigma, mean + RS * sigma);
                //flag ++;
                for( l=0; l<m*n; l++ ) {
                    if ( d[l] > mean + RS * sigma || d[l] < mean - RS * sigma ) {
                        d[l] = v_invalid;
                        //printf( "ok\n" );
                    }
                }
            }
            get_mean_sigma( d, m*n, &v_invalid, &mean, &sigma );
            Bkg_s[ i * w + j ] = mean;
        }
    }

    printf( "flag: %i\n", flag );
    free(d);

    mytimer;
    writelog( 0, "do interpolation ...\n" );

    if ( w != WidthGlobal && h != HeightGlobal ) {

        const gsl_interp2d_type *L = gsl_interp2d_bilinear;
        const gsl_interp2d_type *C = gsl_interp2d_bicubic;
        double *xs, *ys;
        gsl_spline2d *spline;


        switch( All.InterpMethod ) {
            case 0:
                spline = gsl_spline2d_alloc(L, w, h);
                break;
            case 1:
                spline = gsl_spline2d_alloc(C, w, h);
                break;
            default:
                endrun("Failed to determine interpolation method !");
        }

        printf( "Interpolation Method: %s. \n", gsl_spline2d_name( spline ) );
        sprintf( buf2, "%s", gsl_spline2d_name( spline ) );

        gsl_interp_accel *xacc = gsl_interp_accel_alloc();
        gsl_interp_accel *yacc = gsl_interp_accel_alloc();
        xs = malloc( sizeof(double) * w );
        ys = malloc( sizeof(double) * h );
        for( i=0; i<w; i++ ) xs[i] = i * GridN;
        for( i=0; i<h; i++ ) ys[i] = i * GridM;
        gsl_spline2d_init(spline, xs, ys, Bkg_s, w, h);

        for( i=0; i<HeightGlobal; i++ ) {
            for( j=0; j<WidthGlobal; j++ ) {
                //printf( "%i %i\n", i, j );
                Bkg[i * WidthGlobal + j] = gsl_spline2d_eval(spline, j, i, xacc, yacc);
            }
        }

        gsl_spline2d_free(spline);
        gsl_interp_accel_free(xacc);
        gsl_interp_accel_free(yacc);

    }

    mytimer;
    mytimer_end;

    find_vmin_vmax( Bkg, Npixs, NULL, NULL, &vmin, &vmax );
    printf( "[BKG] vmin: %.3e, vmax: %.3e\n", vmin, vmax );

    sprintf( buf, "%s/bkg_%s.hdf5", All.OutputDir, buf2 );
    h5_f = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    h5_dims[0] = HeightGlobal;
    h5_dims[1] = WidthGlobal;
    hdf5_write_data( h5_f, H5T_NATIVE_DOUBLE, h5_dims, 2, "data", DataRaw );

    h5_dims[0] = HeightGlobal;
    h5_dims[1] = WidthGlobal;
    hdf5_write_data( h5_f, H5T_NATIVE_DOUBLE, h5_dims, 2, "bkg", Bkg );

    h5_dims[0] = h;
    h5_dims[1] = w;
    hdf5_write_data( h5_f, H5T_NATIVE_DOUBLE, h5_dims, 2, "bkg_s", Bkg_s );


    H5Fclose( h5_f );
}


void data_cuting() {

    int i, j;
    int index;

    HStartCut = Height * All.CuttingYStart;
    HEndCut = Height * All.CuttingYEnd;
    WStartCut = Width * All.CuttingXStart;
    WEndCut = Width * All.CuttingXEnd;

    if ( HEndCut<=HStartCut || WEndCut<=WStartCut ) {
        endrun( "Invalid cutting parameters." );
    }

    printf( "Height: %i, Width: %i\n", Height, Width );
    printf( "region: (%g, %g), (%g, %g)\n",
                All.CuttingYStart, All.CuttingYEnd, 
                All.CuttingXStart, All.CuttingXEnd );
    printf( "region: (%i, %i), (%i, %i)\n",
            HStartCut, HEndCut, WStartCut, WEndCut );

    for( i=HStartCut, index=0; i<HEndCut; i++ )
        for ( j=WStartCut; j<WEndCut; j++, index++ ) {
            Data[index] = Data[i*Width+j];
            DataRaw[index] = DataRaw[i*Width+j];
        }

    HeightGlobal = Height = HEndCut - HStartCut;
    WidthGlobal = Width = WEndCut - WStartCut;
    NpixsGlobal = Npixs = Height * Width;
}

void normalize() {

    double vmin, vmax, dv;
    int i;

    vmax = -1e100;
    vmin = 1e100;

    for( i=0; i<Npixs; i++ ) {

        if ( isnan( Data[i] ) )
            continue;

        if ( LogNorm ) {
            if ( Data[i] > 0 )
                vmin = ( vmin < Data[i] ) ? vmin : Data[i];
        }
        else {
            vmin = ( vmin < Data[i] ) ? vmin : Data[i];
        } 
        vmax = ( vmax > Data[i] ) ? vmax : Data[i];
    }


    for( i=0; i<Npixs; i++ ) {

        if ( isnan( Data[i] ) )
            Data[i] = vmin;

        if ( LogNorm )
            if ( Data[i] < vmin )
                Data[i] = vmin;
    }

    if ( LogNorm )
        dv =log( vmax / vmin );
    else
        dv = vmax - vmin;

    for ( i=0; i<Npixs; i++ ) {
        if ( LogNorm )
            Data[i] = log( Data[i]/vmin ) / dv;
        else
            Data[i] = ( Data[i] - vmin ) / dv;
    }

}

/*
void ft_clipping( int mode ){

    int i, j, N2;
    FILE *fd;
    double fac, k;

    if ( mode == 0 )
        return;

    fd = fopen( "t.dat", "w" );
    for( i=0; i<Height; i++ ) {
        for( j=0; j<Width; j++ )
            fprintf( fd, "%g ", Data[i*Width+j] );
        fprintf( fd, "\n" );
    }
    fclose( fd );

    fftw_real *img;
    fftw_complex *fft_of_img;
    rfftwnd_plan fft_plan, fft_plan_inv;

    N2 = 2 * ( Width/2+1 );
    img = malloc( sizeof(fftw_real) * Height * N2 ) ;
    fft_of_img = ( fftw_complex* ) img;

    for( i=0; i<Height; i++ )
        for ( j=0; j<Width; j++ )
            //img[i*N2+j] = log(Data[i*Width+j]);
            img[i*N2+j] = Data[i*Width+j];

    fft_plan     = rfftw2d_create_plan( Height, Width, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE );
    fft_plan_inv = rfftw2d_create_plan( Height, Width, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE );

    rfftwnd_one_real_to_complex( fft_plan, img, NULL );

    fac = 0.01;
    fd = fopen( "fft.dat", "w" );
    for( i=0; i<Height; i++ ) {
        for( j=0; j<N2/2; j++ ) {

            if ( i > Height/2 )
                k = sqrt( SQR(i-Height) + SQR(j) );
            else
                k = sqrt( SQR(i) + SQR(j) );
            if ( k < Height/2 * fac )
                fft_of_img[i*N2/2+j].re = fft_of_img[i*N2/2+j].im = 0;

            fprintf( fd, "%g ", 
                    sqrt(
                    SQR(fft_of_img[i*N2/2+j].re)
                   +SQR(fft_of_img[i*N2/2+j].im)
                        ) / Npixs
                    );
        }
        fprintf( fd, "\n" );
    }
    fclose( fd );

    rfftwnd_one_complex_to_real( fft_plan_inv, fft_of_img, NULL );

    rfftwnd_destroy_plan( fft_plan );
    rfftwnd_destroy_plan( fft_plan_inv );

    for( i=0; i<Height; i++ )
        for ( j=0; j<Width; j++ )
            Data[i*Width+j] = img[i*N2+j] / Npixs;

    fd = fopen( "tt.dat", "w" );
    for( i=0; i<Height; i++ ) {
        for( j=0; j<Width; j++ )
            fprintf( fd, "%g ", Data[i*Width+j] );
        fprintf( fd, "\n" );
    }
    fclose( fd );

}
*/

void pre_proc( int mode ) {
    
    put_start(mode);
    double sigma[2], mean[2], rms[2];
    int *flag, xi, yi, p;

    if ( All.OnlyFoF ) {
        data_cuting();
        return;
    }

    if ( mode==0 ) {

        if ( All.DataCutting )
            data_cuting();

        if ( All.SigmaClipping ) {
            RSigma = All.RSigma;
            sigma_clipping(0, 0, 0);
        }

        LogNorm = All.LogNorm;
        normalize();

    }

    if ( mode==1) {
        if ( All.SigmaClipping1 ) {
            flag = malloc( sizeof(int)*Npixs );
            memset( flag, 0, sizeof(int)*Npixs );
            p = lset_Head[CurGroup];
            while( p>=0 ) {
                xi = p % WidthGlobal-XShift;
                yi = p / WidthGlobal-YShift;
                /*
                printf( "%i %i [%i] Npixs: %i, Width: %i, Height: %i\n",
                    xi, yi, p, Npixs, Width, Height );
                */
                if ( xi < 0 || xi >= Width || yi < 0 || yi >= Height  )
                    endrun( "" );

                flag[ yi*Width+xi ] = 1;
                p = lset_Next[p];
            }
            writelog( 1, "[%5i], NPixs: %i, Npixs[inner]: %i, NPixs[outer]: %i\n",
                CurGroup, Npixs, lset_Len[CurGroup], Npixs-lset_Len[CurGroup] );

            printf( "\r[%5i], NPixs: %8i, Npixs[inner]: %8i, NPixs[outer]: %8i",
                CurGroup, Npixs, lset_Len[CurGroup], Npixs-lset_Len[CurGroup] );

            get_mean_sigma_rms( Data, flag, Npixs,  mean, sigma, rms );
            free( flag );
            RSigma = All.RSigma1;
            sigma_clipping( sigma[0], mean[0], 1);
        }
    }
    put_end(mode);
}
