/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"
#include "drfftw.h"

void sigma_clipping( int mode ) {

    int i;
    double sigma, mu, vmin, r;

    if ( mode==0 ) {
        if ( All.SigmaClipping == 0 )
            return;
        else 
            r = All.RSigma;
    }

    if ( mode==1) {
        if ( All.SigmaClipping1 == 0 )
            return;
        else
            r = All.RSigma1;
    }

    vmin = 1e100;
    for ( i=0; i<Npixs; i++ ) {
        if ( isnan( Data[i] ) )
            continue;
        vmin = ( vmin < Data[i] ) ? vmin : Data[i];
    }

    for( i=0, mu=0; i<Npixs; i++ ) {
        if ( isnan( Data[i] ) )
            Data[i] = vmin;

        mu += Data[i];
    }

    mu /= Npixs;

    for( i=0, sigma=0; i<Npixs; i++ ) {
        sigma += SQR( Data[i] - mu );
    }
    sigma = sqrt( sigma / Npixs );
    vmin = mu + r * sigma;  

    for( i=0; i<Npixs; i++ ) {
        if ( Data[i] < vmin )
            Data[i] = vmin; 
    }

}

void data_cuting( int mode ) {

    int i, j;
    int index;

    if ( mode == 1 )
        return;

    put_header( "data cutting" );
    HStartCut = Height * All.CuttingYStart;
    HEndCut = Height * All.CuttingYEnd;
    WStartCut = Width * All.CuttingXStart;
    WEndCut = Width * All.CuttingXEnd;
    //printf( "%i %i %i %i\n", Height, Width, h, w  );

    if ( HEndCut<=HStartCut || WEndCut<=WStartCut ) {
        endrun( "Invalid cutting parameters." );
    }

    writelog( "Height: %i, Width: %i\n", Height, Width );
    writelog( "region: (%g, %g), (%g, %g)\n",
                All.CuttingYStart, All.CuttingYEnd, 
                All.CuttingXStart, All.CuttingXEnd );
    writelog( "region: (%i, %i), (%i, %i)\n",
            HStartCut, HEndCut, WStartCut, WEndCut );

    //output_data( "before_cutting.dat" );

    for( i=HStartCut, index=0; i<HEndCut; i++ )
        for ( j=WStartCut; j<WEndCut; j++, index++ ) {
            Data[index] = Data[i*Width+j];
            DataRaw[index] = DataRaw[i*Width+j];
       }

    HeightGlobal = Height = HEndCut - HStartCut;
    WidthGlobal = Width = WEndCut - WStartCut;
    NpixsGlobal = Npixs = Height * Width;

    put_end();
    //output_data( "after_cutting.dat" );
    //endrun( "test" );
}

void normalize( int mode ) {

    double vmin, vmax, dv, logflag;
    int i;

    vmax = -1e100;
    vmin = 1e100;

    if ( mode == 0 )
        logflag = All.LogNorm;
    if ( mode == 1 )
        logflag = All.LogNorm1;

    for( i=0; i<Npixs; i++ ) {

        if ( isnan( Data[i] ) )
            continue;

        if ( logflag ) {
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

        if ( logflag )
            if ( Data[i] < vmin )
                Data[i] = vmin;
    }

    if ( logflag )
        dv =log( vmax / vmin );
    else
        dv = vmax - vmin;

    for ( i=0; i<Npixs; i++ ) {
        if ( logflag )
            Data[i] = log( Data[i]/vmin ) / dv;
        else
            Data[i] = ( Data[i] - vmin ) / dv;
    }

}

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

void pre_proc( int mode ) {
    
    put_header( "pre proc" );
    if ( All.FTClipping )
        ft_clipping( mode );

    sigma_clipping( mode );

    if ( All.DataCutting )
        data_cuting( mode );

    normalize( mode );
    put_end();
}
