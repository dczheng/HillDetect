/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"
#include "drfftw.h"

void sigma_clipping() {

    long i;
    double sigma, mu, vmin;

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
    vmin = mu + All.RSigma * sigma;  

    for( i=0; i<Npixs; i++ ) {
        if ( Data[i] < vmin )
            Data[i] = vmin; 
    }

}

void normalize() {

    double vmin, vmax, dv;
    long i;

    vmax = -1e100;
    vmin = 1e100;

    for( i=0; i<Npixs; i++ ) {

        if ( isnan( Data[i] ) )
            continue;

        if ( All.LogNorm ) {
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

        if ( All.LogNorm )
            if ( Data[i] < vmin )
                Data[i] = vmin;
    }

    if ( All.LogNorm )
        dv =log( vmax / vmin );
    else
        dv = vmax - vmin;

    for ( i=0; i<Npixs; i++ ) {
        if ( All.LogNorm )
            Data[i] = log( Data[i]/vmin ) / dv;
        else
            Data[i] = ( Data[i] - vmin ) / dv;
    }

}

void ft_clipping(){

    int i, j, N2;
    FILE *fd;
    double fac;
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
            img[i*N2+j] = log(Data[i*Width+j]);

    fft_plan     = rfftw2d_create_plan( Height, Width, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE );
    fft_plan_inv = rfftw2d_create_plan( Height, Width, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE );

    rfftwnd_one_real_to_complex( fft_plan, img, NULL );

    fac = 0.05;
    fd = fopen( "fft.dat", "w" );
    for( i=0; i<Height; i++ ) {
        for( j=0; j<N2/2; j++ ) {
            if ( i < Height*fac && j < N2/2*fac ) 
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
            fprintf( fd, "%g ", exp(Data[i*Width+j]) );
        fprintf( fd, "\n" );
    }
    fclose( fd );

    endrun(0);

    
}

void pre_proc() {
    
    if ( All.FTClipping )
        ft_clipping();

    if ( All.SigmaClipping )
        sigma_clipping();


    normalize();
}
