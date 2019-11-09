/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

void read_fits( char *fits_fn ) {

    int flag;
    double vmin, vmax;
    put_start;
    fitsfile *ffp;
    char comment[ FLEN_COMMENT ];
#define MAXDIMS 8
    int naxis, status=0, bitpix, anynull;
    long naxes[MAXDIMS], na, start_pos[MAXDIMS], i;

    //printf( "read fits data...\n" );
    for( i=0; i<MAXDIMS; i++ )
        start_pos[i] = 1;

    fits_open_file(&ffp, fits_fn, READONLY, &status);
    fits_get_img_dim(ffp, &naxis, &status);
    fits_get_img_size(ffp, MAXDIMS, naxes, &status);
    fits_get_img_type( ffp, &bitpix, &status );

    //printf( "NAXIS: %i\n", naxis );
    //printf( "NAXES: " );
    for( i=0, na=1; i<naxis; i++ ) {
        na *= naxes[i];
        //printf( "%li ", naxes[i] );
    }
    //printf( "\n" );
    //printf( "BITPIX: %i\n", bitpix );
    //printf( "na: %li\n", na );

    //float *tmp;
    //tmp = (float*)malloc(na*sizeof(float));
    Data = (double*)malloc(na*sizeof(double));


    fits_read_pix(ffp, TDOUBLE, start_pos, na, 0,
            Data, &anynull, &status);
    //fits_read_pix(ffp, TFLOAT, start_pos, na, 0,
    //        tmp, &anynull, &status);
    //for( i=0; i<na; i++ ) {
    //    Data[i] = tmp[i];
    //}

    //myfree(tmp);

    fits_read_key( ffp, TINT, "CRPIX1", &CRPIX1, comment, &status );
    fits_read_key( ffp, TINT, "CRPIX2", &CRPIX2, comment, &status );

    fits_read_key( ffp, TDOUBLE, "CRVAL1", &CRVAL1, comment, &status );
    fits_read_key( ffp, TDOUBLE, "CRVAL2", &CRVAL2, comment, &status );

    fits_read_key( ffp, TDOUBLE, "CDELT1", &CDELT1, comment, &status );
    fits_read_key( ffp, TDOUBLE, "CDELT2", &CDELT2, comment, &status );

    //fits_read_key( ffp, TDOUBLE, "CRVAL3", &FREQ, comment, &status );

    printf1( CRPIX1 );
    printf1( CRPIX2 );

    printf2( CRVAL1 );
    printf2( CRVAL2 );

    printf2( CDELT1 );
    printf2( CDELT2 );

    //FREQ /= 1e6;
    //printf2( FREQ );

    fits_close_file(ffp, &status);
    fits_report_error(stderr, status);

    WidthGlobal = Width = naxes[0];
    HeightGlobal = Height = naxes[1];
    NpixsGlobal = Npixs = naxes[0] * naxes[1];

    //print_data( Data, 230, 290, 580, 640, 0 );

    flag = 0;
    vmin = 1e100;
    vmax = -vmin;
    for( i=0; i<Npixs; i++ ) {
        if ( isnan( Data[i] ) ) {
            flag = 1;
            continue;
        }
        vmin = ( Data[i] < vmin ) ? Data[i] : vmin;
        vmax = ( Data[i] > vmax ) ? Data[i] : vmax;
    }

    printf( "vmin: %g, vmax: %g\n", vmin, vmax );

    if ( flag ) {
        printf( "HAVE NAN: set to vmin\n" );
        for( i=0; i<Npixs; i++ )
            if ( isnan( Data[i] ) )
                Data[i] = vmin;
    }

    DataRaw = (double*)malloc(na*sizeof(double));
    for( i=0; i<Npixs; i++ ) {
        DataRaw[i] = Data[i];
    }

    //print_data( DataRaw, 230, 640, 580, 640, 0 );
    put_end;
}

void free_fits(){
    free( Data );
    free( DataRaw );
}
