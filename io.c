/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

void read_fits( char *fits_fn ) {

    put_header( "read fits", 0 );
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

    fits_read_key( ffp, TDOUBLE, "CRVAL3", &FREQ, comment, &status );

    writelog1( CRPIX1 );
    writelog1( CRPIX2 );

    writelog2( CRVAL1 );
    writelog2( CRVAL2 );

    writelog2( CDELT1 );
    writelog2( CDELT2 );

    FREQ /= 1e6;
    writelog2( FREQ );

    fits_close_file(ffp, &status);
    fits_report_error(stderr, status);

    WidthGlobal = Width = naxes[0];
    HeightGlobal = Height = naxes[1];
    NpixsGlobal = Npixs = naxes[0] * naxes[1];

    //print_data( Data, 230, 290, 580, 640, 0 );

    DataRaw = (double*)malloc(na*sizeof(double));
    for( i=0; i<Npixs; i++ )
        DataRaw[i] = Data[i];

    //print_data( DataRaw, 230, 640, 580, 640, 0 );

    put_end();
}

void read_file_names() {

    FILE *fd;
    char buf[MYFILENAME_MAX];
    int i;
    FileNum = 0;
    if ( ThisTask == 0 ) {
        fd = fopen( All.FileNameList, "r" );
        if ( NULL == fd )
            endrun("");
        do {
            fgets( buf, MYFILENAME_MAX, fd );
            FileNum ++;
        } while( !feof( fd ) );
        FileNum -= 1;
        fclose( fd );
    }

    MPI_Bcast( &FileNum, 1, MPI_INT, 0, MPI_COMM_WORLD );

    /*
    //printf( "Task: %i, FileNum: %i\n", ThisTask, FileNum );
    fflush( stdout );
    MPI_Barrier( MPI_COMM_WORLD );
    */

    AllFileNames = malloc( FileNum * MYFILENAME_MAX );

    if ( ThisTask == 0 ) {
        fd = fopen( All.FileNameList, "r" );

        for( i=0; i<FileNum; i++ )
            fgets( AllFileNames + i*MYFILENAME_MAX,  MYFILENAME_MAX, fd );

        fclose( fd );
    }

    MPI_Bcast( AllFileNames, FileNum*MYFILENAME_MAX, MPI_BYTE, 0, MPI_COMM_WORLD );

}

void free_fits(){
    free( Data );
    free( DataRaw );
}

void output_map( FILE *fd, int W, int H, double *D, int *Xs, int *Ys ) {
    int i, j, xmin, ymin, xmax, ymax, ww, hh;
    (void)hh;

    if ( Xs == NULL ) {
        xmin = 0;
        xmax = W;
    }
    else {
        xmin = Xs[0];
        xmax = Xs[1];
    }

    if ( Ys == NULL ) {
        ymin = 0;
        ymax = H;
    }
    else {
        ymin = Ys[0];
        ymax = Ys[1];
    }

    ww = xmax - xmin;
    hh = ymax - ymin;

    fprintf( fd, "%i %i %i %i ", xmin, ymin, xmax, ymax );
    for( i=0; i<ww-4; i++ )
        fprintf( fd, "0 " );
    fprintf( fd, "\n" );

    for( i=ymin; i<ymax; i++ ) {
        for( j=xmin; j<xmax; j++ ) {
            fprintf( fd, "%g ", D[i*W+j] );
        }
        fprintf( fd, "\n" );
    }
}
