/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

void read_fits( char *fits_fn ) {

    fitsfile *ffp;
#define MAXDIMS 8
    int naxis, status=0, i, bitpix, anynull;
    long naxes[MAXDIMS], na, start_pos[MAXDIMS];

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

    Data = (double*)malloc(na*sizeof(double));
    fits_read_pix(ffp, TDOUBLE, start_pos, na, 0,
            Data, &anynull, &status);

    fits_close_file(ffp, &status);
    fits_report_error(stderr, status);

    Width = naxes[0];
    Height = naxes[1];
    Npixs = naxes[0] * naxes[1];

}

void read_file_names() {

    FILE *fd;
    char buf[MYFILENAME_MAX];
    int i;
    FileNum = 0;
    if ( ThisTask == 0 ) {
        fd = fopen( All.FileNameList, "r" );
        if ( NULL == fd )
            endrun( 20190727 );
        while( !feof( fd ) ) {
            fgets( buf, MYFILENAME_MAX, fd );
            FileNum ++;
        }
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
}
