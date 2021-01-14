/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

#define printf1( a ) writelog( 0, "%-20s: %i\n", #a, a )
#define printf2( a ) writelog( 0, "%-20s: %g\n", #a, a )
#define printf3( a, b ) writelog( 0, "%-20s: %s\n", a, b )

void write_fits( char *fn, int w, int h, double *data ) {
    fitsfile *ffp; 
    long naxis[2];
    int status;
    status = 0;
    char buf[300];

    sprintf( buf, "%s/%s", All.OutputDir, fn );
    fn = buf;
    if ( access( fn, 0 ) != -1 ) {
        remove( fn );
    }
    naxis[0] = w;
    naxis[1] = h;
    fits_create_file( &ffp, fn, &status  );
    fits_create_img( ffp, DOUBLE_IMG, 2, naxis, &status );
    fits_write_img( ffp, TDOUBLE, 1, (LONGLONG)(naxis[0]*naxis[1]), data, &status  );
    fits_close_file( ffp, &status  );
    fits_report_error( stderr, status  );
}

void read_fits_int( char *fits_fn, int **data, int *W, int *H ) {
    char buf[20];
    fitsfile *ffp;
#define MAXDIMS 8
    int naxis, status=0, bitpix, anynull;
    long naxes[MAXDIMS], na, start_pos[MAXDIMS], i;

    for( i=0; i<MAXDIMS; i++ )
        start_pos[i] = 1;

    fits_open_file(&ffp, fits_fn, READONLY, &status);
    fits_get_img_dim(ffp, &naxis, &status);
    fits_get_img_size(ffp, MAXDIMS, naxes, &status);
    fits_get_img_type( ffp, &bitpix, &status );

    for( i=0, na=1; i<naxis; i++ ) {
        na *= naxes[i];
    }

    *data = (int*)malloc(na*sizeof(int));
    fits_read_pix(ffp, TINT, start_pos, na, 0,
            *data, &anynull, &status);
    fits_close_file(ffp, &status);
    fits_report_error(stderr, status);

    if ( W != NULL && H != NULL ) {
        *W = naxes[0];
        *H = naxes[1];
    }

    sprintf( buf, "%li x %li", naxes[0], naxes[1] );
    printf3( "Dim", buf );

}
void read_fits_dbl( char *fits_fn, double **data, int *W, int *H ) {
    char buf[20];
    fitsfile *ffp;
#define MAXDIMS 8
    int naxis, status=0, bitpix, anynull;
    long naxes[MAXDIMS], na, start_pos[MAXDIMS], i;

    for( i=0; i<MAXDIMS; i++ )
        start_pos[i] = 1;

    fits_open_file(&ffp, fits_fn, READONLY, &status);
    fits_get_img_dim(ffp, &naxis, &status);
    fits_get_img_size(ffp, MAXDIMS, naxes, &status);
    fits_get_img_type( ffp, &bitpix, &status );

    for( i=0, na=1; i<naxis; i++ ) {
        na *= naxes[i];
    }

    *data = (double*)malloc(na*sizeof(double));
    fits_read_pix(ffp, TDOUBLE, start_pos, na, 0,
            *data, &anynull, &status);
    fits_close_file(ffp, &status);
    fits_report_error(stderr, status);

    if ( W != NULL && H != NULL ) {
        *W = naxes[0];
        *H = naxes[1];
    }

    sprintf( buf, "%li x %li", naxes[0], naxes[1] );
    printf3( "Dim", buf );

}

void read_input_fits( char *fits_fn ) {
    int i;
    fitsfile *ffp;
    int status=0;
    char comment[ FLEN_COMMENT ];

    put_start(0);

    read_fits_dbl( fits_fn, &Data, &Width, &Height );

    fits_open_file(&ffp, fits_fn, READONLY, &status);
    fits_read_key( ffp, TDOUBLE, "CDELT1", &CDELT1, comment, &status  );
    fits_read_key( ffp, TDOUBLE, "CDELT2", &CDELT2, comment, &status  );
    fits_close_file(ffp, &status);
    fits_report_error(stderr, status);

    WidthGlobal = Width;
    HeightGlobal = Height;
    NpixsGlobal = Npixs = Width * Height;

    DataRaw = (double*)malloc(Npixs*sizeof(double));
    for( i=0; i<Npixs; i++ ) {
        DataRaw[i] = Data[i];
    }
    find_vmin_vmax( Data, Npixs, &DataMin, &DataMax );

    VInvalid = DataMin / All.VInvalid;
    writelog( 0, "DataMin: %g, DataMax: %g, VInvalid: %g\n", DataMin, DataMax, VInvalid );
    DataRawMin = DataMin;
    DataRawMax = DataMax;
    writelog( 0, "DataRawMin: %g, DataRawMax: %g\n", DataRawMin, DataRawMax );
    put_end(0);
}

void free_fits(){
    free( Data );
    free( DataRaw );
}
