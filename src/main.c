/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

hid_t h5_map, h5_map_after;

void run_first_finder() {

    put_start(0);
    mytimer_start;

    pre_proc(0);

    lset_Next = malloc( sizeof(int) * Npixs );
    lset_Head = malloc( sizeof(int) * Npixs );
    lset_Len = malloc( sizeof(int) * Npixs );

    lset();

    mytimer_end;
    put_end(0);
}

void open_files_for_second_finder() {

    char buf[120];
    int Xs[2];
    hsize_t h5_dims[2];

    h5_map = hdf5_create( "map.hdf5" );
 
    sprintf( buf, "%s/map_after.hdf5", All.OutputDir );
    h5_map_after = hdf5_create( "map_after.hdf5" );

    Xs[0] = HStartCut;
    Xs[1] = WStartCut;
    h5_dims[0] = 2;
    hdf5_write_attr_nd( h5_map, H5T_NATIVE_INT, h5_dims, 1, "CutStart", Xs );
    hdf5_write_attr_nd( h5_map_after, H5T_NATIVE_INT, h5_dims, 1, "CutStart", Xs );

    Xs[0] = HEndCut;
    Xs[1] = WEndCut;
    h5_dims[0] = 2;
    hdf5_write_attr_nd( h5_map, H5T_NATIVE_INT, h5_dims, 1, "CutEnd", Xs );
    hdf5_write_attr_nd( h5_map_after, H5T_NATIVE_INT, h5_dims, 1, "CutEnd", Xs );

    h5_fof = hdf5_create( "fof_regs.hdf5" );

    hdf5_write_attr_scalar( h5_map, H5T_NATIVE_INT, "Ngroup", &lset_Nreg );
    hdf5_write_attr_scalar( h5_map_after, H5T_NATIVE_INT, "Ngroup", &lset_Nreg );
 
}

void close_files_for_second_finder() {
    hdf5_close( h5_map );
    hdf5_close( h5_map_after );
    hdf5_close( h5_fof );
}

void run_second_finder() {

    char buf[100];
    hid_t h5_g;
    hsize_t h5_dims[2];
    int i, p, j, k, Xs[2],
         xmin, xmax, ymin, ymax,
         x, y, w, h;

    put_start(0);
    Nsource = Ngroup = 0;
    open_files_for_second_finder();

    for( k=0; k<lset_Nreg; k++ ) {

        CurGroup = k;
        xmin = ymin = INT_MAX;
        xmax = ymax = -xmin;
        p = lset_Head[k];

        while(p>=0) {
            x = p % WidthGlobal;
            y = p / WidthGlobal;
            xmin = ( x<xmin ) ? x : xmin;
            xmax = ( x>xmax ) ? x : xmax;
            ymin = ( y<ymin ) ? y : ymin;
            ymax = ( y>ymax ) ? y : ymax;
            p = lset_Next[p];
        }

        ymin -= All.SecondFinderPad;
        if ( ymin<=0 )
            ymin = 0;
        ymax += All.SecondFinderPad; 
        if ( ymax > HeightGlobal-1 )
            ymax = HeightGlobal-1;

        xmin -= All.SecondFinderPad;
        if ( xmin<=0 )
            xmin = 0;
        xmax += All.SecondFinderPad; 
        if ( xmax > WidthGlobal-1 )
            xmax = WidthGlobal-1;

        h = ymax - ymin + 1;
        w = xmax - xmin + 1;

        XShift = xmin;
        YShift = ymin;

        Width = w;
        Height = h;
        Npixs = w*h;

        for( i=ymin; i<=ymax; i++ ) {
            for( j=xmin; j<=xmax; j++ ) {
                Data[(i-ymin)*w+(j-xmin)] = DataRaw[i*WidthGlobal+j];
            }
        }

        sprintf( buf, "Group%i", k );
        h5_g = H5Gcreate( h5_map, buf, 0 );

        Xs[0] = ymin+HStartCut;
        Xs[1] = xmin+WStartCut;
        h5_dims[0] = 2;
        hdf5_write_attr_nd( h5_g, H5T_NATIVE_INT, h5_dims, 1, "CRPIX", Xs );

        h5_dims[0] = Height;
        h5_dims[1] = Width;
        hdf5_write_data( h5_g, H5T_NATIVE_DOUBLE, h5_dims, 2, "map", Data );

        H5Gclose( h5_g );

        pre_proc(1);

        sprintf( buf, "Group%i", k );
        h5_g = H5Gcreate( h5_map_after, buf, 0 );

        Xs[0] = ymin+HStartCut;
        Xs[1] = xmin+WStartCut;
        h5_dims[0] = 2;
        hdf5_write_attr_nd( h5_g, H5T_NATIVE_INT, h5_dims, 1, "CRPIX", Xs );

        h5_dims[0] = Height;
        h5_dims[1] = Width;
        hdf5_write_data( h5_g, H5T_NATIVE_DOUBLE, h5_dims, 2, "map", Data );

        group_finder();
        put_sep(1);
    }
    printf( "\n" );

    hdf5_write_attr_scalar( h5_fof, H5T_NATIVE_INT, "Ngroup", &Ngroup );

    writelog( 0, "Ngroup: %i, Nsource: %i\n", Ngroup, Nsource );

    close_files_for_second_finder();
    put_end(0);
}


void run() {

    put_sep(0);

    read_input_fits( All.FileName );
    writelog( 0, "read sourde data ...\n" );
    read_fits_dbl( All.SrcFileName, &SrcData, NULL, NULL );

    if (  All.BkgFittingPolyOrder ) {
        double bkg_params[5];
        bkg_params[0] = All.BkgFittingPolyOrder;
        bkg_params[1] = All.BkgFittingM;
        bkg_params[2] = All.BkgFittingN;
        bkg_params[3] = All.BkgFittingPadding;
        bkg_params[4] = All.BkgRSigmaBeforeFitting;
        Bkg = background( DataRaw, HeightGlobal, WidthGlobal, 0, bkg_params);

        if ( NULL == Bkg )
            endrun( "failed to estimate background!\n" );

        write_fits( "bkg.fits", WidthGlobal, HeightGlobal, Bkg );
    }

    //background_estimation();
    //noise_estimation();

    return;
    run_first_finder();
    put_sep(0);

    if ( !All.DisableSecondFinder ) {
        run_second_finder();
        put_sep(0);
    }

    free( lset_Next );
    free( lset_Head );
    free( lset_Len );

    free_fits();

}

void run_fof_only() {

    read_input_fits( All.FileName );
    put_sep(0);

    pre_proc( 0 );
    XShift = YShift = SigmaClippingVmin = 0;
    Ngroup = 0;

    h5_fof = hdf5_create( "only_fof_regs.hdf5" );

    group_finder();

    hdf5_write_attr_scalar( h5_fof, H5T_NATIVE_INT, "Ngroup", &Ngroup );

    hdf5_close( h5_fof );

    free_fits();

}

int main( int argc, char **argv ) {

    char buf[120];
    memset( sep_str, '-', SEP_LEN-2 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';

    LOG_FILE =  myfopen( "w", "%s.log", argv[1] );

    mytimer_start;
    read_parameters( argv[1] );
    put_sep(0);

    InputBaseName = basename( All.FileName );
    create_dir( All.OutputDir );
    sprintf( buf, "%s/%s", All.OutputDir, InputBaseName );
    sprintf( All.OutputDir, "%s", buf );
    create_dir(buf);

    if ( All.OnlyFoF )
        run_fof_only();
    else 
        run();

    fclose( LOG_FILE );

    mytimer_end;

    return 0;

}

