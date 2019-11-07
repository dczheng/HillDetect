/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

hid_t h5_map, h5_map_after;

void run_first_finder() {

    put_start;
    mytimer_start;

    pre_proc(0);

    lset_Next = malloc( sizeof(int) * Npixs );
    lset_Head = malloc( sizeof(int) * Npixs );
    lset_Len = malloc( sizeof(int) * Npixs );

    lset();

    mytimer_end;
    put_end;
}

void open_files_for_second_finder() {

    char buf[100];
    int Xs[2];
    hsize_t h5_dims[2];


    sprintf( buf, "%s/map.hdf5", All.OutputDir );
    h5_map = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
 
    sprintf( buf, "%s/map_after.hdf5", All.OutputDir );
    h5_map_after = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

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


    sprintf( buf,"%s/fof_regs.hdf5", All.OutputDir );
    h5_fof = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    hdf5_write_attr_scalar( h5_map, H5T_NATIVE_INT, "Ngroup", &lset_Nreg );
    hdf5_write_attr_scalar( h5_map_after, H5T_NATIVE_INT, "Ngroup", &lset_Nreg );
    hdf5_write_attr_scalar( h5_fof, H5T_NATIVE_INT, "Ngroup", &lset_Nreg );
 
}

void close_files_for_second_finder() {
    H5Fclose( h5_map );
    H5Fclose( h5_map_after );
    H5Fclose( h5_fof );
}

void run_second_finder() {

    char buf[100];
    hid_t h5_g;
    hsize_t h5_dims[2];

    int i, p, j, k, Xs[2],
         xmin, xmax, ymin, ymax, x, y, w, h;
    put_start;

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

        printf( "group %05i [%05i] "
                "region: (%i, %i) - (%i, %i)\n",
                 k, lset_Len[k],
                 xmin + WStartCut,
                 ymin + HStartCut,
                 xmax + WStartCut,
                 ymax + HStartCut );
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
    }

    close_files_for_second_finder();
    put_end;

}


void run() {

    read_fits( All.FileName );
    put_sep;

    run_first_finder();
    put_sep;

    if ( !All.DisableSecondFinder ) {
        run_second_finder();
        put_sep;
    }

    free( lset_Next );
    free( lset_Head );
    free( lset_Len );

    free_fits();

}

int main( int argc, char **argv ) {

    memset( sep_str, '-', SEP_LEN-2 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';

    put_sep;

    read_parameters( argv[1] );
    InputBaseName = basename( All.FileName );
    create_dir( All.OutputDir );
    sprintf( All.OutputDir, "%s/%s", All.OutputDir, InputBaseName );
    create_dir( All.OutputDir );

    run();

    return 0;

}

