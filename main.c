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
    hid_t h5_dsp, h5_attr;
    int Xs[2];
    int h5_ndims;
    hsize_t h5_dims[2];


    sprintf( buf, "%s/map.hdf5", All.OutputDir );
    h5_map = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
 
    sprintf( buf, "%s/map_after.hdf5", All.OutputDir );
    h5_map_after = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    Xs[0] = HStartCut;
    Xs[1] = WStartCut;
    h5_ndims = 1;
    h5_dims[0] = 2;
    h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
    h5_attr = H5Acreate( h5_map, "CutStart", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
    H5Aclose( h5_attr );
    h5_attr = H5Acreate( h5_map_after, "CutStart", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp );

    Xs[0] = HEndCut;
    Xs[1] = WEndCut;
    h5_ndims = 1;
    h5_dims[0] = 2;
    h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
    h5_attr = H5Acreate( h5_map, "CutEnd", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
    H5Aclose( h5_attr );
    h5_attr = H5Acreate( h5_map_after, "CutEnd", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp );

    sprintf( buf,"%s/fof_regs.hdf5", All.OutputDir );
    h5_fof = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
 
    h5_dsp = H5Screate( H5S_SCALAR );
    h5_attr = H5Acreate( h5_fof, "Ngroup", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, &lset_Nreg );
    H5Aclose( h5_attr );
    h5_attr = H5Acreate( h5_map, "Ngroup", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, &lset_Nreg );
    H5Aclose( h5_attr );
    h5_attr = H5Acreate( h5_map_after, "Ngroup", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
    H5Awrite( h5_attr, H5T_NATIVE_INT, &lset_Nreg );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp );
 
}

void close_files_for_second_finder() {
    H5Fclose( h5_map );
    H5Fclose( h5_map_after );
    H5Fclose( h5_fof );
}

void run_second_finder() {

    char buf[100];
    hid_t h5_ds, h5_dsp, h5_group, h5_attr;
    int h5_ndims;
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
        

        h = ymax - ymin;
        w = xmax - xmin;

        XShift = xmin;
        YShift = ymin;

        Width = w;
        Height = h;
        Npixs = w * h;

        for( i=ymin; i<ymax; i++ ) {
            for( j=xmin; j<xmax; j++ ) {
                Data[(i-ymin)*w+(j-xmin)] = DataRaw[i*WidthGlobal+j];
            }
        }

        sprintf( buf, "Group%i", k );
        h5_group = H5Gcreate( h5_map, buf, 0 );

        Xs[0] = ymin+HStartCut;
        Xs[1] = xmin+WStartCut;
        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_attr = H5Acreate( h5_group, "CRPIX", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
        H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp );

        h5_ndims = 2;
        h5_dims[0] = Height;
        h5_dims[1] = Width;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_ds = H5Dcreate( h5_group, "map", H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT );
        H5Dwrite( h5_ds, H5T_NATIVE_DOUBLE, h5_dsp, H5S_ALL, H5P_DEFAULT, Data );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp );
        H5Gclose( h5_group );

        pre_proc(1);

        sprintf( buf, "Group%i", k );
        h5_group = H5Gcreate( h5_map_after, buf, 0 );

        Xs[0] = ymin+HStartCut;
        Xs[1] = xmin+WStartCut;
        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_attr = H5Acreate( h5_group, "CRPIX", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
        H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp );

        h5_ndims = 2;
        h5_dims[0] = Height;
        h5_dims[1] = Width;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_ds = H5Dcreate( h5_group, "map", H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT );
        H5Dwrite( h5_ds, H5T_NATIVE_DOUBLE, h5_dsp, H5S_ALL, H5P_DEFAULT, Data );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp );
        H5Gclose( h5_group );

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

