/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"
//#define RUN_DEBUG

int *next_lset0, *head_lset0, *len_lset0, gn;

void run0() {

    int i, N, j, sigs;
    char buf[MYFILENAME_MAX];

    read_file_names();

    sigs = FileNum / 100;
    for( i=0; i<FileNum; i++ ) {

        if ( sigs != 0 && i % sigs == 0 && ThisTask == 0 )
            printf( "%.2f%%\n", ((double)i)/FileNum * 100 );

        if ( i % NTask == ThisTask ) {
            sprintf( FileName, "%s", AllFileNames+i*MYFILENAME_MAX );

            N = strlen( FileName );
            for( j=0; j<N; j++ )
                if ( FileName[j] == '\n' ) {
                    FileName[j] = '\0';
                    break;
                }

            sprintf( buf, "%s/%s", All.InputDir, FileName );
            //printf( "[%i] %s ...\n", ThisTask, buf );
            //fflush( stdout );

            read_fits( buf );
            pre_proc(0);
            //print_data( DataRaw, 230, 290, 580, 640, 0 );
            //return 0;

            lset(0);

            free_fits();
        }

        //MPI_Barrier( MPI_COMM_WORLD );
    }
    free( AllFileNames );
}

void run_first_finder() {

    char buf[100];
    hid_t h5_ds, h5_dsp, h5_attr;
    int h5_ndims, Xs[2];
    hsize_t h5_dims[2];

    pre_proc(0);
    put_sep;

    put_header( "run first finder", 0 );
    mytimer_start();

    if ( ThisTask == 0 ) {
        LsetErrFd = myfopen( "w", "%s/%s_lset0_err.dat", All.OutputDir, InputBaseName );

        sprintf( buf, "%s/%s_lset0_map.hdf5", All.OutputDir, InputBaseName );
        h5_Lset0Map = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        h5_ndims = 2;
        h5_dims[0] = HeightGlobal;
        h5_dims[1] = WidthGlobal;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_ds = H5Dcreate( h5_Lset0Map, "map", H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT );
        H5Dwrite( h5_ds, H5T_NATIVE_DOUBLE, h5_dsp, H5S_ALL, H5P_DEFAULT, DataRaw );
        H5Dclose( h5_ds );
        H5Sclose( h5_dsp );
        H5Fclose( h5_Lset0Map );


        sprintf( buf, "%s/%s_lset0_edges_regs.hdf5", All.OutputDir, InputBaseName );
        h5_EdgesRegs = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        h5_EdgesGroup = H5Gcreate( h5_EdgesRegs, "Edges", 0 );
        h5_RegsGroup =  H5Gcreate( h5_EdgesRegs, "Regs", 0 );

        h5_ndims = 1;
        h5_dims[0] = 2;
        Xs[0] = HStartCut;
        Xs[1] = WStartCut;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_attr = H5Acreate( h5_EdgesRegs, "CRPIX", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
        H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp );

        Xs[0] = HEndCut-HStartCut;
        Xs[1] = WEndCut-WStartCut;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_attr = H5Acreate( h5_EdgesRegs, "NAXIS", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
        H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
        H5Aclose( h5_attr );
        H5Sclose( h5_dsp );

        sprintf( buf, "%s/%s_lset0_lines.hdf5", All.OutputDir, InputBaseName );
        h5_Lines = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        h5_LinesGroup = H5Gcreate( h5_Lines, "Lines", 0 );


    }

    XShift = 0;
    YShift = 0;
    lset(0);

    writelog( "Nfof: %i\n", Nfof );
    next_lset0 = malloc( sizeof(int) * Npixs );
    head_lset0 = malloc( sizeof(int) * Npixs );
    len_lset0 = malloc( sizeof(int) * Npixs );
    put_sep;

    gn = Nfof;
    memcpy( next_lset0, Next, sizeof(int) * Npixs );
    memcpy( head_lset0, Head, sizeof(int) * Npixs );
    memcpy( len_lset0,  Len, sizeof(int) * Npixs );

    gn = 0;
    while( len_lset0[gn]> All.MinEdgeInSecondFinder )
        gn++;

    group_finder_free();

    do_sync;

    if ( ThisTask == 0 ) {
        fclose( LsetErrFd );
        H5Gclose( h5_LinesGroup );
        H5Fclose( h5_Lines );
        H5Gclose( h5_EdgesGroup );
        H5Fclose( h5_EdgesRegs );

    }
    mytimer_end();

}

void merge_map() {

    int k, task, h5_ndims;
    hsize_t h5_dims[2], h5_maxdims[2];
    char buf[100];
    double *data_buf;
    hid_t h5_ds, h5_dsp, h5_g1, h5_g2, h5_attr, h5_map_f1, h5_map_f2,
        h5_map_after_f1, h5_map_after_f2;
    if ( ThisTask )
        return;

    sprintf( buf, "%s/%s_map1.hdf5", All.OutputDir, InputBaseName );
    h5_map_f1 = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    sprintf( buf, "%s/%s_map1_after.hdf5", All.OutputDir, InputBaseName );
    h5_map_after_f1 = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    data_buf = malloc( sizeof(double) * NpixsGlobal );
    for ( task=0; task<NTask; task++ ) {

        sprintf( buf, "%s/%s_map1_%03i.hdf5", All.OutputDir, InputBaseName, task );
        h5_map_f2 = H5Fopen( buf, H5F_ACC_RDWR, H5P_DEFAULT );

        sprintf( buf, "%s/%s_map1_%03i_after.hdf5", All.OutputDir, InputBaseName, task );
        h5_map_after_f2 = H5Fopen( buf, H5F_ACC_RDWR, H5P_DEFAULT );

        for ( k=0; k<gn; k++ ) {
            if ( k % NTask != task )
                continue;

        /************************map**********************/
            sprintf( buf, "Group%i", k );
            h5_g1 = H5Gcreate( h5_map_f1,  buf, 0 );
            h5_g2 = H5Gopen( h5_map_f2,  buf );

            h5_attr = H5Aopen_name( h5_g2, "REPIXS" );
            H5Aread( h5_attr, H5T_NATIVE_INT, data_buf );
            H5Aclose( h5_attr );

            h5_ndims = 1;
            h5_dims[0] = 2;
            h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
            h5_attr = H5Acreate( h5_g1, "REPIXS", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
            H5Awrite( h5_attr, H5T_NATIVE_INT, data_buf );
            H5Aclose( h5_attr );
            H5Sclose( h5_dsp );

            h5_ds = H5Dopen( h5_g2, "map" );
            H5Dread( h5_ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_buf );
            h5_dsp = H5Dget_space( h5_ds );
            H5Dclose( h5_ds );
            H5Sget_simple_extent_dims( h5_dsp, h5_dims, h5_maxdims );
            h5_ndims = H5Sget_simple_extent_ndims( h5_dsp );
            H5Sclose( h5_dsp );

            h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
            h5_ds = H5Dcreate( h5_g1, "map", H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT );
            H5Dwrite( h5_ds, H5T_NATIVE_DOUBLE, h5_dsp, H5S_ALL, H5P_DEFAULT, data_buf );
            H5Dclose( h5_ds );
            H5Sclose( h5_dsp );

            H5Gclose( h5_g1 );
            H5Gclose( h5_g2 );

        /************************map-after**********************/
            sprintf( buf, "Group%i", k );
            h5_g1 = H5Gcreate( h5_map_after_f1,  buf, 0 );
            h5_g2 = H5Gopen( h5_map_after_f2,  buf );

            h5_attr = H5Aopen_name( h5_g2, "REPIXS" );
            H5Aread( h5_attr, H5T_NATIVE_INT, data_buf );
            H5Aclose( h5_attr );

            h5_ndims = 1;
            h5_dims[0] = 2;
            h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
            h5_attr = H5Acreate( h5_g1, "REPIXS", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
            H5Awrite( h5_attr, H5T_NATIVE_INT, data_buf );
            H5Aclose( h5_attr );
            H5Sclose( h5_dsp );

            h5_ds = H5Dopen( h5_g2, "map" );
            H5Dread( h5_ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_buf );
            h5_dsp = H5Dget_space( h5_ds );
            H5Dclose( h5_ds );
            H5Sget_simple_extent_dims( h5_dsp, h5_dims, h5_maxdims );
            h5_ndims = H5Sget_simple_extent_ndims( h5_dsp );
            H5Sclose( h5_dsp );

            h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
            h5_ds = H5Dcreate( h5_g1, "map", H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT );
            H5Dwrite( h5_ds, H5T_NATIVE_DOUBLE, h5_dsp, H5S_ALL, H5P_DEFAULT, data_buf );
            H5Dclose( h5_ds );
            H5Sclose( h5_dsp );

            H5Gclose( h5_g1 );
            H5Gclose( h5_g2 );



        /***************************************************/

        }

        H5Fclose( h5_map_f2 );
        H5Fclose( h5_map_after_f2 );
        sprintf( buf, "%s/%s_map1_%03i.hdf5", All.OutputDir, InputBaseName, task );
        remove( buf );
        sprintf( buf, "%s/%s_map1_%03i_after.hdf5", All.OutputDir, InputBaseName, task );
        remove( buf );

    }
    H5Fclose( h5_map_f1 );
    H5Fclose( h5_map_after_f1 );
    free( data_buf );

}

void merge_lset1() {

    int k, task, h5_ndims, i, iters;
    hsize_t h5_dims[2], h5_maxdims[2];
    char buf[100];
    double *data_buf;
    hid_t h5_ds, h5_dsp, h5_g1, h5_g2, h5_attr,
        h5_lines_f1, h5_lines_f2;
    if ( ThisTask )
        return;

    sprintf( buf, "%s/%s_lset1_lines.hdf5", All.OutputDir, InputBaseName );
    h5_lines_f1 = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    data_buf = malloc( sizeof(double) * NpixsGlobal );
    for ( task=0; task<NTask; task++ ) {

        sprintf( buf, "%s/%s_lset1_lines_%03i.hdf5", All.OutputDir, InputBaseName, task );
        h5_lines_f2 = H5Fopen( buf, H5F_ACC_RDWR, H5P_DEFAULT );

        for ( k=0; k<gn; k++ ) {
            if ( k % NTask != task )
                continue;

        /************************lines**********************/
            sprintf( buf, "Group%i", k );
            h5_g1 = H5Gcreate( h5_lines_f1,  buf, 0 );
            h5_g2 = H5Gopen( h5_lines_f2,  buf );

            h5_attr = H5Aopen_name( h5_g2, "iters" );
            H5Aread( h5_attr, H5T_NATIVE_INT, &iters );
            H5Aclose( h5_attr );

            if ( iters > 0 ) {

                h5_dsp = H5Screate( H5S_SCALAR );
                h5_attr = H5Acreate( h5_g1, "iters", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
                H5Awrite( h5_attr, H5T_NATIVE_INT, &iters );
                H5Aclose( h5_attr );
                H5Sclose( h5_dsp );

            }

            for ( i=1; i<=iters; i++ ) {

                sprintf( buf, "S1S2-%i", i );

                h5_ndims = 1;
                h5_dims[0] = 2;
                h5_attr = H5Aopen_name( h5_g2, buf );
                H5Aread( h5_attr, H5T_NATIVE_DOUBLE, data_buf );
                H5Aclose( h5_attr );

                h5_ndims = 1;
                h5_dims[0] = 2;
                h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL  );
                h5_attr = H5Acreate( h5_g1, buf, H5T_NATIVE_DOUBLE, h5_dsp, H5P_DEFAULT );
                H5Awrite( h5_attr, H5T_NATIVE_DOUBLE, data_buf );
                H5Aclose( h5_attr );
                H5Sclose( h5_dsp );

                sprintf( buf, "lines-%i", i );

                h5_ds = H5Dopen( h5_g2, buf );
                H5Dread( h5_ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_buf );
                h5_dsp = H5Dget_space( h5_ds );
                H5Dclose( h5_ds );
                H5Sget_simple_extent_dims( h5_dsp, h5_dims, h5_maxdims );
                h5_ndims = H5Sget_simple_extent_ndims( h5_dsp );
                H5Sclose( h5_dsp );

                h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
                h5_ds = H5Dcreate( h5_g1, buf, H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
                H5Dwrite( h5_ds, H5T_NATIVE_INT, h5_dsp, H5S_ALL, H5P_DEFAULT, data_buf );
                H5Dclose( h5_ds );
                H5Sclose( h5_dsp );

            }

            H5Gclose( h5_g1 );
            H5Gclose( h5_g2 );
        /***************************************************/

        }

        H5Fclose( h5_lines_f2 );
        sprintf( buf, "%s/%s_lset1_lines_%03i.hdf5", All.OutputDir, InputBaseName, task );
        remove( buf );

    }
    H5Fclose( h5_lines_f1 );
    free( data_buf );

}

void run_second_finder() {

    char buf[100];
    hid_t h5_ds, h5_dsp, h5_group, h5_attr;
    int h5_ndims;
    hsize_t h5_dims[2];

    int i, p, j, k, Xs[2],
         xmin, xmax, ymin, ymax, x, y, w, h, flag, index;
    put_header( "run second finder", 0 );
    do_sync;

//#define FIXEDSIZE
#ifdef FIXEDSIZE
    int wh;
#endif
#ifdef ADD_PAD
    int lmin;
#endif

#ifdef ADD_PAD
    lmin = 50;
#endif

    for( k=0, index=0; k<gn; k++ ) {
        flag = 0;
        /*
        if ( len_lset0[k] < 30 ) 
            continue;
        */
        p = head_lset0[k];
        while( p>=0 ) {
            x = p % WidthGlobal;
            y = p / WidthGlobal;
            p = next_lset0[p];
            if ( x == 0 || y == 0 ) {
                flag = 1;
                break;
            }
        }
        if ( flag == 0 ) {
            head_lset0[index] = head_lset0[k];
            len_lset0[index] = len_lset0[k]; 
            index++;
        }
    }

    gn = index;
    writelog( "gn: %i\n", gn );

//#define DISABLE_FIRST_FINDER
#ifdef DISABLE_FIRST_FINDER
    gn = 1;
#endif

    if ( NTask != 1 ) {


        sprintf( buf, "%s/%s_map1_%03i.hdf5", All.OutputDir, InputBaseName, ThisTask );
        h5_Lset1Map = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

        sprintf( buf, "%s/%s_map1_%03i_after.hdf5", All.OutputDir, InputBaseName, ThisTask );
        h5_Lset1Map_after = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

        if ( All.Lset1  ) {
            LsetErrFd = myfopen( "w", "%s/%s_lset1_err_%03i.dat",
                All.OutputDir, InputBaseName, ThisTask );
            sprintf( buf,"%s/%s_lset1_lines_%03i.hdf5",
                All.OutputDir, InputBaseName, ThisTask );
            h5_Lines = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

            sprintf( buf,"%s/%s_lset1_edge_%03i.hdf5",
                All.OutputDir, InputBaseName, ThisTask );
            h5_EdgesRegs = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        }
        else {

            sprintf( buf,"%s/%s_fof_regs_%03i.hdf5",
                All.OutputDir, InputBaseName, ThisTask );
            h5_Regs = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        }

    }
    else {

        sprintf( buf, "%s/%s_map1.hdf5", All.OutputDir, InputBaseName );
        h5_Lset1Map = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

        sprintf( buf, "%s/%s_map1_after.hdf5", All.OutputDir, InputBaseName );
        h5_Lset1Map_after = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        
        if ( All.Lset1 ) {
            LsetErrFd = myfopen( "w", "%s/%s_lset1_err.dat",
                All.OutputDir, InputBaseName );
            sprintf( buf,"%s/%s_lset1_edge.hdf5",
                All.OutputDir, InputBaseName );
            h5_EdgesRegs = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

            sprintf( buf,"%s/%s_lset1_lines.hdf5",
                All.OutputDir, InputBaseName );
            h5_Lines = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        }
        else {
            sprintf( buf,"%s/%s_fof_regs.hdf5",
                All.OutputDir, InputBaseName );
            h5_Regs = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        }

    }

            h5_dsp = H5Screate( H5S_SCALAR );
            h5_attr = H5Acreate( h5_Regs, "NGroups", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
            H5Awrite( h5_attr, H5T_NATIVE_INT, &gn );
            H5Aclose( h5_attr );
            H5Sclose( h5_dsp );

            h5_ndims = 1;
            h5_dims[0] = 2;
            Xs[0] = HStartCut;
            Xs[1] = WStartCut;
            h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
            h5_attr = H5Acreate( h5_Regs, "GCRPIX", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
            H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
            H5Aclose( h5_attr );
            H5Sclose( h5_dsp );

            Xs[0] = HEndCut-HStartCut;
            Xs[1] = WEndCut-WStartCut;
            h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
            h5_attr = H5Acreate( h5_Regs, "GNAXIS", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
            H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
            H5Aclose( h5_attr );
            H5Sclose( h5_dsp );

#ifdef RUN_DEBUG
    for( k=0; k<NTask; k++ )
#else
    for( k=0; k<gn; k++ )
#endif
    {
        if ( k % NTask != ThisTask )
            continue;

        CurGroup = k;
        xmin = ymin = INT_MAX;
        xmax = ymax = -xmin;
        p = head_lset0[k];
        flag = 0;
        while(p>=0) {
            x = p % WidthGlobal;
            y = p / WidthGlobal;
            xmin = ( x<xmin ) ? x : xmin;
            xmax = ( x>xmax ) ? x : xmax;
            ymin = ( y<ymin ) ? y : ymin;
            ymax = ( y>ymax ) ? y : ymax;
            p = next_lset0[p];
        }

        if ( NTask > 1 )
        printf( "task: %03i, group %05i [%05i] "
                "region: (%i, %i) - (%i, %i)\n",
                 ThisTask, k,
                 len_lset0[k], xmin, ymin, xmax, ymax );
        
#ifdef DISABLE_FIRST_FINDER
        ymax = HeightGlobal;
        ymin = 0;
        xmax = WidthGlobal;
        xmin = 0;
#endif

        h = ymax - ymin;
        w = xmax - xmin;
#ifdef ADD_PAD
        if ( h < lmin ) {
            ymin = ymin - (lmin-h) / 2;
            h = lmin;
            ymax = ymin + h;
        }
        if ( w < lmin ) {
            xmin = xmin - (lmin-w) / 2;
            w = lmin;
            xmax = xmin + w;
        }
#endif

        XShift = xmin;
        YShift = ymin;

#ifdef FIXEDSIZE
        if ( w > h )
            wh = w;
        else
            wh = h;

        xmax = xmin + wh;
        ymax = ymin + wh;
        Width = wh;
        Height = wh;
        Npixs = wh * wh;
#else
        Width = w;
        Height = h;
        Npixs = w * h;
#endif
        for( i=ymin; i<ymax; i++ ) {
            for( j=xmin; j<xmax; j++ ) {
#ifdef FIXEDSIZE
                Data[(i-ymin)*wh+(j-xmin)] = DataRaw[i*WidthGlobal+j];
#else
                Data[(i-ymin)*w+(j-xmin)] = DataRaw[i*WidthGlobal+j];
#endif
            }
        }

        Xs[0] = ymin;
        Xs[1] = xmin;

        sprintf( buf, "Group%i", k );
        h5_group = H5Gcreate( h5_Lset1Map, buf, 0 );

        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_attr = H5Acreate( h5_group, "REPIXS", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
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
        h5_group = H5Gcreate( h5_Lset1Map_after, buf, 0 );

        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dsp = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_attr = H5Acreate( h5_group, "REPIXS", H5T_NATIVE_INT, h5_dsp, H5P_DEFAULT );
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

        if ( All.Lset1 ) {
            fprintf( LsetErrFd,   "Group: %03i\n", k);
            sprintf( buf, "Group%i", k );
            h5_LinesGroup = H5Gcreate( h5_Lines, buf, 0 );
            sprintf( buf, "Group%i", k );
            h5_EdgesGroup = H5Gcreate( h5_EdgesRegs, buf, 0 );
            lset(1);
            H5Gclose( h5_LinesGroup );
            H5Gclose( h5_EdgesGroup );
            fclose( LsetErrFd );
        }
        else {

            sprintf( buf, "Group%i", k );
            h5_RegsGroup = H5Gcreate( h5_Regs, buf, 0 );
            group_finder();
            H5Gclose( h5_RegsGroup );

        }

    }


    H5Fclose( h5_Lset1Map );
    H5Fclose( h5_Lset1Map_after );

    if ( All.Lset1 ) {
        H5Fclose( h5_Lines );
        H5Fclose( h5_EdgesRegs );
    }
    else {
        H5Fclose( h5_Regs );
    }

    do_sync;
    if ( NTask > 1 ) {
        merge_map();
        if ( All.Lset1 )
            merge_lset1();
    }
    put_end();
}


void run() {

    read_fits( All.FileName );
    put_sep;

    next_lset0 = malloc( sizeof(int) * Npixs );
    head_lset0 = malloc( sizeof(int) * Npixs );
    len_lset0 = malloc( sizeof(int) * Npixs );

    run_first_finder();
    run_second_finder();

    free( next_lset0 );
    free( head_lset0 );
    free( len_lset0 );

    free_fits();

}

void global_init() {
    put_header( "global init", 0 );
    InputBaseName = basename( All.FileName );
    XShift = YShift = 0;
    create_dir( All.OutputDir );
    create_dir( All.PhiDir );
    put_end();

}
void global_free() {
}

void pipeline_test() {
}

int main( int argc, char **argv ) {

    char *bname;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NTask );
    
    bname = basename( argv[1] );
    if ( ThisTask == 0 ) {
        if ( access( "fgext.log", 0 ) == -1 ) {
            printf( "`fgext.log` is created by Task 0\n" );
            if ( mkdir( "fgext.log", 0755 ) == -1 )
                endrun( "failed create directory `fgext.log`.\n"  );
        }
    }
    memset( sep_str, '-', SEP_LEN-2 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';

    do_sync;

    LogFileFd = myfopen( "w", "./fgext.log/%s_%03i.log", bname, ThisTask );
    put_sep;

#ifdef ZDEBUG
    writelog( "Assign `SIGSEGV` to signal hander function.\n" );
    signal( SIGSEGV, signal_hander );
    init_zsig();
    put_sep;
#endif

    read_parameters( argv[1] );
    put_sep;

    if ( All.Lset1 == 0 && NTask > 1 )
        endrun( "paralle is not supported for Lset1=0!" );

    global_init();
    put_sep;

//#define PIPELINE_TEST
#ifdef PIPELINE_TEST
    pipeline_test();
#else
    switch (All.ParalleLevel) {
        case 0:
            run0();
            break;
        case 1:
            run();
            break;
        default:
            printf( "Unsupported paralle level.\n" );
            endrun( "" );
    } 
    put_sep;
#endif

    global_free();
    fclose( LogFileFd );
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;

}

