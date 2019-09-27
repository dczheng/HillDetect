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
            pre_proc(1);
            //print_data( DataRaw, 230, 290, 580, 640, 0 );
            //return 0;

            lset(1);

            free_fits();
        }

        //MPI_Barrier( MPI_COMM_WORLD );
    }
    free( AllFileNames );
}

void run_first_finder() {

    char buf[100];
    hid_t h5_dataset, h5_dataspace;
    int h5_ndims;
    hsize_t h5_dims[2];

    pre_proc(0);
    put_sep;

    put_header( "run first finder", 0 );
    mytimer_start();

    if ( ThisTask == 0 ) {
        sprintf( buf, "%s/%s_lset0_map.hdf5", All.OutputDir, InputBaseName );
        h5_Lset0Map = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        h5_ndims = 2;
        h5_dims[0] = HeightGlobal;
        h5_dims[1] = WidthGlobal;
        h5_dataspace = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_dataset = H5Dcreate( h5_Lset0Map, "map", H5T_NATIVE_DOUBLE, h5_dataspace, H5P_DEFAULT );
        H5Dwrite( h5_dataset, H5T_NATIVE_DOUBLE, h5_dataspace, H5S_ALL, H5P_DEFAULT, DataRaw );
        H5Dclose( h5_dataset );
        H5Sclose( h5_dataspace );
        H5Fclose( h5_Lset0Map );
    }

    if ( ThisTask == 0 ) {
        LsetErrFd = myfopen( "w", "%s/%s_lset0_err.dat", All.OutputDir, InputBaseName );
        EdgesFd = myfopen( "w","%s/%s_lset0_edges.dat", All.OutputDir, InputBaseName );

        sprintf( buf, "%s/%s_lset0_lines.hdf5", All.OutputDir, InputBaseName );
        h5_Lines = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        h5_LinesGroup = H5Gcreate( h5_Lines, "Lines", 0 );

    }

    XShift = 0;
    YShift = 0;
    lset(0);

    writelog( "NfofEdge: %i\n", NfofEdge );
    next_lset0 = malloc( sizeof(int) * Npixs );
    head_lset0 = malloc( sizeof(int) * Npixs );
    len_lset0 = malloc( sizeof(int) * Npixs );
    put_sep;

    gn = NfofEdge;
    memcpy( next_lset0, Next, sizeof(int) * Npixs );
    memcpy( head_lset0, Head, sizeof(int) * Npixs );
    memcpy( len_lset0,  Len, sizeof(int) * Npixs );
    find_region_free();

    do_sync;

    if ( ThisTask == 0 ) {
        fclose( LsetErrFd );
        fclose( EdgesFd );
        H5Gclose( h5_LinesGroup );
        H5Fclose( h5_Lines );

    }
    mytimer_end();

}

void run_second_finder() {

    char buf[100];
    hid_t h5_dataset, h5_dataspace, h5_group, h5_attr;
    int h5_ndims;
    hsize_t h5_dims[2];

    int i, p, j, k, Xs[2],
         xmin, xmax, ymin, ymax, x, y, w, h, flag, index;
    put_header( "run second finder", 0 );

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

    LsetErrFd = myfopen( "w", "%s/%s_lset1_err_%03i.dat",
                All.OutputDir, InputBaseName, ThisTask );
    EdgesFd = myfopen( "w","%s/%s_lset1_edges_%03i.dat",
                All.OutputDir, InputBaseName, ThisTask );
    RegsFd = myfopen( "w","%s/%s_lset1_regs_%03i.dat",
                All.OutputDir, InputBaseName, ThisTask );

    sprintf( buf,"%s/%s_lset1_lines_%03i.hdf5",
                All.OutputDir, InputBaseName, ThisTask );
    h5_Lines = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    sprintf( buf, "%s/%s_lset1_%03i.hdf5", All.OutputDir, InputBaseName, ThisTask );
    h5_Lset1Map = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

#ifdef RUN_DEBUG
    for( k=0; k<NTask; k++ )
#else
    for( k=0; k<gn; k++ )
#endif
    {
        if ( k % NTask != ThisTask )
            continue;
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


        printf( "task: %03i, group %05i [%05i] "
                "region: (%i, %i) - (%i, %i)\n",
                 ThisTask, k,
                 len_lset0[k], xmin, ymin, xmax, ymax );
        
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
        Xs[0] = xmin;
        Xs[1] = ymin;

        fprintf( LsetErrFd,   "Group: %03i\n", k);
        fprintf( EdgesFd,     "Group: %03i\n", k);
        fprintf( RegsFd,      "Group: %03i\n", k);

        sprintf( buf, "Group%i", k );
        h5_group = H5Gcreate( h5_Lset1Map, buf, 0 );

        sprintf( buf, "Group%i", k );
        h5_LinesGroup = H5Gcreate( h5_Lines, buf, 0 );

        h5_ndims = 1;
        h5_dims[0] = 2;
        h5_dataspace = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_attr = H5Acreate( h5_group, "REPIXS", H5T_NATIVE_INT, h5_dataspace, H5P_DEFAULT );
        H5Awrite( h5_attr, H5T_NATIVE_INT, Xs );
        H5Aclose( h5_attr );
        H5Sclose( h5_dataspace );

        h5_ndims = 2;
        h5_dims[0] = Height;
        h5_dims[1] = Width;
        h5_dataspace = H5Screate_simple( h5_ndims, h5_dims, NULL );
        h5_dataset = H5Dcreate( h5_group, "map", H5T_NATIVE_DOUBLE, h5_dataspace, H5P_DEFAULT );
        H5Dwrite( h5_dataset, H5T_NATIVE_DOUBLE, h5_dataspace, H5S_ALL, H5P_DEFAULT, Data );
        H5Dclose( h5_dataset );
        H5Sclose( h5_dataspace );


        pre_proc(1);
        lset(1);

        H5Gclose( h5_group );
        H5Gclose( h5_LinesGroup );

    }

    fclose( LsetErrFd );
    fclose( EdgesFd );
    fclose( RegsFd );
    H5Fclose( h5_Lset1Map );
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

