/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

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

    pre_proc(0);
    put_sep;

    put_header( "run first finder" );
    mytimer_start();
    if ( ThisTask == 0 ) {
        sprintf( buf, "%s/%s_lset0_map.dat", All.OutputDir, InputBaseName );
        output_map( buf, WidthGlobal, HeightGlobal, DataRaw, NULL, NULL );
    }

    if ( ThisTask == 0 ) {
        LsetErrFd = myfopen( "w", "%s/%s_lset0_err.dat", All.OutputDir, InputBaseName );
        LsetLinesFd = myfopen( "w","%s/%s_lset0_lines.dat", All.OutputDir, InputBaseName );
        EdgesFd = myfopen( "w","%s/%s_lset0_edges.dat", All.OutputDir, InputBaseName );
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
        fclose( LsetLinesFd );
        fclose( EdgesFd );
    }
    mytimer_end();

}

void run_second_finder() {

   int i, p, j, k,
         Xs[2], Ys[2],
         xmin, xmax, ymin, ymax, x, y, w, h, flag, index;
    char buf[ MYFILENAME_MAX ];
    put_header( "run second finder" );

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
    LsetLinesFd = myfopen( "w","%s/%s_lset1_lines_%03i.dat",
                All.OutputDir, InputBaseName, ThisTask );
    EdgesFd = myfopen( "w","%s/%s_lset1_edges_%03i.dat",
                All.OutputDir, InputBaseName, ThisTask );
    RegsFd = myfopen( "w","%s/%s_lset1_regs_%03i.dat",
                All.OutputDir, InputBaseName, ThisTask );

#define RUN1_DEBUG
#ifdef RUN1_DEBUG
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
        Xs[1] = xmax;
        Ys[0] = ymin;
        Ys[1] = ymax;

        sprintf( buf, "%s/%s_lset1_raw_%04i.dat",
        All.OutputDir, InputBaseName, k );
        output_map( buf, WidthGlobal,  HeightGlobal, DataRaw, Xs, Ys );

        sprintf( buf, "%s/%s_lset1_before_pre_proc_%04i.dat",
        All.OutputDir, InputBaseName, k );
        output_map( buf, Width,  Height, Data, NULL, NULL );

        pre_proc(1);

        sprintf( buf, "%s/%s_lset1_after_pre_proc_%04i.dat",
        All.OutputDir, InputBaseName, k );
        output_map( buf, Width,  Height, Data, NULL, NULL );

        fprintf( LsetErrFd,   "Group: %03i\n", k);
        fprintf( LsetLinesFd, "Group: %03i\n", k);
        fprintf( EdgesFd,     "Group: %03i\n", k);
        fprintf( RegsFd,      "Group: %03i\n", k);
        lset(1);

    }

    fclose( LsetErrFd );
    fclose( LsetLinesFd );
    fclose( EdgesFd );
    fclose( RegsFd );
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
    put_header( "global init" );
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

