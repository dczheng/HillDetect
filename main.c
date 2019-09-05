/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

char *bname;

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

void run1() {

    int *next_lset0, *head_lset0, *len_lset0, gn,
         i, p, j, k, Width_global, Height_global,
         Xs[2], Ys[2],
         xmin, xmax, ymin, ymax, x, y, w, h, flag, index;
    double t0, t1;
    char buf[ MYFILENAME_MAX ];
//#define FIXEDSIZE
#ifdef FIXEDSIZE
    int wh;
#endif
#ifdef ADD_PAD
    int lmin;
#endif

//#define RUN1_DEBUG
    read_fits( All.FileName );
    bname = basename( All.FileName );
    pre_proc(0);

    put_sep;
    sprintf( OutputPrefix, "%s", bname );
    writelog( "run lset on level 0 ...\n" );
    t0 = second();
    XShift = 0;
    YShift = 0;
    lset(0);
    t1 = second();
    writelog( "time :%.3f sec\n", t1-t0 );
    writelog( "run lset on level 0 ... done.\n" );
    writelog( "NfofEdge: %i\n", NfofEdge );
    put_sep;

    gn = NfofEdge;
    Width_global = Width;
    Height_global = Height;
    next_lset0 = malloc( sizeof(int) * Npixs );
    head_lset0 = malloc( sizeof(int) * Npixs );
    len_lset0 = malloc( sizeof(int) * Npixs );
    memcpy( next_lset0, Next, sizeof(int) * Npixs );
    memcpy( head_lset0, Head, sizeof(int) * Npixs );
    memcpy( len_lset0,  Len, sizeof(int) * Npixs );
    find_region_free();
    do_sync;
#ifdef ADD_PAD
    lmin = 50;
#endif

    for( k=0, index=0; k<gn; k++ ) {
        flag = 0;
        p = head_lset0[k];
        while( p>=0 ) {
            x = p % Width_global;
            y = p / Width_global;
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

    if ( ThisTask == 0 ) {
        sprintf( buf, "%s/%s_map.dat", All.OutputDir, OutputPrefix );
        output_map( buf, Width, Height, DataRaw, NULL, NULL );
    }

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
            x = p % Width_global;
            y = p / Width_global;
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
        sprintf( OutputPrefix, "%s_%05i", bname, k );

        for( i=ymin; i<ymax; i++ ) {
            for( j=xmin; j<xmax; j++ ) {
#ifdef FIXEDSIZE
                Data[(i-ymin)*wh+(j-xmin)] = DataRaw[i*Width_global+j];
#else
                Data[(i-ymin)*w+(j-xmin)] = DataRaw[i*Width_global+j];
#endif
            }
        }
        Xs[0] = xmin;
        Xs[1] = xmax;
        Ys[0] = ymin;
        Ys[1] = ymax;
        sprintf( buf, "%s/%s_map.dat", All.OutputDir, OutputPrefix );
        output_map( buf, Width_global,  Height_global, DataRaw, Xs, Ys );

#ifdef RUN1_DEBUG
        sprintf( buf, "%s/%s_%03i_before.dat", All.OutputDir, OutputPrefix, ThisTask );
        output_map( buf, Width,  Height, Data, NULL, NULL );
#endif
        pre_proc(1);
#ifdef RUN1_DEBUG
        sprintf( buf, "%s/%s_%03i_after.dat", All.OutputDir, OutputPrefix, ThisTask );
        output_map( buf, Width,  Height, Data, NULL, NULL );
#endif
        lset(1);

    }

    free( next_lset0 );
    free( head_lset0 );
    free( len_lset0 );
    free_fits();

}

void global_init() {

    XShift = YShift = 0;
    create_dir( All.OutputDir );
    create_dir( All.PhiDir );

}
void global_free() {
}

int main( int argc, char **argv ) {

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
    read_parameters( argv[1] );
    global_init();

    switch (All.ParalleLevel) {
        case 0:
            run0();
            break;
        case 1:
            run1();
            break;
        default:
            printf( "Unsupported paralle level.\n" );
            endrun( "" );
    } 
    put_sep;

    global_free();
    fclose( LogFileFd );
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;

}

