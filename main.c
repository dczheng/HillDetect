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
            pre_proc();
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

    int *next, *head, *tail, *len, gn,
         i, p, j, k, Width_global,
         xmin, xmax, ymin, ymax, x, y, w, h,
         lmin;

#define RUN1_DEBUG
    read_fits( All.FileName );
    pre_proc();
    All.DataCutting = 0;

    lset(0);
    writelog( "NfofEdge: %i\n", NfofEdge );
    gn = NfofEdge;
    Width_global = Width;
    next = malloc( sizeof(int) * Npixs );
    head = malloc( sizeof(int) * Npixs );
    tail = malloc( sizeof(int) * Npixs );
    len = malloc( sizeof(int) * Npixs );
    memcpy( next, Next, sizeof(int) * Npixs );
    memcpy( head, Head, sizeof(int) * Npixs );
    memcpy( tail, Tail, sizeof(int) * Npixs );
    memcpy( len,  Len, sizeof(int) * Npixs );
    find_region_free();
    do_sync;
    lmin = 50;

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
        p = head[k];
        while(p>0) {
            x = p % Width_global;
            y = p / Width_global;
            xmin = ( x<xmin ) ? x : xmin;
            xmax = ( x>xmax ) ? x : xmax;
            ymin = ( y<ymin ) ? y : ymin;
            ymax = ( y>ymax ) ? y : ymax;
            p = next[p];
        }
        printf( "task: %03i, group %05i [%05i] "
                "region: (%i, %i) - (%i, %i)\n",
                 ThisTask, k,
                 len[k], xmin, ymin, xmax, ymax );
        
        h = ymax - ymin;
        w = xmax - xmin;
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

    }

    free( next );
    free( head );
    free( tail );
    free( len );
    free_fits();

}

void global_init() {

    if ( ThisTask == 0 ) {
        if ( access( "fgext.log", 0 ) == -1 ) {
            printf( "`fgext.log` is created by Task 0\n" );
            if ( mkdir( "fgext.log", 0755 ) == -1 )
                endrun( "failed create directory `fgext.log`.\n"  );
        }
    }
    do_sync;

    LogFileFd = myfopen( "w", "./fgext.log/%s.log", bname );
    memset( sep_str, '-', SEP_LEN-2 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';
    XShift = YShift = 0;

}
void global_free() {
    fclose( LogFileFd );
}

int main( int argc, char **argv ) {

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NTask );

    bname = basename( argv[1] );
    global_init();
    read_parameters( argv[1] );

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

    global_free();
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;

}

