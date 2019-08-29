/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

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

            lset();

            free_fits();
        }

        //MPI_Barrier( MPI_COMM_WORLD );
    }
    free( AllFileNames );
}

void run1() {

    read_fits( All.FileName );
    pre_proc();
    free_fits();

}

int main( int argc, char **argv ) {

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NTask );

    memset( sep_str, '-', SEP_LEN-2 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';

    read_parameters( argv[1] );
    put_sep;

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


    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;

}

