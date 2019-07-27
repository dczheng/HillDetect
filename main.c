/*
    dczheng
    created 2019-07-28
*/

#include "allvars.h"

void init_sep_str() {
    memset( sep_str, '-', SEP_LEN-2 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';
}

void global_free() {
    free( Phi );
    free( AllFileNames );
}

void global_init() {
}

int main( int argc, char **argv ) {

    int i, N, j, sigs;
    char buf[MYFILENAME_MAX];
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NTask );

    init_sep_str();

    read_parameters( argv[1] );
    put_sep;

    global_init();

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

            if ( All.SigmaClipping )
                sigma_clipping();

            normalize();
            lset_init();
            lset();
            lset_free();

            free_fits();
        }

        //MPI_Barrier( MPI_COMM_WORLD );
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;

}
