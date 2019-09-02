#include "allvars.h"

double t0, t1;

void create_dir0( char *s ) {

    if ( ThisTask == 0 ){
        if ( access( s, 0 ) == -1 ){
            writelog( "`%s` is created by Task 0\n", s );

            if ( mkdir( s, 0755) == -1 )
                endrun0( "failed create directory %s.\n", s );
        }
    }

    MPI_Barrier( MPI_COMM_WORLD );

}

double second() {
    return ( (double) clock() / CLOCKS_PER_SEC );
}

