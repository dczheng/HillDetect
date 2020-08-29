#include "allvars.h"

double t0, t1;

void create_dir0( char *s ) {

    if ( access( s, 0 ) == -1 ){
        printf( "create `%s`\n", s );

        if ( mkdir( s, 0755) == -1 )
            endrun0( "failed create directory %s.\n", s );
    }

}

double second() {
    return ( (double) clock() / CLOCKS_PER_SEC );
}

