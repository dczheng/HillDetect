/*
    dczheng
    created 2019-07-28
*/

#define SQR( X )  ( (X)*(X) )
#define CUBE( X )  ( SQR(X)*(X) )

#define put_sep  { \
    if ( ThisTask == 0 ) \
        printf( "%s", sep_str  ); \
        fflush( stdout );\
}

#define endrun( a ) {\
    fprintf( stderr, "%i\n", a );\
    MPI_Abort( MPI_COMM_WORLD, 0 );\
    exit(0);\
}
