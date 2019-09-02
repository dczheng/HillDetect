#define SQR(X)  ( (X)*(X) )
#define writelog( fmt, ... ) { \
    fprintf( LogFileFd, fmt, ##__VA_ARGS__ ); \
    fflush( LogFileFd ); \
    if ( ThisTask == 0 ) { \
        printf( fmt, ##__VA_ARGS__ ); \
    }\
}

#define myfopen( opt, fmt, ... ) ({\
    char buf[100];\
    FILE *fd;\
    do {\
        sprintf( buf, fmt, ##__VA_ARGS__ );\
        fd = fopen( buf, opt );\
        if ( NULL == fd ) {\
            printf( "can not open `%s`.\n", buf );\
            endrun("");\
        }\
    }while(0);\
    fd;\
})

#define create_dir( fmt, ... ) {\
    char buf[100];\
    sprintf( buf, fmt, ##__VA_ARGS__ );\
    create_dir0( buf );\
}

#define endrun0( fmt, ... ) {\
    fprintf( stderr, fmt, ##__VA_ARGS__ ); \
    fprintf( stderr, "END IN: ( %s %s %i )\n" , __FILE__, __FUNCTION__, __LINE__ ); \
    MPI_Abort( MPI_COMM_WORLD, 0 ); \
    exit( 0 ); \
}

#define endrun( s ) {\
    endrun0( "error info: %s\n", s ); \
}

/*
#define endrun( i ) {\
    endrun0( "error level: %i\n", i ); \
}
*/


#define do_sync MPI_Barrier( MPI_COMM_WORLD );

#define put_sep writelog( sep_str );

#define mytimer_start() \
    double timer1, timer2;\
    (void) timer1;\
    (void) timer2;\
writelog( "[Timer Start in `%s`]\n", __FUNCTION__ ); \
    timer1 = second(); \
    timer2 = timer1;

#define mytimer() \
    writelog( "[Time in `%s`]: %g sec\n", __FUNCTION__, second() - timer2 ); \
    timer2 = second();

#define mytimer_end() \
    writelog( "[Total Time in `%s`]: [%g sec]\n", __FUNCTION__, second() - timer1 ); \

