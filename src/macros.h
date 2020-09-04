#define SQR(X)  ( (X)*(X) )

#define myfopen( opt, _fmt, ... ) ({\
    char fopenbuf[120];\
    FILE *fd;\
    do {\
        sprintf( fopenbuf, _fmt, ##__VA_ARGS__ );\
        fd = fopen( fopenbuf, opt );\
        if ( NULL == fd ) {\
            printf( "can not open `%s`.\n", fopenbuf );\
            endrun("");\
        }\
    }while(0);\
    fd;\
})

#define create_dir( _fmt, ... ) {\
    char _buf[150];\
    sprintf( _buf, _fmt, ##__VA_ARGS__ );\
    create_dir0( _buf );\
}

#define endrun0( _fmt, ... ) {\
    fprintf( stderr, _fmt, ##__VA_ARGS__ ); \
    fprintf( stderr, "END IN: ( %s %s %i )\n" , __FILE__, __FUNCTION__, __LINE__ ); \
    exit( 0 ); \
}

#define writelog( _flag, _fmt, ... ) {\
    fprintf( LOG_FILE, _fmt, ##__VA_ARGS__ );\
    fflush( LOG_FILE );\
    if ( !(_flag) ) {\
        printf( _fmt, ##__VA_ARGS__ );\
        fflush( stdout );\
    }\
}

#define endrun( s ) {\
    endrun0( "error info: %s\n", s ); \
}

#define put_sep( _flag )  {\
    writelog( _flag, sep_str );\
}

#define mytimer_start \
    double timer1, timer2;\
    (void) timer1;\
    (void) timer2;\
    writelog( 0, "[Timer Start in `%s`]\n", __FUNCTION__ ); \
    timer1 = second(); \
    timer2 = timer1;

#define mytimer \
    writelog( 0, "[Time in `%s`]: %g sec\n", __FUNCTION__, second() - timer2 ); \
    timer2 = second();

#define mytimer_end \
    printf( "[Total Time in `%s`]: [%g sec]\n", __FUNCTION__, second() - timer1 ); \

#define put_start( _flag ) {\
    writelog( _flag, ">>> start `%s`\n", __FUNCTION__ );\
}

#define put_end( _flag ) {\
    writelog( _flag, ">>> end `%s`\n", __FUNCTION__ );\
}
