/*
    dczheng
    created 2019-08-14
*/

#include "allvars.h"

// a implementation of 2D fof.
// need global variable: long *Next, *Head, *Tail, *Len;

//#define FOF_SINGLE_DEBUG
void fof_single_finder( int *map, int W, int H, long p0 ) {

    int  i, j, x0, y0;
    long p1, l, s, ss, p_next, N;

    if ( map[p0] == 0 )
        return;

    N = W * H;
    x0 = p0 % W;
    y0 = p0 / W;
    p_next = -1;

#ifdef FOF_SINGLE_DEBUG
    printf( "start\n" );
    fflush( stdout );
#endif
    for( i=-1; i<2; i++ )
        for( j=-1; j<2; j++ ) {
            p1 = (y0+i) * W + (x0+j);
#ifdef FOF_SINGLE_DEBUG
            printf( "p0: %li, p1:  %li, NPixs: %li (%i*%i)\n",
            p0, p1, N, W, H );
            fflush( stdout );
#endif
            if ( p1 == p0     ||
                 p1 < 0       ||
                 p1 > N-1 ||
                 map[p1]  == 0 )
                continue;

#ifdef FOF_SINGLE_DEBUG
            printf( "1\n" );
            fflush( stdout );
#endif
            if ( Head[p1] == Head[p0] )
                continue;

            if ( Len[Head[p1]] > Len[Head[p0]] ) {
                l = p1;
                s = p0;
             }
             else {
                l = p0;
                s = p1;
             }

#ifdef FOF_SINGLE_DEBUG
            printf( "2\n" );
            fflush( stdout );
#endif
            Next[Tail[Head[l]]] = Head[s];
            Tail[Head[l]] = Tail[Head[s]];
            Len[Head[l]] += Len[Head[s]];
            Len[Head[s]] = 1;

            ss = Head[s];
            do{
                /*
                printf( "p0: %li, p1: %li l: %li, s: %li, ss: %li\n"
                    "Head[l]: %li, Head[s]: %li\n",
                p0, p1, l, s, ss, Head[l], Head[s] );
                */
                Head[ss] = Head[l];
            } while( (ss=Next[ss]) >= 0 );

#ifdef FOF_SINGLE_DEBUG
            printf( "3\n" );
            fflush( stdout );
#endif
            p_next = p1;

    }

#ifdef FOF_SINGLE_DEBUG
    printf( "end\n\n" );
    fflush( stdout );
#endif

    if ( p_next >= 0 )
        fof_single_finder( map, W, H, p_next );

}

void fof_init( long N ) {

    long i;
    Next = malloc( sizeof(long) * N );
    Head = malloc( sizeof(long) * N );
    Len = malloc( sizeof(long) * N );
    Tail = malloc( sizeof(long) * N );
    for( i=0; i<N; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i] = 1;
    }
}

void fof_reset( long N ) {

    long i;
    for( i=0; i<N; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i] = 1;
    }
}

void fof_free() {

    free( Next );
    free( Head );
    free( Tail );
    free( Len );

}

int fof_compare_len( const void *a, const void *b ) {
    return ( (*(long*)a) < *((long*)b) ) ? 1 : -1;
}

void fof_sort ( long N ) {

    long *tmp, p;
    tmp = malloc( sizeof(long) * N * 3 );
    for ( p=0; p<N; p++ ) {
        tmp[ 3*p ] = Len[p];
        tmp[ 3*p+1 ] = Head[p];
        tmp[ 3*p+2 ] = Tail[p];
    }
    
    qsort( tmp, N, sizeof(long)*3, fof_compare_len );

    for ( p=0; p<N; p++ ) {
        Len[p] = tmp[ 3*p ];
        Head[p] = tmp[ 3*p+1 ];
        Tail[p] = tmp[ 3*p+2 ];
    }
    free( tmp );
}

void fof( int *map, int W, int H ) {
    long p;
    for( p=0; p<W*H; p++ ) {
        fof_single_finder( map, W, H , p );
    }
    fof_sort( W*H );
}
