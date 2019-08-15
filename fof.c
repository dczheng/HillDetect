/*
    dczheng
    created 2019-08-14
*/

#include "fof.h"

// a implementation of 2D fof.

long *Next, *Head, *Tail, *Len, Width, Height, FoFN;
int *FoFMap;
//#define FOF_SINGLE_DEBUG
void fof_single_finder( long p0 ) {

    int  i, j, x0, y0;
    long p1, l, s, ss, p_next;

    if ( FoFMap[p0] == 0 )
        return;

    x0 = p0 % Width;
    y0 = p0 / Width;
    p_next = -1;

#ifdef FOF_SINGLE_DEBUG
    printf( "start\n" );
    fflush( stdout );
#endif
    for( i=-1; i<2; i++ )
        for( j=-1; j<2; j++ ) {
            p1 = (y0+i) * Width + (x0+j);
#ifdef FOF_SINGLE_DEBUG
            printf( "p0: %li, p1:  %li, NPixs: %li (%i*%i)\n",
            p0, p1, FoFN, Width, Height );
            fflush( stdout );
#endif
            if ( p1 == p0     ||
                 p1 < 0       ||
                 p1 > FoFN-1 ||
                 FoFMap[p1]  == 0 )
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
        fof_single_finder( p_next );

}

void fof_init( int *map, int W, int H ) {

    long i;

    Width = W;
    Height = H;
    FoFN = W * H;
    FoFMap = map;

    Next = malloc( sizeof(long) * FoFN );
    Head = malloc( sizeof(long) * FoFN );
    Len = malloc( sizeof(long) * FoFN );
    Tail = malloc( sizeof(long) * FoFN );

    for( i=0; i<FoFN; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i] = 1;
    }
}

void fof_reset() {

    long i;
    for( i=0; i<FoFN; i++ ) {
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

void fof_sort () {

    long *tmp, p;
    tmp = malloc( sizeof(long) * FoFN * 3 );
    for ( p=0; p<FoFN; p++ ) {
        tmp[ 3*p ] = Len[p];
        tmp[ 3*p+1 ] = Head[p];
        tmp[ 3*p+2 ] = Tail[p];
    }
    
    qsort( tmp, FoFN, sizeof(long)*3, fof_compare_len );

    for ( p=0; p<FoFN; p++ ) {
        Len[p] = tmp[ 3*p ];
        Head[p] = tmp[ 3*p+1 ];
        Tail[p] = tmp[ 3*p+2 ];
    }
    free( tmp );
}

void fof() {
    long p;
    for( p=0; p<FoFN; p++ ) {
        fof_single_finder(  p );
    }
    fof_sort();
}
