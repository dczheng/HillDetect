/*
    dczheng
    created 2019-07-28

    Reference
    ---------
    http://dx.doi.org/10.5201/ipol.2012.g-cv
*/

#include "allvars.h"

int *edgex, *edgey;
long edgen;
FILE *fd_lines;

void get_c1_c2( double *c1, double *c2 ) {

    long i;
    double s1, s2, n1, n2;
    n1 = n2 = 0;
    s1 = s2 = 0;
    for( i=0; i<Npixs; i++ ) {
        if ( Phi[i] > 0 ) {
            n1 ++;
            s1 += Data[i];
        }
        else {
            n2 ++;
            s2 += Data[i];
        }
    }

    *c1 = n1 ? ( s1/n1 ) : 0;
    *c2 = n2 ? ( s2/n2 ) : 0;

}

void get_s1_s2( double *s1, double *s2 ) {

    long i;
    *s1 = *s2 = 0;
    for( i=0; i<Npixs; i++ ) {
        if ( Phi[i] > 0 ) {
            *s1 += Data[i];
        }
        else {
            *s2 += Data[i];
        }
    }

}

void save_phi( int iter ) {

    char buf[110];
    FILE *fd;
    int i, j;
    long index;

    sprintf( buf, "%s/%s_phi_%i", All.PhiDir, FileName, iter );
    fd = fopen( buf, "w" );

    for( i=0, index=0; i<Height; i++ ) {
        for( j=0; j<Width; j++, index++ ) {
            fprintf( fd, "%g ", Phi[index] );
        }
        fprintf( fd, "\n" );
    }

    fclose( fd );
}

void save_line( int iter ) {

    int i,j;
    long index;
    double s1, s2;
    get_s1_s2( &s1, &s2 );

    for( i=0, index=0, edgen=0; i<Height; i++ ) {
        for( j=0; j<Width; j++, index++) {
            if ( Phi[index] > 0 &&
                 (
                  (j>0 && Phi[index-1]<0) ||
                  (j+1<Width && Phi[index+1]<0) ||
                  (i>0 && Phi[index-Width]<0) ||
                  (i+1<Height && Phi[index+Width]<0)
                  )
               ) {
                edgex[edgen] = j;
                edgey[edgen] = i;
                edgen++;
            }

        }
    }

    //printf( "edgen: %li\n", edgen );

#define myprintf( xy ) { \
    fprintf( fd_lines, "%i %li %g %g ", iter, edgen, s1, s2 ); \
    for( index=0; index<edgen; index++ ) { \
        fprintf( fd_lines, "%i ", xy[index] ); \
    } \
    fprintf( fd_lines, "\n" ); \
}

    myprintf( edgex );
    myprintf( edgey );

#undef myprintf

}

void init_phi() {

    int i, j;
    long index;

    Phi = malloc( Width * Height * sizeof(double) );

    for( i=0, index=0; i<Height; i++ )
        for ( j=0; j<Width; j++, index++ ) {
            Phi[index]= sin(i*M_PI/5.0) * sin(j*M_PI/5.0);
        }

/*
    FILE *fd;
    char buf[110];
    sprintf( buf, "%s-phi0.dat", fits_fn );
    fd = fopen( buf, "w" );
    for( i=0, index=0; i<height; i++ ) {
        for ( j=0; j<width; j++, index++ ) {
            fprintf( fd, "%g ", phi[index] );
        }
        fprintf( fd, "\n" );
    }
*/

}

void lset_init() {

     char buf[MYFILENAME_MAX];

     edgex = malloc( sizeof(int) * Npixs );
     edgey = malloc( sizeof(int) * Npixs );

     sprintf( buf, "%s/%s.out", All.OutputDir, FileName );
     //printf( "%s\n", buf );
     fd_lines = fopen( buf, "w" );
     init_phi();

}

void lset_free() {

    free( edgex );
    free( edgey );
    free( Phi );
    fclose( fd_lines );
}

#define DIVIDE_EPS 1e-16
void lset() {

    /*
      This function is copied from `http://www.ipol.im/pub/art/2012/g-cv/chanvese_20120715.tar.gz`
    */

    const long NumPixels = Width * Height;
    const long NumEl = NumPixels;
    const double *fPtr;
    double PhiDiffNorm, PhiDiff;
    double *PhiPtr;
    double  PhiLast, Delta, PhiX, PhiY, IDivU, IDivD, IDivL, IDivR;
    double Dist1, Dist2, PhiTol, dt, c1, c2;
    int Iter, i, j;
    int iu, id, il, ir;

    PhiTol = All.Tol;
    PhiDiffNorm = (All.Tol > 0) ? All.Tol*1000 : 1000;
    dt = All.TimeStep;

    lset_init();

    get_c1_c2( &c1, &c2 );

    for(Iter = 1; Iter <= All.MaxIters; Iter++)
    {
        PhiPtr = Phi;
        fPtr = Data;
        PhiDiffNorm = 0;

        for(j = 0; j < Height; j++)
        {
            iu = (j == 0) ? 0 : -Width;
            id = (j == Height - 1) ? 0 : Width;

            for(i = 0; i < Width; i++, PhiPtr++, fPtr++)
            {
                il = (i == 0) ? 0 : -1;
                ir = (i == Width - 1) ? 0 : 1;

                Delta = dt/(M_PI*(1 + PhiPtr[0]*PhiPtr[0]));
                PhiX = PhiPtr[ir] - PhiPtr[0];
                PhiY = (PhiPtr[id] - PhiPtr[iu])/2;
                IDivR = 1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY);
                PhiX = PhiPtr[0] - PhiPtr[il];
                IDivL = 1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY);
                PhiX = (PhiPtr[ir] - PhiPtr[il])/2;
                PhiY =  PhiPtr[id] - PhiPtr[0];
                IDivD = 1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY);
                PhiY = PhiPtr[0] - PhiPtr[iu];
                IDivU = 1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY);

                Dist1 = fPtr[0] - c1;
                Dist2 = fPtr[0] - c2;
                Dist1 *= Dist1;
                Dist2 *= Dist2;

                /* Semi-implicit update of phi at the current point */
                PhiLast = PhiPtr[0];
                PhiPtr[0] = (PhiPtr[0] + Delta*(
                        All.Mu*(PhiPtr[ir]*IDivR + PhiPtr[il]*IDivL
                            + PhiPtr[id]*IDivD + PhiPtr[iu]*IDivU)
                        - All.Nu - All.Lambda1*Dist1 + All.Lambda2*Dist2) ) /
                    (1 + Delta*All.Mu*(IDivR + IDivL + IDivD + IDivU));
                PhiDiff = (PhiPtr[0] - PhiLast);
                PhiDiffNorm += PhiDiff * PhiDiff;
            }
        }

        PhiDiffNorm = sqrt(PhiDiffNorm/NumEl);
        get_c1_c2( &c1, &c2 );

/*
        printf( "[iter: %i]  Delta: %e,  c1: %e, c2: %e\n",
                Iter, PhiDiffNorm, c1, c2 );
*/

        if(Iter >= 2 && PhiDiffNorm <= PhiTol)
            break;

        save_line( Iter );

        if ( All.IsSavePhi )
            save_phi(Iter);

    }

    lset_free();    

}

