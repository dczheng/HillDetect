/*
    dczheng
    created 2019-07-28

    Reference
    ---------
    http://dx.doi.org/10.5201/ipol.2012.g-cv
*/

#include "allvars.h"

FILE *fdlf;

void get_c1_c2( double *c1, double *c2 ) {

    int i;
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

    int i;
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
    int index;

    sprintf( buf, "%s/%s_phi_%i", All.PhiDir, OutputPrefix, iter );
    fd = fopen( buf, "w" );

    for( i=0, index=0; i<Height; i++ ) {
        for( j=0; j<Width; j++, index++ ) {
            fprintf( fd, "%g ", Phi[index] );
        }
        fprintf( fd, "\n" );
    }

    fclose( fd );
}

void lset_find_line() {
    int i,j;
    int index;
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

}

void save_line( int iter ) {

    //printf( "edgen: %i\n", edgen );
    int index;
    double s1, s2;
    get_s1_s2( &s1, &s2 );

#define myprintf( xy, dis ) { \
    fprintf( fdlf, "%i %i %g %g ", iter, edgen, s1, s2 ); \
    for( index=0; index<edgen; index++ ) { \
        fprintf( fdlf, "%i ", xy[index] + dis ); \
    } \
    fprintf( fdlf, "\n" ); \
}

    myprintf( edgex, XShift );
    myprintf( edgey, YShift );

#undef myprintf

}

void init_phi() {

    int i, j;
    int index;

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

#define DIVIDE_EPS 1e-16
void lset( int mode ) {

    /*
      This function is copied from `http://www.ipol.im/pub/art/2012/g-cv/chanvese_20120715.tar.gz`
    */

    const int NumPixels = Width * Height;
    const int NumEl = NumPixels;
    const double *fPtr;
    double PhiDiffNorm, PhiDiff;
    double *PhiPtr;
    double  PhiLast, Delta, PhiX, PhiY, IDivU, IDivD, IDivL, IDivR;
    double Dist1, Dist2, PhiTol, dt, c1, c2;
    int Iter, i, j;
    int iu, id, il, ir;
    FILE *fd;

    if ( mode == 0 ) {
        PhiTol = All.Tol;
        PhiDiffNorm = (All.Tol > 0) ? All.Tol*1000 : 1000;
    }
    else {
        PhiTol = All.Tol1;
        PhiDiffNorm = (All.Tol1 > 0) ? All.Tol1*1000 : 1000;
    }
    dt = All.TimeStep;

    edgex = malloc( sizeof(int) * Npixs );
    edgey = malloc( sizeof(int) * Npixs );

    init_phi();
    find_region_init();

    fd = myfopen( "w", "%s/%s_lset_err.dat", All.OutputDir, OutputPrefix );

    if ( mode == 1 ) {
        fdlf = myfopen( "w", "%s/%s_lset_lines.dat", 
            All.OutputDir, OutputPrefix );
    }
    if ( mode == 0 ) {
        if ( ThisTask == 0 )
            fdlf = myfopen( "w", "%s/%s_lset_lines.dat", 
                All.OutputDir, OutputPrefix );
    }
    

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

        fprintf( fd, "[%i]  Delta: %e\nc1: %e, c2: %e\n",
                            Iter, PhiDiffNorm, c1, c2 );
        lset_find_line();

        if ( mode == 1 ) {
            save_line( Iter );
            if ( All.IsSavePhi )
                save_phi(Iter);
            //find_region( Iter, mode );
        }

        if ( mode == 0 ) {
            if ( ThisTask == 0 ) {
                save_line( Iter );
                if ( All.IsSavePhi )
                    save_phi(Iter);
            }
        }

        if(Iter >= 2 && PhiDiffNorm <= PhiTol)
            break;

    }

    if ( mode == 0 )
       find_region( Iter, mode );

    free( edgex );
    free( edgey );
    free( Phi );

    if ( mode == 1 ) {
//        find_region_free();
        fclose( fdlf );
    }
    if ( mode == 0 ) {
        if ( ThisTask == 0 )
            fclose( fdlf );
    }

    fclose( fd );

}

