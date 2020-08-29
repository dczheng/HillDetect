/*
    dczheng
    created 2019-07-28

    Reference
    ---------
    http://dx.doi.org/10.5201/ipol.2012.g-cv
*/

#include "allvars.h"

hid_t f_lset_h5;

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

/*
void save_phi( int iter ) {

    FILE *fd;
    int i, j;
    int index;

    fd = myfopen( "w", "%s/%s_phi_%i", All.PhiDir, InputBaseName, iter );

    for( i=0, index=0; i<Height; i++ ) {
        for( j=0; j<Width; j++, index++ ) {
            fprintf( fd, "%g ", Phi[index] );
        }
        fprintf( fd, "\n" );
    }

    fclose( fd );
}
*/

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
    int i, *edges;
    char buf[20];
    double s[2];
    hid_t h5_dsp,h5_attr, h5_g;
    hsize_t h5_dims[2]; 


    get_s1_s2( &s[0], &s[1] );

    h5_dsp = H5Screate( H5S_SCALAR );

    if ( iter == 1 ) {
        h5_attr = H5Acreate( 
                f_lset_h5, "iters", H5T_NATIVE_INT,
                        h5_dsp, H5P_DEFAULT );
    }
    else {
        h5_attr = H5Aopen( f_lset_h5, "iters", H5P_DEFAULT );
    }
    H5Awrite( h5_attr, H5T_NATIVE_INT, &iter );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp );


    sprintf( buf, "iter%i", iter );
    h5_g = H5Gcreate( f_lset_h5, buf, 0 );

    h5_dims[0] = 2;
    hdf5_write_attr_nd( h5_g, H5T_NATIVE_DOUBLE, h5_dims, 1, "S1S2", s );

    h5_dims[0] = 2;
    if  ( edgen == 0 ) {
        h5_dims[1] = 1;
        edges = malloc( sizeof(int) * 2 );
        edges[0] = edges[1] = 0;;
    }
    else {
        h5_dims[1] = edgen;
        edges = malloc( sizeof(int) * edgen * 2 );
        for ( i=0; i<edgen; i++ ) {
            edges[i] = edgey[i] + HStartCut;
            edges[i+edgen] = edgex[i] + WStartCut;
        }

    }

    hdf5_write_data( h5_g, H5T_NATIVE_INT, h5_dims, 2, "line", edges );

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
    fclose(fd);
*/

}

void open_files_for_lset() {

    char buf[120];

    sprintf( buf, "%s/lset_lines.hdf5", All.OutputDir );
    f_lset_h5 = H5Fcreate( buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT  );

}

void close_files_for_lset() {
    H5Fclose( f_lset_h5  );
}

#define DIVIDE_EPS 1e-16

void lset() {
    /*
      This function is copied from `http://www.ipol.im/pub/art/2012/g-cv/chanvese_20120715.tar.gz`
    */

    put_start(0);
    const int NumPixels = Width * Height;
    const int NumEl = NumPixels;
    const double *fPtr;
    double PhiDiffNorm, PhiDiff;
    double *PhiPtr;
    double  PhiLast, Delta, PhiX, PhiY, IDivU, IDivD, IDivL, IDivR;
    double Dist1, Dist2, PhiTol, dt, c1, c2;
    int Iter, i, j;
    int iu, id, il, ir;
    double Mu, Nu, Tol, Lambda1, Lambda2, TimeStep;
    int MaxIters;

    Mu = All.Mu;
    Nu = All.Nu;
    Tol = All.Tol;
    Lambda1 = All.Lambda1;
    Lambda2 = All.Lambda2;
    MaxIters = All.MaxIters;
    TimeStep = All.TimeStep;

    PhiTol = Tol;
    dt = TimeStep;
    PhiDiffNorm = (Tol > 0) ? Tol*1000 : 1000;

    edgex = malloc( sizeof(int) * Npixs );
    edgey = malloc( sizeof(int) * Npixs );

    open_files_for_lset();
    init_phi();

    get_c1_c2( &c1, &c2 );
    for(Iter = 1; Iter <= MaxIters; Iter++)
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
                        Mu*(PhiPtr[ir]*IDivR + PhiPtr[il]*IDivL
                            + PhiPtr[id]*IDivD + PhiPtr[iu]*IDivU)
                        - Nu - Lambda1*Dist1 + Lambda2*Dist2) ) /
                    (1 + Delta*Mu*(IDivR + IDivL + IDivD + IDivU));
                PhiDiff = (PhiPtr[0] - PhiLast);
                PhiDiffNorm += PhiDiff * PhiDiff;
            }
        }

        PhiDiffNorm = sqrt(PhiDiffNorm/NumEl);
        get_c1_c2( &c1, &c2 );

        lset_find_line();

        fprintf( LOG_FILE, "[%5i]  Delta: %e, c1: %e, c2: %e\n",
                  Iter, PhiDiffNorm, c1, c2 );
        printf( "\r[%5i]  Delta: %e, c1: %e, c2: %e",
                  Iter, PhiDiffNorm, c1, c2 );
        fflush( stdout );
       fflush( LOG_FILE );
        save_line( Iter );

        if(Iter >= 2 && PhiDiffNorm <= PhiTol)
            break;

    }
    printf( "\n" ); 
    close_files_for_lset();

    lset_group_finder();

    free( edgex );
    free( edgey );
    free( Phi );
    put_end(0);
}
