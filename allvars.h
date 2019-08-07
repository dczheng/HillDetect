/*
    dczheng
    created 2019-07-28
*/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"
#include "mpi.h"
#include "signal.h"
#include "limits.h"
#include "fitsio.h"
#include "macros.h"
#include "protos.h"

#define MYFILENAME_MAX 100
#define SEP_LEN  50

extern char sep_str[SEP_LEN];
extern double *Data, *Phi;
extern int Width, Height, ThisTask, NTask, FileNum;
extern long Npixs;
extern char FileName[ MYFILENAME_MAX ], *AllFileNames;

extern struct global_parameters_struct {
    char InputDir[ MYFILENAME_MAX ],
         OutputDir[ MYFILENAME_MAX ],
         PhiDir[ MYFILENAME_MAX ],
         FileNameList[ MYFILENAME_MAX ];

    int  LogNorm, MaxIters, IsSavePhi, SigmaClipping, FTClipping;
    double  Mu, Nu, Tol, Lambda1, Lambda2, TimeStep, RSigma; 

}All;
