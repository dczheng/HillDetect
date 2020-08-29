/*
    dczheng
    created 2019-07-28
*/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "unistd.h"
#include "math.h"
#include "time.h"
#include "dirent.h"
#include "mpi.h"
#include "signal.h"
#include "limits.h"
#include "sys/stat.h"
#include "macros.h"
#include "libgen.h"
#include "protos.h"
#include "fitsio.h"
#include "hdf5.h"

#define MYFILENAME_MAX 100
#define SEP_LEN  50

typedef struct GlobalParams {
    char OutputDir[ MYFILENAME_MAX ],
         FileName[ MYFILENAME_MAX ];

    int  LogNorm, MaxIters,
         SigmaClipping,
         SigmaClipping1,
         DataCutting, PeakCenterFlag,
         LsetPixMin, SecondFinderPixMin, OnlyFoFPixMin,
         DisableSecondFinder, SecondFinderPad, OnlyFoF;
    double  Mu, Nu, Tol, Lambda1, Lambda2, TimeStep,
            RSigma, RSigma1,
            CuttingXStart, CuttingXEnd,
            CuttingYStart, CuttingYEnd, Beam; 

}GlobalParams;

//extern_start
extern char sep_str[SEP_LEN];
extern double *Data, *Phi, *DataRaw, CRVAL1, CRVAL2, CDELT1, CDELT2,
              SigmaClippingVmin;
extern int 
            Width, Height, WidthGlobal, HeightGlobal, NpixsGlobal,
            FileNum,
            *edgex, *edgey,
            Npixs, edgen, XShift, YShift, CRPIX1, CRPIX2,
            HStartCut, HEndCut, WStartCut, WEndCut, CurGroup,
            *lset_Next, *lset_Head, *lset_Len, lset_Nreg, Nsource, Ngroup;

extern char FileName[ MYFILENAME_MAX ],
            *AllFileNames, *InputBaseName;
extern GlobalParams All;
hid_t h5_fof;
FILE *LOG_FILE;
//extern_end

extern int *Next, *Head, *Len, *Tail;
