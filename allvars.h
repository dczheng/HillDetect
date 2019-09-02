/*
    dczheng
    created 2019-07-28
*/

#include "hdf5.h"
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

#define MYFILENAME_MAX 100
#define SEP_LEN  50

typedef struct GlobalParams {
    char InputDir[ MYFILENAME_MAX ],
         OutputDir[ MYFILENAME_MAX ],
         PhiDir[ MYFILENAME_MAX ],
         FileName[ MYFILENAME_MAX ],
         FileNameList[ MYFILENAME_MAX ];

    int  LogNorm, MaxIters, IsSavePhi, SigmaClipping, FTClipping,
         DataCutting, ParalleLevel;
    double  Mu, Nu, Tol, Lambda1, Lambda2, TimeStep, RSigma, 
            CuttingXStart, CuttingXEnd,
            CuttingYStart, CuttingYEnd ; 

}GlobalParams;

//extern_start
extern char sep_str[SEP_LEN];
extern double *Data, *Phi, *DataRaw;
extern int 
            Width, Height, ThisTask, NTask, FileNum,
            NfofEdge, NfofRegion,
            *edgex, *edgey,
            Npixs, edgen, XShift, YShift;

extern char FileName[ MYFILENAME_MAX ], *AllFileNames;
extern GlobalParams All;
extern FILE *LogFileFd, *LsetLinesFd;
//extern_end

extern int *Next, *Head, *Len, *Tail;

