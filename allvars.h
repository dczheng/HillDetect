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
    char InputDir[ MYFILENAME_MAX ],
         OutputDir[ MYFILENAME_MAX ],
         PhiDir[ MYFILENAME_MAX ],
         FileName[ MYFILENAME_MAX ],
         FileNameList[ MYFILENAME_MAX ];

    int  LogNorm, LogNorm1, MaxIters, MaxIters1, IsSavePhi,
         SigmaClipping,
         SigmaClipping1,
         FTClipping,
         DataCutting, ParalleLevel;
    double  Mu, Nu, Tol, Lambda1, Lambda2, TimeStep,
            Mu1, Nu1, Tol1, Lambda11, Lambda21, TimeStep1,
    RSigma, RSigma1,
            CuttingXStart, CuttingXEnd,
            CuttingYStart, CuttingYEnd ; 

}GlobalParams;

//extern_start
extern char sep_str[SEP_LEN];
extern double *Data, *Phi, *DataRaw, CRVAL1, CRVAL2, CDELT1, CDELT2,
              FREQ;
extern int 
            Width, Height, WidthGlobal, HeightGlobal, NpixsGlobal,
            ThisTask, NTask, FileNum,
            NfofEdge, NfofRegion,
            *edgex, *edgey,
            Npixs, edgen, XShift, YShift, CRPIX1, CRPIX2,
            HStartCut, HEndCut, WStartCut, WEndCut;

extern char FileName[ MYFILENAME_MAX ],
            *AllFileNames, *InputBaseName;
extern GlobalParams All;
extern FILE *LogFileFd, *LsetErrFd,
            *RegsFd, *EdgesFd;
//extern_end

hid_t   h5_Regs, h5_Lines, h5_LinesGroup, h5_Edges, h5_Lset0Map, h5_Lset1Map;

extern int *Next, *Head, *Len, *Tail;

#define ZDEBUG

#include "signal_hander.h"

