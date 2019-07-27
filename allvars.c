/*
    dczheng
    created 2019-07-28
*/


#include "allvars.h"

char sep_str[SEP_LEN], FileName[ MYFILENAME_MAX ], *AllFileNames;
double *Data, *Phi;
int Width, Height, ThisTask, NTask, FileNum;
long Npixs;
struct global_parameters_struct All;
