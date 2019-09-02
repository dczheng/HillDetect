#include "stdio.h"
#include "stdlib.h"
extern int *Next, *Head, *Tail, *Len;

void fof_init( int *map, int W, int H );
void fof_free();
void fof_reset();
void fof();
