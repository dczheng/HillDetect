GSL_INCL   = 
GSL_LIBS   = -lm 

FITS_INCL  =
FITS_OPTS  = 
FITS_LIBS  = -lcfitsio

FFTW_INCL  = 
FFTW_LIBS  = -ldrfftw -ldfftw

PYTHON ?= $(shell which python3)

LIBS       = $(GSL_LIBS) $(FITS_LIBS) $(FFTW_LIBS) -lm
INCL       = $(GSL_INCL) $(FITS_INCL) $(FFTW_INCL)

OPTS       = $(FITS_OPTS) -Wall #-O2 #-O3
DEBUG     ?=
CC         =  mpicc 

EXEC       =  ./bin/fgext
SRCS       = allvars.c fof.c io.c lset.c main.c pre_proc.c read_params.c debug.c
MY_INCL    =  allvars.h  protos.h add_params.h
OBJS       =  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -o $(EXEC)

%.o:%.c $(MY_INCL) Makefile
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@

allvars.c: allvars.h preprocessor.py
	$(shell ${PYTHON} ./preprocessor.py)

add_params.h: allvars.h preprocessor.py

protos.h: $(SRCS) preprocessor.py


.PHONY: clean

clean:
	-rm  $(EXEC)  $(OBJS) allvars.c add_params.h protos.h
