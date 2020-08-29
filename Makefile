FITS_INCL  =
FITS_LIBS  = -lcfitsio

HDF5_INCL  =
HDF5_OPTS  = -DH5_USE_16_API
HDF5_LIBS  = -lhdf5

FFTW_INCL  = 
FFTW_LIBS  =

PYTHON    ?= $(shell which python3)

LIBS       = $(GSL_LIBS) $(FITS_LIBS) $(FFTW_LIBS) $(HDF5_LIBS) -lm
INCL       = $(GSL_INCL) $(FITS_INCL) $(FFTW_INCL) $(HDF5_INCL) -I./src

OPTS       = $(FITS_OPTS)  $(HDF5_OPTS) -Wall #-O2 #-O3
DEBUG     ?=
CC         =  gcc 

EXEC       =  HillDetect

TEMP_FILES   = ./src/add_params.h ./src/protos.h ./src/allvars.c

SRCS       =  $(wildcard ./src/*.c) ./src/allvars.c 
MY_INCL    =  $(wildcard ./src/*.h) ./src/add_params.h ./src/protos.h
OBJS       =  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -o $(EXEC)

%.o:%.c $(MY_INCL)
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@

$(TEMP_FILES): ./src/preprocessor.py Makefile
	@./src/preprocessor.py

	
.PHONY: clean

clean:
	-rm  ./src/*.o $(TEMP_FILES)
	-rm  $(EXEC)

test:
	@echo $(OBJS)
	@echo $(MY_INCL)
