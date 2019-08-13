LIBS       =  -lcfitsio -lm -ldrfftw -ldfftw
INCL       = 
OPTS       = -Wall #-O2 #-O3
DEBUG     ?=
CC         =  mpicc

EXEC       =  ./bin/fgext
SRCS       =  $(wildcard *.c)
MY_INCL    =  $(wildcard *.h)
OBJS       =  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -o $(EXEC)

%.o:%.c $(MY_INCL) Makefile
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@
test: $(EXEC)
	./bin/fgext ./example.in

.PHONY: clean

clean:
	-rm  $(EXEC)  $(OBJS)
