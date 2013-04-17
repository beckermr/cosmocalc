#options
CC          = gcc
OPTIMIZE    = -g -O3 -Wall
GSLI        = -I/opt/local/include
GSLL        = -L/opt/local/lib
EXTRACFLAGS = -I./cosmocalc/c
EXTRACLIB   = 

CLINK=$(CC)
CFLAGS=$(OPTIMIZE) $(GSLI) $(EXTRACFLAGS)
CLIB=$(EXTRACLIB) $(GSLL) -lm -lgsl -lgslcblas

#set it all up
EXEC = computecosmo
all: $(EXEC)

CSOURCES = main.c
COBJECTS = $(CSOURCES:.c=.o)

CCALCSOURCES=$(FFTLOG) global.c distances.c linear_powspec.c transfer_function.c \
	growth_function.c hubble.c peakheight.c mass_bias_functions.c \
	linear_corrfunc.c utils.c nonlinear_powspec.c nonlinear_corrfunc.c weaklens.c \
	haloprofs.c
CCALCOBJECTS1=$(CCALCSOURCES:.c=.o)
CCALCOBJECTS=$(addprefix ./cosmocalc/c/,$(CCALCOBJECTS1))

$(EXEC): $(COBJECTS) $(CCALCOBJECTS)
	$(CLINK) -o $@ $(COBJECTS) $(CLIB) $(CCALCOBJECTS)

$(COBJECTS): $(CSOURCES) Makefile $(CCALCOBJECTS)

.PHONY : clean
clean: 
	rm -f *.o

.PHONY : spotless
spotless: 
	rm -f *.o $(EXEC) *~

