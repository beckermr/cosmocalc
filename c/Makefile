#compile time options
#OPTS += -DFFTLOG #define to compile fft log function
#OPTS += -DDEBUG_IO #define for some debuggin I/O - none implemented so far
#OPTS += -DDEBUG -DDEBUG_LEVEL=0 #leave undefined for no debugging - 0,1, and 2 give progressively more output to stderr

#select your computer
#COMP="orange"
#COMP="mandor-icc"
#COMP="mandor-gcc"

CC          = gcc
OPTIMIZE    = -g -O3 -Wall
GSLI        = -I/opt/local/include
GSLL        = -L/opt/local/lib

#edit/add to match your machine
ifeq ($(COMP),"orange")
CC          = mpicc
OPTIMIZE    = -Wall -O3 -wd2259 -wd981
GSLI        =  -I/afs/slac.stanford.edu/g/ki/software/gsl/include 
GSLL        =  -L/afs/slac.stanford.edu/g/ki/software/gsl/lib
FFTWI       =  -I/afs/slac.stanford.edu/g/ki/software/fftw/3.2.2/include
FFTWL       =  -L/afs/slac.stanford.edu/g/ki/software/fftw/3.2.2/lib 
EXTRACFLAGS =
EXTRACLIB   =  -lgsl -lgslcblas -lfftw3 -lfftw3f
endif

ifeq ($(COMP),"mandor-icc")
CC          =  icc 
OPTIMIZE    =  -g -O3
GSLI        =  -I/home/beckermr/include
GSLL        =  -L/home/beckermr/lib -lgsl -lgslcblas
FFTWI       =  -I/home/beckermr/include
FFTWL       =  -L/home/beckermr/lib -lfftw3 -lfftw3f
EXTRACFLAGS =  -Wall -wd981 #-wd1419 -wd810
EXTRACLIB   =  
endif

ifeq ($(COMP),"mandor-gcc")
CC          =  gcc
OPTIMIZE    =  -g -O3 #-Werror
GSLI        =  -I/home/beckermr/include
GSLL        =  -L/home/beckermr/lib -lgsl -lgslcblas
FFTWI       =  -I/home/beckermr/include
FFTWL       =  -L/home/beckermr/lib -lfftw3 -lfftw3f
EXTRACFLAGS =  -ansi -std=c99 -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
	-Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline= 
EXTRACLIB   =  
endif

CLINK=$(CC)
CFLAGS=$(OPTIMIZE) $(GSLI) $(FFTWI) $(FITSI) $(EXTRACFLAGS) $(OPTS)
CLIB=$(EXTRACLIB) $(GSLL) $(FFTWL) $(FITSL) $(GSLL) -lm -lgsl -lgslcblas

#set it all up
EXEC = computecosmo
all: $(EXEC)
lib: staticlib
staticlib: libcosmocalc.a
sharedlib: libcosmocalc.so

ifneq (FFTLOG,$(findstring FFTLOG,$(CFLAGS)))
FFTLOG=
else
FFTLOG=fftlog.o
endif

CCALCSOURCES = $(FFTLOG) distances.c linear_powspec.c transfer_function.c \
	growth_function.c hubble.c peakheight.c mass_bias_functions.c \
	linear_corrfunc.c utils.c nonlinear_powspec.c nonlinear_corrfunc.c weaklens.c
CCALCOBJECTS=$(CCALCSOURCES:.c=.o)

CSOURCES = main.c
COBJECTS = $(CSOURCES:.c=.o)

$(EXEC): $(COBJECTS) libcosmocalc.a
	$(CLINK) $(GSLL) $(COBJECTS) libcosmocalc.a -o $@ $(CLIB) 

libcosmocalc.so: $(CCALCSOURCES) cosmocalc.h Makefile
	rm -f $(CCALCOBJECTS)
	$(CC) $(CFLAGS) -fPIC -c $(CCALCSOURCES) 
	$(CLINK) $(CFLAGS) $(CLIB) -shared -o libcosmocalc.so $(CCALCOBJECTS)

libcosmocalc.a: $(CCALCSOURCES) cosmocalc.h Makefile
	$(CC) $(CFLAGS) -c $(CCALCSOURCES) 
	ar rcs libcosmocalc.a $(CCALCOBJECTS)

$(COBJECTS): $(CSOURCES) $(CCALCSOURCES) cosmocalc.h Makefile
$(CCALCOBJECTS): $(CCALCSOURCES) cosmocalc.h Makefile

.PHONY : clean
clean: 
	rm -f *.o

.PHONY : spotless
spotless: 
	rm -f *.o $(EXEC) libcosmocalc.a libcosmocalc.so *~

