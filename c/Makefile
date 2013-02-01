#compile time options
#OPTS += -DFFTLOG #define to compile fft log function
#OPTS += -DDEBUG_IO #define for some debuggin I/O - none implemented so far
#OPTS += -DDEBUG -DDEBUG_LEVEL=0 #leave undefined for no debugging - 0,1, and 2 give progressively more output to stderr

CC          = gcc
OPTIMIZE    = -g -O3 -Wall
GSLI        = -I/opt/local/include
GSLL        = -L/opt/local/lib
FFTWI       =  
FFTWL       =  
EXTRACFLAGS =
EXTRACLIB   =  

CLINK=$(CC)
CFLAGS=$(OPTIMIZE) $(GSLI) $(FFTWI) $(FITSI) $(EXTRACFLAGS) $(OPTS)
CLIB=$(EXTRACLIB) $(GSLL) $(FFTWL) $(FITSL) -lm -lgsl -lgslcblas

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

CCALCSOURCES = $(FFTLOG) global.c distances.c linear_powspec.c transfer_function.c \
	growth_function.c hubble.c peakheight.c mass_bias_functions.c \
	linear_corrfunc.c utils.c nonlinear_powspec.c nonlinear_corrfunc.c weaklens.c
CCALCOBJECTS=$(CCALCSOURCES:.c=.o)

CSOURCES = main.c
COBJECTS = $(CSOURCES:.c=.o)

$(EXEC): $(COBJECTS) $(CCALCOBJECTS)
	$(CLINK) -o $@ $(COBJECTS) $(CLIB) $(CCALCOBJECTS)

libcosmocalc.so: $(CCALCSOURCES) cosmocalc.h Makefile
	rm -f $(CCALCOBJECTS)
	$(CC) $(CFLAGS) -fPIC -c $(CCALCSOURCES) 
	$(CLINK) $(CFLAGS) $(CLIB) -shared -o libcosmocalc.so $(CCALCOBJECTS)

libcosmocalc.a: $(CCALCSOURCES) cosmocalc.h Makefile
	$(CC) $(CFLAGS) -fPIC -c $(CCALCSOURCES) 
	ar rcs libcosmocalc.a $(CCALCOBJECTS)

$(COBJECTS): $(CSOURCES) $(CCALCSOURCES) cosmocalc.h Makefile
$(CCALCOBJECTS): $(CCALCSOURCES) cosmocalc.h Makefile

.PHONY : clean
clean: 
	rm -f *.o

.PHONY : spotless
spotless: 
	rm -f *.o $(EXEC) libcosmocalc.a libcosmocalc.so *~

