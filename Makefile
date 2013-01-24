#vars
CXX         = g++ -fopenmp
GSLI        = -I/opt/local/include 
CXXFLAGS    = -Wall -O3 ${GSLI}
GSLL        = -L/opt/local/lib
LDFLAGS     = ${GSLL} -lm -lgsl -lgslcblas 

EXEC = computecosmo

all: $(EXEC)
sharedlib: libcosmocalc.so
staticlib: libcosmocalc.a

CXXSOURCES = main.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)

$(EXEC): $(CXXOBJECTS)
	$(CXX) -o $@ $(CXXOBJECTS) $(LDFLAGS) -L./ -lcosmocalc

$(CXXOBJECTS): $(CXXSOURCES) libcosmocalc.so Makefile

CCALCSOURCES = cosmocalc_init.cpp cosmocalc_distances.cpp cosmocalc_growth.cpp cosmocalc_transfer_function.cpp \
	cosmocalc_linear_powspec.cpp
CCALCOBJECTS = $(CCALCSOURCES:.cpp=.o)

libcosmocalc.so: $(CCALCSOURCES) cosmocalc.h cosmocalc_assert.h Makefile
	$(CXX) $(CXXFLAGS) -fPIC -c $(CCALCSOURCES) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o libcosmocalc.so $(CCALCOBJECTS)

libcosmocalc.a: $(CCALCSOURCES) cosmocalc.h cosmocalc_assert.h Makefile
	$(CXX) $(CXXFLAGS) -c $(CCALCSOURCES) 
	ar rcs libcosmocalc.a $(CCALCOBJECTS)

.PHONY : clean
clean: 
	rm -f *.o 

.PHONY : spotless
spotless: 
	rm -f *.o *.a *.so $(EXEC) *~

