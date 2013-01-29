#vars
CXX         = g++ -fopenmp
GSLI        = -I/opt/local/include
GSLL        = -L/opt/local/lib
CXXFLAGS    = -Wall -g -O3 ${GSLI}
LDFLAGS     = ${GSLL} -lm -lgsl -lgslcblas

EXEC = computecosmo

all: $(EXEC)
lib: sharedlib staticlib
sharedlib: libcosmocalc.so
staticlib: libcosmocalc.a

CXXSOURCES = main.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)

CCALCHEADERS = cosmocalc.h cosmocalc_assert.h flrw_distances.h w0wacosmo.h eh98_transfunct.h linpowspec.h halofit.h
CCALCSOURCES = flrw_distances.cpp w0wacosmo.cpp eh98_transfunct.cpp linpowspec.cpp halofit.cpp
CCALCOBJECTS = $(CCALCSOURCES:.cpp=.o)

$(EXEC): $(CXXOBJECTS) $(CCALCOBJECTS)
	$(CXX) $(CXXFLAGS) -c $(CCALCSOURCES)
	$(CXX) $(GSLL) $(CXXOBJECTS) $(CCALCOBJECTS) -o $@ $(LDFLAGS) 

libcosmocalc.so: $(CCALCSOURCES) $(CCALCHEADERS) Makefile
	rm -f $(CCALCOBJECTS)
	$(CXX) $(CXXFLAGS) -fPIC -c $(CCALCSOURCES) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o libcosmocalc.so $(CCALCOBJECTS)

libcosmocalc.a: $(CCALCSOURCES) $(CCALCHEADERS) Makefile
	$(CXX) $(CXXFLAGS) -c $(CCALCSOURCES) 
	ar rcs libcosmocalc.a $(CCALCOBJECTS)

$(CXXOBJECTS): $(CXXSOURCES) $(CCALCSOURCES) $(CCALCHEADERS) Makefile
$(CCALCOBJECTS): $(CCALCSOURCES) $(CCALCHEADERS) Makefile


.PHONY : clean
clean: 
	rm -f *.o 

.PHONY : spotless
spotless: 
	rm -f *.o *.a *.so $(EXEC) *~
