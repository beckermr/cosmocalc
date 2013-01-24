#vars
CXX         = g++ -g
GSLI        = -I/opt/local/include
CXXFLAGS    = -Wall -O3 ${GSLI}
GSLL        = -L/opt/local/lib
LDFLAGS     = -lm -lgsl -lgslcblas

EXEC = computecosmo

all: $(EXEC)
lib: sharedlib staticlib
sharedlib: libcosmocalc.so
staticlib: libcosmocalc.a

CXXSOURCES = main.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)

CCALCSOURCES = cosmocalc_init.cpp cosmocalc_distances.cpp cosmocalc_growth.cpp cosmocalc_transfer_function.cpp \
	cosmocalc_linear_powspec.cpp cosmocalc_nonlinear_powspec.cpp
CCALCOBJECTS = $(CCALCSOURCES:.cpp=.o)

$(EXEC): $(CXXOBJECTS) $(CCALCOBJECTS)
	$(CXX) $(CXXFLAGS) -c $(CCALCSOURCES)
	$(CXX) $(GSLL) $(CXXOBJECTS) $(CCALCOBJECTS) -o $@ $(LDFLAGS) 

libcosmocalc.so: $(CCALCSOURCES) cosmocalc.h cosmocalc_assert.h Makefile
	rm -f $(CCALCOBJECTS)
	$(CXX) $(CXXFLAGS) -fPIC -c $(CCALCSOURCES) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o libcosmocalc.so $(CCALCOBJECTS)

libcosmocalc.a: $(CCALCSOURCES) cosmocalc.h cosmocalc_assert.h Makefile
	$(CXX) $(CXXFLAGS) -c $(CCALCSOURCES) 
	ar rcs libcosmocalc.a $(CCALCOBJECTS)

$(CXXOBJECTS): $(CXXSOURCES) $(CCALCSOURCES) cosmocalc.h Makefile
$(CCALCOBJECTS): $(CCALCSOURCES) cosmocalc.h Makefile


.PHONY : clean
clean: 
	rm -f *.o 

.PHONY : spotless
spotless: 
	rm -f *.o *.a *.so $(EXEC) *~

