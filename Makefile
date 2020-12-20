CPPFLAGS= -g -std=c++11 
SOURCES=src/arED.cpp \
		ext/essaMEM/sparseSA.cpp  ext/essaMEM/sssort_compact.cc \
		ext/edlib/edlib.cpp

all:
	$(CXX) $(CPPFLAGS) -I ext/ -o redit $(SOURCES) -lz

clean:
	rm -f redit
