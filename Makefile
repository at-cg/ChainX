CPPFLAGS= -g -std=c++11 
SOURCES1=src/redit.cpp \
				ext/essaMEM/sparseSA.cpp  ext/essaMEM/sssort_compact.cc

SOURCES2=src/redit_edlib.cpp \
				 ext/essaMEM/sparseSA.cpp  ext/essaMEM/sssort_compact.cc ext/edlib/edlib.cpp

all:
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o redit $(SOURCES1) -lz
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o redit_edlib $(SOURCES2) -lz

clean:
	rm -f redit redit_edlib
