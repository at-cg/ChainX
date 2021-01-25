CPPFLAGS= -DNDEBUG -std=c++11 -O3
SOURCES1=src/chainx.cpp \
				 ext/essaMEM/sparseSA.cpp  ext/essaMEM/sssort_compact.cc

SOURCES2=src/edlib_wrapper.cpp \
				 ext/edlib/edlib.cpp

SOURCES3=src/printanchors.cpp \
				 ext/essaMEM/sparseSA.cpp  ext/essaMEM/sssort_compact.cc

all:
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o chainX $(SOURCES1) -lz
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o edlib_wrapper $(SOURCES2) -lz
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o printanchors $(SOURCES3) -lz

clean:
	rm -f chainX edlib_wrapper printanchors
