CPPFLAGS= -DNDEBUG -std=c++11 -O3
export CC=$(CXX)
SOURCES1=src/chainx.cpp \
				 ext/essaMEM/sparseSA.cpp  ext/essaMEM/sssort_compact.cc

SOURCES2=src/edlib_wrapper.cpp \
				 ext/edlib/edlib.cpp

SOURCES3=src/printanchors.cpp \
				 ext/essaMEM/sparseSA.cpp  ext/essaMEM/sssort_compact.cc

SOURCES4=src/chainx-mininimizer.cpp

all:
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o chainX $(SOURCES1) -lz
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o edlib_wrapper $(SOURCES2) -lz
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o printanchors $(SOURCES3) -lz
	+$(MAKE) -C ext/minimap2-2.24
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o chainX-mininimizer $(SOURCES4) ext/minimap2-2.24/libminimap2.a -lz -lm -lpthread

clean:
	+$(MAKE) -C ext/minimap2-2.24 clean
	rm -f chainX edlib_wrapper printanchors chainX-mininimizer
