# define additional include directories
INCS = -IMatrixLib -IRxNLib -IVectorLib -I.. -I.
LIBDIRS = -LMatrixLib -LRxNLib -LVectorLib
LIBS   = -lMatrix -lRxN -lVector

all:    *.o
	g++ $(LIBDIRS) $(LIBS) -o Bluebird *.o RxNLib/libRxN.a VectorLib/libVector.a MatrixLib/libMatrix.a
	#gcc -static $(LIBDIRS) $(LIBS)
	rm -f *.o

*.o:    *.cpp *.h
	g++ -fpermissive -w -Wno-deprecated -c -O2 -DCLK_TCK=20000000 $(INCS) *.cpp
	#gcc -static -w1 -wd858 -wd1125 -wd470 -Wno-deprecated -c -O2 -DUSE_MPI_STUB -DCLK_TCK=20000000 $(INCS) *.cpp
	#gcc -O2 $(INCS) *.cpp
clean:
	rm -f *.o
