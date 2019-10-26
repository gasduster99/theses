#define additional include directories
INCS = -I../MatrixLib -I../RxNLib -I../VectorLib -I.. -I.
  
all:	*.o
	ar r libRxN.a *.o
	ar t libRxN.a
	rm -f *.o

*.o:	*.cpp *.h
	## ricpc -static -w1 -wd858 -wd1125 -wd470 -Wno-deprecated -c -O2 -DUSE_MPI_STUB -DCLK_TCK=20000000 $(INCS) *.cpp
	g++ -fpermissive -w -c -O2 -DCLK_TCK=20000000 $(INCS) *.cpp

clean:
	rm -f *.o
