#define additional include directories
INCS = -I../MatrixLib -I../RxNLib -I../VectorLib -I.. -I.

all:    *.o
	ar r libMatrix.a *.o
	ar t libMatrix.a
	rm -f *.o

*.o:    *.cpp *.h
	#icpc -static -w1 -wd858 -wd1125 -wd470 -c -O2 -DUSE_MPI_STUB $(INCS) *.cpp
	g++ -w -c -O2 $(INCS) *.cpp
clean:
	rm -f *.o
