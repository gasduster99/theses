
all:
	-@ echo ""
	-@ echo "Must specify a target"
	-@ echo "   MW_SER : MuddyWaters serial"
	-@ echo "   MW_MPI : MuddyWaters MPI"
	-@ echo "   CW_SER : Clearwater serial"
	-@ echo "   CW_MPI : Clearwater MPI"
	-@ echo "   CW_ISO : Clearwater IsoFit"
	-@ echo "   MW_ISO : MuddyWaters IsoFit"
	-@ echo "   W32_ISO : Windows IsoFit"
	-@ echo "   W32_OST : Windows Ostrich"
	-@ echo ""

W32_ISO:	*.o_W32_ISO
	icpc -static -o IsoFit *.o
	rm -f *.o


W32_OST:	*.o_W32_OST
	icpc -static -o Ostrich *.o
	rm -f *.o

CW_ISO:	*.o_ISO
	icpc -static -o IsoFit *.o
	rm -f *.o

MW_ISO:	*.o_ISO
	icpc -static -o IsoFit *.o
	rm -f *.o

MW_MPI:	*.o_MPI
	mpiicpc -o Ostrich *.o
	rm -f *.o
	mv -f Ostrich OstrichMPI

CW_MPI:	*.o_MPI
	icpc -static -o Ostrich *.o
	rm -f *.o

MW_SER:	*.o_SER
	g++ -O2 -o Ostrich *.o
	rm -f *.o
	mv -f Ostrich OstrichSerial

CW_SER:	*.o_SER
	icpc -static -o Ostrich *.o
	rm -f *.o

*.o_MPI:	*.cpp *.c *.h
	mpiicpc -c -O2 -DUSE_MPI *.cpp
	mpiicpc -c -O2 -DUSE_MPI *.c

*.o_SER:	*.cpp *.c *.h
	## icpc -static -w1 -wd858 -wd470 -c -O2  -DUSE_MPI_STUB *.cpp
	g++ -c -O2  -DUSE_MPI_STUB *.cpp
	## icpc -static -w1 -wd858 -wd470 -c -O2  -DUSE_MPI_STUB *.c
	g++ -c -O2  -DUSE_MPI_STUB *.c

*.o_ISO:	*.cpp *.c *.h
	icpc -static -w1 -wd858 -wd470 -c -O2  -DISOFIT_BUILD -DUSE_MPI_STUB *.cpp
	icpc -static -w1 -wd858 -wd470 -c -O2  -DISOFIT_BUILD -DUSE_MPI_STUB *.c

*.o_W32_ISO:	*.cpp *.c *.h
	icpc -static -w -c -O2  -DISOFIT_BUILD -DUSE_MPI_STUB *.cpp
	icpc -static -w -c -O2  -DISOFIT_BUILD -DUSE_MPI_STUB *.c

*.o_W32_OST:	*.cpp *.c *.h
	icpc -static -w -c -O2  -DUSE_MPI_STUB *.cpp
	icpc -static -w -c -O2  -DUSE_MPI_STUB *.c

clean:
	rm -f *.o
