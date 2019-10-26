#!/bin/csh

cd Bluebird
cd RxNLib
make
cd ../VectorLib
make
#-02 => -O2
cd ../MatrixLib
make
cd ..
make
cd ..
mkdir executables
mv Bluebird/Bluebird executables/

cd Ostrich
make MW_SER
mv OstrichSerial ../executables/Ostrich
