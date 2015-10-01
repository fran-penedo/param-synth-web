#!/bin/bash

DREAL=dReal-3.15.09.01-linux-shared-libs.tar.gz
DREALURL=https://github.com/dreal/dreal3/releases/download/v3.15.09.01/$DREAL
NUSMV=http://nusmv.fbk.eu/distrib/NuSMV-2.5.4.tar.gz

mkdir lib

# NuSMV
wget -O - $NUSMV | tar -xz --strip 1 -C lib
cd lib/cudd*
make -f Makefile_64bit
cd ../nusmv
./configure
make
cd ../..

#dReal
wget -O - $DREALURL | tar -xz -C lib
mv lib/$DREAL lib/dReal

