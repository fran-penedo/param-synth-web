#!/bin/bash

# NuSMV
wget -O - http://nusmv.fbk.eu/distrib/NuSMV-2.5.4.tar.gz | tar -xz --strip 1 -C lib
cd lib/cudd*
make -f Makefile_64bit
cd ../nusmv
./configure
make
cd ../..

#dReal
cd lib/dreal
./build-dreal.sh
