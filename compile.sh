#!/bin/bash

rm -rf build
mkdir -p build
cd build

echo $HOSTNAME | grep "rc.usf.edu"

if [ $? == 0 ]; then
  module purge
  module load compilers/gcc/5.1.0
  module load compilers/intel/2015_cluster_xe
  export CC=icc
  export CXX=icpc
fi

echo $HOSTNAME | grep "bridges.psc.edu"

if [ $? == 0 ]; then
  export CC=icc
  export CXX=icpc
fi

echo $HOSTNAME | grep ".sdsc.edu"

if [ $? == 0 ]; then
  module load cmake/3.9.1
  module load gnu/4.9.2
  module load intel/2015.2.164

  export CC=icc
  export CXX=icpc
fi

if [ "$1" = "debug" ]; then
  cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DOPENCL=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=Debug -Wno-dev ../
elif [ "$1" = "cuda" ]; then
  cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DOPENCL=OFF -DCUDA=ON -DCMAKE_BUILD_TYPE=release -Wno-dev ../
elif [ "$1" = "mpi" ]; then
  cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=ON -DOPENCL=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=release -Wno-dev ../
else
  cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DOPENCL=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=Release -Wno-dev ../
fi

make
