#!/bin/bash

if [ "$1" = "vdw" ] || [ "$1" == "all" ] || [ "$1" == "cuda" ]; then
    git submodule update --init --recursive
    cd lapack/
    mkdir -p build
    cd build/
    cmake -DCMAKE_INSTALL_LIBDIR=`pwd`/../install_dir ..
    cmake --build . -j --target install
    cd ../..
fi

rm -rf build
mkdir -p build
cd build

echo $HOSTNAME | grep "rc.usf.edu"

if [ $? == 0 ]; then
    module purge
    module load compilers/gcc/5.1.0
    module load compilers/intel/2015_cluster_xe
    if [ "$1" = "cuda" ] || [ "$1" == "all" ]; then
        module load apps/cuda/8.0
    fi
    export CC=icc
    export CXX=icpc
fi

echo $HOSTNAME | grep "bridges.psc.edu"

if [ $? == 0 ]; then
    module purge
    module load gcc/5.3.0
    if [ "$1" = "cuda" ] || [ "$1" == "all" ]; then
        module load cuda/8.0
    fi
    export CC=gcc
    export CXX=g++
fi

echo $HOSTNAME | grep ".sdsc.edu"

if [ $? == 0 ]; then
    module purge
    module purge
    module load cmake/3.9.1
    module load gnu/4.9.2
    if [ "$1" = "cuda" ] || [ "$1" == "all" ]; then
        module load cuda/8.0
    fi
    export CC=gcc
    export CXX=g++
fi

if [ "$1" = "debug" ]; then
    cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=Debug -Wno-dev ../
elif [ "$1" = "cuda" ]; then
    cmake -DQM_ROTATION=ON -DVDW=ON -DMPI=OFF -DCUDA=ON -DCMAKE_BUILD_TYPE=release -Wno-dev ../
elif [ "$1" = "mpi" ]; then
    cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=ON -DCUDA=OFF -DCMAKE_BUILD_TYPE=release -Wno-dev ../
elif [ "$1" = "vdw" ]; then
    cmake -DQM_ROTATION=ON -DVDW=ON -DMPI=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=release -Wno-dev ../
elif [ "$1" = "all" ]; then
    cmake -DQM_ROTATION=ON -DVDW=ON -DMPI=ON -DCUDA=ON -DCMAKE_BUILD_TYPE=release -Wno-dev ../
else
    cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=Release -Wno-dev ../
fi

make -j
