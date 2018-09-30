#!/bin/bash

rm -rf build
mkdir -p build
cd build

if [ "$1" = "debug" ]; then
  cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DOPENCL=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=Debug -Wno-dev ../
else
  cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DOPENCL=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=Release -Wno-dev ../
fi

make
