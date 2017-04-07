#!/bin/bash


mkdir build

cd build

cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DOPENCL=OFF -DCUDA=OFF ../

make
