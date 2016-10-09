#!/usr/bin/env bash

./tools/build.sh -t ndi
rm -rf snapshots/vtk snapshots/detector snapshots/zaxis
rm -f *.log
mkdir -p snapshots/vtk snapshots/detector snapshots/zaxis
./build/ndi
gnuplot tools/gnuplot/gpltmp \
&& eog res.png

