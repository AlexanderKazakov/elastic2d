#!/usr/bin/env bash

rm -rf snapshots/vtk snapshots/detector snapshots/zaxis
rm -f *.log
mkdir -p snapshots/vtk snapshots/detector snapshots/zaxis
./build/ndi

#mkdir -p snapshots/zaxis/all
#
#for i in {0..360..2}
#do
#	first=$(printf "snapshots/zaxis/mesh0core00snap%04d.txt" "$i")
#	second=$(printf "snapshots/zaxis/mesh1core00snap%04d.txt" "$i")
#	common=$(printf "snapshots/zaxis/all/snap%04d.txt" "$i")
#	cat $first $second > $common
#done

gnuplot ./tools/gnuplot/ndi && \
eog ./snapshots/detector/
