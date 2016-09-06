#!/usr/bin/env bash

./tools/run.sh -n 1 -t tmp

mkdir -p snapshots/zaxis/all

for i in {0..360..2}
do
	first=$(printf "snapshots/zaxis/mesh0core00snap%04d.txt" "$i")
	second=$(printf "snapshots/zaxis/mesh1core00snap%04d.txt" "$i")
	common=$(printf "snapshots/zaxis/all/snap%04d.txt" "$i")
	cat $first $second > $common
done

gnuplot ./tools/gnuplot/ndi && \
eog ./snapshots/detector/
