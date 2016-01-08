#!/usr/bin/env bash

rm -f time.png speedup.png
cpus=`cat /proc/cpuinfo | grep processor | wc -l`

for np in $(seq 1 $((2*$cpus)))
do
    rm -f snaps/*.vtk;
    echo "Running with" $np "processes";
    /usr/bin/time -f "$np %e" -o speedup.txt -a mpirun -np $np ./build/gcm_exe;
done

gnuplot tools/plotter
rm speedup.txt
eog time.png
