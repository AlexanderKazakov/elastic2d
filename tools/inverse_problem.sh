#!/usr/bin/env bash

usage() {
    echo "Usage: $0
    [-n] number_of_processes"
    1>&2; exit 1;
}

statement() {
# $1 - number of layers
    rm -rf snapshots/1dseismo
    rm -rf snapshots/vtk
    rm -f *.log
    mkdir -p snapshots/1dseismo
    mkdir -p snapshots/vtk
    echo "Start ./build/gcm_inverse_problem $1 with $np processes ..."
    mpirun -np $np ./build/gcm_inverse_problem $1
    
    gnuplot tools/gnuplot-1d-binary.txt #&& eog snapshots/1dseismo/core00statement0000.bin.png

    mv snapshots $1
    mv $1 saved_snaps/$1 
}

np=`cat /proc/cpuinfo | grep processor | wc -l`

while getopts ":n:" o; do
    case "${o}" in
        n)
            np=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

rm -rf saved_snaps
mkdir saved_snaps

statement 2
statement 3
statement 5
statement 8
