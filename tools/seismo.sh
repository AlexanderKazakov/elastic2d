#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-n] number_of_processes"
    1>&2; exit 1; }

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

rm -rf snapshots/2dseismo
rm -f *.log
mkdir -p snapshots/2dseismo
echo "Start ./build/gcm_seismo with $np processes ..."
mpirun -np $np ./build/gcm_seismo

gnuplot tools/gnuplot-seismo.txt && eog snapshots/2dseismo/core00statement0000.bin.png


