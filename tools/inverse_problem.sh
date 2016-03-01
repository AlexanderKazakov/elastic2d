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

rm -rf snapshots/1dseismo
rm -f *.log
mkdir -p snapshots/1dseismo
echo "Start ./build/gcm_inverse_problem with $np processes ..."
mpirun -np $np ./build/gcm_inverse_problem

gnuplot tools/gnuplot-1d-binary.txt && eog snapshots/1dseismo/seismo.png

