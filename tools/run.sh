#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-n] number_of_processes
    [-p] (to open paraview after calculation - only for single-process mode)"
    1>&2; exit 1; }

np=`cat /proc/cpuinfo | grep processor | wc -l`
run_paraview=0

while getopts ":n:p" o; do
    case "${o}" in
        n)
            np=${OPTARG}
            ;;
        p)
            run_paraview=1
            ;;
        *)
            usage
            ;;
    esac
done

rm -f snaps/*
echo "running calculation with $np processes"
mpirun -np $np ./build/gcm

if (($np == 1 && run_paraview)); then
    echo "running paraview"
    paraview --data=snaps/core00_snapshot..vtk
fi