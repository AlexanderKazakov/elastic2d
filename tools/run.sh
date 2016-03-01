#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-n] number_of_processes
    [-p] (to open Paraview after calculation)"
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

rm -rf snapshots/vtk
rm -f *.log
mkdir -p snapshots/vtk
echo "Start ./build/gcm_exe with $np processes ..."
mpirun -np $np ./build/gcm_exe

if ((run_paraview)); then
    echo "Run Paraview ..."
    if (($np == 1)); then
        if (ls snapshots/vtk | grep vts); then
            paraview --data=snapshots/vtk/core00step..vts
        else
            paraview --data=snapshots/vtk/core00step..vtu
        fi
    else
        paraview
    fi
fi
