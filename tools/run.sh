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

rm -rf snaps
rm -f *.log
mkdir snaps
echo "Start ./build/gcm_exe with $np processes ..."
mpirun -np $np ./build/gcm_exe

if ((run_paraview)); then
    echo "Run Paraview ..."
    if (($np == 1)); then
        if (ls snaps | grep vts); then
            paraview --data=snaps/core00meshstep..vts
        else
            paraview --data=snaps/core00meshstep..vtu
        fi
    else
        paraview
    fi
fi
