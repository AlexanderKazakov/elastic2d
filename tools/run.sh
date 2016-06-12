#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-n] number_of_processes
    [-t] task_name
    [-p] (to open Paraview after calculation)"
    1>&2; exit 1; }

np=`cat /proc/cpuinfo | grep processor | wc -l`
task="cgal2d"
run_paraview=0

while getopts ":n:t:p" o; do
    case "${o}" in
        n)
            np=${OPTARG}
            ;;
        t)
            task=${OPTARG}
            ;;
        p)
            run_paraview=1
            ;;
        *)
            usage
            ;;
    esac
done

rm -rf snapshots/vtk snapshots/detector snapshots/zaxis
rm -f *.log
mkdir -p snapshots/vtk snapshots/detector snapshots/zaxis
echo "Start ./build/gcm_exe with $np processes ..."
mpirun -np $np ./build/gcm_exe --task ${task}

if ((run_paraview)); then
    echo "Run Paraview ..."
    if (($np == 1)); then
        if (ls snapshots/vtk | grep vts); then
            paraview --data=snapshots/vtk/core00snap..vts
        else
            paraview --data=snapshots/vtk/core00snap..vtu
        fi
    else
        paraview
    fi
fi
