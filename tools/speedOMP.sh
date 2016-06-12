#!/usr/bin/env bash

statement() {
# $1 - number of threads
    mkdir -p snapshots/zaxis
    mkdir -p snapshots/detector
    mkdir -p snapshots/vtk
    export OMP_NUM_THREADS=$1
    echo "Start ./build/gcm_exe with $1 threads ..."
    (time mpirun -np 1 ./build/gcm_exe -t cgalani) 2>&1 | tee -a tmp 
}


statement 12
statement 11
statement 10
statement 9
statement 8
statement 7
statement 6
statement 5
statement 4
statement 3
statement 2
statement 1


statement 12
statement 11
statement 10
statement 9
statement 8
statement 7
statement 6
statement 5
statement 4
statement 3
statement 2
statement 1


statement 12
statement 11
statement 10
statement 9
statement 8
statement 7
statement 6
statement 5
statement 4
statement 3
statement 2
statement 1
