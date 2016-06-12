#!/usr/bin/env bash

statement() {
# $1 - number of cores
    mkdir -p snapshots/zaxis
    mkdir -p snapshots/detector
    mkdir -p snapshots/vtk
    echo "Start ./build/gcm_exe with $1 processes ..."
    (/usr/bin/time --format "%e" mpirun -np $1 ./build/gcm_exe -t layers) 2>&1 | tee -a tmp 
}


statement 8
statement 7
statement 6
statement 5
statement 4
statement 3
statement 2
statement 1

statement 8
statement 7
statement 6
statement 5
statement 4
statement 3
statement 2
statement 1
