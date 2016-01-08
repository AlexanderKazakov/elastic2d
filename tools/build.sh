#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-t] (to compile tests)
    [-d] (to compile in Debug mode)
    [-c] (to clean all before compiling)"
    1>&2; exit 1; }

np=`cat /proc/cpuinfo | grep processor | wc -l`

rm -f snaps/*
cmake_line="cmake .."

while getopts ":dc" o; do
    case "${o}" in
        c)
            rm -rf build
            mkdir build
            ;;
        d)
            cmake_line="$cmake_line -DCMAKE_BUILD_TYPE=Debug"
            ;;
        *)
            usage
            ;;
    esac
done

cd build
eval $cmake_line
make -j$np
cd ..