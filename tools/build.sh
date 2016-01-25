#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-o] (enable additional optimization)
    [-d] (to compile in Debug mode)
    [-c] (to clean all before compiling)
    [-p] (enable profiling)
    [-v] (show compiler output)"
    1>&2; exit 1; }

np=`cat /proc/cpuinfo | grep processor | wc -l`

rm -f snaps/*
cmake_line="cmake .."

while getopts ":dcopv" option; do
    case "${option}" in
        c)
            rm -rf build
            mkdir build
            ;;
        d)
            cmake_line="$cmake_line -DCMAKE_BUILD_TYPE=Debug"
            ;;
        o)
            cmake_line="$cmake_line -DADDITIONAL_OPTIMIZE=ON"
            ;;
        p)
            cmake_line="$cmake_line -DPROFILE=ON"
            ;;
        v)
            cmake_line="$cmake_line -DVERBOSE_MAKE=ON"
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