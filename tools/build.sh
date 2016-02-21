#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-o] (enable additional optimization)
    [-d] (to compile in Debug mode)
    [-c] (to clean all before compiling)
    [-v] (show make output)
    [-g] (add -g flag)"
    1>&2; exit 1; }

np=`cat /proc/cpuinfo | grep processor | wc -l`

rm -f snaps/*
cmake_line="cmake .."

while getopts ":dcopvg" option; do
    case "${option}" in
        c)
            rm -rf build CMakeCache.txt CMakeFiles/ cmake_install.cmake
            mkdir build
            ;;
        d)
            cmake_line="$cmake_line -DCMAKE_BUILD_TYPE=Debug"
            ;;
        o)
            cmake_line="$cmake_line -DADDITIONAL_OPTIMIZE=ON"
            ;;
        v)
            cmake_line="$cmake_line -DVERBOSE_MAKE=ON"
            ;;
        g)
            cmake_line="$cmake_line -DDEBUGG=ON"
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
