#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-t] (to compile tests)
    [-d] (to compile in Debug mode)
    [-c] (to clean all before compiling)"
    1>&2; exit 1; }

rm -f snaps/*
cmake_line="cmake .."

while getopts ":tdc" o; do
    case "${o}" in
        c)
            rm -rf build
            mkdir build
            ;;
        t)
            cmake_line="$cmake_line -DCOMPILE_TESTS=ON"
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
make
cd ..