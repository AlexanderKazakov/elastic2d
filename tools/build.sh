#!/usr/bin/env bash

usage() { echo "Usage: $0
    [-c] (to clean all before compiling and before the very first build)
    [-d] (to compile in Debug mode)
    [-r] (to compile in Release mode)
    [-o] (enable additional optimization)
    [-v] (show make output)"
    1>&2; exit 1; }


np=`cat /proc/cpuinfo | grep processor | wc -l`

rm -f snaps/*
cmake_line="cmake .."
build_type="RelWithDebInfo"

while getopts ":cdrov" option; do
    case "${option}" in
        c)
            rm -rf build CMakeCache.txt CMakeFiles/ cmake_install.cmake
            mkdir build
            ;;
        d)
            build_type="Debug"
            ;;
        r)
            build_type="Release"
            ;;
        o)
            cmake_line="$cmake_line -DADDITIONAL_UNSAFE_OPTIMIZE=ON"
            ;;
        v)
            cmake_line="$cmake_line -DVERBOSE_MAKE=ON"
            ;;
        *)
            usage
            ;;
    esac
done

cmake_line="$cmake_line -DCMAKE_BUILD_TYPE=$build_type"

cd build
eval $cmake_line
make -j$np
cd ..
