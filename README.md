## About

Two-dimensional simulation of wave processes in elastic body. Rectangular area, structured grid. Grid-characteristic method.

## Build

base build
```
mkdir build
cd build
cmake ..
make
```

build with google tests, gprof profiling and in Debug mode (by default, all this options are disabled)
```
mkdir build
cd build
cmake .. -DCOMPILE_TESTS=ON -DPROFILE=ON -DCMAKE_BUILD_TYPE=Debug
make
```

## Run

```
mpirun -np <NUMBER OF PROCESSES> ./elastic2d
```

## Run tests

```
mpirun -np <NUMBER OF PROCESSES> ./elastic2d_mpi_tests
mpirun -np 1 ./elastic2d_tests
```

