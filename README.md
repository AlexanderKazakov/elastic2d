## About

1-2-3-dimensional simulation of wave and associated processes

## Build

```
./tools/build.sh [-c]
```
```
-c  clear all before build
```

## Run

```
./tools/run.sh [-n] [-p]
```
```
-n  number of processes (default is number of cores on the machine)
-p  run Paraview after calculation (only in one-process mode)
```

## Run tests

On MPI connections
```
mpirun -np <NUMBER OF PROCESSES> ./build/gcm_mpi_tests
```

On other functionality
```
./build/gcm_tests
```
