## About

1-2-3-dimensional simulation of wave and associated processes

## Build

```
./tools/build.sh -h
```

## Run

```
./tools/run.sh -h
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
