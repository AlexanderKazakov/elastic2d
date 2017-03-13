#!/usr/bin/env bash

while [ -e /proc/10089 ]; do
	echo "waiting ..."
	sleep 60
done
echo "finished!"

mv snapshots/vtk/ ~/work/vtk/vtk 
./tools/build.sh
./tools/run.sh -t titan

