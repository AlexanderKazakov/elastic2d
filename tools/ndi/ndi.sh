#!/usr/bin/env bash

savefoldername="tools/ndi/sim"
thickness=$1
filenameforreport=$thickness
if [[ -n "$4" ]]; then filenameforreport=$4; fi
echo $filenameforreport
savefilenamepng="$savefoldername/png/$filenameforreport.png"
savefilenametxt="$savefoldername/txt/$filenameforreport.txt"
mkdir -p  tools/ndi/sim/png tools/ndi/sim/txt

nodesPerMm=$2
numSteps=$[75*$nodesPerMm]
filename=$(printf "snapshots/detector/mesh1core00snap%04d.txt" "$numSteps")

recalculate=$3
if [ -z $recalculate ]; then
	recalculate=0;
fi

if [ $recalculate -eq 1 ]; then
	./tools/build.sh -t ndi
	rm -rf snapshots/vtk snapshots/detector snapshots/zaxis
	rm -f *.log
	mkdir -p snapshots/vtk snapshots/detector snapshots/zaxis
	./build/ndi $thickness $nodesPerMm
	cp $filename $savefilenametxt
fi

gnuplot -persist <<-EOFMarker
set grid
set title "G1-G2: $thickness"
set xlabel "time, us"
set ylabel "sensor value"
unset key

#plot '$savefilenametxt' using ((\$1)*1e+6):(abs(\$2)*1e+8) w l lw 2 lt -1
set xrange [6.5 : 16.5]

#plot '$savefilenametxt' using 0:(abs(\$2)*1e+8) w l lw 2 lt -1
#set xrange [GPVAL_DATA_X_MIN+0.1*(GPVAL_DATA_X_MAX-GPVAL_DATA_X_MIN) : ]

set yrange [0.5 :]
#plot '$savefilenametxt' using ((\$1)*1e+6):(abs(\$2)*1e+8) w l lw 2 lt -1
#plot '$savefilenametxt' using 0:(abs(\$2)*1e+8) w l lw 2 lt -1


set terminal pngcairo enhanced
set output '$savefilenamepng'
plot '$savefilenametxt' using ((\$1)*1e+6):(abs(\$2)*1e+8) w l lw 2 lt -1
#plot '$savefilenametxt' using 0:(abs(\$2)*1e+8) w l lw 2 lt -1
EOFMarker


