#!/usr/bin/env bash

rm -rf tools/ndi/sim/png
# nodesPerMm = 10 for appropriate testing quality, = 20 for report quality
nodesPerMm=$1
recalculate=$2
if [ -z $recalculate ]; then
        recalculate=0;
fi
if [ $recalculate -eq 1 ]; then rm -rf tools/ndi/sim/txt; fi

declare -a depths=(
# depths from 2016:
#	"0.00"	"5.03"	"5.14"	"5.10"	"5.21"	"6.45"	"9999"	"6.48"	"5.47"	"2.29"	"2.96"	"4.20"	"1.74"	"9999	"2.41"	"3.19"	"3.69"
# depths from 2017:
    "0.00" "3.19" "3.83" "4.60" "3.19" "3.20" "4.44" "4.63" "9999" "1.51" "1.85" "1.71" "1.54" "1.79" "9999"
)

counter=0
for depth in "${depths[@]}"
do
	echo "Doing $counter'th for G1-G2: $depth"
	resfilename="$counter.depth=$depth"
	./tools/ndi/ndi.sh $depth $nodesPerMm $recalculate $resfilename
	counter=$[1 + $counter]
done

cd tools/ndi/report
pdflatex report.tex
cd -
