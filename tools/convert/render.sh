#!/usr/bin/env bash

# append three images from three folders
# and create video from result


for i in {0..70}; do
	onename=$(printf "detector/%02d.png" "$i") 
	another=$(printf "zaxis/%02d.png" "$i")
	yetanother=$(printf "3d/%02d.png" "$i")
	firstname=$(printf "first%02d.png" "$i")
	resname=$(printf "%02d.png" "$i")
	convert +append $onename $another $firstname
	convert -append $firstname $yetanother $resname
	echo $resname
	rm $firstname
done

avconv -framerate 7/1 -i %02d.png -vf "scale=1464:1366" res.mp4
