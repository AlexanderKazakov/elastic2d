#!/usr/bin/env bash

# append images from folders
# and create video from result

rm -f ndi/*.png ndi/*.mp4

for i in {0..180}; do
	
	velocity=$(printf "ndi/3d/velocity/step.%04d.png" "$i")
	pressure=$(printf "ndi/3d/pressure/step.%04d.png" "$i")
	slice=$(printf "ndi/zaxis/snap%04d.txt.png" "$((2*i))")
	detector=$(printf "ndi/detector/mesh1core00snap%04d.txt.png" "$((2*i))")
	
	one=$(printf "ndi/one%04d.png" "$i")
	two=$(printf "ndi/two%04d.png" "$i")
	res=$(printf "ndi/res%04d.png" "$i")
	convert +append $slice $detector $one
	convert +append $velocity $pressure $two
	convert -append $one $two $res
	
	echo $res
	rm $one $two
done

avconv -framerate 7/1 -i ndi/res%04d.png -vf "scale=1464:1366" ndi/res.mp4
