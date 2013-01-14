#!/bin/bash

for file in energy.*; do
	echo "======================= $file ======================="
	sed 1,1d $file | awk '{ printf("%5d: %14.6f %14.6f %14.6f %14.6f\n", $9, $2/$9, $3/$9, $4/$9, $5/$9); }'
done
