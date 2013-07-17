#!/bin/bash

[ $# -ne 1 ] && echo usage: $0 COORDINATE_FILE && exit 1				# Require coordinate file for processing

while read line; do
	[ ! "`echo $line | cut -f 1 -d " "`" == "ATOM" ] && echo $line && continue	# Only interested in lines containing coordinates
	xyz=($(echo $line | awk '{print $7" "$8" "$9}'))				# Extract %11.6f x-, y- and z-coords
	x=$(echo ${xyz[0]}); y=$(echo ${xyz[1]}); z=$(echo ${xyz[2]})			# Could have done this all in 1 line with awk
	x1=$(echo ${x:0:-3}); y1=$(echo ${y:0:-3}); z1=$(echo ${z:0:-3})		# Strip off the last 3 digits, essentially converting %11.6f to %8.3f
	x2=$(echo -n $x | tail -c 3); y2=$(echo -n $y | tail -c 3); z2=$(echo -n $z | tail -c 3) # Keep those last 3 digits
	occ=$(echo $x2$y2)
	tFac=$z2

	# Now print a PDB file with the last 3 digits of x-, y- and z- as column data that VMD will read in
	printf "%s%8.3f%8.3f%8.3f   %d%d%d\n" "$(echo "$line" | awk '{print $0}' | head -c 30)" "$x1" "$y1" "$z1" "$x2" "$y2" "$z2"
done < $1

