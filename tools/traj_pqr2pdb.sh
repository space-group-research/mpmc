#!/bin/bash

#       Copyright (C) 2014 Christian Cioce, University of South Florida
#       This is free software. You may redistribute copies of it under the terms of
#       the GNU General Public License.

# NOTE: This is a painstakingly slow script for longer trajectories. It really should 
#       be rewritten for efficiency (Python might be a good choice), but for now it's
#       what it is.

[ $# -ne 1 ] && echo usage: $0 COORDINATE_FILE && exit 1				# Require coordinate file for processing

base=$(echo $1 | rev | cut -d. -f2- | rev)						# Strip the .pqr extension, preserving the filename
nCycles=$(grep -c END $1)								# Get the total number of iterations for output statistics
((i=0))

while read line; do
	[ "$line" == "ENDMDL" ] && ((i++)) && echo -ne "Cycle $i/$nCycles\r"			 # When we reach the end of a cycle, increment the counter and refresh the status output
	[ ! "`echo $line | cut -f 1 -d " "`" == "ATOM" ] && echo $line >> $base.pdb && continue	 # Only interested in lines containing coordinates
	xyz=($(echo $line | awk '{print $7" "$8" "$9}'))					 # Extract %11.6f x-, y- and z-coords
	x=$(echo ${xyz[0]}); y=$(echo ${xyz[1]}); z=$(echo ${xyz[2]})				 # Could have done this all in 1 line with awk
	x1=$(echo ${x:0:-3}); y1=$(echo ${y:0:-3}); z1=$(echo ${z:0:-3})			 # Strip off the last 3 digits, essentially converting %11.6f to %8.3f
	x2=$(echo -n $x | tail -c 3); y2=$(echo -n $y | tail -c 3); z2=$(echo -n $z | tail -c 3) # Keep those last 3 digits
	occ=$(echo $x2$y2)
	tFac=$z2

	# Now print a PDB file with the last 3 digits of x-, y- and z- as column data that VMD will read in
	printf "%s%8.3f%8.3f%8.3f   %s%s%s\n" "$(echo "$line" | awk '{print $0}' | head -c 30)" "$x1" "$y1" "$z1" "$x2" "$y2" "$z2" >> $base.pdb
done < $1

exit 0

