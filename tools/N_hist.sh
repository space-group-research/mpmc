#!/bin/bash
#
# Spit out a histogram of N for a gcmc trajectory
#
# Jon Belof
# Space Research Group
# Department of Chemistry
# University of South Florida


# usage
if [ $# -ne 1 ];
then
        echo usage: $0 "[trajectoryfile]"
        exit 1
fi

# check the cmd line argument
file=$1
if [ ! -e $file ];
then
        echo "$0: couldn't access the file $file"
        exit 1
fi

for i in `seq 0 500`
do
	bin=`cat $file | grep "moveable_molecules" | awk '{print $3}' | grep -x "moveable_molecules=$i" | wc -l`
	echo "$i $bin"
done

