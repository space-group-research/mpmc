#!/bin/bash
#
#	Copyright (C) 2014 Christian Cioce, University of South Florida
#	This is free software. You may redistribute copies of it under the terms of
#	the GNU General Public License.
#

# TODO: Add functional support for main loops

# Error Check Input File for Correct Starting Conditions
[[ $# -lt 1 ]] && echo $0 usage: "[mpmc suft_fit input file]" && exit 1

# Declare function to draw sample from a uniform distribution over the half-open interval [low, high).
# *** NOTE: Requires python & numpy be installed on your system ***
getRand () {
npRand=`python << EOF
import numpy as np
print np.random.uniform(-0.5,0.5)
EOF`
}

extract () {
	awk '{print $1, $2}' ./fit > EOE.fit.dat
	awk '{print $1, $3}' ./fit > PAR.fit.dat
	awk '{print $1, $4}' ./fit > S.fit.dat
	awk '{print $1, $5}' ./fit > T.fit.dat
	awk '{print $1, $6}' ./fit > X.fit.dat
}

mpmc="mpmc"
sched=`grep "fit_schedule" $1 | awk '{print $2}'`
maxE=`grep "fit_start_temp" $1 | awk '{print $2}'`

[[ $sched != 0.999 ]] && sed -i "s/$sched/0.999/g" $1		# Force default schedule at 0.999
[[ $maxE -ne 50000 ]] && sed -i "s/$maxE/50000/g" $1		# Force default start temp at 50000 K

initPQR=`grep "pqr_input" $1 | awk '{print $2}'`
[[ `grep ^ATOM $initPQR | wc -l` -ne 10 ]] && echo "Input model $initPQR is not a 5-site model" && exit 1	# TODO: I really should make this entirely general, but whatev for now

if [ "$2" == "random" ]; then
	initEps=(`grep ^ATOM $initPQR | head -5 | awk '{print $13}' | uniq | grep -wv 0.00000`)
	initSig=(`grep ^ATOM $initPQR | head -5 | awk '{print $14}' | uniq | grep -wv 0.00000`)

	newEps=()
	newSig=()
	for i in "${initEps[@]}"; do getRand; nEVal=`echo "$i - ($i * $npRand)" | bc -l`; nEVal=`printf "%.5f" "$nEVal"`; newEps=("${newEps[@]}" "$nEVal"); done
	for i in "${initSig[@]}"; do getRand; nSVal=`echo "$i - ($i * $npRand)" | bc -l`; nSVal=`printf "%.5f" "$nSVal"`; newSig=("${newSig[@]}" "$nSVal"); done

	cp $initPQR $initPQR.0
	i=0
	for j in "${initEps[@]}"; do
		sed -i "s/$j/${newEps[$i]}/g" $initPQR
		((i++))
	done

	i=0
	for j in "${initSig[@]}"; do
		sed -i "s/$j/${newSig[$i]}/g" $initPQR
		((i++))
	done
fi

#################################
# ********** PART  I ********** #
#################################

count=1
while :
do
	echo "Starting Cycle: $count"
	~/$mpmc/build/mpmc $1 >& ./runlog.log

	if [ -e ./fit_geometry.pqr ]; then
		resDir="./$count-Results-999"
		mkdir $resDir
		[[ $? -ne 0 ]] && echo "ERROR: Cannot mkdir $resDir" && exit 1

		# Move Necessary Files
		mv ./runlog.log $resDir
		cp ./fit_geometry.pqr $resDir

		# Descend Into Directory
		cd $resDir

		# Perform Extraction
		tail -104 ./runlog.log | head -101 > ./fit
		extract

		# Ascend Back To The Main Directory
		cd ..

		# Prepare For Next Cycle
		mv ./fit_geometry.pqr ./$initPQR
		((count++))
	else
		# MPMC has found the "best" fits
		echo "fit_geometry.pqr not found..."
		break
	fi
done

#################################
# ********** PART II ********** #
#################################

echo; echo "Changing schedule (0.999 --> 0.9999) and start_temp (50,000 --> 25,000)"; echo

# Now that fitting @ start_temp= 50,000 K & schedule = 0.999 has been exhausted, half the start_temp & cool at a slower rate
sed -i "s/$sched/0.9999/g" $1
sed -i "s/$maxE/25000/g" $1

while :
do
        echo "Starting Cycle: $count"
        ~/$mpmc/build/mpmc $1 >& ./runlog.log

        if [ -e ./fit_geometry.pqr ]; then
		resDir="./$count-Results-9999"
                mkdir $resDir
		[[ $? -ne 0 ]] && echo "ERROR: Cannot mkdir $resDir" && exit 1

                # Move Necessary Files
                mv ./runlog.log $resDir
                cp ./fit_geometry.pqr $resDir

                # Descend Into Directory
                cd $resDir

                # Perform Extraction
                tail -104 ./runlog.log | head -101 > ./fit
                extract

                # Ascend Back To The Main Directory
                cd ..

                # Prepare For Next Cycle
                mv ./fit_geometry.pqr ./$initPQR
                ((count++))
        else
                # MPMC has found the "best" fits
                echo "fit_geometry.pqr not found...beginning to exhaustively reattempt"
                break
        fi
done

# Retry up to 2 times before quitting
cc=1

until [ $cc -eq 3 ]; do
        echo "Starting Cycle: $count        cc: $cc / 2"
        ~/$mpmc/build/mpmc $1 >& ./runlog.log

        if [ -e ./fit_geometry.pqr ]; then
		resDir="./$count-Results-9999"
                mkdir $resDir

                # Move Necessary Files
                mv ./runlog.log $resDir
                cp ./fit_geometry.pqr $resDir

                # Descend Into Directory
                cd $resDir

                # Perform Extraction
                tail -104 ./runlog.log | head -101 > ./fit
                extract

                # Ascend Back To The Main Directory
                cd ..

                # Prepare For Next Cycle
                mv ./fit_geometry.pqr ./$initPQR
                ((count++))

                # Reset the Ticker to 1
                cc=1
        else
                # Update Counter
                echo "fit_geometry.pqr not found...updating the counter"
                ((cc++))
        fi
done

echo
echo "Found the best possible fits after exhaustion!"

