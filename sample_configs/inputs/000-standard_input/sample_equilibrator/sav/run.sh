#!/bin/bash

# ======================
#Space group 2016
# D. Franz
# Example MPMC equilibration script.
# ======================

# set energy diff. threshold
threshold="100.0" # energy units in K

# check for already existing final_energy.dat from a previous run.                                            
if [ -f ./final_energy.dat ]; then                                                                            
        previous_energy=$(cat final_energy.dat)     
	#echo 'Previous energy found.'                                                          
else
	previous_energy=0
	#echo 'No previous energy found. Setting to zero.'
fi  

# run MPMC. You should change directory to executable as needed.
myExe='/disk2/home/dfranz/mpmc/build/mpmc'
$myExe *.inp > runlog.log # change "| tee" to "> " to avoid seeing MPMC output.

# retrieve "final" energy from runlog after finished and write new file.
cat runlog.log | grep "potential energy" | tail -1 | awk '{print $5}' > final_energy.dat
final_energy=$(cat final_energy.dat)

#echo "final E: "$final_energy
#echo "previous E: "$previous_energy

# get energy difference magnitude
diff=$(echo "scale=6; sqrt(($final_energy - $previous_energy)^2)" | bc -l)

#echo "energy difference: "$diff
echo 'Energy: '$final_energy' ; Threshold: '$threshold' ; Diff: '$diff 


# check for equilibration.
if (( $(echo "$diff > $threshold" | bc -l) )); then
	mv *restart.pqr input.pqr
	bash run.sh
else
	echo 'Reached equilibration. Energy: '$final_energy' ; Threshold: '$threshold' ; Diff: '$diff                         
        exit $?
fi

