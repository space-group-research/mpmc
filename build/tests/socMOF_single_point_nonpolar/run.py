from subprocess import check_output
import os

check = -70447.13152
energy = float(check_output('../../mpmc test.input | grep "potential energy" | grep -o [-.0-9]*',shell=True))

if (abs(energy-check)<=0.001):
	print os.getcwd().split(os.sep)[-1] + ": TEST PASSED"
else:
	print os.getcwd().split(os.sep)[-1] + ": TEST FAILED"
