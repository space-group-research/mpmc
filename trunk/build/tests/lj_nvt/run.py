from subprocess import check_output
import os

check = -245
energy = float(check_output('../../mpmc test.input | grep "potential energy" | grep -o "= [-.0-9]*" | grep -o [-.0-9]* | tail -1',shell=True))

if (abs(energy-check)<=30):
	print os.getcwd().split(os.sep)[-1] + ": TEST PASSED"
else:
	print os.getcwd().split(os.sep)[-1] + ": TEST FAILED"
