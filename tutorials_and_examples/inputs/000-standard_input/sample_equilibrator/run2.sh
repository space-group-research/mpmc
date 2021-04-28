#!/bin/bash

myExe='/disk2/home/dfranz/mpmc/build/mpmc'

$myExe *.inp > runlog.log &
tail -f runlog.log | grep "poten" | awk {'print $5 $7'}
