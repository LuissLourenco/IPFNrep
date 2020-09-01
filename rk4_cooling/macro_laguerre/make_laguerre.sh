#!/bin/bash

# NEEDS TO BE RUN WITH bash macro_laguerre.sh

g++ -o2 macro_laguerre.cpp -lm -o macro_laguerre `root-config --cflags --glibs`

PX=-2000
KDAMP=1.18E-8

PROCESS=0
#for l in {1..3}			# RUN FOR RAFAEL
for l in {0..0} 		# RUN FOR LUIS
do
	gnome-terminal -- bash -ic "./macro_laguerre $PROCESS $l 0 $PX $KDAMP;	./macro_laguerre $PROCESS $l 1 $PX $KDAMP; pushover 'Estágio' 'Process $PROCESS has finished'; rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "
	let "PROCESS= $PROCESS + 1"
	gnome-terminal -- bash -ic "./macro_laguerre $PROCESS $l 2 $PX $KDAMP;	./macro_laguerre $PROCESS $l 3 $PX $KDAMP; pushover 'Estágio' 'Process $PROCESS has finished'; rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "
	let "PROCESS= $PROCESS + 1"
done


