#!/bin/bash

g++ -o2 macro_laguerre.cpp -lm -o macro_laguerre `root-config --cflags --glibs`

PX=100

PROCESS=0
for l in {1..3}			# RUN FOR RAFAEL
#for l in {0..0} 		# RUN FOR LUIS
do
	gnome-terminal -- bash -ic "./macro_laguerre $PROCESS $l 0 $PX;	./macro_laguerre $PROCESS $l 1 $PX;	pushover 'Estágio' 'Process $PROCESS has finished'; rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "
	let "PROCESS= $PROCESS + 1"
	gnome-terminal -- bash -ic "./macro_laguerre $PROCESS $l 2 $PX;	./macro_laguerre $PROCESS $l 3 $PX; pushover 'Estágio' 'Process $PROCESS has finished'; rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "
	let "PROCESS= $PROCESS + 1"
done
