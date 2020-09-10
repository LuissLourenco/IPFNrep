#!/bin/bash

# NEEDS TO BE RUN WITH bash macro_laguerre.sh

g++ -o2 macro_laguerre.cpp -lm -o macro_laguerre `root-config --cflags --glibs`

PX=-2000
KDAMP=1.18E-8

#for l in {1..3}			# RUN FOR RAFAEL
#for l in {0..0} 		# RUN FOR LUIS
#do
#	gnome-terminal -- bash -ic "./macro_laguerre $PROCESS $l 0 $PX $KDAMP;	./macro_laguerre $PROCESS $l 1 $PX $KDAMP; pushover 'Estágio' 'Process $PROCESS has finished'; rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "
#	let "PROCESS= $PROCESS + 1"
#	gnome-terminal -- bash -ic "./macro_laguerre $PROCESS $l 2 $PX $KDAMP;	./macro_laguerre $PROCESS $l 3 $PX $KDAMP; pushover 'Estágio' 'Process $PROCESS has finished'; rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "
#	let "PROCESS= $PROCESS + 1"
#done

#PROCESS=0
#gnome-terminal -- bash -ic "\
#./macro_laguerre $PROCESS 3 3 $PX $KDAMP;\
#./macro_laguerre $PROCESS 1 1 $PX $KDAMP;\
#pushover 'Estágio' 'Process $PROCESS has finished';\
#rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "

PROCESS=0
gnome-terminal -- bash -ic "\
./macro_laguerre $PROCESS 0 0 -2000 0;\
./macro_laguerre $PROCESS 0 1 -2000 0;\
./macro_laguerre $PROCESS 0 2 -2000 0;\
rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "

PROCESS=1
gnome-terminal -- bash -ic "\
./macro_laguerre $PROCESS 1 0 -2000 0;\
./macro_laguerre $PROCESS 1 1 -2000 0;\
./macro_laguerre $PROCESS 1 2 -2000 0;\
rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "

<<COMMENT1

PROCESS=2
gnome-terminal -- bash -ic "\
./macro_laguerre $PROCESS 2 0 -2000 0;\
./macro_laguerre $PROCESS 2 1 -2000 0;\
./macro_laguerre $PROCESS 2 2 -2000 0;\
rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "

PROCESS=3
gnome-terminal -- bash -ic "\
./macro_laguerre $PROCESS 0 0 -2000 1.18E-8;\
./macro_laguerre $PROCESS 0 1 -2000 1.18E-8;\
./macro_laguerre $PROCESS 0 2 -2000 1.18E-8;\
rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "

PROCESS=4
gnome-terminal -- bash -ic "\
./macro_laguerre $PROCESS 1 0 -2000 1.18E-8;\
./macro_laguerre $PROCESS 1 1 -2000 1.18E-8;\
./macro_laguerre $PROCESS 1 2 -2000 1.18E-8;\
rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "

PROCESS=5
gnome-terminal -- bash -ic "\
./macro_laguerre $PROCESS 2 0 -2000 1.18E-8;\
./macro_laguerre $PROCESS 2 1 -2000 1.18E-8;\
./macro_laguerre $PROCESS 2 2 -2000 1.18E-8;\
rm InputToBatch$PROCESS.txt; rm Out$PROCESS.txt; "

COMMENT1
