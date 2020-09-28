#!/bin/bash

rm Outputs/*

if [ $# -ne 5 ] ; then
	exit 1
fi

a0=$1
r0=$2
p0=$3
p=$4
l=$5
ns=10
mkdir 'Out_a0_'$a0'_r0_'$r0'_p0_'$p0'_pl_'$p$l
./run_multiple $a0 $r0 $p0 $p $l $ns
mv Outputs/* 'Out_a0_'$a0'_r0_'$r0'_p0_'$p0'_pl_'$p$l
