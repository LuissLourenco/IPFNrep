
rm Outputs/*


if [ $# -eq 0 ] ; then
	./run_multiple
else
	a0=$1
	r0=$2
	p0=$3
	p=$4
	l=$5
	ns=20
	./run_multiple $a0 $r0 $p0 $p $l $ns
fi