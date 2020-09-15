rm Outputs/*



a0=30
p0=0
ns=20
phi0=0




p=2
l=0
for phi0 in 0 45 90
do
	echo 'pl='$p$l' phi='$phi0
	mkdir 'Out_a0_'$a0'_phi0_'$phi0'_p0_'$p0'_pl_'$p$l
	./run_multiple $a0 $phi0 $p0 $p $l $ns
	mv Outputs/* 'Out_a0_'$a0'_phi0_'$phi0'_p0_'$p0'_pl_'$p$l
done


<<comment
p=0
l=0
for r0 in 1 3
do
	echo 'pl='$p$l' r0='$r0
	mkdir 'Out_a0_'$a0'_r0_'$r0'_p0_'$p0'_pl_'$p$l
	./run_multiple $a0 $r0 $p0 $p $l $ns
	mv Outputs/* 'Out_a0_'$a0'_r0_'$r0'_p0_'$p0'_pl_'$p$l
done
comment