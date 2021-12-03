name=C2_J1
pinoffset=8
ncpu=4

gmx mdrun -deffnm prod_$name -px prod_${name}_pullx.xvg -pf prod_${name}_pullf.xvg -cpi prod_$name.cpt -v -quiet -pin on -pinoffset $pinoffset -ntomp $ncpu -ntmpi 1

