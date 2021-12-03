source ~/.bashrc
#path=/sansom/s157/bioc1642/Desktop/prj_Notch/2/PCPS_21/CG2AT_2021-09-23_12-51-13/FINAL #J1, pdb 4cbz
#name=C2_J1
path=/sansom/s157/bioc1642/Desktop/prj_Notch/3/PCPS_14/CG2AT_2021-10-11_13-01-00/FINAL #J2apo, pdb 5mw7
name=C2_J2apo
pdb=final_cg2at_aligned.pdb
top=topol_final.top

pinoffset=8
ncpu=4

cd $path

#cp mdp files
min=/sansom/s157/bioc1642/Desktop/prj_Notch/AT/min.mdp
cp $min .
nvt=/sansom/s157/bioc1642/Desktop/prj_Notch/AT/nvt.mdp
cp $nvt .
npt=/sansom/s157/bioc1642/Desktop/prj_Notch/AT/npt.mdp
cp $npt .
prod=/sansom/s157/bioc1642/Desktop/prj_Notch/AT/prod.mdp
cp $prod .

## make index file
cat <<EOF>> index_input
13 | 14 | 20 | 21
name 24 LIP
q
EOF
gmx make_ndx -f final_cg2at_aligned.pdb -quiet < index_input

## nvt equilibration without minimisation
gmx grompp -f nvt.mdp -c $pdb -r $pdb -p $top -o nvt_$name.tpr -maxwarn 1 -quiet -n index.ndx
gmx mdrun -deffnm nvt_$name -v -quiet -pin on -pinoffset $pinoffset -ntomp $ncpu -ntmpi 1 -nsteps 5000

## check temperature
echo 16 0 | gmx energy -f nvt_$name.edr -o temperature.xvg -quiet
xmgrace temperature.xvg

## npt equilibration
gmx grompp -f npt.mdp -c nvt_$name.gro -r nvt_$name.gro -p $top -n index.ndx -o npt_$name.tpr -maxwarn 1 -quiet 
gmx mdrun -deffnm npt_$name -v -quiet -pin on -pinoffset $pinoffset -ntomp $ncpu -ntmpi 1 -nsteps 20000

## check pressure
echo 17 0 | gmx energy -f npt_$name.edr -o pressure.xvg -quiet
xmgrace pressure.xvg

## check density
echo 23 0 | gmx energy -f npt_$name.edr -o density.xvg -quiet
xmgrace density.xvg

## production run (with continuation)
gmx grompp -f prod.mdp -c npt_$name.gro -r npt_$name.gro -p $top -n index.ndx -o prod_$name.tpr -maxwarn 1 -quiet 
gmx mdrun -deffnm prod_$name -px prod_${name}_pullx.xvg -pf prod_${name}_pullf.xvg -cpi prod_$name.cpt -v -quiet -pin on -pinoffset $pinoffset -ntomp $ncpu -ntmpi 1 -nsteps 5000

## continue produciton run
cp ../../../../continue.sh .
sed -i -e "s/C2_J1/$name/g" continue.sh # change name in continue.sh file
nohup bash continue.sh > nohup_continue.out &

## check status
tail -f nohup_continue.out
htop
nvidia-smi -l

## center protein and fix periodic boundary conditions (optinally reduce traj with 10-fold with -dt 100 (1 frame per 100 ps instead of 1 per 10 ps)
#echo 1 0 | gmx trjconv -f prod_$name.xtc -s prod_$name.tpr -n index.ndx -pbc mol -ur compact -center -o prod_${name}_cent.xtc -quiet
echo 1 0 | gmx trjconv -f prod_$name.xtc -s prod_$name.tpr -n index.ndx -pbc mol -ur compact -center -dt 100 -o prod_${name}_cent.xtc -quiet

## clean up
rm \#*
	
cd ../../../..
