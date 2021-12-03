#!/bin/sh
#

# --------------------------------------------------------------
# 
# This bash script runs simulations with C2 domains from notch ligand DLL4
##
# C2 + POPC:GM1:GM2(90:5:5) Bilayer
#
# C2 + POPC Bilayer
#
# C2 + POPC:DOPS (75:25) Bilayer
#
# CG MARTINI 2.1
# 
# Andreas Haahr Larsen
#
# March 2019
# 
# --------------------------------------------------------------

## ensure modules and other settings from batchrc file
source ~/.bashrc

## set parameters

###################### expect user input from here ##################################

# select protein (1=DLL4, 2=Jagged1, 3=Jagged2)
protein_folder=3

# set lipids to loop over (1=PC:GM1:GM2, 2=PC:PS:PIP2, 3=PC:PS, 4=PC, 5=PC:PS:GM1)
lip_first=3
lip_last=3

# set number of repetitions to run for each lip
rep_first=1
rep_last=24

# activate/deactivate overall modules of script
PREPARE=0
MDRUN=0
ANALYSIS=1

# submodules of ANALYSIS (only relevant if ANALYSIS=1)
ANALYSIS_CENT=0
ANALYSIS_DIST=0
ANALYSIS_DIST_COM=0
ANALYSIS_DENS=0
ANALYSIS_RZZ=1
ANALYSIS_RMSD=0

# parallelisation with mdrun
pinoffset=0
noCPUs=4
noGPUs=1

########################## to here ###################################################

# set parameters depending on selected protein
if [ $protein_folder -eq 1 ]
then
    protein_name=DLL4
    pdb_name=5mvx_t
    lip_for_ref=1
    key_frame_rep=0
    key_frame=0
    NEUTRAL_PROTEIN=1
elif [ $protein_folder -eq 2 ]
then
    protein_name=Jagged1
    pdb_name=4cbz_c2
    lip_for_ref=1
    key_frame_rep=0
    key_frame=0
    NEUTRAL_PROTEIN=1
elif [ $protein_folder -eq 3 ]
then
    protein_name=Jagged2apo
    pdb_name=5mw7_c2
    lip_for_ref=3
    key_frame_rep=0
    key_frame=0
    NEUTRAL_PROTEIN=1
fi

# define paths to scripts 
py2=/sansom/s157/bioc1642/anaconda2/bin/python2.7
py3=/sansom/s157/bioc1642/anaconda3/bin/python3.7
martinize=/sansom/s157/bioc1642/Desktop/Scripts/martinize_GROMACS_2018_plumed.py # edited line 1851 to have "/gromacs/top/" instead of "/top/"
insane=/sansom/s157/bioc1642/Desktop/Scripts/insane.py
cg2at=/sansom/s157/bioc1642/Desktop/Scripts/cg2at/cg2at-gmx2019-martini2.pl

# define paths to mdp files
min=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/minimization.mdp
eq=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/equilibration.mdp
prod=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/production.mdp

# define path to pdb file
protein=/sansom/s157/bioc1642/Desktop/prj_C2/Structures/${pdb_name}_clean.pdb # truncated, missing residues added

# define path to topology files
ffdir=/sansom/s157/bioc1642/Desktop/Scripts/martini2

# plumed files
plumed_file=/sansom/s157/bioc1642/Desktop/Scripts/plumed/dist_ang_restraint_${protein_name}.dat

# simulation time in ps to start collecting density maps
time=1996000

## create and go to new protein working directory
mkdir -p $protein_folder
cd $protein_folder

## loop over lipid compositions
for j in $(seq $lip_first $lip_last)
do

  # folder name
  if [ $j -eq 1 ]
  then
    Folderprefix=PCGM1GM3
  elif [ $j -eq 2 ]
  then
    Folderprefix=PCPSPIP2
  elif [ $j -eq 3 ]
  then
    Folderprefix=PCPS
  elif [ $j -eq 4 ]
  then
    Folderprefix=PC
  elif [ $j -eq 5 ]
  then
    Folderprefix=PCPSGM1  
  fi

  ## loop over replicas
  for i in $(seq $rep_first $rep_last)
  do
    
    ############ GENERATE FOLDERS ETC #######################################################
    
    ## create (if not already existing) and go to working directory
    folder=${Folderprefix}_$i
    mkdir -p $folder
    cd $folder
    
    ############ PREPARE TO RUN SIMULATION ###################################################
    if [ $PREPARE -eq 1 ]
    then
      ## rotate protein 
      gmx insert-molecules -ci $protein -o ${protein_name}_rot.pdb -nmol 1 -box 6 6 6 -rot xyz -quiet

      ## center protein
      gmx editconf -f ${protein_name}_rot.pdb -o ${protein_name}_rot_cent.pdb -c -quiet

      ## CG protein  with martinize 2.2
      $py3 $martinize -f ${protein_name}_rot_cent.pdb -x ${protein_name}_CG.pdb -o topol.top -v -ff martini22 -elastic -dssp dssp

      ## Build bilayer + protein with insane with 10% antifreeze
      box_xy=7
      box_z=18
      dist=4.4 # four times cutoff length (see mpd file)
      if [ $j -eq 1 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:90 -l DPG1:5 -l DPG3:5 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 2 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:80 -l POPS:15 -l POP2:5 -sol W:9 -sol WF:1 -salt 0 -dm $dist 
      elif [ $j -eq 3 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:75 -l DOPS:25 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 4 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 5 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:70 -l DOPS:25 -l DPG1:5 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      fi

      ## change typology file for system
      sed -i -e 's/#include "martini.itp"//g' topol.top
      cat << EOF > topol.add
#include "$ffdir/martini_v2.2.itp"
#include "$ffdir/martini_v2.0_ions.itp"
#include "$ffdir/martini_v2.0_lipids_all_201506.itp"
#include "Protein.itp"

#ifdef POSRES
#include "posre.itp"
#endif
EOF
      cat topol.add topol.top > tmp
      rm topol.add
      mv tmp topol.top

      ## minimization
      gmx grompp -f $min -c ${protein_name}_Bilayer.gro -p topol.top -o min.tpr -quiet
      gmx mdrun -deffnm min -quiet -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset

      ## make index file
      if [ $j -eq 1 ]
      then
        cat << EOF > index.input
16 | 17 | 18 | 19 | 23 | 24 | 25
name 26 SOL
13 | 14 | 15 | 20 | 21 | 22
name 27 LIP
a BB
14 | 21
name 29 GM1
15 | 22
name 30 GM3
q
EOF
      elif [ $j -eq 2 ]
      then
        cat << EOF > index.input
16 | 17 | 18 | 23 | 24 | 25
name 26 SOL
13 | 14 | 15 | 20 | 21 | 22
name 27 LIP
a BB
14 | 21
name 29 PS
15 | 22
name 30 PIP2
q
EOF
      elif [ $j -eq 3 ]
      then 
        cat << EOF > index.input
15 | 16 | 17 | 21 | 22 | 23
name 24 SOL
13 | 14 | 19 | 20
name 25 LIP
a BB
14 | 20
name 27 PS
q
EOF
      elif [ $j -eq 4 -a $NEUTRAL_PROTEIN -eq 1 ]
      then
        cat << EOF > index.input
14 | 15
name 16 SOL
13
name 17 LIP
a BB
q
EOF
      elif [ $j -eq 4 -a $NEUTRAL_PROTEIN -eq 0 ]
      then
        cat << EOF > index.input
14 | 15 | 16 | 19 | 20 | 21
name 22 SOL
13
name 23 LIP
a BB
q
EOF
      elif [ $j -eq 5 ]
      then
        cat << EOF > index.input
16 | 17 | 18 | 23 | 24 | 25
name 26 SOL
13 | 14 | 15 | 20 | 21 | 22
name 27 LIP
a BB
14 | 21
name 29 PS
15 | 22
name 30 GM1
q
EOF
      fi
      gmx make_ndx -f ${protein_name}_Bilayer.gro -quiet < index.input
      rm index.input

      ## see overview of index file 
      # gmx make_ndx -f ${protein_name}_Bilayer.gro -n index.ndx -quiet
  
      ## generate position restraint for protein
      echo 1 | gmx genrestr -f ${protein_name}_CG.pdb -fc 1000 1000 1000 -o posre.itp -quiet

      ## equilibrate, protein restrained, wall restraint
      cp $eq eq_dt.mdp
      if [ $j -le 2 ] || [ $j -eq 5 ] 
      then
        sed -i -e 's/dt                       = 0.03/dt                       = 0.02/g' eq_dt.mdp # lower dt from 0.03 to 0.02      
        sed -i -e 's/nsteps                   = 300000/nsteps                   = 500000/g' eq_dt.mdp
        sed -i -e 's/define                   = -DSTRONG_POSRES   ; Prevent protein from moving too much/define                   = -DSTRONG_POSRES -DFLEXIBLE   ; Prevent protein from moving too much and change some constraints in PIP to bonds/g' eq_dt.mdp # change some constrants in pip to bonds
        gmx grompp -f eq_dt.mdp -c min.gro -r min.gro -p topol.top -o eq.tpr -quiet -n index.ndx
      else
        gmx grompp -f eq_dt.mdp -c min.gro -r min.gro -p topol.top -o eq.tpr -quiet -n index.ndx -maxwarn 1
      fi
      gmx mdrun -deffnm eq -v -quiet -plumed $plumed_file -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset
      mv COLVAR COLVAR_equilibration

      ## prepare production run
      cp $prod production_short.mdp
      sed -i -e 's/nsteps                   = 86000000/nsteps                   = 57333333/g' production_short.mdp # lower production time from 3 to 2 us
      gmx grompp -f production_short.mdp -c eq.gro -r eq.gro -p topol.top -o md.tpr -quiet -n index.ndx
    
    # end PREPARE if statement
    fi
    
    ############ RUN SIMULATION ##############################################################
    if [ $MDRUN -eq 1 ]
    then
    	gmx mdrun -deffnm md -v -quiet -plumed $plumed_file -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset
    fi
    
    ############ ANALYSIS #####################################################################
    if [ $ANALYSIS -eq 1 ]
    then

      ## center protein in box      
      if [ $ANALYSIS_CENT -eq 1 ]
      then
        echo 1 0 | gmx trjconv -s md.tpr -f md.xtc -o md_cent.xtc -pbc mol -center -ur compact -quiet
      fi
      
      ## calculate min and max distances from protein (reference) to membrane (selection)
      if [ $ANALYSIS_DIST -eq 1 ]
      then
        # global distances
        gmx pairdist -f md.xtc -n index.ndx -tu ns -type min -ref 1 -sel LIP -o dist_min.xvg -quiet
        gmx pairdist -f md.xtc -n index.ndx -tu ns -type max -ref 1 -sel LIP -o dist_max.xvg -quiet
      
        # residue distances (between each residue and closest lipid, format: column1: time, column2:res1-lip-dist, columnN+1=resN-lip-dist)
        gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel LIP -o dist_res_min.xvg -quiet
        gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type max -ref 1 -refgrouping res -sel LIP -o dist_res_max.xvg -quiet
	gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel GM1 -o dist_res_min_GM1.xvg -quiet
	gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel GM3 -o dist_res_min_GM3.xvg -quiet
      fi
      
      ## calculate COM distance from protein to membrane
      if [ $ANALYSIS_DIST_COM -eq 1 ]
      then
	gmx distance -f md_cent.xtc -s md.tpr -n index.ndx -select "com of group LIP plus com of group Protein" -tu ns -oall dist_com_cent -quiet
      fi

      ## calculate 2D density maps in xy plane 
      if [ $ANALYSIS_DENS -eq 1 ]
      then
        echo 1  | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PROT -quiet
        echo 13 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PC -quiet
        if [ $j -eq 1 ]
        then
          echo 14 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PS -quiet
          echo 15 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PIP2 -quiet
	  echo 16 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PIP3 -quiet
        fi
        if [ $j -eq 2 ]
        then
          echo 14 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PS -quiet
	  echo 15 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PIP2 -quiet
        fi
        if [ $j -eq 3 ]
        then
          echo 14 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PS -quiet
        fi
      fi
      
      ## make reference structure (tpr file) for rotmat 
      if [ $key_frame -eq 0 ] 
      then
        #using the last frame of the first rep
        if [ $j -eq $lip_for_ref -a $i -eq 0 ] 
        then
          gmx grompp -f $prod -c md.gro -r md.gro -p topol.top -o md_ref.tpr -quiet -n index.ndx
          mv md_ref.tpr ../.
        fi
      else  
        # using the frame with closest bound residue (key frame)
        # generate key frame index file
        if [ $j -eq $lip_for_ref -a $i -eq $key_frame_rep ] 
        then
          cat << EOF > frame_index.ndx
[ frames ]

$key_frame

EOF
          # extract key frame from trajectory
          echo 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -o key_frame.gro -fr frame_index.ndx -quiet
	  
	  # generate tpr file from key frame
          gmx grompp -f $prod -c key_frame.gro -r key_frame.gro -p topol.top -o md_ref.tpr -quiet -n index.ndx

	  # center, for vizualization
	  echo 1 0 | gmx trjconv -s md_ref.tpr -f key_frame.gro -o key_frame_cent.gro -pbc mol -center -ur compact -quiet
	
	  # copy tpr file to parent directory
          mv md_ref.tpr ../.

        fi
      fi

      ## calculate rotation matrix of all frames with respect to reference structure (../md_ref.tpr)
      if [ $ANALYSIS_RZZ -eq 1 ]
      then
	echo 1 | gmx trjconv -f md.xtc -s md.tpr -o md_prot.xtc -quiet
	echo 1 | gmx trjconv -f md.gro -s md.tpr -o md_prot.gro -quiet # for visualization in vmd
        echo 1 | gmx rotmat -s ../md_ref.tpr -f md_prot.xtc -fitxy -o Rzz.xvg -quiet
      fi

      ## calculate RMSD
      if [ $ANALYSIS_RMSD -eq 1 ]
      then
        BB=$(($sel + 1))
        echo $BB $BB | gmx rms -f eq.xtc -s eq.tpr -tu ns -quiet -o rmsd_eq.xvg -n index.ndx
        echo $BB $BB | gmx rms -f md.xtc -s md.tpr -tu ns -quiet -o rmsd_md.xvg -n index.ndx
      fi

      ## make pdb of system
      echo 0 | gmx trjconv -f md.gro -s md.tpr -o md.pdb -quiet

    # end ANALYSIS if statement
    fi

    ############ FINISH #####################################################################

    ## clean up
    rm \#*

    ## navigate back to protein working directory
    cd ..

  #end loop over replicas
  done

#end loop over lipid compositions
done

## navigate back to parent directory
cd ..
