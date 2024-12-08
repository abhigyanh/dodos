#!/bin/bash

# First use 'gmx trjconv' to make all molecules whole
echo "System" | gmx trjconv -f traj.trr -s example.tpr -o traj_whole.trr -pbc whole 

# Then use 'gmx traj' to dump coordinates and velocities in numpy-readable format.
mkdir ./dodos
gmx traj -f traj_whole.trr -s example.tpr -n not_name_MW.ndx -ov ./dodos/veloc.xvg -ox ./dodos/pos.xvg -ob ./dodos/box.xvg

