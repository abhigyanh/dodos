#!/bin/bash

rm -v confout.gro ener.edr md.log state.cpt traj_whole.trr traj.trr
cd dodos
	rm -v *.png
	rm -v DoS.txt entropy.txt log.txt
	rm -v box.xvg pos.xvg veloc.xvg
cd ..
