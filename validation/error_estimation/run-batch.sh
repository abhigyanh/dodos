#!/bin/bash

for n in {1..10}
do
    cd run-$n
        gmx mdrun -v -s example.tpr
        bash dump-from-trajectory.sh
        cd dodos
            dodos --temperature 298 --threads $(nproc)
            rm pos.xvg veloc.xvg box.xvg
        cd ..
    cd ..
done
