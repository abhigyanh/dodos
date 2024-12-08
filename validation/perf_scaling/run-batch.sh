#!/bin/bash

for n in {1..32}
do (
    cd $n-nt
        gmx mdrun -v -s example.tpr
        bash dump-from-trajectory.sh
        cd dodos
            dodos --temperature 298 --threads $(($(nproc)/2))
        cd ..
    cd ..
) &
done
wait