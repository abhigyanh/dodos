#!/bin/bash

echo "# num_processes velocity_decomposition_time(sec)"

for n in {1..32}; do
    logfile=nt-$n/dodos/log.txt
    perf_line=$(grep "ms per atom" $logfile)

    echo "$n $perf_line"
done