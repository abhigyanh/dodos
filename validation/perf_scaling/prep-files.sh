#!/bin/bash

for n in {1..32}
do
    mkdir nt-$n
    cp -v -r run-files/* nt-$n/
done