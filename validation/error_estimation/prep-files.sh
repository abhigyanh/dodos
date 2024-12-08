#!/bin/bash

for n in {1..10}
do
    mkdir run-$n
    cp -v -r run-files/* run-$n/
done