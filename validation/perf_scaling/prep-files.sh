#!/bin/bash

for n in {1..32}
do
    mkdir $n-threads
    cp -v -r run-files/* $n-threads/
done