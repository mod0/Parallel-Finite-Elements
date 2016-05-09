#!/bin/bash

source ./shared.sh

gridSize=99
overlaps=(17 18 19 20 21 22)
numThreads=2

build -fopenmp

for o in "${overlaps[@]}"
do
    echo "---${o}-Cell Overlap---"
    timeRun ${numThreads} ${gridSize} ${o}
done
