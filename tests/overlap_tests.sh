#!/bin/bash

source ./shared.sh

gridSize=199
overlaps=(9 13 20)
numThreads=5

build -fopenmp

for o in "${overlaps[@]}"
do
    echo "---${o}-Cell Overlap---"
    timeRun ${numThreads} ${gridSize} ${o}
done
