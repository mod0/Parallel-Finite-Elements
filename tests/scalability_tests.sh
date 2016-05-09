#!/bin/bash

source ./shared.sh

gridSizes=(63 127 255 511)
overlapsPerGridSize=(1 1 1 1)
numThreads=(2 4 8 16 32)

numGridSizes=${#gridSizes}


for (( i=0; i<${numGridSizes}; i++ ));
do
    curGridSize=${gridSizes[$i]}

    echo "-----------------------------------------"
    echo "Starting tests with grid size ${curGridSize}"
    echo "-----------------------------------------"

    build
    echo "---Serial---"
    timeRun 1 ${gridSizes[$i]} 0

    build -fopenmp
    echo "---1 Thead---"
    timeRun 1 ${gridSizes[$i]} 0

    for t in "${numThreads[@]}"
    do
        echo "---${t} Theads---"
        timeRun ${t} ${curGridSize} ${overlapsPerGridSize[i]}
    done
done
