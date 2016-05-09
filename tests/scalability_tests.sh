#!/bin/bash

source ./shared.sh

gridSizes=(63 127 255 511)
numThreads=(2 4 8 16 32)

numGridSizes=${#gridSizes}


for gs in "${gridSizes[@]}"
do

    echo "-----------------------------------------"
    echo "Starting tests with grid size ${gs}"
    echo "-----------------------------------------"

    build
    echo "---Serial---"
    timeRun 1 ${gs} 0

    build -fopenmp
    echo "---1 Thead---"
    timeRun 1 ${gs} 0

    for t in "${numThreads[@]}"
    do
        ((overlap = 1 + gs / (2 * t)))
        echo "---${t} Theads with Overlap of ${overlap} Cells---"
        timeRun ${t} ${gs} ${overlap}
    done
done
