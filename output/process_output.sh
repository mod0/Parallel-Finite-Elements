#!/bin/bash

# This script combines the solution data for different subdomains in one and
# prepends a header
shopt -s extglob

# Glob for the raw subdomain files
readonly RAW_DATA=subdomain+([0-9]).csv.+([0-9])

# Glob for concatenated subdomain files
readonly CAT_DATA=data.csv.+([0-9])

# Header to prepend to concatenated files
readonly HEADER="x, y, Temperature"

# remove the old files
rm -f ${CAT_DATA}

for f in ${RAW_DATA}; do
    timestep=$(echo ${f} | grep -o "[0-9]\+$")
    cat $f >> data.csv.${timestep}
done

# Prepend the header
sed -i "1i${HEADER}\n" ${CAT_DATA}

# Move the raw data into a directory called raw
rm -rf raw
mkdir raw
mv ${RAW_DATA} raw
