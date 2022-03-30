#!/bin/bash 
if [ -f results_${1}_surface.csv ]; then
  rm results_${1}_surface.csv
fi
touch results_${1}_surface.csv
echo 'structure_file,enthalpy,time' >> results_${1}_surface.csv
STRUC_DIR=out
for file in $(ls $STRUC_DIR)
do
./surface.out $STRUC_DIR/$file >> results_${1}_surface.csv 
done
