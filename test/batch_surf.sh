#!/bin/bash 
if [ -f results_${1}_surf_${2}.csv ]; then
  rm results_${1}_surf_${2}.csv
fi
touch results_${1}_surf_${2}.csv
echo 'structure_file,enthalpy,time' >> results_${1}_surf_${2}.csv
STRUC_DIR=out
for file in $(ls $STRUC_DIR)
do
./surf.out $STRUC_DIR/$file force_field_mixing_rules.def 298.0 12.0 $1 Xe $2 >> results_${1}_surf_${2}.csv 
done
