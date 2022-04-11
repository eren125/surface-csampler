#!/bin/bash 
if [ -f results_${1}.csv ]; then
  rm results_${1}.csv
fi
touch results_${1}.csv
echo 'structure_file,enthalpy,henry,time' >> results_${1}.csv
STRUC_DIR=out
for file in $(ls $STRUC_DIR)
do
./a.out $STRUC_DIR/$file force_field_mixing_rules.def 298.0 12.0 $1 Xe >> results_${1}.csv 
done
