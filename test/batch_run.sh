#!/bin/bash 
if [ -f results.txt ]; then
  rm results.txt
fi
touch results.txt
sed -i '1s/^/structure_name,enthalpy,time(s)\n/' results.txt
STRUC_DIR=out
for file in $(ls $STRUC_DIR)
do
./a.out $STRUC_DIR/$file force_field_mixing_rules.def 298.0 12.0 100 >> results.txt 
done
