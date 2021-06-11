#!/bin/bash 
STRUC_DIR=/home/emmanuel/Simulations/share/raspa/structures/cif/
for file in $(ls $STRUC_DIR)
do
bash convert_cif.sh $STRUC_DIR/$file out/
done
