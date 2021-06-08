#!/bin/bash
filename=`basename $1`
atomsk $1 $2/$filename cif
sed -i '1s/^/\n/' $2/$filename
sed -i '1s/^/data_\n/' $2/$filename
