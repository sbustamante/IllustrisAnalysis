#!/bin/bash

FILE=$1
NGAS=$2
NDM=$3
NOTH=$4

echo "=========================================================================="
echo "Shifting FOF file"
echo "=========================================================================="

#Adding gas particles to FOF catalog
if [ $NGAS -ne 0 ]; then
    echo $(($NDM+$NGAS+$NOTH)) > fof.grp.temp
    for (( i=1; i<=$NGAS; i++ )); do
	echo "0" >> fof.grp.temp
    done
    
    tail -n +2 fof.grp >> fof.grp.temp
    
    for (( i=1; i<=$NOTH; i++ )); do
	echo "0" >> fof.grp.temp
    done
    
    mv fof.grp.temp $FILE

fi 
