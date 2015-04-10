#!/bin/bash

#Simulation
SIMULATION="SPH_20_064_11"
#Snap
SNAP=${SIMULATION}"_102"
#Snapbase
SNAPBASE="$HOME/data/Simulations/$SIMULATION/$SNAP"
#Number of snapshots
SNAPFILES=1
#Output filename
OUTPUT="temp.tmp"
#Type of particles  [ 0-all      1-gas      2-DM ]
TYPE=2
#Linking lenght
LINKING=0.2
#Minimum number of considered particles for FOF
MINIMPART=30
#Only dark matter in gadget 1 format   #[ 0-False	1-True ]
ONLYDM=1
#Shift gas particles 	#[ 0-False	1-True ]
SHIFTGAS=1

#Axis for cut
AXIS=0
#Coordinate of slide
SLIDE=10000
#Thickness of slide
DX=1000
#Jump in sampling
SAMPLING=1


#Running Cutter [ 1 ]
if [ $1 -eq 1 ]; then
    rm ./Cutter.out
    make Cutter 
    #Running code
    time ./Cutter.out $SNAPBASE $SNAPFILES $AXIS $SLIDE $DX $SAMPLING $OUTPUT $TYPE

#Running Ascii [ 2 ]
elif [ $1 -eq 2 ]; then
    rm ./Ascii.out
    make Ascii
    #Running code
    time ./Ascii.out $SNAPBASE $SNAPFILES $SAMPLING $OUTPUT $TYPE
    
#Running Halos [ 3 ]
elif [ $1 -eq 3 ]; then
    rm ./Halos.out
    make Halos
    #Running code
    time ./Halos.out $SNAPBASE $SNAPFILES $OUTPUT $TYPE $LINKING $MINIMPART $SHIFTGAS

#Running DomainIdentifier [ 4 ]
elif [ $1 -eq 4 ]; then
    #Recompliling
    rm ./g2tog1.out 
    make g2tog1
    cd ./domain_identifier/
    rm ./Domain_identifier.x
    make Domain_identifier
    cd -
    
    #Converting G2 format to G1
    time ./g2tog1.out $SNAPBASE $SNAPFILES $ONLYDM
    
    #Calculating Catalog of Halos
    bash run.sh 3
    cp fof.grp "$SNAPBASE.FullSnap.gad1.grp"
    
    #Running DomainIdentifier
    printf "\n\n\n===================================================="
    printf "\nDOMAIN IDENTIFIER\n"
    printf "====================================================\n\n\n"
    time mpiexec -n 1 ./domain_identifier/Domain_identifier.x "$SNAPBASE.FullSnap.gad1" ./domain_identifier/parameters.dat 1
    #Deleting tmp file
    rm "$SNAPBASE.FullSnap.gad1"
    mv "fof.grp" "../data/$SIMULATION/$SNAP.members"
    #Creating copies of produced data
    cp "$SNAPBASE.FullSnap.gad1.parts.rescue" "../data/$SIMULATION/$SNAP.domain"
    rm $SNAPBASE.FullSnap.gad1.*
    cp "Halo_Catalog.dat" "../data/$SIMULATION/$SNAP.halo_catalog"
    rm "Halo_Catalog.dat"
    cp "Mass_function_fof.dat" "../data/$SIMULATION/$SNAP.mass_function_fof"
    rm "Mass_function_fof.dat"
    #Creating Ascii copy of all particles
    time ./Ascii.out $SNAPBASE $SNAPFILES $SAMPLING $OUTPUT 0
    mv "temp.tmp" "../data/$SIMULATION/part.dat"
    
fi