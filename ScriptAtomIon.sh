#!/bin/bash

BOUNDSTATES=3

# ––– DATA DIRECTORY NAME ––– 
AIDIR=${BOUNDSTATES}bs_1200xPoints_12channels_test

# ––– NAME OF PARAMETER FILE ––– 
PARFILE=AtomIon1D_${BOUNDSTATES}bs.par


# —————————————————————————————————————————————————————
# things below this line should not need to be modified

NORMAL=$(tput sgr0)
RED=$(tput setaf 1)
BLUE=$(tput setaf 4)
GREEN=$(tput setaf 2)
# colors for printed text


if [ ! -d "../data_local/$AIDIR" ]
then
    mkdir ../data_local/$AIDIR
else 
    printf "\nSubdirectory %s already exists.\nOverwrite? [Type number, then press ENTER]:\n" "${BLUE}$AIDIR${NORMAL}"
    select yn in "Yes" "No"; 
    do 
      case $yn in
          Yes ) echo "${BLUE}Overwriting...${NORMAL}"; 
                sleep 1; 
                break;;
          No  ) echo "${RED}Change name of AtomIon subdirectory.${NORMAL}"; 
                exit 1;;
      esac
    done
fi
# creates subdirectory if it does not exist
# asks you to overwrite if already exists.

check=0
i=1
while [ $check == 0 ]
do
    CHECKDIR=./workdir${i}
    if [ -d "$CHECKDIR" ]
    then
        i=$((i+1))
    else
        check=1
    fi
done

WORKDIR=workdir${i}
mkdir ./$WORKDIR
echo "Working directory ${BLUE}$i${NORMAL} created;"

make
# compiles AtomIon1D.f90

if [ $? -ne 0 ]
then
    echo "${RED}AtomIon1D does not compile.${NORMAL}"
    exit 1
fi

cp ./AtomIon1D.x ./$WORKDIR
cp ./Legendre.dat ./$WORKDIR

cp ./input_Atomion/$PARFILE ./$WORKDIR

cd ./$WORKDIR
pwd

sleep 1

./AtomIon1D.x < $PARFILE
# runs AtomIon1D.x

if [ $? -ne 0 ]
then
    echo "${RED}Some problems occurred while running AtomIon1D.${NORMAL}"
else
    echo "${GREEN}Radial code complete;${NORMAL}"
fi

mv $PARFILE ../../data_local/$AIDIR/
mv DiabaticEnergies.dat ../../data_local/$AIDIR/
mv AdiabaticEnergies.dat ../../data_local/$AIDIR/
mv VQmat.dat ../../data_local/$AIDIR/
mv Pmat.dat ../../data_local/$AIDIR/
# mv QuickEffV.dat ../../data_local/$AIDIR/
# mv QuickPMat.dat ../../data_local/$AIDIR/
# mv QuickQtil.dat ../../data_local/$AIDIR/
mv QTilmat.dat ../../data_local/$AIDIR/
mv ParAndThresh.dat ../../data_local/$AIDIR/
# moves all output files to data directory


if [ $? -ne 0 ]
then
    echo "${RED}Some problem occurred while moving the files.${NORMAL}"
    exit 1
else 
    echo "Output moved subdirectory ${BLUE}$AIDIR${NORMAL} in ${BLUE}data_local${NORMAL};"
    cd ..
    rm -rf $WORKDIR
    echo "Working directory deleted;" 
fi

